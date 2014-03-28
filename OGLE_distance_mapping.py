from scipy import *
from scipy import stats
from scipy import interpolate
from scipy import mgrid
from pylab import *
from matplotlib.pylab import *
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
from matplotlib.ticker import MaxNLocator
import matplotlib.cm as cm
import sys
from mayavi import mlab

rcParams['axes.unicode_minus'] = False
matplotlib.rc('font', family="serif")
matplotlib.rc('font', serif="Times New Roman")
matplotlib.rc('xtick', labelsize=12)
matplotlib.rc('ytick', labelsize=12)
matplotlib.rc('axes', labelsize=14)



ogle_data = []

ogle_rrl_cat = file("ogle_III_LMC_RRLs.txt", "r")
# ID	                Field	 StarID	RA	        Decl	    Type	I	    V	    P_1	        dP_1	    T0_1	    A_1	    R21_1	phi21_1	R31_1	phi31_1	P_2	    dP_2	T0_2	A_2	    R21_2	phi21_2	R31_2	phi31_2	ID_OGLE_II	ID_MACHO	ID_GCVS	ID_OTHER	Remarks
for line in ogle_rrl_cat:
    if line[0] != "#":
        star_id = line.split()[0]
        star_ra = float(line.split()[3])*15
        star_dec = float(line.split()[4])
        star_type = line.split()[5]
        star_I = float(line.split()[6])
        star_V = float(line.split()[7])
        star_per = float(line.split()[8])
        ogle_data.append((star_ra, star_dec, star_id, star_type, star_V, star_I, star_per))
ogle_rrl_cat.close()

ogle_datatype = dtype([('ra', 'float'),
                     ('dec', 'float'),
                     ('ogle_id', '|S20'),
                     ('type', '|S4'),
                     ('V_mag', 'float'),
                     ('I_mag', 'float'),
                     ('period', 'float')])
ogle_data = array(ogle_data, dtype=ogle_datatype)



vmag_mask = ogle_data["V_mag"]>0
imag_mask = ogle_data["I_mag"]>0

P_0 = 0.52853966619770265

rrl_fundamental_periods = where(ogle_data["type"]=="RRc", 10**(log10(ogle_data["period"])+0.127), ogle_data["period"])

rrl_log_normalized_fundamental_periods = log10(rrl_fundamental_periods/P_0)



V_calibrated_M_0 = 0.4319
V_calibrated_M_0_err = 0.0184
I_calibrated_M_0 = 0.1065
I_calibrated_M_0_err = 0.0380


def fit_line(x, y):
    slope, intercept, r, prob2, see = stats.linregress(x, y)
    mx = x.mean()
    sx2 = ((x-mx)**2).sum()
    sd_intercept = see * sqrt(1./len(x) + mx*mx/sx2)
    sd_slope = see * sqrt(1./sx2)
    res = (intercept + slope*x) - y

    slope, intercept, r, prob2, see = stats.linregress(x[abs(res)<1*res.std()], y[abs(res)<1*res.std()])
    mx = x[abs(res)<1*res.std()].mean()
    sx2 = ((x[abs(res)<1*res.std()]-mx)**2).sum()
    sd_intercept = see * sqrt(1./len(x) + mx*mx/sx2)
    sd_slope = see * sqrt(1./sx2)

    return intercept, sd_intercept, slope, sd_slope, abs(res)<1*res.std()



V_M_0, V_M_0_err, V_alpha, V_alpha_err, V_fit_mask = fit_line(rrl_log_normalized_fundamental_periods[vmag_mask], ogle_data["V_mag"][vmag_mask])
I_M_0, I_M_0_err, I_alpha, I_alpha_err, I_fit_mask = fit_line(rrl_log_normalized_fundamental_periods[imag_mask], ogle_data["I_mag"][imag_mask])

V_mu = V_M_0 - V_calibrated_M_0
V_mu_err = hypot(V_calibrated_M_0_err, V_M_0_err)

I_mu = I_M_0 - I_calibrated_M_0
I_mu_err = hypot(I_calibrated_M_0_err, I_M_0_err)


def scatter2fract_dist(scatter_mag):
    upper = ((10**((scatter_mag+5)/5))/10)-1
    lower = 1-((10**((-scatter_mag+5)/5))/10)
    return (upper + lower)/2

def compute_distance(mu, mu_err, A=0, A_err=0):
    effective_mu = mu-A
    effective_mu_err = (mu_err**2 + A_err**2)**0.5
    dist = 10**(effective_mu/5. + 1)
    dist_fract_err = scatter2fract_dist(effective_mu_err)
    dist_err = dist_fract_err * dist
    return dist, dist_err

V_dist, V_dist_err =  compute_distance(V_mu, V_mu_err)
I_dist, I_dist_err =  compute_distance(I_mu, I_mu_err)




V_distance_mus = ogle_data["V_mag"][vmag_mask][V_fit_mask] - (V_calibrated_M_0 + V_alpha*rrl_log_normalized_fundamental_periods[vmag_mask][V_fit_mask])

V_distance_mu_errs = hypot(V_calibrated_M_0_err, V_alpha_err*rrl_log_normalized_fundamental_periods[vmag_mask][V_fit_mask])
V_distance_pc = []
V_distance_pc_err = []
for n in range(len(V_distance_mus)):
    dist, dist_err = compute_distance(V_distance_mus[n], V_distance_mu_errs[n])
    V_distance_pc.append(dist)
    V_distance_pc_err.append(dist_err)
V_ra = ogle_data["ra"][vmag_mask][V_fit_mask]
V_dec = ogle_data["dec"][vmag_mask][V_fit_mask]
V_distance_pc = array(V_distance_pc)/1000.
V_distance_pc_err = array(V_distance_pc_err)/1000.

I_distance_mus = ogle_data["I_mag"][imag_mask][I_fit_mask] - (I_calibrated_M_0 + I_alpha*rrl_log_normalized_fundamental_periods[imag_mask][I_fit_mask])

I_distance_mu_errs = hypot(I_calibrated_M_0_err, I_alpha_err*rrl_log_normalized_fundamental_periods[imag_mask][I_fit_mask])
I_distance_pc = []
I_distance_pc_err = []
for n in range(len(I_distance_mus)):
    dist, dist_err = compute_distance(I_distance_mus[n], I_distance_mu_errs[n])
    I_distance_pc.append(dist)
    I_distance_pc_err.append(dist_err)
I_ra = ogle_data["ra"][imag_mask][I_fit_mask]
I_dec = ogle_data["dec"][imag_mask][I_fit_mask]
I_distance_pc = array(I_distance_pc)/1000.
I_distance_pc_err = array(I_distance_pc_err)/1000.





sky_pixel_size = 0.2    # deg
distance_pixel_size = sin(radians(sky_pixel_size))*I_dist/1000.

dec_grid = arange(-73.0, -66.0, sky_pixel_size)
ra_factor = cos(radians((dec_grid.max() + dec_grid.min())/2.0))

ra_grid = arange(65.0, 95.0, sky_pixel_size/ra_factor)[::-1]

distance_grid = arange(I_distance_pc.min()-distance_pixel_size/2.0, I_distance_pc.max()+distance_pixel_size/2.0, distance_pixel_size)

density_array = zeros((len(distance_grid)-1, len(dec_grid)-1, len(ra_grid)-1)) # distance, dec, ra

density_array_3d = zeros(( len(ra_grid)-1, len(dec_grid)-1, len(distance_grid)-1 )) # ra, dec, distance

for r in range(len(ra_grid)-1):
    ra_high_bound = ra_grid[r]
    ra_low_bound = ra_grid[r+1]
    ra_mask = (I_ra>ra_low_bound) & (I_ra<ra_high_bound)
    for d in range(len(dec_grid)-1):
        dec_low_bound = dec_grid[d]
        dec_high_bound = dec_grid[d+1]
        ra_dec_mask = (ra_mask) & (I_dec>dec_low_bound) & (I_dec<dec_high_bound)
        for m in range(len(density_array)-1):
            dist_low_bound = distance_grid[m]
            dist_high_bound = distance_grid[m+1]
            stars_in_box = (stats.norm.cdf(dist_high_bound, I_distance_pc[ra_dec_mask], I_distance_pc_err[ra_dec_mask]) - stats.norm.cdf(dist_low_bound, I_distance_pc[ra_dec_mask], I_distance_pc_err[ra_dec_mask])).sum()
            density_array[m][d][r] = stars_in_box
            density_array_3d[r][d][m] = stars_in_box

max_density = density_array.flatten().max()
colorbar_saturation_value = max_density*1.0

"""
grid_ra, grid_dec, grid_dist = mgrid[95.0:65.0:52j, -73.0:-66.0:34j, 47.357210344325154:65.250851213998004:92j]

mlab.contour3d(grid_ra, grid_dec, grid_dist, density_array_3d, contours=[0.1], opacity=0.2,
    extent=(ra_grid[-1], ra_grid[0], dec_grid[0], dec_grid[-1], distance_grid[0], distance_grid[-1]) )
"""

sys.exit()





mgrid[95.0:65.0:53j, -73.0:-66.0:35j, 47.357210344325154:65.250851213998004:93]




interpolate.griddata(points, values, xi, method='linear', fill_value=0.0)





fig = plt.figure(figsize = (6.6, 6.0))
ax1 = subplot(111)
# extent=(left, right, bottom, top)
lmc_map = ax1.imshow(density_array.sum(axis=0), origin="lower", interpolation="nearest", 
    extent=(ra_grid.max(), ra_grid.min(), dec_grid.min(), dec_grid.max()),
    aspect=1.0/ra_factor, cmap=cm.gist_yarg, vmin=0, vmax=density_array.sum(axis=0).flatten().max())


cbaxes = fig.add_axes([0.05, 0.10, 0.90, 0.05]) 
cb = plt.colorbar(lmc_map, cax=cbaxes, orientation="horizontal")  
cbaxes.set_xlabel(r"RRL Density ($%.1f^\circ$ by $%.1f^\circ$ by $%.1f$ kpc)" % (sky_pixel_size, sky_pixel_size, distance_pixel_size))

# This code draws major and minor tick lines. Major ticks get number labels.
majorLocator_x = MultipleLocator(5)
minorLocator_x = MultipleLocator(1)
ax1.xaxis.set_major_locator(majorLocator_x)
ax1.xaxis.set_minor_locator(minorLocator_x)

majorLocator_y1 = MultipleLocator(2)
minorLocator_y1 = MultipleLocator(.2)
ax1.yaxis.set_major_locator(majorLocator_y1)
ax1.yaxis.set_minor_locator(minorLocator_y1)

ax1.set_ylabel(r"Declination (deg)", labelpad=2)
ax1.set_xlabel(r"Right Ascension (deg)")

# pos =         [left, bottom, width, height]
ax1.set_position([0.08, 0.15, 0.90, 0.90])

canvas = FigureCanvas(fig)
canvas.print_figure("lmc_maps/ogle_map_sum.jpg", dpi=300)
close("all")



for slice in range(len(distance_grid)-1):
    slice_dist = (distance_grid[slice] + distance_grid[slice+1])/2.0

    fig = plt.figure(figsize = (6.6, 6.0))
    ax1 = subplot(111)
    # extent=(left, right, bottom, top)
    lmc_map = ax1.imshow(density_array[slice], origin="lower", interpolation="nearest", 
        extent=(ra_grid.max(), ra_grid.min(), dec_grid.min(), dec_grid.max()),
        aspect=1.0/ra_factor, cmap=cm.gist_yarg, vmin=0, vmax=colorbar_saturation_value)


    cbaxes = fig.add_axes([0.05, 0.10, 0.90, 0.05]) 
    cb = plt.colorbar(lmc_map, cax=cbaxes, orientation="horizontal")  
    cbaxes.set_xlabel(r"RRL Density ($%.1f^\circ$ by $%.1f^\circ$ by $%.1f$ kpc)" % (sky_pixel_size, sky_pixel_size, distance_pixel_size))

    # This code draws major and minor tick lines. Major ticks get number labels.
    majorLocator_x = MultipleLocator(5)
    minorLocator_x = MultipleLocator(1)
    ax1.xaxis.set_major_locator(majorLocator_x)
    ax1.xaxis.set_minor_locator(minorLocator_x)

    majorLocator_y1 = MultipleLocator(2)
    minorLocator_y1 = MultipleLocator(.2)
    ax1.yaxis.set_major_locator(majorLocator_y1)
    ax1.yaxis.set_minor_locator(minorLocator_y1)

    ax1.set_ylabel(r"Declination (deg)", labelpad=2)
    ax1.set_xlabel(r"Right Ascension (deg)")

    ax1.text(0.05, 0.875, "%.1f kpc" % (slice_dist), horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes, fontsize=12)

    # pos =         [left, bottom, width, height]
    ax1.set_position([0.08, 0.15, 0.90, 0.90])

    canvas = FigureCanvas(fig)
    canvas.print_figure("lmc_maps/ogle_map_%03d.jpg" % slice, dpi=300)
    close("all")

