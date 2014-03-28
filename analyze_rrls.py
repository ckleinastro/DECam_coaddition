from scipy import loadtxt, dtype, log10, where, polyfit, stats, linspace, hypot
import socket
import sys

from pylab import *
from matplotlib.pylab import *
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
from matplotlib.ticker import MaxNLocator

from numpy import invert

rcParams['axes.unicode_minus'] = False
matplotlib.rc('font', family="serif")
matplotlib.rc('font', serif="Times New Roman")
matplotlib.rc('xtick', labelsize=12)
matplotlib.rc('ytick', labelsize=12)
matplotlib.rc('axes', labelsize=14)


####    Declare the query_field and query_chip.
# if len(sys.argv) != 3:
#     print "Improper usage. As arguments, supply the field number and chip ID."
#     print "Ex: coadd_images.py 100 C13"
#     sys.exit()
# query_field = sys.argv[1]
# query_chip = sys.argv[2]

# query_field = "110"
# query_chip = "C21"



""" Local laptop testing and development modification. """
if socket.gethostname() == "Christopher-Kleins-MacBook-Pro.local" or           \
    socket.gethostname()[:8] == "airbears":
    coadd_directory = "/Volumes/Extra_HDD/Research/DECam_Data/coadds/chip21_coadds_and_data"
else:
    coadd_directory = "/global/scratch2/sd/cenko/DES/Jan2014/coadds/"
""" End of local laptop testing and development modification. """

chip_num_list = range(1, 61)
chip_num_list.append(62)

rrl_datatype = dtype([('ra', 'float'),
                     ('dec', 'float'),
                     ('ogle_id', '|S20'),
                     ('type', '|S4'),
                     ('period', 'float'),
                     ('V_mag', 'float'),
                     ('I_mag', 'float'),
                     ('decam_z_mag', 'float'),
                     ('decam_z_magerr', 'float'),
                     ('sdss_z_mag', 'float'),
                     ('sdss_z_magerr', 'float')])
rrl_data_flag = False
for chip_num in chip_num_list:
    chip_str = "C%02d" % (chip_num,)
    for field in range(30):
        field_str = "1" + "%02d" % (field,)
        rrl_datafile = coadd_directory + "/" + field_str + "-" + chip_str      \
                                                     + "_psf_coadd_rrl_data.txt"
        if not rrl_data_flag:
            try:
                rrl_data = loadtxt(rrl_datafile, dtype=rrl_datatype)
                rrl_data_flag = True
            except:
                continue
        else:
            try:
                rrl_data = append(rrl_data, loadtxt(rrl_datafile, dtype=rrl_datatype))
            except:
                continue




vmag_mask = rrl_data["V_mag"]>0
imag_mask = rrl_data["I_mag"]>0

P_0 = 0.52853966619770265

rrl_fundamental_periods = where(rrl_data["type"]=="RRc",  10**(log10(rrl_data["period"])+0.127), rrl_data["period"])

rrl_log_normalized_fundamental_periods = log10(rrl_fundamental_periods/P_0)

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


V_calibrated_M_0 = 0.4319
V_calibrated_M_0_err = 0.0184
I_calibrated_M_0 = 0.1065
I_calibrated_M_0_err = 0.0380
z_calibrated_M_0 = 0.5406
z_calibrated_M_0_err = 0.0539

V_M_0, V_M_0_err, V_alpha, V_alpha_err, V_fit_mask = fit_line(rrl_log_normalized_fundamental_periods[vmag_mask], rrl_data["V_mag"][vmag_mask])
I_M_0, I_M_0_err, I_alpha, I_alpha_err, I_fit_mask = fit_line(rrl_log_normalized_fundamental_periods[imag_mask], rrl_data["I_mag"][imag_mask])
z_M_0, z_M_0_err, z_alpha, z_alpha_err, z_fit_mask = fit_line(rrl_log_normalized_fundamental_periods, rrl_data["sdss_z_mag"])

V_mu = V_M_0 - V_calibrated_M_0
V_mu_err = hypot(V_calibrated_M_0_err, V_M_0_err)

I_mu = I_M_0 - I_calibrated_M_0
I_mu_err = hypot(I_calibrated_M_0_err, I_M_0_err)

z_mu = z_M_0 - z_calibrated_M_0
z_mu_err = hypot(z_calibrated_M_0_err, z_M_0_err)

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
z_dist, z_dist_err = compute_distance(z_mu, z_mu_err)


V_plot_offset = -4.0
I_plot_offset = -1.6

fig = plt.figure(figsize = (3.3, 2.5))
ax1 = subplot(111)

xgrid = linspace(-0.4, 0.4, 1000)
logp_xgrid = log10(P_0*(10**xgrid))

ax1.errorbar(log10(rrl_fundamental_periods[vmag_mask][invert(V_fit_mask)]), rrl_data["V_mag"][vmag_mask][invert(V_fit_mask)]+V_plot_offset,
    linestyle="none", marker="x", ms=3, color="blue", alpha=0.3)

ax1.errorbar(log10(rrl_fundamental_periods[imag_mask][invert(I_fit_mask)]), rrl_data["I_mag"][imag_mask][invert(I_fit_mask)]+I_plot_offset,
    linestyle="none", marker="x", ms=3, color="green", alpha=0.3)

ax1.errorbar(log10(rrl_fundamental_periods[invert(z_fit_mask)]), rrl_data["sdss_z_mag"][invert(z_fit_mask)], rrl_data["sdss_z_magerr"][invert(z_fit_mask)],
    linestyle="none", marker="x", ms=3, color="red", alpha=0.3)

ax1.errorbar(log10(rrl_fundamental_periods[vmag_mask][V_fit_mask]), rrl_data["V_mag"][vmag_mask][V_fit_mask]+V_plot_offset,
    linestyle="none", marker="s", ms=3, color="blue", label=r"$V-%.1f$" % (-1*V_plot_offset), alpha=0.3)

ax1.errorbar(log10(rrl_fundamental_periods[imag_mask][I_fit_mask]), rrl_data["I_mag"][imag_mask][I_fit_mask]+I_plot_offset,
    linestyle="none", marker="s", ms=3, color="green", label=r"$I-%.1f$" % (-1*I_plot_offset), alpha=0.3)

ax1.errorbar(log10(rrl_fundamental_periods[z_fit_mask]), rrl_data["sdss_z_mag"][z_fit_mask], rrl_data["sdss_z_magerr"][z_fit_mask],
    linestyle="none", marker="s", ms=3, color="red", label=r"$z_{\rm SDSS}$", alpha=0.3)


ax1.plot(logp_xgrid, V_M_0 + V_alpha*xgrid + V_plot_offset, color="k")
ax1.plot(logp_xgrid, I_M_0 + I_alpha*xgrid + I_plot_offset, color="k")
ax1.plot(logp_xgrid, z_M_0 + z_alpha*xgrid, color="k")


ax1.legend(loc="upper left", ncol=1, fontsize=6, numpoints=1, handletextpad=0.05, columnspacing=0.2)

# This code draws major and minor tick lines. Major ticks get number labels.
majorLocator_x = MultipleLocator(0.1)
minorLocator_x = MultipleLocator(0.05)
ax1.xaxis.set_major_locator(majorLocator_x)
ax1.xaxis.set_minor_locator(minorLocator_x)

majorLocator_y1 = MultipleLocator(2)
minorLocator_y1 = MultipleLocator(1)
ax1.yaxis.set_major_locator(majorLocator_y1)
ax1.yaxis.set_minor_locator(minorLocator_y1)

ax1.set_xlabel(r"$\log_{10}\left({\rm Period [d]}\right)$", labelpad=2)
ax1.set_ylabel(r"$m$ $+$ offset")

ax1.text(0.03, 0.64, r"$V_{P=%.3f}=%.3f$" % (P_0, V_M_0), horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes, fontsize=8)

ax1.text(0.03, 0.38, r"$I_{P=%.3f}=%.3f$" % (P_0, I_M_0), horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes, fontsize=8)

ax1.text(0.03, 0.10, r"$z_{P=%.3f}=%.3f$" % (P_0, z_M_0), horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes, fontsize=8)

ax1.set_ylim(21.0, 13.8)
# ax1.set_xlim(-0.620, 0.011)
ax1.set_xlim(-0.5, -0.08)

# pos =         [left, bottom, width, height]
ax1.set_position([0.18, 0.190, 0.79, 0.79])

canvas = FigureCanvas(fig)
canvas.print_figure("plr.pdf", dpi=300)
close("all")

