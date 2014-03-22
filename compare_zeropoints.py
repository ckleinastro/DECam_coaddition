from pylab import *
from matplotlib.pylab import *
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
from matplotlib.ticker import MaxNLocator
from numpy import array, loadtxt, where

import sys

coadds_dir = "/Volumes/Extra_HDD/Research/DECam_Data/coadds/chip21_coadds_and_data/"

rcParams['axes.unicode_minus'] = False
matplotlib.rc('font', family="serif")
matplotlib.rc('font', serif="Times New Roman")
matplotlib.rc('xtick', labelsize=12)
matplotlib.rc('ytick', labelsize=12)
matplotlib.rc('axes', labelsize=14)

standard_field_data = loadtxt("zp_outliers_cal.dat")

zpt = standard_field_data[:,0][standard_field_data[:,8]==1]
ciz = standard_field_data[:,1][standard_field_data[:,8]==1]
cam = standard_field_data[:,2][standard_field_data[:,8]==1]
cdx = standard_field_data[:,3][standard_field_data[:,8]==1]
cdy = standard_field_data[:,4][standard_field_data[:,8]==1]
rms = standard_field_data[:,5][standard_field_data[:,8]==1]
std_dates = standard_field_data[:,7][standard_field_data[:,8]==1].astype('|S8')


relative_zp_filename = "102-C21_relZP_list.txt"

relative_zp_type=dtype([('image_path', '|S150'),
                 ('airmass', 'float'),
                 ('date-obs', '|S26'),
                 ('psf_zp', 'float'),
                 ('psf_zp_err', 'float'),
                 ('aper_zp', 'float'),
                 ('aper_zp_err', 'float')])

relative_zp_data = loadtxt(coadds_dir + relative_zp_filename, dtype=relative_zp_type)

rel_dates = []
rel_airmasses = []
rel_psf_zps = []
rel_psf_zp_errs = []
rel_aper_zps = []
rel_aper_zp_errs = []
image_paths = []
for rel_entry in relative_zp_data:
    filepath = rel_entry["image_path"]
    night_str = filepath.split("done/")[1].split("/C")[0]
    if night_str in std_dates:
        image_paths.append(filepath)
        rel_dates.append(night_str)
        rel_airmasses.append(rel_entry["airmass"])
        rel_psf_zps.append(rel_entry["psf_zp"])
        rel_psf_zp_errs.append(rel_entry["psf_zp_err"])
        rel_aper_zps.append(rel_entry["aper_zp"])
        rel_aper_zp_errs.append(rel_entry["aper_zp_err"])
        
rel_dates = array(rel_dates)
rel_airmasses = array(rel_airmasses)
rel_psf_zps = array(rel_psf_zps)
rel_psf_zp_errs = array(rel_psf_zp_errs)
rel_aper_zps = array(rel_aper_zps)
rel_aper_zp_errs = array(rel_aper_zp_errs)


rel_compared_data = []
night_in_nov = []
night_strs = []
filepath_strs = []

std_zpt_anchor = 0

for n in range(len(std_dates)):
    night_str = std_dates[n]
    if night_str in rel_dates:
        rel_indices = where(rel_dates==night_str)[0]
        for m in rel_indices:
            rel_compared_data.append([zpt[n], ciz[n], cam[n], cdx[n], cdy[n], rms[n],
                                      rel_airmasses[m], 
                                      rel_psf_zps[m], rel_psf_zp_errs[m], 
                                      rel_aper_zps[m], rel_aper_zp_errs[m]])
            night_strs.append(night_str)
            filepath_strs.append(image_paths[m])
            night_in_nov.append(int(night_str[-2:]))
            if rel_psf_zps[m] == 0 and rel_aper_zps[m] == 0:
                std_zpt_anchor = zpt[n]
rel_compared_data = array(rel_compared_data)

"""
0   zpt
1   ciz
2   cam
3   cdx
4   cdy
5   rms
6   airmass
7   psf_zp
8   psf_zp_err
9   aper_zp
10  aper_zp_err
"""


translated_std_zpts = rel_compared_data[:,0] - std_zpt_anchor


fig = plt.figure(figsize = (3.3, 2.5))
ax1 = subplot(111)

# ax1.errorbar(night_in_nov, translated_std_zpts, rel_compared_data[:,5], linestyle="none", marker="o", ms=3, color="black")

ax1.errorbar(night_in_nov, rel_compared_data[:,7]-translated_std_zpts, rel_compared_data[:,8], linestyle="none", marker="o", ms=3, color="red", label=r"${\rm ZP}_{\rm rel-psf} - {\rm ZP}_{\rm std}$")

# ax1.errorbar(night_in_nov, rel_compared_data[:,9]-translated_std_zpts, rel_compared_data[:,10], linestyle="none", marker="o", ms=3, color="blue", label=r"${\rm ZP}_{\rm rel-aper} - {\rm ZP}_{\rm std}$")

# This code draws major and minor tick lines. Major ticks get number labels.
majorLocator_x = MultipleLocator(2)
minorLocator_x = MultipleLocator(1)
ax1.xaxis.set_major_locator(majorLocator_x)
ax1.xaxis.set_minor_locator(minorLocator_x)

majorLocator_y1 = MultipleLocator(0.1)
minorLocator_y1 = MultipleLocator(0.05)
ax1.yaxis.set_major_locator(majorLocator_y1)
ax1.yaxis.set_minor_locator(minorLocator_y1)

ax1.set_xlabel(r"Night of Observation")
ax1.set_ylabel(r"Relative Zeropoint") 

ax1.hlines(0, 0, 19, colors="k", linestyles="--")
ax1.legend(loc="lower right", fontsize=8, numpoints=1)
ax1.set_xlim(0, 19)

# pos =         [left, bottom, width, height]
ax1.set_position([0.17, 0.180, 0.79, 0.80])

canvas = FigureCanvas(fig)
canvas.print_figure("zeropoint_comparison.pdf" , dpi=300)
close("all")