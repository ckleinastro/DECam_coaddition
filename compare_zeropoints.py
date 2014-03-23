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

# Strip out the standard data for 20121109, this night never fits well, probably wasn't photometric
standard_field_data = standard_field_data[standard_field_data[:,7]!=20121109.0]

zpt = standard_field_data[:,0][standard_field_data[:,8]==1]
ciz = standard_field_data[:,1][standard_field_data[:,8]==1]
cam = standard_field_data[:,2][standard_field_data[:,8]==1]
cdx = standard_field_data[:,3][standard_field_data[:,8]==1]
cdy = standard_field_data[:,4][standard_field_data[:,8]==1]
rms = standard_field_data[:,5][standard_field_data[:,8]==1]
std_dates = standard_field_data[:,7][standard_field_data[:,8]==1].astype('|S8')

for field in range(30):
    field_str = "1" + "%02d" % (field,)
    relative_zp_filename = field_str + "-C21_relZP_list.txt"

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

    translation_constant = -1*median(rel_compared_data[:,7] - rel_compared_data[:,0])

    translated_std_zpts = rel_compared_data[:,0] - translation_constant


    fig = plt.figure(figsize = (3.3, 4.5))
    ax1 = subplot(211)

    ax2 = subplot(212)

    # ax1.errorbar(night_in_nov, translated_std_zpts, rel_compared_data[:,5], linestyle="none", marker="o", ms=3, color="black")

    ax1.errorbar(night_in_nov, rel_compared_data[:,7]-translated_std_zpts, rel_compared_data[:,8], linestyle="none", marker="o", ms=3, color="red", label=r"${\rm ZP}_{\rm rel-psf} - {\rm ZP}_{\rm std}$")

    # ax1.errorbar(night_in_nov, rel_compared_data[:,9]-translated_std_zpts, rel_compared_data[:,10], linestyle="none", marker="o", ms=3, color="blue", label=r"${\rm ZP}_{\rm rel-aper} - {\rm ZP}_{\rm std}$")

    # This code draws major and minor tick lines. Major ticks get number labels.
    majorLocator_x = MultipleLocator(2)
    minorLocator_x = MultipleLocator(1)
    ax1.xaxis.set_major_locator(majorLocator_x)
    ax1.xaxis.set_minor_locator(minorLocator_x)

    majorLocator_y1 = MultipleLocator(0.05)
    minorLocator_y1 = MultipleLocator(0.025)
    ax1.yaxis.set_major_locator(majorLocator_y1)
    ax1.yaxis.set_minor_locator(minorLocator_y1)

    # ax1.set_xlabel(r"Night of Observation")
    ax1.set_ylabel(r"${\rm ZP}_{\rm rel-psf} - {\rm ZP}_{\rm std} + %.3f$" % translation_constant, labelpad=-2) 

    ax1.hlines(0, 0, 19, colors="k", linestyles="--")
    # ax1.legend(loc="lower right", fontsize=8, numpoints=1)
    ax1.set_xlim(0, 19)
    ax1.set_ylim((rel_compared_data[:,7]-translated_std_zpts).min()-0.06, (rel_compared_data[:,7]-translated_std_zpts).max()+0.06)

    ax2.plot(night_in_nov, rel_compared_data[:,6], marker="s", ms=3, linestyle="none", color="green")

    ax2.set_xlabel(r"Night of Observation")
    ax2.set_ylabel(r"Airmass") 

    ax2.set_ylim(rel_compared_data[:,6].min()-0.015, rel_compared_data[:,6].max()+0.015)

    ax2.set_xlim(0, 19)

    ax2.xaxis.set_major_locator(majorLocator_x)
    ax2.xaxis.set_minor_locator(minorLocator_x)

    majorLocator_y2 = MultipleLocator(0.04)
    minorLocator_y2 = MultipleLocator(0.02)
    ax2.yaxis.set_major_locator(majorLocator_y2)
    ax2.yaxis.set_minor_locator(minorLocator_y2)

    # pos =         [left, bottom, width, height]
    ax1.set_position([0.18, 0.380, 0.76, 0.60])
    ax2.set_position([0.18, 0.120, 0.76, 0.20])

    canvas = FigureCanvas(fig)
    canvas.print_figure(relative_zp_filename.split(".txt")[0] + "_zeropoint_comparison.pdf" , dpi=300)
    close("all")