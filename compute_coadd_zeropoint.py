import sys
from scipy import loadtxt, median, array, dtype, where, hypot, cos
import pickle
import socket
from os import system
from operator import itemgetter


####    Declare the query_field and query_chip.
if len(sys.argv) != 3:
    print "Improper usage. As arguments, supply the field number and chip ID."
    print "Ex: coadd_images.py 100 C13"
    sys.exit()
query_field = sys.argv[1]
query_chip = sys.argv[2]

# query_field = "100"
# query_chip = "C21"

""" Local laptop testing and development modification. """
if socket.gethostname() == "Christopher-Kleins-MacBook-Pro.local" or           \
    socket.gethostname()[:8] == "airbears":
    coadd_directory = "/Volumes/Extra_HDD/Research/DECam_Data/coadds/chip21_coadds_and_data"
    #### Define the paths to the sex and swarp bins on my laptop.
    SEXCMD = "sex"
    SWARPCMD = "swarp"
else:
    coadd_directory = "/global/scratch2/sd/cenko/DES/Jan2014/coadds/"
    #### Define the paths to the sex and swarp bins on carver.
    SEXCMD = "/project/projectdirs/dessn/local/carver/bin/sex"
    SWARPCMD = "/project/projectdirs/dessn/local/carver/bin/swarp"
""" End of local laptop testing and development modification. """


####    Create a list of the paths on NERSC to all the reduced images.
####    Most of this work was done in creating the "imagefile_dict.p" 
####    pickled dictionary file.
imagefile_dict = pickle.load( open( "imagefile_dict.p", "rb" ) )
# The field and chip will change with each run. 45 fields and 61 chips total 
# means that this "coadd_images.py" script will run 2745 times.
image_path_list = imagefile_dict[query_field][query_chip]


""" Local laptop testing and development modification. """
if socket.gethostname() == "Christopher-Kleins-MacBook-Pro.local" or           \
    socket.gethostname()[:8] == "airbears":
    for n in range(len(image_path_list)):
        image_path_list[n] = "/Volumes/Extra_HDD/Research/DECam_Data" +        \
            image_path_list[n]
""" End of local laptop testing and development modification. """



# Define median Absolute Deviation clipping for input array of numbers.
def mad_clipping_mask(input_data, sigma_clip_level):
    medval = median(input_data)
    sigma = 1.4826 * median(abs(medval - input_data))
    high_sigma_clip_limit = medval + sigma_clip_level * sigma
    low_sigma_clip_limit = medval - sigma_clip_level * sigma
    return (input_data>(low_sigma_clip_limit)) & (input_data<(high_sigma_clip_limit))

def rrl_decam_z_to_sloan_z(decam_z, decam_z_err, period, color_coefficient):
    sdss_z = decam_z + color_coefficient*(-0.1827+0.2301*period)
    # This equation is derived from earlier PLR work, from the code that fits the SEDs
    sdss_z_err = hypot(decam_z_err, color_coefficient*(0.14))
    return sdss_z, sdss_z_err


"""
Zero-point (zpt)
Color term (ciz)
Airmass term (cam)
Delta x term (cdx)
Delta y term (cdy)
RMS for fit (well, robust scatter actually)
Number of stars rejected in clipping
Date
Photometric Flat (1 = Good, 0 = bad)
"""


# Read in the absolute photometric zeropoints (provided by Brad, using Adam's
# code to solve for the zeropoints from standard field observations) 
zpt = []
ciz = []
cam = []
cdx = []
cdy = []
std_dates = []
standard_field_datafile = file("/Users/cklein/Desktop/Ongoing_Research/DECam_Science_Verification/des_pcal/des_pcal.%s.dat" % query_chip[1:], "r")
for line in standard_field_datafile:
    if line.split()[-1] == "1":
        if line.split()[7] not in ["20121109", "20121110"]:
            zpt.append(float(line.split()[0]) - 2.94022814764)
            # subtract 2.5 * log10(15.0) = 2.940 because the standard field
            # observations were 15 sec and the target observations were 1 sec.
            ciz.append(float(line.split()[1]))
            cam.append(float(line.split()[2]))
            cdx.append(float(line.split()[3]))
            cdy.append(float(line.split()[4]))
            std_dates.append(line.split()[7])
zpt = array(zpt)
ciz = array(ciz)
cam = array(cam)
cdx = array(cdx)
cdy = array(cdy)

# Read in the relative zeropoint data produced during the image coaddition 
# process (coadd_images.py)
relative_zp_filename = coadd_directory + "/"                                   \
                            + query_field + "-" + query_chip + "_relZP_list.txt"
relative_zp_type=dtype([('image_path', '|S150'),
                 ('airmass', 'float'),
                 ('date-obs', '|S26'),
                 ('psf_zp', 'float'),
                 ('psf_zp_err', 'float'),
                 ('aper_zp', 'float'),
                 ('aper_zp_err', 'float')])
relative_zp_data = loadtxt(relative_zp_filename, dtype=relative_zp_type)

# We want to compute a zp translation offset from the photometric zeropoints 
# (from the fit to the standard star fields) to the field-chip images. Some 
# field-chip images were observed more than once a night (or not at all on some
# nights) and so we need to collate the relative zp data with each photometric
# night for which we have a fitted absolute zeropoint.
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


# This list of filepaths is the list field-chip images which are considered to
# be taken during photometric nights.
filepath_strs = []
rel_compared_data = []
for n in range(len(std_dates)):
    night_str = std_dates[n]
    if night_str in rel_dates:
        rel_indices = where(rel_dates==night_str)[0]
        for m in rel_indices:
            rel_compared_data.append([zpt[n],
                                      rel_airmasses[m], 
                                      rel_psf_zps[m], rel_psf_zp_errs[m], 
                                      rel_aper_zps[m], rel_aper_zp_errs[m],
                                      ciz[n], cam[n], cdx[n], cdy[n]])
            filepath_strs.append(image_paths[m])
rel_compared_data = array(rel_compared_data)

"""
0   zpt
1   airmass
2   psf_zp
3   psf_zp_err
4   aper_zp
5   aper_zp_err
6   ciz
7   cam
8   cdx
9   cdy
"""



translation_constant = -1*median(rel_compared_data[:,2] - rel_compared_data[:,0])
translated_std_zpts = rel_compared_data[:,0] - translation_constant

error_in_median_zp = (1/((1/rel_compared_data[:,3][rel_compared_data[:,3]!=0])**2).sum())**0.5
rms_about_median_zp = (((rel_compared_data[:,2]-translated_std_zpts)**2).mean())**0.5


# Read in the coordinates of the calibrator stars
calibrator_coords = []
zp_calibrators_region_filename = coadd_directory + "/"                         \
                 + query_field + "-" + query_chip + "_candidate_calibrators.reg"
calibrator_regionfile = file(zp_calibrators_region_filename, "r")
calibrator_regionfile.readline()
calibrator_regionfile.readline()
for line in calibrator_regionfile:
    calibrator_ra = float(line.split("(")[1].split(",")[0])
    calibrator_dec = float(line.split(",")[1].split(",")[0])
    calibrator_coords.append((calibrator_ra, calibrator_dec))
calibrator_regionfile.close()
calibrator_coords = array(calibrator_coords)

absolute_photometry = {}
for cal_coord in calibrator_coords:
    absolute_photometry[tuple(cal_coord.tolist())] = []

for n in range(len(filepath_strs)):
    image_path = filepath_strs[n]
    absolute_zeropoint = translation_constant - rel_compared_data[:,2][n]
    # calibrated_photometry = instrumental_mag + absolute_zeropoint
    relative_photometry = loadtxt(image_path.replace(".fits",                  \
                                                            ".calibrators.txt"))
    for m in range(len(calibrator_coords)):
        cal_ra = calibrator_coords[m][0]
        cal_dec = calibrator_coords[m][1]

        distance_array = 3600*hypot(cos(cal_dec*0.01745329252)*(cal_ra-relative_photometry[:,0]),           \
                                               cal_dec-relative_photometry[:,1])
        # if this detected star is within 2 arcsec of a identified calibrator,
        # then write it out
        if distance_array.min() < 2:
            rel_row = relative_photometry[                                     \
                            where(distance_array == distance_array.min())[0][0]]
            absolute_photometry[tuple(calibrator_coords[m].tolist())].append(  \
                    (n, rel_row[2] + absolute_zeropoint,                       \
                    hypot(rel_row[3], rms_about_median_zp)))





absolute_calibrator_photometry = []
for cal_coord in calibrator_coords:
    light_curve = array(absolute_photometry[tuple(cal_coord.tolist())])
    # remove outliers using mad_clipping_mask
    light_curve = light_curve[mad_clipping_mask(light_curve[:,1], 3)]
    mean_mag = light_curve[:,1].mean()
    error_about_mean = (1.0/(((1.0/light_curve[:,2])**2).sum()))**0.5
    absolute_calibrator_photometry.append((cal_coord[0], cal_coord[1],         \
                                                    mean_mag, error_about_mean))
absolute_calibrator_photometry = array(absolute_calibrator_photometry)


image_path = coadd_directory + "/" + query_field + "-" + query_chip            \
                                                             + "_psf_coadd.fits"
weightmap_image = coadd_directory + "/" + query_field + "-" + query_chip       \
                                                      + "_psf_coadd.weight.fits"

system(SEXCMD + " " + image_path + " -c decam_rrl_detect.sex" +            
    " -WEIGHT_IMAGE " + weightmap_image +                                  
    " -CATALOG_NAME " + image_path.replace(".fits", ".ldac"))
system("psfex " + image_path.replace(".fits", ".ldac") +                   
    " -c decam_rrl.psfex")
system("rm " + image_path.replace(".fits", ".ldac"))
system(SEXCMD + " " + image_path + " -c decam_coadd_rrl_psf.sex" +               
    " -CATALOG_NAME " + image_path.replace(".fits", ".sexcat") +           
    " -WEIGHT_IMAGE " + weightmap_image +                                  
    " -PSF_NAME " + image_path.replace(".fits", ".psf"))
# We write out a regions file from the sexcat for each image for the 
# non-flagged, good detections.
source_extractor_data = loadtxt(image_path.replace(".fits", ".sexcat"))
# Keep only the non-flagged, good detections.
#       source_extractor_data[:,29] is the FLAGS column
#       source_extractor_data[:,5] is the FLUX_PSF column
sex_ra = source_extractor_data[:,3][(source_extractor_data[:,29]<=3) &     \
                                    (source_extractor_data[:,5]!=0)]
sex_dec = source_extractor_data[:,4][(source_extractor_data[:,29]<=3) &    \
                                     (source_extractor_data[:,5]!=0)]
sex_mag_psf = source_extractor_data[:,21][(source_extractor_data[:,29]<=3) \
                                        & (source_extractor_data[:,5]!=0)]
sex_magerr_psf = source_extractor_data[:,22][(source_extractor_data[:,29]<=3) \
                                        & (source_extractor_data[:,5]!=0)]


coadd_zp_list = []
coadd_zp_err_list = []
for calibrator in absolute_calibrator_photometry:
    cal_ra = calibrator[0]
    cal_dec = calibrator[1]
    cal_mag = calibrator[2]
    cal_magerr = calibrator[3]
    distance_array = 3600*hypot(                                               \
                cos(cal_dec*0.01745329252)*(sex_ra-cal_ra), sex_dec-cal_dec)
    # if this detected star is within 2 arcsec of a identified calibrator,
    # then write it out
    if distance_array.min() < 2:
        sex_index = where(distance_array == distance_array.min())[0][0]
        # calibrated_photometry = instrumental_mag + absolute_zeropoint
        # absolute_zeropoint = calibrated_photometry - instrumental_mag
        coadd_zp_list.append(cal_mag - sex_mag_psf[sex_index])
        coadd_zp_err_list.append(hypot(cal_magerr, sex_magerr_psf[sex_index]))
coadd_zp_array = array(coadd_zp_list)
coadd_zp_err_array = array(coadd_zp_err_list)

zp_mask = mad_clipping_mask(coadd_zp_array, 3)
coadd_zp_array = coadd_zp_array[zp_mask]
coadd_zp_err_array = coadd_zp_err_array[zp_mask]

coadd_zp = coadd_zp_array.mean()
coadd_zp_std = coadd_zp_array.std()
coadd_zp_mean_err = (1.0/(((1.0/coadd_zp_err_array)**2).sum()))**0.5

zeropoint_file = file(image_path.replace(".fits", ".zerpoint.txt"), "w")
zeropoint_file.write("%.5f\t%.5f\n" % (coadd_zp, coadd_zp_mean_err))
zeropoint_file.close()




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




rrl_data = []

mag_diffs = []
mag_sigmas = []
comp_coords = []
mag_err_ratios = []

psf_mags = []
psf_magerrs = []

cal_psf_mags = []
cal_psf_magerrs = []

output_rrl_data = []

output_region_file = file(image_path.replace(".fits", ".reg"), "w")
output_region_file.write("global color=blue dashlist=8 3 width=1 " + 
    "font='helvetica 12 bold' select=1 highlite=1 dash=0 fixed=0 edit=1 " + 
    "move=1 delete=1 include=1 source=1\nfk5\n")
for n in range(len(sex_ra)):
    star_ra = sex_ra[n]
    star_dec = sex_dec[n]
    star_psfmag_z = sex_mag_psf[n] + coadd_zp
    star_psfmagerr_z = hypot(sex_magerr_psf[n], coadd_zp_mean_err)
    # star_psfmagerr_z = hypot(sex_magerr_psf[n], coadd_zp_std)
    output_region_file.write(('''circle(%f,%f,1") # color=blue ''' + 
        "width=1 text={MAG_PSF = %.3f +/- %.3f}\n") % 
        (star_ra, star_dec, star_psfmag_z, star_psfmagerr_z))

    distance_array = 3600*hypot(cos(star_dec*0.01745329252)*                   \
                            (star_ra-absolute_calibrator_photometry[:,0]),     \
                            star_dec-absolute_calibrator_photometry[:,1])
    if distance_array.min() < 2.0:
        cal_index = where(distance_array == distance_array.min())[0][0]
        cal_ra = absolute_calibrator_photometry[cal_index][0]
        cal_dec = absolute_calibrator_photometry[cal_index][1]
        cal_mag = absolute_calibrator_photometry[cal_index][2]
        cal_magerr = absolute_calibrator_photometry[cal_index][3]
        output_region_file.write(('''circle(%f,%f,0.5") # color=green ''' + 
            "width=1\n") % (cal_ra, cal_dec))
        output_region_file.write(("#text(%f,%f) text={CAL_MAG = %.3f +/- %.3f} color=blue\n") % 
            (cal_ra, cal_dec + 2.0/3600, cal_mag, cal_magerr))
        
        mag_err_ratios.append(cal_magerr/star_psfmagerr_z)
        comp_coords.append((cal_ra, cal_dec))
        mag_diffs.append(star_psfmag_z-cal_mag)
        mag_sigmas.append((star_psfmag_z-cal_mag)/hypot(star_psfmagerr_z, cal_magerr))
        
        cal_psf_mags.append(star_psfmag_z)
        cal_psf_magerrs.append(star_psfmagerr_z)
    else:
        rrl_distance_array = 3600*hypot(cos(star_dec*0.01745329252)*           \
                            (star_ra-ogle_data["ra"]), star_dec-ogle_data["dec"])
        if rrl_distance_array.min() < 0.6:
            rrl_index = where(rrl_distance_array == rrl_distance_array.min())[0][0]
            rrl_data = ogle_data[rrl_index]
            rrl_ra = rrl_data["ra"]
            rrl_dec = rrl_data["dec"]
            ogle_id = rrl_data["ogle_id"]
            rrl_type = rrl_data["type"]
            rrl_V = rrl_data["V_mag"]
            rrl_I = rrl_data["I_mag"]
            rrl_per = rrl_data["period"]
            output_rrl_data.append((rrl_ra, rrl_dec, ogle_id, rrl_per, rrl_V, rrl_I, star_psfmag_z, star_psfmagerr_z, rrl_type, rrl_distance_array.min()))
        else:
            psf_mags.append(star_psfmag_z)
            psf_magerrs.append(star_psfmagerr_z)
output_region_file.close()

psf_mags = array(psf_mags)
psf_magerrs = array(psf_magerrs)
cal_psf_mags = array(cal_psf_mags)
cal_psf_magerrs = array(cal_psf_magerrs)

mag_diffs = array(mag_diffs)
mag_sigmas = array(mag_sigmas)
mag_err_ratios = array(mag_err_ratios)

mag_diffs = mag_diffs[mad_clipping_mask(mag_sigmas, 3)]
mag_sigmas = mag_sigmas[mad_clipping_mask(mag_sigmas, 3)]

# print "percent within 1 sigma is %.2f%%" % (100*((mag_sigmas<1) & (mag_sigmas>-1)).sum()/float(len(mag_sigmas)))


rrl_psf_mags = []
rrl_psf_magerrs = []

if len(output_rrl_data) == 0:
    print "No RRLs detected in Field %s Chip %s, exiting." % (query_field, query_chip)
    sys.exit()

set_of_ogle_ids = set(array(output_rrl_data)[:,2])


output_region_file = file(image_path.replace(".fits", ".rrl.reg"), "w")
output_region_file.write("global color=blue dashlist=8 3 width=1 " + 
    "font='helvetica 12 bold' select=1 highlite=1 dash=0 fixed=0 edit=1 " + 
    "move=1 delete=1 include=1 source=1\nfk5\n")


output_rrl_datafile = file(image_path.replace(".fits", "_rrl_data.txt"), "w")
unique_list_of_ogle_ids = list(set_of_ogle_ids)
for ogle_id in unique_list_of_ogle_ids:
    rrl_duplicates = []
    for rrl_row in output_rrl_data:
        if ogle_id == rrl_row[2]:
            rrl_duplicates.append(rrl_row)
    sorted_rrl_duplicates = sorted(rrl_duplicates, key=itemgetter(9))
    
    rrl_row = sorted_rrl_duplicates[0]
    rrl_psf_mags.append(rrl_row[6])
    rrl_psf_magerrs.append(rrl_row[7])
    output_region_file.write(('''circle(%f,%f,0.75") # color=cyan ''' + 
        "width=1\n") % (rrl_row[0], rrl_row[1]))
    output_region_file.write(("#text(%f,%f) text={RRL_P = %.3f (%s)} color=cyan\n") % 
        (rrl_row[0], rrl_row[1]+0.15/3600, rrl_row[3], rrl_row[8]))
    output_region_file.write(("#text(%f,%f) text={z = %.3f +/- %.3f} color=cyan\n") % 
        (rrl_row[0], rrl_row[1]-0.15/3600, rrl_row[6], rrl_row[7]))
    
    period = rrl_row[3]
    if rrl_row[8] == "RRc":
        fundamental_period = 10**(log10(period)+0.127)
    else:
        fundamental_period = period
    
    decam_z_mag = rrl_row[6]
    decam_z_magerr = rrl_row[7]
    
    sdss_z, sdss_z_err = rrl_decam_z_to_sloan_z(decam_z_mag, decam_z_magerr, fundamental_period, ciz.mean())
    
    if rrl_row[8] == "RRab":
        output_rrl_datafile.write("%.7f\t%.7f\t%s\t%s\t%.5f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" % (rrl_row[0], rrl_row[1], rrl_row[2], rrl_row[8], rrl_row[3], rrl_row[4], rrl_row[5], rrl_row[6], rrl_row[7], sdss_z, sdss_z_err))
    else:
        output_rrl_datafile.write("%.7f\t%.7f\t%s\t%s\t\t%.5f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" % (rrl_row[0], rrl_row[1], rrl_row[2], rrl_row[8], rrl_row[3], rrl_row[4], rrl_row[5], rrl_row[6], rrl_row[7], sdss_z, sdss_z_err))

output_rrl_datafile.close()
output_region_file.close()

rrl_psf_mags = array(rrl_psf_mags)
rrl_psf_magerrs = array(rrl_psf_magerrs)

from pylab import *
from matplotlib.pylab import *
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
from matplotlib.ticker import MaxNLocator

rcParams['axes.unicode_minus'] = False
matplotlib.rc('font', family="serif")
matplotlib.rc('font', serif="Times New Roman")
matplotlib.rc('xtick', labelsize=12)
matplotlib.rc('ytick', labelsize=12)
matplotlib.rc('axes', labelsize=14)

fig = plt.figure(figsize = (3.3, 2.5))
ax1 = subplot(111)

ax1.hlines(0, 13.8, 25, color="k", linestyle="--", linewidth=0.5, alpha=0.3)

ax1.scatter(psf_mags, psf_magerrs, marker="o", s=3, color="blue", linewidths=0, alpha=0.3, label="non-calibrators (%i)" % len(psf_mags))
ax1.scatter(cal_psf_mags, cal_psf_magerrs, marker="o", s=3, color="red", linewidths=0, alpha=0.3, label="calibrators (%i)" % len(cal_psf_mags))
ax1.scatter(rrl_psf_mags, rrl_psf_magerrs, marker="o", s=3, color="cyan", linewidths=0, alpha=0.7, label="RRL stars (%i)" % len(rrl_psf_mags))
ax1.scatter(absolute_calibrator_photometry[:,2], absolute_calibrator_photometry[:,3], marker="o", s=3, color="green", linewidths=0, alpha=0.3, label="calibrators from light curve")

ax1.legend(loc="upper left", fontsize=6, scatterpoints=1)

# This code draws major and minor tick lines. Major ticks get number labels.
majorLocator_x = MultipleLocator(1.0)
minorLocator_x = MultipleLocator(0.5)
ax1.xaxis.set_major_locator(majorLocator_x)
ax1.xaxis.set_minor_locator(minorLocator_x)

majorLocator_y1 = MultipleLocator(0.02)
minorLocator_y1 = MultipleLocator(0.01)
ax1.yaxis.set_major_locator(majorLocator_y1)
ax1.yaxis.set_minor_locator(minorLocator_y1)

ax1.set_xlabel(r"$z$ (DECam filter)")
ax1.set_ylabel(r"$\sigma_z$") 

ax1.set_ylim(-0.004, 0.11)
ax1.set_xlim(10.85, 20.3)

# pos =         [left, bottom, width, height]
ax1.set_position([0.18, 0.190, 0.79, 0.79])

canvas = FigureCanvas(fig)
canvas.print_figure(image_path.replace(".fits", "_magerr.pdf"), dpi=300)
close("all")
