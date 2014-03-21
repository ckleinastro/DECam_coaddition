import ephem
from os import system

query_field = "100"
query_chip = "C13"

coadd_dir = "/Users/cklein/Desktop/Ongoing_Research/DECam_Science_Verification/create_relative_coadds/"

system("xy2sky -d " + coadd_dir + query_field + "-" + query_chip + "_psf_coadd.fits 2048 1024 > field_coords.txt")
# 82.39570 -64.99591 J2000 2048.000 1024.000 
coords_file = file("field_coords.txt", "r")
line = coords_file.readline()
coords_file.close()
target_ra = float(line.split()[0])
target_dec = float(line.split()[1])


search_radius = 1000    # arcseconds

ogle_data = []
ogle_seps = []

output_region_file = file(coadd_dir + query_field + "-" + query_chip + "_ogle_rrls.reg", "w")
output_region_file.write("global color=blue dashlist=8 3 width=1 " + 
    "font='helvetica 12 bold' select=1 highlite=1 dash=0 fixed=0 edit=1 " + 
    "move=1 delete=1 include=1 source=1\nfk5\n")
ogle_rrl_cat = file("ogle_III_LMC_RRLs.txt", "r")
# ID	                Field	 StarID	RA	        Decl	    Type	I	    V	    P_1	        dP_1	    T0_1	    A_1	    R21_1	phi21_1	R31_1	phi31_1	P_2	    dP_2	T0_2	A_2	    R21_2	phi21_2	R31_2	phi31_2	ID_OGLE_II	ID_MACHO	ID_GCVS	ID_OTHER	Remarks
for line in ogle_rrl_cat:
    if line[0] != "#":
        star_ra = float(line.split()[3])*15
        star_dec = float(line.split()[4])
        target_separation_arcsec = 206264.806247*(float(ephem.separation(
        (target_ra * 0.01745329252, target_dec * 0.01745329252), 
        (star_ra * 0.01745329252, star_dec * 0.01745329252))))
        ogle_seps.append(target_separation_arcsec)
        if target_separation_arcsec < search_radius:
            star_id = line.split()[0]
            star_type = line.split()[5]
            star_I = float(line.split()[6])
            star_per = float(line.split()[8])
            ogle_data.append([star_id, star_type, star_per, star_I, star_ra, star_dec])

            output_region_file.write(('''circle(%f,%f,1") # color=green ''' + 
                "width=1 text={%s, P=%.3f, I=%.3f}\n") % (star_ra, star_dec, star_type, 
                star_per, star_I))
ogle_rrl_cat.close()
output_region_file.close()


