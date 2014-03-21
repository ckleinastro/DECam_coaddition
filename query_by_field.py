import pickle
from os import system

imagefile_dict = pickle.load( open( "imagefile_dict.p", "rb" ) )

local_top_dir = "/Volumes/Extra_HDD/Research/DECam_Data/global/scratch2/sd/cenko/DES/Jan2014/proc/done/"

query_field = "107"

for q in range(30):
    query_field = "1" + "%02d" % (q,)

    query_chip = "C21"

    image_path_list = imagefile_dict[query_field][query_chip]
    output_file = file("cp_images_script_%s_%s.txt" % (query_field, query_chip), "w")

    output_file.write("mkdir " + local_top_dir + "20121101/" + query_chip + "\n")
    output_file.write("mkdir " + local_top_dir + "20121101/" + query_chip + "/zband" + "\n")

    output_file.write("mkdir " + local_top_dir + "20121102/" + query_chip + "\n")
    output_file.write("mkdir " + local_top_dir + "20121102/" + query_chip + "/zband" + "\n")

    output_file.write("mkdir " + local_top_dir + "20121103/" + query_chip + "\n")
    output_file.write("mkdir " + local_top_dir + "20121103/" + query_chip + "/zband" + "\n")

    output_file.write("mkdir " + local_top_dir + "20121104/" + query_chip + "\n")
    output_file.write("mkdir " + local_top_dir + "20121104/" + query_chip + "/zband" + "\n")

    output_file.write("mkdir " + local_top_dir + "20121105/" + query_chip + "\n")
    output_file.write("mkdir " + local_top_dir + "20121105/" + query_chip + "/zband" + "\n")

    output_file.write("mkdir " + local_top_dir + "20121106/" + query_chip + "\n")
    output_file.write("mkdir " + local_top_dir + "20121106/" + query_chip + "/zband" + "\n")

    output_file.write("mkdir " + local_top_dir + "20121107/" + query_chip + "\n")
    output_file.write("mkdir " + local_top_dir + "20121107/" + query_chip + "/zband" + "\n")

    output_file.write("mkdir " + local_top_dir + "20121108/" + query_chip + "\n")
    output_file.write("mkdir " + local_top_dir + "20121108/" + query_chip + "/zband" + "\n")

    output_file.write("mkdir " + local_top_dir + "20121109/" + query_chip + "\n")
    output_file.write("mkdir " + local_top_dir + "20121109/" + query_chip + "/zband" + "\n")

    output_file.write("mkdir " + local_top_dir + "20121110/" + query_chip + "\n")
    output_file.write("mkdir " + local_top_dir + "20121110/" + query_chip + "/zband" + "\n")

    output_file.write("mkdir " + local_top_dir + "20121111/" + query_chip + "\n")
    output_file.write("mkdir " + local_top_dir + "20121111/" + query_chip + "/zband" + "\n")

    output_file.write("mkdir " + local_top_dir + "20121112/" + query_chip + "\n")
    output_file.write("mkdir " + local_top_dir + "20121112/" + query_chip + "/zband" + "\n")

    output_file.write("mkdir " + local_top_dir + "20121113/" + query_chip + "\n")
    output_file.write("mkdir " + local_top_dir + "20121113/" + query_chip + "/zband" + "\n")

    output_file.write("mkdir " + local_top_dir + "20121114/" + query_chip + "\n")
    output_file.write("mkdir " + local_top_dir + "20121114/" + query_chip + "/zband" + "\n")

    output_file.write("mkdir " + local_top_dir + "20121116/" + query_chip + "\n")
    output_file.write("mkdir " + local_top_dir + "20121116/" + query_chip + "/zband" + "\n")

    output_file.write("mkdir " + local_top_dir + "20121118/" + query_chip + "\n")
    output_file.write("mkdir " + local_top_dir + "20121118/" + query_chip + "/zband" + "\n")



    for image_path in image_path_list:
        output_file.write("scp cklein@carver.nersc.gov:" + image_path + " " + 
        local_top_dir + image_path.split("done/")[1] + "\n")
    
    image_path = image_path_list[0]
    output_file.write("scp cklein@carver.nersc.gov:" +     
        image_path[:image_path.rfind("/")] + "/mask.fits " + 
        local_top_dir + "/" + image_path.split("/")[9] + "/" + 
        image_path.split("/")[10] + "/zband/mask.fits\n")

    output_file.write("scp cklein@carver.nersc.gov:" + image_path[:image_path.rfind("/")] + "/flat.fits " + 
        local_top_dir + "/" + image_path.split("/")[9] + "/" + 
        image_path.split("/")[10] + "/zband/flat.fits\n")
    
    
    output_file.close()
