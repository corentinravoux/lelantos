import os
from lslyatomo import tomography


#### Options #####

merge_output_maps = True
rebin_map = False
create_distance_map = True


#### Input parameters #####

tomography_path = os.path.join(os.getcwd(),"tomography_data_folder")
map_name = os.path.join(os.getcwd(),"tomography","map_reconstructed.bin")
property_file = os.path.join(os.getcwd(),"pixel","property_file.pickle")


#### Merging parameters #####

launching_file_name = os.path.join(os.getcwd(),"pixel","DR16_launch_data.pickle")


#### Rebin parameters #####

map_rebin_name = os.path.join(os.getcwd(),"tomography","map_reconstructed_rebin.bin")
new_shape = (16,14,150)
operation = "mean"
rebin_prop_name = os.path.join(tomography_path,"property_file_rebin.pickle")


#### Dist map parameters #####

name_pixel_base = os.path.join(os.getcwd(),"pixel","pixel_stripe82_DR16.bin")
name_dist_map  = os.path.join(tomography_path,"dist_map_to_los.bin")
nb_process = 1
radius_local = 50



if __name__ =="__main__":

    if(merge_output_maps):tomography.create_merged_map(tomography_path,launching_file_name,map_name,property_file)
    if(rebin_map):tomography.rebin_map(map_name,property_file,new_shape,map_rebin_name,rebin_prop_name,operation=operation)
    if(create_distance_map):tomography.create_distance_map(name_dist_map,name_pixel_base,property_file,nb_process=nb_process,radius_local=radius_local)
