import os
pwd = os.getcwd()
from lslyatomo import tomography


#### Options #####

merge_output_maps = False
rebin_map = True
create_distance_map = True

#### Input parameters #####

name_pixel_base = "pixel_stripe82_DR16.bin"
shape_map = None
MapName = "./Python_Results/map_reconstructed.bin"
mapsize = None #(1589.3355138492882, 180.36550604086563, 834.6641917663637)


#### Merging parameters #####

number_chunks = (6,1)
shape_sub_map = (20,20,20)
overlaping = 10.0

#### Rebin parameters #####

name_out = "./Python_Results/map_reconstructed_rebin.bin"
new_shape = (17,15,10)
operation = "mean"

#### Dist map parameters #####

name_dist_map  = "dist_map_to_los.bin"
nb_process = 1 


property_file = "property_file.pickle"



launching_file_name = "data_launch_dachshund.pickle"
overlaping = 10.0
MapName = "./Python_Results/map_reconstructed.bin"


if __name__ =="__main__":

    map_merged = tomography.create_merged_map(launching_file_name,MapName,property_file)
    tomography.rebin_map(MapName,property_file,new_shape,name_out)
    dist_map = tomography.create_distance_map(name_dist_map,name_pixel_base,property_file,nb_process=1,radius_local=50)
