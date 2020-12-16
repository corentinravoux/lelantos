import os
pwd = os.getcwd()
from lsstomo import threedview

#### General parameters ####

return_map = False
return_los = False
return_voids = True



#### Map parameters ####

convert_map_to_vtk = False
map_name = "map_reconstructed.bin"
shape_map = (1464, 180, 834)
size_map = (1587.5079899286238, 180.58889824100757, 834.0181249831885)
rebin_param = 2
map_name_out = "map_reconstructed.vtk"



#### Map parameters ####

pixel_name = "pixel_stripe82_DR16.bin"
pixels_out = "pixel_for_3D.txt"



#### Void/Proto-cluster parameters ####

void_catalog = "Dictionary_Voids_SPHERICAL_0.14threshold_0.12average_7rmin.pickle"
void_out_name = "voids_for_3D.txt"
type_file = "PICKLE"
move_axis = None


    

if __name__ == '__main__' :
    view = threedview.View_3D_Tomography(pwd)
    if(return_map):
        view.treat_map_for_3d(convert_map_to_vtk,map_name,shape_map,size_map,map_name_out,rebin_param=rebin_param)
    if(return_los):
        view.PixelTo3D(pixel_name,pixels_out)
    if(return_voids):
        view.voids_to_3D(void_catalog,void_out_name,type_file,moveaxis=move_axis)
