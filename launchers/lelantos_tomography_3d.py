import os
from lelantos import tomography,voidfinder

#### General parameters ####

convert_map_to_vtk = True
mask_map = True
return_pixel = True
return_void = True
return_qso = True

tomography_3d_path = os.path.join(os.getcwd(),"tomography_3d")


#### Map parameters ####

map_name = os.path.join(os.getcwd(),"tomography","map_reconstructed.bin")
property_file = os.path.join(os.getcwd(),"pixel","property_file.pickle")
distance_map = os.path.join(os.getcwd(),"tomography","dist_map_to_los.bin")
criteria_distance_mask = 20

vtk_name_map = os.path.join(tomography_3d_path,"map_reconstructed")
masked_name_map = os.path.join(tomography_3d_path,"masked_map_reconstructed.bin")

#### Pixel parameters ####

pixel_name = os.path.join(os.getcwd(),"pixel","pixel_stripe82_DR16.bin")
pixel_out_name = os.path.join(tomography_3d_path,"pixel_for_3D.txt")



#### Void parameters ####

void_catalog = os.path.join(os.getcwd(),"void","Catalog_Voids_SPHERICAL_0.14threshold_0.12average_4rmin_ITERATION_deletion.fits")
void_out_name = os.path.join(tomography_3d_path,"void_for_3D.txt")
move_axis_void = None


#### QSO parameters ####

qso_catalog = os.path.join(os.getcwd(),"pixel","qso_catalog.fits")
qso_out_name = os.path.join(tomography_3d_path,"qso_for_3D.txt")
move_axis_qso = None


if __name__ == '__main__' :
    os.makedirs(tomography_3d_path,exist_ok=True)
    if(convert_map_to_vtk):
        tomography.convert_to_vtk(map_name,property_file,vtk_name_map)
    if(mask_map):
        tomography.mask_map_to_3d(map_name,property_file,masked_name_map,distance_map,criteria_distance_mask)
    if(return_pixel):
        tomography.pixel_to_3d(pixel_name,pixel_out_name)
    if(return_qso):
        voidfinder.qso_to_3d(qso_catalog,qso_out_name,moveaxis=move_axis_qso)
    if(return_void):
        voidfinder.void_to_3d(void_catalog,void_out_name,moveaxis=move_axis_void)
