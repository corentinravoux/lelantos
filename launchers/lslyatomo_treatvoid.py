import os
pwd = os.getcwd()
from lslyatomo import tomographic_objects

void_name = "Dictionary_Voids_WATERSHED_0.11threshold_1.5dist_clusters_4rmin_CLUSTERS_deletion.fits"
void_name2 = "Dictionary_Voids_SPHERICAL_0.14threshold_0.12average_4rmin_ITERATION_deletion.fits"

pixel_name = "pixel_stripe82_DR16.bin"

qso_name = "qso_catalog.fits"


method_cut =["DIST"]
cut_crossing_param = 2
cut_radius = (7,30)
distance_map_name = "dist_map_to_los.bin"
distance_map_prop = "property_file.pickle"
distance_map_param = 20
distance_map_percent = 0.95
cut_border_prop = "property_file.pickle"


if __name__ == '__main__' :
    void = tomographic_objects.VoidCatalog.init_from_fits(void_name)

    void.create_crossing_criteria(pixel_name)
    void.compute_filling_factor(property_name=distance_map_prop)
    
    # qso = void.get_crossing_qso(qso_name)
    
    void_cut = tomographic_objects.VoidCatalog.init_from_fits(void_name) 
    void_cut.cut_catalog_void(method_cut,cut_crossing_param=cut_crossing_param,pixel_name=pixel_name,cut_radius=cut_radius,distance_map_name=distance_map_name,distance_map_prop=distance_map_prop,distance_map_param=distance_map_param,distance_map_percent=distance_map_percent,cut_border_prop=cut_border_prop)
    
    
    void_merged = tomographic_objects.VoidCatalog.init_by_merging([void_name,void_name2],name="merged")