import os
pwd = os.getcwd()
from lslyatomo import tomography



#### Options ####

stack_voids = True
stack_clusters = True
stack_qso = True
stack_hsc = False


#### General options ####

shapeMap = (1464, 180, 834)
PixelName = "pixel_stripe82_DR16.bin"
MapName = "map_reconstructed.bin"
mapsize = (1587.5079899286238, 180.58889824100757, 834.0181249831885)
part_of_stripe = 1





#### Stack clusters options ####

deltamin_cluster,deltamax_cluster = -0.2,0.0
clusters = "Dictionary_Clusters_SPHERICAL_-0.14threshold_-0.12average_7rmin.pickle"
size_stack_cluster = 20
normalized_by_radius_cluster = None
name_out_cluster = "clusters_{}".format(part_of_stripe)
if normalized_by_radius_cluster is not None :
    name_out_void = name_out_void + "_renorm{}".format(normalized_by_radius_cluster)


#### Stack quasars options ####

deltamin_qso,deltamax_qso = -0.04,0.0
qso = "DataQSOposition.pickle"
size_stack_qso= 60
name_out_qso = "quasars_{}".format(part_of_stripe)




#### Stack HSC options ####

hsc_name_file = "pdr2_test_photoz.fits"
hsc_type_file = "FITS"
hsc_zmin , hsc_zmax = 2.1,3.2
hsc_confidence = 0.0
hsc_std = 5
hsc_mag_max = 23.0
hsc_omega_m = 0.3089
deltamin_hsc,deltamax_hsc = -0.05,0.05
size_stack_hsc = 20
name_out_hsc = "hsc_{}".format(part_of_stripe)


#### Stack voids options ####

name_map = "Python_Results/map_reconstructed.bin"

property_file_map = "property_file.pickle"
property_file_stack = "property_file_stack_void.pickle"
voids = "Dictionary_Voids_WATERSHED_0.11threshold_1.5dist_clusters_4rmin_CLUSTERS_deletion.fits" #"todo"
type_catalog ="void"
size_stack_void = 20
name_out_void = "voids.bin"


normalized_by_radius_void = None # To implement



if __name__ == '__main__' :
    
    tomography.stack_map(name_out_void,voids,type_catalog,name_map,property_file_map,property_file_stack,size_stack_void)


    # Treat = tomography.TreatClamato(pwd,MapName,shapeMap,PixelName)
    # if stack_voids : Treat.create_and_plot_stack_voids_or_cluster(mapsize,voids,size_stack_void,deltamin_void,deltamax_void,name_out_void,normalized=normalized_by_radius_void,number=None,coordinate_correction=coordinate_correction)
    # if stack_clusters : Treat.create_and_plot_stack_voids_or_cluster(mapsize,clusters,size_stack_cluster,deltamin_cluster,deltamax_cluster,name_out_cluster,normalized=normalized_by_radius_cluster,number=None)
    # if stack_qso : Treat.create_and_plot_stack_quasars(mapsize,qso,size_stack_qso,deltamin_qso,deltamax_qso,name_out_qso)
    # if stack_hsc : Treat.create_and_plot_stack_hsc(mapsize,hsc_name_file,hsc_type_file,hsc_zmin,hsc_zmax,hsc_confidence,hsc_std,hsc_omega_m,hsc_mag_max,size_stack_hsc,deltamin_hsc,deltamax_hsc,name_out_hsc)