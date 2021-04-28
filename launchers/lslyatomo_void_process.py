import os
from lslyatomo import voidfinder

compute_additional_stats = True
cut_catalog = False
merge_catalog = False
create_fake_qso_catalog = True

### General factor parameters ###

catalog_name = os.path.join(os.getcwd(),"void","Catalog_Voids_SPHERICAL_0.14threshold_0.12average_4rmin_ITERATION_deletion.fits")
catalog_name = os.path.join(os.getcwd(),"void","Catalog_Voids_WATERSHED_0.11threshold_1.5dist_clusters_4rmin_CLUSTERS_deletion.fits")

### Additional stats parameters ###

pixel_name = os.path.join(os.getcwd(),"pixel","pixel_stripe82_DR16.bin")


### Cut parameters ###

method_cut ="ALL"
cut_crossing_param = 2
cut_radius = (7,30)
distance_map_name = os.path.join(os.getcwd(),"tomography","dist_map_to_los.bin")
distance_map_prop = os.path.join(os.getcwd(),"pixel","property_file.pickle")
distance_map_param = 20
distance_map_percent = 0.95
cut_border_prop = os.path.join(os.getcwd(),"pixel","property_file.pickle")

print(os.getcwd())


### Merge parameters ###


list_catalog_name = [catalog_name,
                     os.path.join(os.getcwd(),"void","Catalog_Voids_WATERSHED_0.11threshold_1.5dist_clusters_4rmin_CLUSTERS_deletion.fits")]
merged_catalog_name = "WATERSHED_SPHERICAL_merged_catalog.fits"


### QSO like parameters ###



catalog_name_qso_like = catalog_name


if __name__ == '__main__' :
    pwd = os.path.join(os.getcwd(),"void")
    if(not(os.path.isdir(pwd))):os.mkdir(pwd)


    if(compute_additional_stats):
        voidfinder.compute_additional_stats(catalog_name,pixel_name)
    if(cut_catalog):
        voidfinder.cut_catalog(pwd,catalog_name,method_cut,
                               cut_crossing_param=cut_crossing_param,
                               cut_radius=cut_radius,
                               distance_map_name=distance_map_name,
                               distance_map_prop=distance_map_prop,
                               distance_map_param=distance_map_param,
                               distance_map_percent=distance_map_percent,
                               cut_border_prop=cut_border_prop)
    if(merge_catalog):
        voidfinder.create_merged_catalog(pwd,list_catalog_name,merged_catalog_name)

    if(create_fake_qso_catalog):
        voidfinder.create_qso_like_catalog(catalog_name_qso_like)
