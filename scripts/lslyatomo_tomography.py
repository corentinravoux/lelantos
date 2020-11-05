import os
pwd = os.getcwd()
from lslyatomo import tomography




#### Options ####

plot_map = True
plot_centered_maps = False
plot_signal_to_noise = False
plot_delta_histogram = False
plot_delta_histogram_comparison = False
plot_integrated_map = False


#### General options ####

pixel_name = None #"test_pixels.bin"
map_name = "map_reconstructed.bin"
shape_map =  (6264, 180, 1914)
map_size = (6244.895769918219, 195.77463172987967, 1911.1878071122878)

#### Plot map options ####

plot_ra_dec_map = False
space = 10
deltamin = -1.0
name_plot_map = "map_full_angle"
deltamax = 0.5
direction_plot = "x"
print_voids=None #"Dictionary_Voids_SPHERICAL_0.14threshold_0.12average_7rmin_cutcrossing1_cutdistance_20Mpc_0.95percent.pickle" # None to not plot them
print_hsc=None #"Dict_pdr2_test_photoz.fits.pickle" # None to not plot them
no_LOS = True
middle= 85  # None to print all slices
minimal_void_crossing = None #1 # None to not cut
Om = 0.3147 
max_list_name = None #"list_of_maximums_of_data_cube.pickle"
redshift_axis= None #(Om,max_list_name) # None to not plot it
name_dist_los = None #"dist_map_to_los.bin" # None to not cut
criteria_outside_LOS_range = None #20
compute_redshift_out_LOS = None # 0.25 # None to not compute redshift
rotate = True
hd = True
sub_plotting = None #(1/4,1,1)


#### Plot centered map options ####

catalog_centered="Dictionary_Voids_SPHERICAL_0.14threshold_0.12average_7rmin.pickle" # None to not plot them
catalog_centered="Dictionary_Clusters_WATERSHED_threshold_-0.272814605average_7rmin_2_cutcrossing6.pickle" # None to not plot them
name_plot_centered = "Stripe82_void_centered"
nb_plot = 5
radius_centered = 30


#### Plot histogram options ####

nb_bins = 50
normalization = True
name_dist_los_histo = "dist_map_to_los.bin"
criteria_outside_LOS_range_histo = 20
name_out_histo = "histogram_deltas"
log_scale = False


#### Plot histogram comparison options ####

legend_comparison = ["Stripe 82 DR16","Saclay mocks"]
map_name_to_compare = "/local/home/cravoux/Documents/Tomography/Saclay_mocks/deltas_merged/eBOSS/noise/45_new_cosmo/1A/1450_exptime/v4.7/RSD/delta1/map_reconstructed.bin"
name_dist_los_comparison = "/local/home/cravoux/Documents/Tomography/Saclay_mocks/deltas_merged/eBOSS/noise/45_new_cosmo/1A/1450_exptime/v4.7/RSD/delta1/dist_map_to_los.bin"
name_out_histo_comparison = "histogram_deltas_comparison"
log_scale_comparison = False


#### Plot integrated map options ####

zmax = 2.5
name_out_integrated = "Stripe82"
vmin_integrate = -0.1
vmax_integrate = 0.1


#### Plot SNR map ####

name_out_snr = "signal_to_noise_ratio_Stripe82"







### new

name_map = "./Python_Results/map_reconstructed.bin"
prop = "property_file.pickle"
pixel_name = "pixel_stripe82_DR16.bin"
qso = "qso_catalog.fits"
void = "Dictionary_Voids_WATERSHED_0.11threshold_1.5dist_clusters_4rmin_CLUSTERS_deletion.fits" #"todo"
distance_mask = "dist_map_to_los.bin"
criteria_distance_mask = 20
hd=True
direction = "y"
space = 10
deltamin = -1.0
deltamax = 0.5
center_mpc = 50
name_plot_map = "test"


galaxy = None # to test
minimal_void_crossing = None # to test

rotate = False # To debug
redshift_axis=False # To debug


additional_args = {"linewidth_pixel" : 0.2,"marker_pixel":"_"}


#### Plot integrated map options ####

zmax = 100
name_out_integrated = "Stripe82"
vmin_integrate = -0.1
vmax_integrate = 0.1



#### Plot centered map options ####

catalog_centered= void
name_plot_centered = "Stripe82_void_centered"
nb_plot = 5
radius_centered = 50


#### Plot histogram options ####

nb_bins = 50
normalization = True
criteria_distance_mask_histo = 60
name_out_histo = "histogram_deltas"
log_scale = False
gauss_fit = True


if __name__ == '__main__' :
    Treat = tomography.TomographyPlot(pwd,map_name=name_map,map_shape=None,pixel_name=pixel_name,property_file=prop,**additional_args)

    Treat.plot(name_plot_map,direction,space,center_mpc,deltamin,deltamax,qso=qso,void=void,galaxy=galaxy,distance_mask = distance_mask,criteria_distance_mask = criteria_distance_mask,rotate = rotate,hd=hd,minimal_void_crossing = minimal_void_crossing,redshift_axis=redshift_axis)

    # Treat.plot_all_slice(name_plot_map,direction,space,deltamin,deltamax,qso=qso,void=void,galaxy=galaxy,distance_mask = distance_mask,criteria_distance_mask = criteria_distance_mask,rotate = rotate,hd=hd,minimal_void_crossing = minimal_void_crossing,redshift_axis=redshift_axis)

    # from lslyatomo import tomographic_objects
    # map_tomo = tomographic_objects.TomographicMap.init_classic(name=map_name,property_file=prop)

    # Treat.compute_integrate_image(zmax,name_out_integrated,vmin_integrate,vmax_integrate,rotate=rotate,hd=hd)


    # if(plot_centered_maps):
    Treat.plot_catalog_centered_maps(direction,name_plot_centered,space,deltamin,deltamax,catalog_centered,nb_plot,radius_centered,qso=qso,rotate = rotate, hd=hd)


    # if(plot_delta_histogram):
    Treat.plot_delta_histogram(name_dist_los_histo,nb_bins,gauss_fit=gauss_fit,norm=normalization,distance_mask=distance_mask,criteria_distance_mask=criteria_distance_mask_histo,log_scale=log_scale)


    # if(plot_delta_histogram):
    #     Treat.plotHistogramDeltas(nb_bins,name_out_histo,norm = normalization,mask_distance_name=name_dist_los_histo,criteria_outside_LOS_range=criteria_outside_LOS_range_histo,log_scale=log_scale)

    # if(plot_delta_histogram_comparison):
    #     Treat.plot_histogram_delta_comparison(map_name,map_name_to_compare,nb_bins,legend_comparison,name_out_histo_comparison,log_scale=log_scale_comparison,mask1=name_dist_los,mask2=name_dist_los_comparison,param_cut=criteria_outside_LOS_range_histo)

    # if(plot_signal_to_noise):
    #     Treat.plotSNR(pixel_name,name_out_snr)
    
    
    
    
    
    
    
    
    
    
    
    
    
