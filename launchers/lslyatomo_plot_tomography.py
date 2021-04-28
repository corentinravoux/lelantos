import os
from lslyatomo import tomography




#### Options ####

plot_map = True
plot_centered_maps = False
plot_delta_histogram = False
plot_delta_histogram_comparison = False
plot_integrated_map = False




#### General options ####

name_map = os.path.join(os.getcwd(),"tomography","map_reconstructed.bin")
prop = os.path.join(os.getcwd(),"pixel","property_file.pickle")
pixel_name = os.path.join(os.getcwd(),"pixel","pixel_stripe82_DR16.bin")
qso = os.path.join(os.getcwd(),"pixel","qso_catalog.fits")
void = os.path.join(os.getcwd(),"void","Catalog_Voids_SPHERICAL_0.12threshold_0.14average_4rmin_ITERATION_deletion.fits")
minimal_void_crossing = None
distance_mask = os.path.join(os.getcwd(),"tomography","dist_map_to_los.bin")
criteria_distance_mask = 20
galaxy = None

additional_args = {"map_delta_min" : -1.0,
                   "map_delta_max" : 0.5,
                   "map_color" : 'jet_r',
                   "map_xlim_min" : 0,
                   "map_interpolation" : 'bilinear',
                   "map_dpi" : 200,

                   "pixel_marker_size" : 2,
                   "pixel_marker_edge_size" : 0.1,
                   "pixel_marker":".",
                   "pixel_marker_color":"k",

                   "pixel_bis_on" :False,
                   "pixel_bis_marker_size" :2,
                   "pixel_bis_marker_edge_size" : 0.1,
                   "pixel_bis_marker":".",
                   "pixel_bis_grey" : "0.4",
                   "pixel_bis_transparency" : 0.4,

                   "qso_marker_size" : 8,
                   "qso_marker":"*",
                   "qso_marker_color":"k",
                   "qso_marker_edge_size":1,

                   "qso_bis_on" :False,
                   "qso_bis_marker_size" : 8,
                   "qso_bis_marker":"*",
                   "qso_bis_marker_color":"k",
                   "qso_bis_marker_edge_size": 0.5,

                   "void_marker_color":"r",
                   "void_bis_marker_color":"k",

                   "color_bar_fraction" : 0.02,

                   "position_redshift_axe":"other",
                   "outward_redshift_axe":50}

#### Plot map options ####

center_mpc = 10 #"all"
direction = "y"
space = 12
name_plot_map = "DR16"
redshift_axis= True
rotate = False
cut_plot= None #(1/2,1,1)



#### Plot histogram options ####

nb_bins = 50
normalization = True
criteria_distance_mask_histo = 60
name_histo = "histogram_deltas"
log_scale = False
gauss_fit = True


#### Plot centered map options ####

catalog_centered= void
name_plot_centered = "Stripe82_void_centered"
nb_plot = 2
radius_centered = 50




#### Plot integrated map options ####

zmin = 100
zmax = 200
name_out_integrated = "Stripe82"
vmin_integrate = -0.1
vmax_integrate = 0.1




if __name__ == '__main__' :
    pwd = os.path.join(os.getcwd(),"plot_tomography")
    if(not(os.path.isdir(pwd))):os.mkdir(pwd)
    Treat = tomography.TomographyPlot(pwd,map_name=name_map,pixel_name=pixel_name,
                                      property_file=prop,**additional_args)
    if(plot_map):
        Treat.plot(name_plot_map,direction,space,center_mpc,qso=qso,
                   void=void,galaxy=galaxy,distance_mask = distance_mask,
                   criteria_distance_mask = criteria_distance_mask,rotate = rotate,
                   minimal_void_crossing = minimal_void_crossing,
                   redshift_axis=redshift_axis,cut_plot=cut_plot)

    if(plot_delta_histogram):
        Treat.plot_delta_histogram(name_histo,nb_bins,
                                   gauss_fit=gauss_fit,norm=normalization,
                                   distance_mask=distance_mask,
                                   criteria_distance_mask=criteria_distance_mask_histo,
                                   log_scale=log_scale)

    if(plot_integrated_map):
        Treat.plot_integrate_image(zmin,zmax,name_out_integrated,vmin_integrate,
                                      vmax_integrate,rotate=rotate)

    if(plot_centered_maps):
        Treat.plot_catalog_centered_maps(direction,name_plot_centered,space,
                                         catalog_centered,
                                         nb_plot,radius_centered,qso=qso,
                                         rotate = rotate)
