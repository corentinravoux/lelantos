#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Date:

Author: Corentin Ravoux

Description:
"""

########################## MODULE IMPORTATION ###############################


import numpy as np
import os
pwd = os.getcwd()
from lslyatomo import cosmology
from lslyatomo import tomographic_objects
import pickle

################################# INPUT ######################################

delta_path = os.path.join(pwd,"data","delta")
Omega_m = 0.3147
z_cut_min=2.1
z_cut_max=3.2
dec_cut_min=0.4
dec_cut_max=1.2
ra_cut_min=-1.0
ra_cut_max=0.75
sigma_min=0.1
sigma_max= None
coordinate_transform = "middle"
software = "dachshund"
mode = "parallel"
repeat = False
return_qso_catalog = "qso_catalog.fits"
return_dla_catalog = None
dla_catalog = None
return_sky_catalogs = True
rebin = None
plot_delta_properties = False
number_chunks = (2,2)
shape_sub_map = (20,20,300)
overlaping = 13
sigmaf = 0.01
shape_map = None
lperp = 13
lpar = 13
name_pixel_properties = "pixel_stripe82_DR16.bin"
name_map_properties = "pixel_stripe82_DR16.bin"
properties = {"sigma_f":sigmaf , "lperp":lperp ,"lpar":lpar, "shape" : shape_map,
              "name_pixel":name_pixel_properties,"name_map":name_map_properties}





################################# OUTPUT #####################################


launcher_name = os.path.join(pwd,"data","pixel","DR16_launch_data")
property_file_name = os.path.join(pwd,"data","pixel","property_file.pickle")
qso_catalog_cartesian_name = os.path.join(pwd,"data","pixel","qso_catalog.fits")
qso_catalog_sky_name = os.path.join(pwd,"data","pixel","qso_catalog.fits_sky_coordinates")
name_pixel = os.path.join(pwd,"data","pixel","pixel_stripe82_DR16.bin")
name_map = os.path.join(pwd,"data","pixel","map_stripe82_DR16.bin")









def test_transformation():
    delta_converter = cosmology.DeltaConverter(pwd,
                                               Omega_m,
                                               delta_path,
                                               coordinate_transform,
                                               plot_delta_properties,
                                               software,
                                               return_qso_catalog=return_qso_catalog,
                                               return_dla_catalog=return_dla_catalog,
                                               dla_catalog=dla_catalog,
                                               return_sky_catalogs=return_sky_catalogs,
                                               repeat=repeat)
    (cartesian_deltas,
    cartesian_qso_catalog,
    cartesian_dla_catalog,
    sky_deltas,
    sky_qso_catalog,
    sky_dla_catalog,
    properties_map_pixels) = delta_converter.transform_delta_to_pixel_file(rebin=rebin,
                                                                           sigma_min=sigma_min,
                                                                           sigma_max=sigma_max,
                                                                           z_cut_min=z_cut_min,
                                                                           z_cut_max=z_cut_max,
                                                                           dec_cut_min=dec_cut_min,
                                                                           dec_cut_max=dec_cut_max,
                                                                           ra_cut_min=ra_cut_min,
                                                                           ra_cut_max=ra_cut_max)

    pixel = tomographic_objects.Pixel.init_from_property_files(property_file_name,name=name_pixel)
    pixel.read()
    np.testing.assert_equal(cartesian_deltas,pixel.pixel_array)


    qso_catalog_cartesian = tomographic_objects.QSOCatalog.init_from_fits(qso_catalog_cartesian_name)
    np.testing.assert_equal(cartesian_qso_catalog[:,0:3],qso_catalog_cartesian.coord)
    np.testing.assert_equal(cartesian_qso_catalog[:,3],qso_catalog_cartesian.primary_key)


    qso_catalog_sky = tomographic_objects.QSOCatalog.init_from_fits(qso_catalog_sky_name)
    np.testing.assert_equal(sky_qso_catalog[:,0:3],qso_catalog_sky.coord)
    np.testing.assert_equal(sky_qso_catalog[:,3],qso_catalog_sky.primary_key)






def test_parallel_launch():
    delta_converter = cosmology.DeltaConverter(pwd,Omega_m,delta_path,
                                               coordinate_transform,
                                               plot_delta_properties,
                                               software,
                                               return_qso_catalog=return_qso_catalog,
                                               return_dla_catalog=return_dla_catalog,
                                               dla_catalog=dla_catalog,
                                               return_sky_catalogs=return_sky_catalogs,
                                               repeat=repeat)
    (cartesian_deltas,
    cartesian_qso_catalog,
    cartesian_dla_catalog,
    sky_deltas,sky_qso_catalog,
    sky_dla_catalog,
    properties_map_pixels) = delta_converter.transform_delta_to_pixel_file(rebin=rebin,
                                                                           sigma_min=sigma_min,
                                                                           sigma_max=sigma_max,
                                                                           z_cut_min=z_cut_min,
                                                                           z_cut_max=z_cut_max,
                                                                           dec_cut_min=dec_cut_min,
                                                                           dec_cut_max=dec_cut_max,
                                                                           ra_cut_min=ra_cut_min,
                                                                           ra_cut_max=ra_cut_max)


    (parallel_launcher_params,
     filename,
     chunks,
     shape) = delta_converter.create_parallel_input(properties,
                                                    cartesian_deltas,
                                                    number_chunks,
                                                    overlaping,
                                                    shape_sub_map)

    for key in list(chunks.keys()):
        if key != 'overlaping' :
            pixel = tomographic_objects.Pixel(name="{}_{}".format(name_pixel,key))
            pixel.read()
            np.testing.assert_equal(chunks[key]["coord"],pixel.pixel_array)

    data_launch = pickle.load(open(f"{launcher_name}.pickle","rb"))
    np.testing.assert_equal(data_launch[0],filename)
    for i in range(len(parallel_launcher_params)):
        np.testing.assert_equal(data_launch[1][i],parallel_launcher_params[i])
    np.testing.assert_equal(data_launch[2],number_chunks)
    np.testing.assert_equal(data_launch[3],overlaping)


    property_file_to_test = delta_converter.create_dachshund_map_pixel_property_file(property_file_name,
                                                                                     cartesian_deltas,
                                                                                     sky_deltas,
                                                                                     shape,
                                                                                     properties_map_pixels)
    property_file = tomographic_objects.MapPixelProperty(name=property_file_name)
    property_file.read()
    np.testing.assert_equal(property_file_to_test.name,property_file.name)
    np.testing.assert_equal(property_file_to_test.size,property_file.size)
    np.testing.assert_equal(property_file_to_test.shape,property_file.shape)
    np.testing.assert_equal(property_file_to_test.boundary_cartesian_coord,property_file.boundary_cartesian_coord)
    np.testing.assert_equal(property_file_to_test.boundary_sky_coord,property_file.boundary_sky_coord)
    np.testing.assert_equal(property_file_to_test.coordinate_transform,property_file.coordinate_transform)
    np.testing.assert_equal(property_file_to_test.Omega_m,property_file.Omega_m)



if __name__ == "__main__":
    test_transformation()
    test_parallel_launch()
