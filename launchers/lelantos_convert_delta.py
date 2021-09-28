import os
from lelantos import cosmology


delta_path = os.path.join(os.getcwd(),"delta")
pixel_path = os.path.join(os.getcwd(),"pixel")

Omega_m = 0.3147
z_cut_min=2.1
z_cut_max=3.2
dec_cut_min=52.4
dec_cut_max=54.2
ra_cut_min=-147.5
ra_cut_max=-144.5
sigma_min=0.1
sigma_max= None



coordinate_transform = "middle"
software = "dachshund"
mode = "parallel"
nameout = "DR16"
property_file_name = "property_file.pickle"

repeat = False
return_qso_catalog = "qso_catalog.fits"
return_dla_catalog = None
dla_catalog = None
return_sky_catalogs = True
rebin = None


plot_delta_properties = True



# DASCHUND properties

name_pixel = "pixel_stripe82_DR16.bin"
name_map = "map_stripe82_DR16.bin"
shape_map = (32, 28, 300)
sigmaf = 0.01
lperp = 13
lpar = 13
properties = {"shape" :shape_map,
              "sigma_f":sigmaf ,
              "lperp":lperp ,
              "lpar":lpar,
              "name_pixel":name_pixel,
              "name_map":name_map}


number_chunks = (2,2)
shape_sub_map = (20,20,300)
overlaping = 13




if __name__ =="__main__":
    os.makedirs(os.path.join(os.getcwd(),"pixel"),exist_ok=True)

    delta_converter = cosmology.DeltaConverter(pixel_path,
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
    delta_converter.transform_delta(mode,
                                    nameout,
                                    properties,
                                    property_file_name,
                                    rebin=rebin,
                                    sigma_min=sigma_min,
                                    sigma_max=sigma_max,
                                    z_cut_min=z_cut_min,
                                    z_cut_max=z_cut_max,
                                    dec_cut_min=dec_cut_min,
                                    dec_cut_max=dec_cut_max,
                                    ra_cut_min=ra_cut_min,
                                    ra_cut_max=ra_cut_max,
                                    number_chunks=number_chunks,
                                    overlaping=overlaping,
                                    shape_sub_map=shape_sub_map)
