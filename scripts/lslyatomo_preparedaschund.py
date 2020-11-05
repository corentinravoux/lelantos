
import os
pwd = os.getcwd()
from lslyatomo import cosmology

#### Options ####


#create_daschund_inputs = True
#parallel_launch = False
#compute_dperp = False
#plot_histo_mean_distance = False
#plot_comparaison_density_dperp = False


delta_path = "./deltas/"

Omega_m = 0.3147
z_cut_min=2.1
z_cut_max=3.2
dec_cut_min=None
dec_cut_max=None
ra_cut_min=None
ra_cut_max=None
sigma_min=0.1
sigma_max= None

coordinate_transform = "middle"
software = "dachshund"
repeat = False
plot_delta_properties = True # to debug
return_qso_catalog = "qso_catalog.fits"
return_dla_catalog = None #"dla_catalog.fits"
dla_catalog = None #"/local/home/cravoux/Documents/Tomography/DLA/DR16/DR16v2_confcut90_nhcut203.fits"
return_sky_catalogs = True
rebin = None #None
shuffle = None #None

# DASCHUND properties

name_pixel = "pixel_stripe82_DR16.bin"
name_map = "map_stripe82_DR16.bin"
shape_map = (6246,180,1914)
sigmaf = 0.01
lperp = 13
lpar = 13
properties = {"shape" :shape_map , "sigma_f":sigmaf , "lperp":lperp ,"lpar":lpar,"name_pixel":name_pixel,"name_map":name_map}
nameout = "input_test"


number_chunks = (2,2)
shape_sub_map = (20,20,20)
overlaping = 13

mode = "parallel"


delta_path_ra_dec = "/local/home/cravoux/Documents/Crosscorr_void_lya/eBOSS_DR16/deltas/deltas_eBOSS_DR16/"

center_ra_to_zero = True

plot_name = "test"


if __name__ =="__main__":

    treat = cosmology.DeltaConverter(pwd,Omega_m,delta_path,coordinate_transform,plot_delta_properties,software,return_qso_catalog=return_qso_catalog,return_dla_catalog=return_dla_catalog,dla_catalog=dla_catalog,return_sky_catalogs=return_sky_catalogs,repeat=repeat)
    out = treat.transform_delta(mode,nameout,properties,rebin=rebin,shuffle=shuffle,sigma_min=sigma_min,sigma_max=sigma_max,z_cut_min=z_cut_min,z_cut_max=z_cut_max,dec_cut_min=dec_cut_min,dec_cut_max=dec_cut_max,ra_cut_min=ra_cut_min,ra_cut_max=ra_cut_max,number_chunks=number_chunks,overlaping=overlaping,shape_sub_map=shape_sub_map)




#    treat = cosmology.DeltaAnalyzer(pwd,delta_path,center_ra=center_ra_to_zero,z_cut_min=z_cut_min,z_cut_max=z_cut_max,dec_cut_min=dec_cut_min,dec_cut_max=dec_cut_max,ra_cut_min=ra_cut_min,ra_cut_max=ra_cut_max,degree=True)
#    (ra,dec,z,zqso,ids,sigmas,deltas) = treat.get_ra_dec()
#    treat.analyze_deltas(plot_name,plot_ra_dec=True,plot_histo_ra_dec=True,plot_density_ra_dec=True,plot_histo_delta=True,plot_snr=True,plot_binned_delta=True)

    
    




