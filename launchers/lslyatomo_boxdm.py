import os
pwd = os.getcwd()
from lslyatomo import boxdm


# Input DM map properties

box_dir = ["/global/cscratch1/sd/ravouxco/mocks_cross_corr_void_lya/mock_0/chunk_1/boxes","/global/cscratch1/sd/ravouxco/mocks_cross_corr_void_lya/mock_0/chunk_2/boxes"]
box_shape =[(2560,2560,1536),(2560,2560,1536)] # along X (RA), Y (DEC) and Z (z) directions
size_cell = 2.19  # Mpc.h-1
box_bound = [(-50,0,-2,2),(0,50,-2,2)] # RA min, RA max, DEC min, DEC max where the box is used
master_file = "/global/cscratch1/sd/ravouxco/mocks_cross_corr_void_lya/mock_0/output/master.fits"



# Output map properties
interpolation_method="LINEAR"

map_property_file="/global/homes/r/ravouxco/1_Documents/Crosscorr/Saclay_mocks/tomo/flux/FLUX_RSD_DLA_DLAm_Noise_C_CF_Met/pixel/property_file.pickle"
name="/global/cscratch1/sd/ravouxco/mocks_cross_corr_void_lya_output/box_dm/delta_m_smooth5.bin"
# name="/global/cscratch1/sd/ravouxco/mocks_cross_corr_void_lya_output/box_dm/delta_m_growth_smooth5.bin"
rsd_box=False
growth_multiplication=False
matter_field=True
gaussian_smoothing=5


if __name__ == "__main__":
    box = boxdm.BoxExtractor(pwd,box_dir,box_shape,size_cell,box_bound,master_file,interpolation_method=interpolation_method)
    box.create_box(map_property_file,name,rsd_box=rsd_box,growth_multiplication=growth_multiplication,matter_field=matter_field,gaussian_smoothing=gaussian_smoothing)
