import os
from lslyatomo import tomography




#### Stack voids options ####

map_name = os.path.join(os.getcwd(),"tomography","map_reconstructed.bin")

property_file = os.path.join(os.getcwd(),"pixel","property_file.pickle")
property_file_stack = os.path.join(os.getcwd(),"stack","property_file_stack_void.pickle")

catalog_name = os.path.join(os.getcwd(),"void",
                            "Catalog_Voids_SPHERICAL_0.12threshold_0.14average_4rmin_ITERATION_deletion.fits")
type_catalog ="void"
size_stack = 20
shape_stack = None
name_stack = os.path.join(os.getcwd(),"stack","stack_voids.bin")

normalized = None
coordinate_convert = "full_angle"
interpolation_method = "LINEAR"



if __name__ == '__main__' :
    pwd = os.path.join(os.getcwd(),"stack")
    if(not(os.path.isdir(pwd))):os.mkdir(pwd)
    stack = tomography.TomographyStack(pwd,map_name,catalog_name,type_catalog,
                                       property_file_stack,size_stack,name_stack,
                                       shape_stack=shape_stack,
                                       property_file=property_file,
                                       coordinate_convert=coordinate_convert,
                                       interpolation_method=interpolation_method,
                                       normalized=normalized)
    stack.stack()
