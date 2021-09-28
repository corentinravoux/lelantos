import os
from lelantos import tomography




#### Stack voids options ####

map_name = os.path.join(os.getcwd(),"tomography","map_reconstructed.bin")

stack_void_path = os.path.join(os.getcwd(),"stack_void")

property_file = os.path.join(os.getcwd(),"pixel","property_file.pickle")
property_file_stack = os.path.join(stack_void_path,"property_file_stack_void.pickle")

catalog_name = os.path.join(os.getcwd(),"void",
                            "Catalog_Voids_SPHERICAL_0.14threshold_0.12average_4rmin_ITERATION_deletion.fits")
type_catalog ="void"
size_stack = 20
shape_stack = None
name_stack = os.path.join(stack_void_path,"stack_voids.bin")

normalized = None
coordinate_convert = None
interpolation_method = "LINEAR"



if __name__ == '__main__' :
    os.makedirs(stack_void_path,exist_ok=True)
    stack = tomography.TomographyStack(map_name,
                                       catalog_name,
                                       type_catalog,
                                       property_file_stack,
                                       size_stack,
                                       name_stack,
                                       shape_stack=shape_stack,
                                       property_file=property_file,
                                       coordinate_convert=coordinate_convert,
                                       interpolation_method=interpolation_method,
                                       normalized=normalized)
    stack.stack()
