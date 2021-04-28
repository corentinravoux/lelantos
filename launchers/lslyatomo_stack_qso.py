import os
from lslyatomo import tomography




#### Stack voids options ####

map_name = os.path.join(os.getcwd(),"tomography","map_reconstructed.bin")

property_file = os.path.join(os.getcwd(),"pixel","property_file.pickle")
property_file_stack = os.path.join(os.getcwd(),"stack","property_file_stack_qso.pickle")

catalog_name = os.path.join(os.getcwd(),"pixel","qso_catalog.fits")
type_catalog ="qso"
size_stack = 20
shape_stack = None
name_stack = os.path.join(os.getcwd(),"stack","stack_qso.bin")

normalized = None
coordinate_convert = None 
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