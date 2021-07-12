import os
from lslyatomo import tomography



stack_void_path = os.path.join(os.getcwd(),"stack_void")
plot_stack_void_path = os.path.join(os.getcwd(),"plot_stack_void")

property_file_stack = os.path.join(stack_void_path,"property_file_stack_void.pickle")
name_stack = os.path.join(stack_void_path,"stack_voids.bin")
name_plot = os.path.join(plot_stack_void_path,"stack_voids")
ellipticity = True
pixel_file_qso_distance = None
kwargs = {"levels" : 3}




if __name__ == '__main__' :
    os.makedirs(plot_stack_void_path,exist_ok=True)
    tomography.TomographyStack.plot_stack(plot_stack_void_path,
                                          name_stack,
                                          property_file_stack,
                                          name_plot,
                                          ellipticity=ellipticity,
                                          pixel_file_qso_distance=pixel_file_qso_distance,
                                          **kwargs)
