import os
from lslyatomo import tomography




property_file_stack = os.path.join(os.getcwd(),"stack","property_file_stack_void.pickle")
name_stack = os.path.join(os.getcwd(),"stack","stack_voids.bin")
name_plot = os.path.join(os.getcwd(),"plot_stack","stack_voids")
deltamin,deltamax = None,None
ellipticity = True
los_quasar = None
kwargs = {"levels" : 3}




if __name__ == '__main__' :
    pwd = os.path.join(os.getcwd(),"plot_stack")
    if(not(os.path.isdir(pwd))):os.mkdir(pwd) 
    tomography.TomographyStack.plot_stack(name_stack,property_file_stack,name_plot,
                                          deltamin,deltamax,ellipticity=ellipticity,
                                          los_quasar=los_quasar,**kwargs)
