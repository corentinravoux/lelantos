import os
from lslyatomo import tomography



stack_qso_path = os.path.join(os.getcwd(),"stack_qso")
plot_stack_qso_path = os.path.join(os.getcwd(),"plot_stack_qso")

property_file_stack = os.path.join(stack_qso_path,"property_file_stack_qso.pickle")
name_stack = os.path.join(stack_qso_path,"stack_qso.bin")
name_plot = os.path.join(plot_stack_qso_path,"stack_qso")
deltamin,deltamax = None,None
ellipticity = False
los_quasar =os.path.join(os.getcwd(),"pixel","pixel_stripe82_DR16.bin")
kwargs = {"levels" : 3}




if __name__ == '__main__' :
    os.makedirs(plot_stack_qso_path,exist_ok=True)
    tomography.TomographyStack.plot_stack(plot_stack_qso_path,
                                          name_stack,
                                          property_file_stack,
                                          name_plot,
                                          ellipticity=ellipticity,
                                          los_quasar=los_quasar,
                                          **kwargs)
