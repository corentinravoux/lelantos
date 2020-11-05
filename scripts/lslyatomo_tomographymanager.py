import os
pwd = "./"
from lslyatomo import task_manager





#### DACHSHUND parameters #####

name_pixel = "pixel_stripe82_DR16.bin"


#### parallel DACHSHUND parameters #####

launch_file = "data_launch_dachshund.pickle"

software = "dachshund"
machine = "pc"




if __name__ =="__main__":
    manager = task_manager.TomographyManager(pwd,software,machine,name_pixel,launch_file)
    manager.launch_all()
    manager.copy()

