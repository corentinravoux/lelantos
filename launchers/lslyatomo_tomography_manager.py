import os
from lslyatomo import task_manager





#### DACHSHUND parameters #####
name_pixel = os.path.join(os.getcwd(),"pixel","pixel_stripe82_DR16.bin")



#### parallel DACHSHUND parameters #####

launch_file = os.path.join(os.getcwd(),"pixel","DR16_launch_data.pickle")

software = "dachshund"
machine = "bash"

copy_folder = "/local/home/cravoux/Documents/Tomography_lslyatomo/Test_lslyatomo/test/test2/tomography"

kwargs = {"n": 1}

if __name__ =="__main__":
    pwd = "./tomography"
    if(not(os.path.isdir(pwd))):os.mkdir(pwd)

    manager = task_manager.TomographyManager(pwd,software,machine,
                                             name_pixel,launch_file,**kwargs)
    manager.launch_all()
    manager.copy()
    # manager.remove_tmp()
    # manager.displace_tomography_folder(copy_folder)
