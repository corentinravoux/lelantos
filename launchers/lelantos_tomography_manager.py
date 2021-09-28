import os
from lelantos import task_manager

tomography_path = os.path.join(os.getcwd(),"tomography_data_folder")

#### DACHSHUND parameters #####
name_pixel = os.path.join(os.getcwd(),"pixel","pixel_stripe82_DR16.bin")



#### parallel DACHSHUND parameters #####

launch_file = os.path.join(os.getcwd(),"pixel","DR16_launch_data.pickle")
symlink_folder = "./tomography"

software = "dachshund"
machine = "bash"


kwargs = {"n": 1}

if __name__ =="__main__":
    os.makedirs(tomography_path,exist_ok=True)

    manager = task_manager.TomographyManager(tomography_path,
                                             software,
                                             machine,
                                             name_pixel,
                                             launch_file,
                                             symlink_folder = symlink_folder,
                                             **kwargs)
    manager.launch_all()
    manager.copy()
    manager.remove_tmp()
