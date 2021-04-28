import os
from lslyatomo import voidfinder



launch_finder = True
find_cluster = False
restart = True
number_core = 1
map_name = os.path.join(os.getcwd(),"tomography","map_reconstructed.bin")
property_file = os.path.join(os.getcwd(),"pixel","property_file.pickle")

minimal_radius = 4
maximal_radius = 70


method_finder = "SPHERICAL" #"WATERSHED" #"SPHERICAL"


if method_finder == "SPHERICAL":
    threshold = 0.14
    average = 0.12
    radius_step = 1
    delete_option = "ITERATION"
    split_map = (2,2)
    split_overlap = 20
    params_void_finder = {"method":method_finder,"threshold":threshold,"average":average,"minimal_radius":minimal_radius,"maximal_radius":maximal_radius,"radius_step":radius_step}



if method_finder == "WATERSHED":
    threshold = 0.11
    dist_clusters = 1.5
    delete_option = "CLUSTERS"
    split_map = (2,3)
    split_overlap = 20
    params_void_finder = {"method":method_finder,"threshold":threshold,"dist_clusters":dist_clusters,"minimal_radius":minimal_radius,"maximal_radius":maximal_radius}



if __name__ == '__main__' :
    pwd = os.path.join(os.getcwd(),"void")
    if(not(os.path.isdir(pwd))):os.mkdir(pwd)
    void_finder = voidfinder.VoidFinder(pwd,map_name,params_void_finder,map_property_file=property_file,number_core=number_core,find_cluster=find_cluster,split_map=split_map,split_overlap=split_overlap,delete_option=delete_option,restart=restart)
    if launch_finder : map_chunks = void_finder.find_voids()
    void_finder.log.close()
