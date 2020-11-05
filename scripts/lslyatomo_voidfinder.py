import os
pwd = os.getcwd()
from lslyatomo import voidfinder

#### Options ####




#### Map parameters ####



#### Finder parameters ####




method_finder = "SPHERICAL"
threshold = 0.14
average = 0.12
radius_step = 1
minimal_radius = 4
maximal_radius = 70
params_void_finder = {"method":method_finder,"threshold":threshold,"average":average,"minimal_radius":minimal_radius,"maximal_radius":maximal_radius,"radius_step":radius_step}
split_map = None #(2,2)
delete_option = "ITERATION"# or "NONE"

# method_finder = "WATERSHED"
# threshold = 0.11
# dist_clusters = 1.5
# minimal_radius = 4
# maximal_radius = 70
# params_void_finder = {"method":method_finder,"threshold":threshold,"dist_clusters":dist_clusters,"minimal_radius":minimal_radius,"maximal_radius":maximal_radius}
# split_map = None #(2,2)
# delete_option = "CLUSTERS"   #"CLUSTERS" "ITERATION" or "NONE"




launch_finder = True
find_cluster = False
number_core = 2
restart = False
map_name = "./Python_Results/map_reconstructed.bin"
property_file = "property_file.pickle"
split_overlap = 10

if __name__ == '__main__' :
    void_finder = voidfinder.VoidFinder(pwd,map_name,params_void_finder,property_file=property_file,number_core=number_core,restart=restart,find_cluster=find_cluster,split_map=split_map,split_overlap=split_overlap,delete_option=delete_option)
    if launch_finder : void_finder.find_voids()

