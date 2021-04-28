import os
from lslyatomo import cosmology


delta_path = os.path.join(os.getcwd(),"delta")
other_delta_path = os.path.join(os.getcwd(),"delta_xcorr")

center_ra = True
z_cut_min=2.1
z_cut_max=3.2
dec_cut_min=0.4
dec_cut_max=1.2
ra_cut_min=-1.0
ra_cut_max=0.75
degree = True

name_out = "new_deltas"
other_name_out = "new_deltas"
n_cut = 11


number_cut = 2

random_density_parameter = 10
number_repeat = 1


Omega_m = 0.3147
coordinate_transform = "middle"
density_names = ["density.pickle"]
separation_names = ["dperp.pickle"]
iterative_selection_parameters = None #{"density_names":density_names,"separation_names":separation_names,"Om":Omega_m,"coordinate_transform":coordinate_transform}

if __name__ == '__main__' :
    pwd = os.path.join(os.getcwd(),"delta_cut")
    if(not(os.path.isdir(pwd))):os.mkdir(pwd)

    delta_modifier = cosmology.DeltaModifier(pwd,delta_path,center_ra=center_ra,
                                             z_cut_min=z_cut_min,z_cut_max=z_cut_max,
                                             dec_cut_min=dec_cut_min,dec_cut_max=dec_cut_max,
                                             ra_cut_min=ra_cut_min,ra_cut_max=ra_cut_max)
    #
    # delta_modifier.subsample_deltas(name_out,number_cut,
    #                               random_density_parameter=random_density_parameter,
    #                               number_repeat=number_repeat,
    #                               iterative_selection_parameters=iterative_selection_parameters)


    # delta_modifier.shuffle_deltas_cut_z(name_out,n_cut,z_cut_min,z_cut_max,other_delta_path=None,other_name_out=None)
    delta_modifier.shuffle_deltas(name_out,other_delta_path=other_delta_path,other_name_out=other_name_out)
