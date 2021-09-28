import os
from lelantos import cosmology


delta_path = os.path.join(os.getcwd(),"delta")
delta_subsample_path = os.path.join(os.getcwd(),"delta_subsample")
name_out = "delta"


center_ra = True
z_cut_min=2.1
z_cut_max=3.2
dec_cut_min=52.4
dec_cut_max=54.2
ra_cut_min=-147.5
ra_cut_max=-144.5
degree = True



number_cut = 2
random_density_parameter = 10
number_repeat = 1
Omega_m = 0.3147
coordinate_transform = "middle"
density_names = ["density.pickle"]
separation_names = ["dperp.pickle"]
iterative_selection_parameters = None #{"density_names":density_names,
                                      # "separation_names":separation_names,
                                      # "Om":Omega_m,
                                      # "coordinate_transform":coordinate_transform}

if __name__ == '__main__' :
    os.makedirs(delta_subsample_path,exist_ok=True)

    delta_modifier = cosmology.DeltaModifier(delta_subsample_path,
                                             delta_path)

    delta_modifier.subsample_deltas(name_out,
                                    number_cut,
                                    center_ra=center_ra,
                                    z_cut_min=z_cut_min,
                                    z_cut_max=z_cut_max,
                                    dec_cut_min=dec_cut_min,
                                    dec_cut_max=dec_cut_max,
                                    ra_cut_min=ra_cut_min,
                                    ra_cut_max=ra_cut_max,
                                    random_density_parameter=random_density_parameter,
                                    number_repeat=number_repeat,
                                    iterative_selection_parameters=iterative_selection_parameters)
