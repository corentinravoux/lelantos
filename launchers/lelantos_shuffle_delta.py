import os
from lelantos import cosmology


delta_path = os.path.join(os.getcwd(),"delta")
other_delta_path = os.path.join(os.getcwd(),"delta_xcorr")

delta_path_out = os.path.join(os.getcwd(),"delta_shuffle")
other_delta_path_out = os.path.join(os.getcwd(),"delta_xcorr_shuffle")

zmin=2.1
zmax=3.2

n_cut = 11



if __name__ == '__main__' :
    os.makedirs(delta_path_out,exist_ok=True)
    os.makedirs(other_delta_path_out,exist_ok=True)

    delta_modifier = cosmology.DeltaModifier(delta_path_out,
                                             delta_path)

    delta_modifier.shuffle_deltas_cut_z(n_cut,
                                        zmin,
                                        zmax,
                                        other_delta_path=other_delta_path,
                                        other_path_out=other_delta_path_out)
