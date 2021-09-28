
import os
from lelantos import cosmology
import numpy as np

plot_delta_path = os.path.join(os.getcwd(),"plot_delta")
delta_path = os.path.join(os.getcwd(),"delta")
center_ra_to_zero = True
name = "DR16"
degree = True

z_cut_min=2.1
z_cut_max=3.2
dec_cut_min=52.4
dec_cut_max=54.2
ra_cut_min=-147.5
ra_cut_max=-144.5


value_names = ["snr","delta","sigma","redshift"]
plot_histo = True
plot_mean_z_dependence = True
plot_z_dependence = False
plot_comparison = True
ra_dec_plots = True



comparison = [os.path.join(os.getcwd(),"delta_shuffle")]
comparison_legend = ["DR16","DR16 Shuffle"]
name_comparison = "DR16_shuffle_comparison"
print_stats = True

plot_args = {"sigma_bins" : 36,
             "sigma_z_bins" : 50,
             "sigma_min": 0.0,
             "sigma_max" : 2.0,
             "sigma_outlier_insensitive" : False,

             "delta_bins" : 100,
             "delta_z_bins" : 20,
             "delta_min_lim": -2.0,
             "delta_max_lim" : 2.0,
             "delta_linestyle" : "None",
             "delta_lambda_rest" : True,
             "delta_lambda_obs" : False,
             "delta_color" : 'b',

             "ra_bins" : 20,

             "snr_norm" : True,
             "snr_bins" : 100,
             "snr_z_bins": 30,
             "snr_multiplicative_coef" : 1/np.sqrt(0.8),


}

plot_args_comparison = {key: value for key, value in plot_args.items()}

plot_args_comparison.update({"delta_alpha" : 0.5,
                             "sigma_alpha" : 0.5,
                             "redshift_alpha" : 0.5,
})


if __name__ =="__main__":

    os.makedirs(plot_delta_path,exist_ok=True)

    treat = cosmology.DeltaAnalyzer(plot_delta_path,
                                    delta_path,
                                    center_ra=center_ra_to_zero,
                                    z_cut_min=z_cut_min,
                                    z_cut_max=z_cut_max,
                                    dec_cut_min=dec_cut_min,
                                    dec_cut_max=dec_cut_max,
                                    ra_cut_min=ra_cut_min,
                                    ra_cut_max=ra_cut_max,
                                    degree=degree)


    treat.plot(value_names,
               name,
               histo=plot_histo,
               mean_z_dependence=plot_mean_z_dependence,
               z_dependence=plot_z_dependence,
               ra_dec_plots=ra_dec_plots,
               **plot_args)

    if(plot_comparison):
        treat.plot(value_names,
                   name_comparison,
                   comparison=comparison,
                   comparison_legend=comparison_legend,
                   histo=plot_histo,
                   mean_z_dependence=plot_mean_z_dependence,
                   z_dependence=plot_z_dependence,
                   print_stats=print_stats,
                   **plot_args_comparison)
