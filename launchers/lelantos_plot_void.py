import os
from lelantos import voidfinder

plot_void_path = os.path.join(os.getcwd(),"plot_void")

value_names = ["radius","redshift","central_value","los_distance"]
plot_histo = True
plot_mean_z_dependence = True
plot_z_dependence = True
plot_comparison = False

threshold = 0.14
average = 0.12
cat = f"Catalog_Voids_SPHERICAL_{threshold}threshold_{average}average_4rmin_ITERATION_deletion.fits"
name = f"Test_{threshold}_{average}"

void_catalog = os.path.join(os.getcwd(),"void",cat)


# Comparison

comparison = [os.path.join(os.getcwd(),"void","Catalog_Voids_WATERSHED_0.11threshold_1.5dist_clusters_4rmin_CLUSTERS_deletion.fits")]
comparison_legend = ["SO","WS"]
name_comparison = "Test_SO_vs_WS"

# Plot args

plot_args = {"central_value_bins" : 50,
             "central_value_z_bins" : 20,
             "central_value_min" : None,
             "central_value_max" : None,

             "radius_bins" : 36,
             "radius_z_bins" : 20,
             "radius_min": 5,
             "radius_max" : 40,

             "redshift_bins" : 50,
             "redshift_min": None,
             "redshifts_max" : None,
             }

plot_args_comparison = {key: value for key, value in plot_args.items()}

plot_args_comparison.update({"central_value_alpha" : 0.5,
                             "radius_alpha" : 0.5,
                             "redshift_alpha" : 0.5,
})



if __name__ == '__main__' :
    os.makedirs(plot_void_path,exist_ok=True)
    plot = voidfinder.PlotVoid(plot_void_path,void_catalog)

    plot.plot(value_names,name,
              histo=plot_histo,
              mean_z_dependence=plot_mean_z_dependence,
              z_dependence=plot_z_dependence,
              **plot_args)

    if(plot_comparison):
        plot.plot(value_names,name_comparison,
                  comparison=comparison,
                  comparison_legend=comparison_legend,
                  histo=plot_histo,
                  mean_z_dependence=plot_mean_z_dependence,
                  z_dependence=plot_z_dependence,
                  **plot_args_comparison)
