import os
from lslyatomo import cosmology


name_pixel = os.path.join(os.getcwd(),"pixel","pixel_stripe82_DR16.bin")
property_file = os.path.join(os.getcwd(),"pixel","property_file.pickle")
plot_pixel_path = os.path.join(os.getcwd(),"plot_pixel")

compute_histo=True
compute_mean_distance_density=True

histo_zpar = 2.5
name_histo = os.path.join(plot_pixel_path,"histogram")

coupled_plot = True
name_dperp= os.path.join(plot_pixel_path,"density_mean_distance")


if __name__ =="__main__":

    os.makedirs(plot_pixel_path, exist_ok=True)

    plot = cosmology.PixelAnalizer(pixel=name_pixel,
                                   property_file=property_file)

    plot.analyze_pixels(compute_histo,
                        compute_mean_distance_density,
                        histo_zpar=histo_zpar,
                        name_histo=name_histo,
                        name_dperp=name_dperp,
                        coupled_plot=coupled_plot)
