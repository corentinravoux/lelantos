#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Corentin Ravoux

Description : Classes to return a Dachshund input based on data that can be
treated via picca software.
"""

#############################################################################
#############################################################################
########################## MODULE IMPORTATION ###############################
#############################################################################
#############################################################################


import os
import configparser
import ast
from lelantos import cosmology, task_manager, tomography, voidfinder


#############################################################################
#############################################################################
############################### FUNCTIONS ###################################
#############################################################################
#############################################################################


def parse_int_tuple(input):
    if input == "None":
        return None
    else:
        return tuple(int(k.strip()) for k in input.strip().split(","))


def parse_float_tuple(input):
    if input == "None":
        return None
    else:
        return tuple(float(k.strip()) for k in input.strip().split(","))


def parse_str_tuple(input):
    if input == "None":
        return None
    else:
        return tuple(str(k.strip()) for k in input.strip().split(","))


def parse_dict(input):
    if input == "None":
        return None
    else:
        acceptable_string = input.replace("'", '"')
        return ast.literal_eval(acceptable_string)


def parse_float(input):
    if input == "None":
        return None
    else:
        return float(input)


def parse_int(input):
    if input == "None":
        return None
    else:
        return int(input)


def parse_string(input):
    if input == "None":
        return None
    else:
        return str(input)


def main(input_file):
    config = configparser.ConfigParser(
        allow_no_value=True,
        converters={
            "str": parse_string,
            "int": parse_int,
            "float": parse_float,
            "tupleint": parse_int_tuple,
            "tuplefloat": parse_float_tuple,
            "tuplestr": parse_str_tuple,
            "dict": parse_dict,
        },
    )
    config.optionxform = lambda option: option
    config.read(input_file)

    main_config = config["main"]
    main_path = os.path.abspath(main_config["path"])
    os.makedirs(main_path, exist_ok=True)

    delta_transform_config = config["delta transform"]
    delta_config = config["delta convert"]
    software_config = config["tomography software"]
    tomography_config = config["tomography launching"]
    tomography_process_config = config["tomography process"]
    void_finder_config = config["void finder"]
    void_process_config = config["void process"]
    stack_void_config = config["void stack"]
    delta_plot_config = config["delta plot"]
    void_plot_config = config["void plot"]
    tomography_plot_config = config["tomography plot"]
    stack_void_plot_config = config["stack void plot"]

    delta_path = os.path.abspath(main_config.getstr("delta_path"))

    if main_config.getboolean("transform_delta"):
        delta_transform_path = delta_transform_config.getstr("delta_path_out")
        transform_delta(
            delta_path,
            delta_transform_path,
            delta_transform_config,
            delta_config,
            main_config,
        )
        delta_path = os.path.abspath(delta_transform_path)

    plot_delta_path = os.path.join(
        main_path, delta_plot_config.getstr("plot_delta_path")
    )
    if main_config.getboolean("plot_delta"):
        plot_delta(
            plot_delta_path, delta_path, delta_plot_config, delta_config, main_config
        )

    pixel_path = os.path.join(main_path, delta_config.getstr("pixel_path"))
    if main_config.getboolean("convert_delta"):
        convert_delta(
            pixel_path, delta_path, software_config, delta_config, main_config
        )

    tomo_abs_path = os.path.abspath(
        os.path.join(main_path, tomography_config.getstr("tomo_path"))
    )
    tomo_path = os.path.relpath(tomo_abs_path, os.getcwd())
    if main_config.getboolean("launch_tomography"):
        launch_tomography(
            tomo_path, pixel_path, software_config, tomography_config, main_config
        )

    if main_config.getboolean("process_tomography"):
        process_tomography(
            tomo_abs_path,
            pixel_path,
            tomography_process_config,
            delta_config,
            software_config,
            main_config,
        )

    void_path = os.path.abspath(
        os.path.join(main_path, void_finder_config.getstr("void_path"))
    )
    void_catalog_name_default = None
    if main_config.getboolean("find_void"):
        void_catalog_name_default = find_void(
            void_path,
            tomo_abs_path,
            pixel_path,
            void_finder_config,
            tomography_process_config,
            delta_config,
        )

    if main_config.getboolean("process_void"):
        void_catalog_name_default = process_void(
            void_path,
            tomo_abs_path,
            pixel_path,
            void_process_config,
            delta_config,
            tomography_process_config,
            software_config,
            void_catalog_name_default,
        )

    stack_void_path = os.path.abspath(
        os.path.join(main_path, stack_void_config.getstr("stack_void_path"))
    )
    if main_config.getboolean("stack_void"):
        (stack_void_name_default, property_stack_void_name_default) = stack_void(
            stack_void_path,
            tomo_abs_path,
            void_path,
            pixel_path,
            stack_void_config,
            tomography_process_config,
            delta_config,
            void_catalog_name_default,
        )

    plot_void_path = os.path.join(main_path, void_plot_config.getstr("plot_void_path"))
    if main_config.getboolean("plot_void"):
        plot_void(
            plot_void_path,
            void_path,
            void_plot_config,
            main_config,
            void_catalog_name_default,
        )

    plot_tomography_path = os.path.join(
        main_path, tomography_plot_config.getstr("plot_tomography_path")
    )
    if main_config.getboolean("plot_tomography"):
        plot_tomography(
            plot_tomography_path,
            tomo_abs_path,
            pixel_path,
            void_path,
            tomography_plot_config,
            tomography_process_config,
            software_config,
            delta_config,
            main_config,
            void_catalog_name_default,
        )

    plot_stack_void_path = os.path.join(
        main_path, stack_void_plot_config.getstr("plot_stack_void_path")
    )
    if main_config.getboolean("plot_stack_void"):
        plot_stack_void(
            plot_stack_void_path,
            stack_void_path,
            stack_void_plot_config,
            stack_void_name_default,
            property_stack_void_name_default,
        )


def transform_delta(
    delta_path, delta_transform_path, delta_transform_config, delta_config, main_config
):
    os.makedirs(delta_transform_path, exist_ok=True)
    if delta_transform_config.getstr("other_delta_path_out") is not None:
        os.makedirs(
            delta_transform_config.getstr("other_delta_path_out"), exist_ok=True
        )
    delta_modifier = cosmology.DeltaModifier(
        delta_transform_path, main_config.getstr("delta_path")
    )

    delta_modifier.shuffle_deltas_cut_z(
        delta_transform_config.getint("n_cut"),
        delta_config.getfloat("z_cut_min"),
        delta_config.getfloat("z_cut_max"),
        other_delta_path=delta_transform_config.getstr("other_delta_path"),
        other_path_out=delta_transform_config.getstr("other_delta_path_out"),
        seed=delta_transform_config.getint("seed"),
    )


def convert_delta(pixel_path, delta_path, software_config, delta_config, main_config):
    os.makedirs(pixel_path, exist_ok=True)
    properties = {
        "sigma_f": software_config.getfloat("sigmaf"),
        "lperp": software_config.getfloat("lperp"),
        "lpar": software_config.getfloat("lpar"),
        "name_pixel": software_config.getstr("name_pixel"),
        "name_map": software_config.getstr("name_map"),
    }

    delta_converter = cosmology.DeltaConverter(
        pixel_path,
        delta_config.getfloat("Omega_m"),
        delta_path,
        delta_config.getstr("coordinate_transform"),
        delta_config.getboolean("plot_delta_properties"),
        software_config.getstr("software"),
        return_qso_catalog=delta_config.getstr("return_qso_catalog"),
        return_dla_catalog=delta_config.getstr("return_dla_catalog"),
        dla_catalog=delta_config.getstr("dla_catalog"),
        return_sky_catalogs=delta_config.getboolean("return_sky_catalogs"),
        repeat=delta_config.getstr("repeat"),
        center_ra=delta_config.getboolean("center_ra"),
    )

    delta_converter.transform_delta(
        delta_config.getstr("mode"),
        f"{main_config.getstr('name')}",
        properties,
        delta_config.getstr("property_file_name"),
        rebin=delta_config.getstr("rebin"),
        sigma_min=delta_config.getfloat("sigma_min"),
        sigma_max=delta_config.getfloat("sigma_max"),
        z_cut_min=delta_config.getfloat("z_cut_min"),
        z_cut_max=delta_config.getfloat("z_cut_max"),
        dec_cut_min=delta_config.getfloat("dec_cut_min"),
        dec_cut_max=delta_config.getfloat("dec_cut_max"),
        ra_cut_min=delta_config.getfloat("ra_cut_min"),
        ra_cut_max=delta_config.getfloat("ra_cut_max"),
        number_chunks=delta_config.gettupleint("number_chunks"),
        overlaping=delta_config.getfloat("overlaping"),
        shape_sub_map=delta_config.gettupleint("shape_sub_map"),
    )


def plot_delta(
    plot_delta_path, delta_path, delta_plot_config, delta_config, main_config
):
    os.makedirs(plot_delta_path, exist_ok=True)
    delta_plotter = cosmology.DeltaAnalyzer(
        plot_delta_path,
        delta_path,
        center_ra=delta_plot_config.getboolean("center_ra"),
        z_cut_min=delta_config.getfloat("z_cut_min"),
        z_cut_max=delta_config.getfloat("z_cut_max"),
        dec_cut_min=delta_config.getfloat("dec_cut_min"),
        dec_cut_max=delta_config.getfloat("dec_cut_max"),
        ra_cut_min=delta_config.getfloat("ra_cut_min"),
        ra_cut_max=delta_config.getfloat("ra_cut_max"),
        degree=delta_plot_config.getboolean("degree"),
    )

    delta_plotter.plot(
        list(delta_plot_config.gettuplestr("value_names")),
        main_config.getstr("name"),
        histo=delta_plot_config.getboolean("plot_histo"),
        mean_z_dependence=delta_plot_config.getboolean("plot_mean_z_dependence"),
        z_dependence=delta_plot_config.getboolean("plot_z_dependence"),
        ra_dec_plots=delta_plot_config.getboolean("plot_ra_dec"),
        **delta_plot_config.getdict("plot_args"),
    )

    if delta_plot_config.getboolean("plot_comparison"):
        delta_plotter.plot(
            list(delta_plot_config.gettuplestr("value_names")),
            f"{main_config.getstr('name')}_comparison",
            comparison=delta_plot_config.getstr("comparison"),
            comparison_legend=list(delta_plot_config.gettuplestr("comparison_legend")),
            histo=delta_plot_config.getboolean("plot_histo"),
            mean_z_dependence=delta_plot_config.getboolean("plot_mean_z_dependence"),
            z_dependence=delta_plot_config.getboolean("plot_z_dependence"),
            print_stats=delta_plot_config.getboolean("print_stats"),
            **delta_plot_config.getdict("plot_args"),
        )


def launch_tomography(
    tomo_path, pixel_path, software_config, tomography_config, main_config
):
    os.makedirs(tomo_path, exist_ok=True)
    manager = task_manager.TomographyManager(
        tomo_path,
        software_config.getstr("software"),
        tomography_config.getstr("machine"),
        os.path.join(pixel_path, software_config.getstr("name_pixel")),
        os.path.join(pixel_path, f"{main_config.getstr('name')}_launch_data.pickle"),
        symlink_folder=tomography_config.getstr("symlink_folder"),
        **tomography_config.getdict("dict_launch"),
    )
    manager.launch_all()
    manager.copy()
    manager.remove_tmp()


def process_tomography(
    tomo_abs_path,
    pixel_path,
    tomography_process_config,
    delta_config,
    software_config,
    main_config,
):
    if tomography_process_config.getboolean("merge_output_maps"):
        tomography.create_merged_map(
            tomo_abs_path,
            os.path.join(
                pixel_path, f"{main_config.getstr('name')}_launch_data.pickle"
            ),
            os.path.join(tomo_abs_path, tomography_process_config.getstr("map_name")),
            os.path.join(pixel_path, delta_config.getstr("property_file_name")),
        )
    if tomography_process_config.getboolean("create_distance_map"):
        tomography.create_distance_map(
            os.path.join(
                tomo_abs_path, tomography_process_config.getstr("name_dist_map")
            ),
            os.path.join(pixel_path, software_config.getstr("name_pixel")),
            os.path.join(pixel_path, delta_config.getstr("property_file_name")),
            nb_process=tomography_process_config.getint("number_process"),
            radius_local=tomography_process_config.getfloat("radius_local"),
        )


def find_void(
    void_path,
    tomo_abs_path,
    pixel_path,
    void_finder_config,
    tomography_process_config,
    delta_config,
):
    os.makedirs(void_path, exist_ok=True)
    params_void_finder = {
        "method": void_finder_config.getstr("method_finder"),
        "threshold": void_finder_config.getfloat("threshold"),
        "average": void_finder_config.getfloat("average"),
        "minimal_radius": void_finder_config.getfloat("minimal_radius"),
        "maximal_radius": void_finder_config.getfloat("maximal_radius"),
        "radius_step": void_finder_config.getfloat("radius_step"),
        "dist_clusters": void_finder_config.getfloat("dist_clusters"),
    }

    void_finder = voidfinder.VoidFinder(
        void_path,
        os.path.join(tomo_abs_path, tomography_process_config.getstr("map_name")),
        params_void_finder,
        map_property_file=os.path.join(
            pixel_path, delta_config.getstr("property_file_name")
        ),
        number_core=void_finder_config.getint("number_process"),
        find_cluster=void_finder_config.getboolean("find_cluster"),
        split_map=void_finder_config.gettupleint("split_map"),
        split_overlap=void_finder_config.gettupleint("split_overlap"),
        delete_option=void_finder_config.getstr("delete_option"),
        restart=void_finder_config.getboolean("restart"),
    )
    void_catalog_name_default = void_finder.find_voids()
    void_finder.log.close()
    return void_catalog_name_default


def process_void(
    void_path,
    tomo_abs_path,
    pixel_path,
    void_process_config,
    delta_config,
    tomography_process_config,
    software_config,
    void_catalog_name_default,
):
    if void_process_config.getstr("void_catalog") is not None:
        void_catalog_name_process = os.path.join(
            void_path, void_process_config.getstr("void_catalog")
        )
    else:
        void_catalog_name_process = void_catalog_name_default
    if void_process_config.getboolean("compute_void_stats"):
        voidfinder.compute_additional_stats(
            void_catalog_name_process,
            os.path.join(pixel_path, software_config.getstr("name_pixel")),
        )
    if void_process_config.getboolean("cut_void_catalog"):
        void_catalog_name_cut = voidfinder.cut_catalog(
            void_path,
            void_catalog_name_process,
            void_process_config.gettuplestr("method_cut"),
            cut_crossing_param=void_process_config.getfloat("cut_crossing_param"),
            cut_radius=void_process_config.gettuplefloat("cut_radius"),
            distance_map_name=os.path.join(
                tomo_abs_path, tomography_process_config.getstr("name_dist_map")
            ),
            distance_map_prop=os.path.join(
                pixel_path, delta_config.getstr("property_file_name")
            ),
            distance_map_param=void_process_config.getfloat("distance_map_param"),
            distance_map_percent=void_process_config.getfloat("distance_map_percent"),
        )
        void_catalog_name_process = void_catalog_name_cut
    if void_process_config.getboolean("create_xcorr_catalog"):
        voidfinder.create_qso_like_catalog(void_catalog_name_process)
    void_catalog_name_default = void_catalog_name_process
    return void_catalog_name_default


def stack_void(
    stack_void_path,
    tomo_abs_path,
    void_path,
    pixel_path,
    stack_void_config,
    tomography_process_config,
    delta_config,
    void_catalog_name_default,
):
    os.makedirs(stack_void_path, exist_ok=True)
    if stack_void_config.getstr("void_catalog") is not None:
        void_catalog_name = os.path.join(
            void_path, stack_void_config.getstr("void_catalog")
        )
    else:
        void_catalog_name = void_catalog_name_default

    stack = tomography.TomographyStack(
        os.path.join(tomo_abs_path, tomography_process_config.getstr("map_name")),
        void_catalog_name,
        stack_void_config.getstr("type_catalog"),
        os.path.join(stack_void_path, stack_void_config.getstr("property_file_stack")),
        stack_void_config.getfloat("size_stack"),
        os.path.join(stack_void_path, stack_void_config.getstr("name_stack")),
        shape_stack=stack_void_config.gettupleint("shape_stack"),
        property_file=os.path.join(
            pixel_path, delta_config.getstr("property_file_name")
        ),
        coordinate_convert=stack_void_config.getstr("coordinate_convert"),
        interpolation_method=stack_void_config.getstr("interpolation_method"),
        normalized=stack_void_config.getboolean("normalized"),
    )
    stack.stack()
    return (
        os.path.join(stack_void_path, stack_void_config.getstr("name_stack")),
        os.path.join(stack_void_path, stack_void_config.getstr("property_file_stack")),
    )


def plot_void(
    plot_void_path, void_path, void_plot_config, main_config, void_catalog_name_default
):
    if void_plot_config.getstr("void_catalog") is not None:
        void_catalog_name_plot = os.path.join(
            void_path, void_plot_config.getstr("void_catalog")
        )
    else:
        void_catalog_name_plot = void_catalog_name_default

    os.makedirs(plot_void_path, exist_ok=True)
    plot = voidfinder.PlotVoid(plot_void_path, void_catalog_name_plot)

    plot.plot(
        list(void_plot_config.gettuplestr("value_names")),
        main_config.getstr("name"),
        histo=void_plot_config.getboolean("plot_histo"),
        mean_z_dependence=void_plot_config.getboolean("plot_mean_z_dependence"),
        z_dependence=void_plot_config.getboolean("plot_z_dependence"),
        **void_plot_config.getdict("plot_args"),
    )

    if void_plot_config.getboolean("plot_comparison"):
        plot.plot(
            list(void_plot_config.gettuplestr("value_names")),
            f"{main_config.getstr('name')}_comparison",
            comparison=void_plot_config.getstr("comparison"),
            comparison_legend=list(void_plot_config.gettuplestr("comparison_legend")),
            histo=void_plot_config.getboolean("plot_histo"),
            mean_z_dependence=void_plot_config.getboolean("plot_mean_z_dependence"),
            z_dependence=void_plot_config.getboolean("plot_z_dependence"),
            **void_plot_config.getdict("plot_args"),
        )


def plot_tomography(
    plot_tomography_path,
    tomo_abs_path,
    pixel_path,
    void_path,
    tomography_plot_config,
    tomography_process_config,
    software_config,
    delta_config,
    main_config,
    void_catalog_name_default,
):
    os.makedirs(plot_tomography_path, exist_ok=True)
    Treat = tomography.TomographyPlot(
        plot_tomography_path,
        map_name=os.path.join(
            tomo_abs_path, tomography_process_config.getstr("map_name")
        ),
        pixel_name=os.path.join(pixel_path, software_config.getstr("name_pixel")),
        property_file=os.path.join(
            pixel_path, delta_config.getstr("property_file_name")
        ),
        **tomography_plot_config.getdict("plot_args"),
    )
    if tomography_plot_config.getboolean("plot_map"):
        if tomography_plot_config.getstr("void_catalog") is not None:
            void_catalog_name_plot_tomo = os.path.join(
                void_path, tomography_plot_config.getstr("void_catalog")
            )
        else:
            void_catalog_name_plot_tomo = void_catalog_name_default
        try:
            center_mpc = tomography_plot_config.getfloat("center_mpc")
        except:
            center_mpc = tomography_plot_config.getstr("center_mpc")
        Treat.plot(
            main_config.getstr("name"),
            tomography_plot_config.getstr("direction"),
            tomography_plot_config.getfloat("space"),
            center_mpc,
            qso=os.path.join(pixel_path, delta_config.getstr("return_qso_catalog")),
            void=void_catalog_name_plot_tomo,
            galaxy=tomography_plot_config.getstr("galaxy"),
            distance_mask=os.path.join(
                tomo_abs_path, tomography_process_config.getstr("name_dist_map")
            ),
            criteria_distance_mask=tomography_plot_config.getfloat(
                "criteria_distance_mask"
            ),
            rotate=tomography_plot_config.getboolean("rotate"),
            minimal_void_crossing=tomography_plot_config.getfloat(
                "minimal_void_crossing"
            ),
            redshift_axis=tomography_plot_config.getboolean("redshift_axis"),
            cut_plot=tomography_plot_config.gettuplefloat("cut_plot"),
        )

    if tomography_plot_config.getboolean("plot_delta_histogram"):
        Treat.plot_delta_histogram(
            f"{main_config.getstr('name')}_histogram_deltas",
            tomography_plot_config.getint("nb_bins"),
            gauss_fit=tomography_plot_config.getboolean("gauss_fit"),
            norm=tomography_plot_config.getboolean("normalization"),
            distance_mask=os.path.join(
                tomo_abs_path, tomography_process_config.getstr("name_dist_map")
            ),
            criteria_distance_mask=tomography_plot_config.getfloat(
                "criteria_distance_mask_histo"
            ),
            log_scale=tomography_plot_config.getboolean("log_scale"),
        )

    if tomography_plot_config.getboolean("plot_delta_histogram_comparison"):
        Treat.plot_delta_histogram_comparison(
            f"{main_config.getstr('name')}_histogram_deltas",
            tomography_plot_config.getstr("map_comparison"),
            tomography_plot_config.getint("nb_bins"),
            list(tomography_plot_config.gettuplestr("legend_comparison")),
            gauss_fit=tomography_plot_config.getboolean("gauss_fit"),
            norm=tomography_plot_config.getboolean("normalization"),
            distance_mask=os.path.join(
                tomo_abs_path, tomography_process_config.getstr("name_dist_map")
            ),
            distance_mask2=tomography_plot_config.getstr("distance_mask_comparison"),
            criteria_distance_mask=tomography_plot_config.getfloat(
                "criteria_distance_mask_histo"
            ),
            log_scale=tomography_plot_config.getboolean("log_scale"),
        )

    if tomography_plot_config.getboolean("plot_integrated_map"):
        Treat.plot_integrate_image(
            tomography_plot_config.getfloat("zmin_integrated_map"),
            tomography_plot_config.getfloat("zmax_integrated_map"),
            f"{main_config.getstr('name')}_integrated_map",
            rotate=tomography_plot_config.getboolean("rotate"),
        )

    if tomography_plot_config.getboolean("plot_centered_maps"):
        Treat.plot_catalog_centered_maps(
            tomography_plot_config.getstr("direction"),
            f"{main_config.getstr('name')}_centered",
            tomography_plot_config.getfloat("space"),
            os.path.join(void_path, tomography_plot_config.getstr("catalog_centered")),
            tomography_plot_config.getint("nb_plot"),
            tomography_plot_config.getfloat("radius_centered"),
            qso=os.path.join(pixel_path, delta_config.getstr("return_qso_catalog")),
            rotate=tomography_plot_config.getboolean("rotate"),
        )


def plot_stack_void(
    plot_stack_void_path,
    stack_void_path,
    stack_void_plot_config,
    stack_void_name_default,
    property_stack_void_name_default,
):
    os.makedirs(plot_stack_void_path, exist_ok=True)

    if stack_void_plot_config.getstr("stack_name") is not None:
        name_stack = os.path.join(
            stack_void_path, stack_void_plot_config.getstr("stack_name")
        )
    else:
        name_stack = stack_void_name_default

    if stack_void_plot_config.getstr("property_stack_name") is not None:
        property_file_stack = os.path.join(
            stack_void_path, stack_void_plot_config.getstr("stack_name")
        )
    else:
        property_file_stack = property_stack_void_name_default

    tomography.TomographyStack.plot_stack(
        plot_stack_void_path,
        name_stack,
        property_file_stack,
        os.path.join(plot_stack_void_path, stack_void_plot_config.getstr("name_plot")),
        ellipticity=stack_void_plot_config.getboolean("ellipticity"),
        pixel_file_qso_distance=stack_void_plot_config.getstr(
            "pixel_file_qso_distance"
        ),
        **stack_void_plot_config.getdict("plot_args"),
    )
