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


import ast
import configparser
import os

from lelantos import cosmology, task_manager, tomography, utils, voidfinder

#############################################################################
#############################################################################
############################### FUNCTIONS ###################################
#############################################################################
#############################################################################


def parse_int_tuple(input):
    """Parse a config string into a tuple of ints.

    Args:
        input (str): Comma-separated integers (e.g. ``"2,2,300"``) or the
            literal ``"None"``.

    Returns:
        tuple[int] | None: Tuple of ints, or ``None`` if ``input == "None"``.
    """
    if input == "None":
        return None
    else:
        return tuple(int(k.strip()) for k in input.strip().split(","))


def parse_float_tuple(input):
    """Parse a config string into a tuple of floats.

    Args:
        input (str): Comma-separated floats (e.g. ``"7,50"``) or ``"None"``.

    Returns:
        tuple[float] | None: Tuple of floats, or ``None`` if ``input == "None"``.
    """
    if input == "None":
        return None
    else:
        return tuple(float(k.strip()) for k in input.strip().split(","))


def parse_str_tuple(input):
    """Parse a config string into a tuple of stripped strings.

    Args:
        input (str): Comma-separated tokens (e.g. ``"snr , delta"``) or ``"None"``.

    Returns:
        tuple[str] | None: Tuple of strings, or ``None`` if ``input == "None"``.
    """
    if input == "None":
        return None
    else:
        return tuple(str(k.strip()) for k in input.strip().split(","))


def parse_dict(input):
    """Parse a config string into a Python dict.

    Single quotes are normalised to double quotes and the value is evaluated
    with :func:`ast.literal_eval`.

    Args:
        input (str): Python-literal dict (e.g. ``{"n": 8}``) or ``"None"``.

    Returns:
        dict | None: Parsed dict, or ``None`` if ``input == "None"``.
    """
    if input == "None":
        return None
    else:
        acceptable_string = input.replace("'", '"')
        return ast.literal_eval(acceptable_string)


def parse_float(input):
    """Parse a config string into a float.

    Args:
        input (str): A float literal or ``"None"``.

    Returns:
        float | None: The float, or ``None`` if ``input == "None"``.
    """
    if input == "None":
        return None
    else:
        return float(input)


def parse_int(input):
    """Parse a config string into an int.

    Args:
        input (str): An integer literal or ``"None"``.

    Returns:
        int | None: The int, or ``None`` if ``input == "None"``.
    """
    if input == "None":
        return None
    else:
        return int(input)


def parse_string(input):
    """Parse a config string, mapping the literal ``"None"`` to ``None``.

    Args:
        input (str): Any string, or ``"None"``.

    Returns:
        str | None: The string, or ``None`` if ``input == "None"``.
    """
    if input == "None":
        return None
    else:
        return str(input)


def main(input_file):
    """Run the full lelantos tomography pipeline from a config file.

    Reads the ``.ini`` file, then executes each enabled stage in order:
    delta transform, delta plot, delta convert, tomography launch, tomography
    process, void find, void process, void stack, and the void/tomography/stack
    plots. The set of stages is controlled by the boolean gates in the
    ``[main]`` section. See ``scripts/interface_explanatory.ini`` for a full
    description of every option.

    Args:
        input_file (str): Path to the interface ``.ini`` configuration file.
    """
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
    """Shuffle and redshift-cut the delta files (``[delta transform]`` stage).

    Wraps :class:`lelantos.cosmology.DeltaModifier` to shuffle the deltas within
    ``n_cut`` redshift chunks, optionally in lock-step with a second
    (cross-correlation) delta set.

    Args:
        delta_path (str): Input delta directory.
        delta_transform_path (str): Output directory for the shuffled deltas.
        delta_transform_config: ``[delta transform]`` config section.
        delta_config: ``[delta convert]`` config section (for the z window).
        main_config: ``[main]`` config section.
    """
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
    """Convert delta files to solver pixel input (``[delta convert]`` stage).

    Wraps :class:`lelantos.cosmology.DeltaConverter` to apply the footprint /
    redshift / sigma cuts, perform the sky->cartesian transform and write the
    binary pixel file(s) and property file consumed by the tomographic solver.

    Args:
        pixel_path (str): Output directory for pixel/binary products.
        delta_path (str): Input delta directory.
        software_config: ``[tomography software]`` config section.
        delta_config: ``[delta convert]`` config section.
        main_config: ``[main]`` config section.
    """
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
    """Produce delta diagnostic plots (``[delta plot]`` stage).

    Wraps :class:`lelantos.cosmology.DeltaAnalyzer` to plot histograms, sky
    scatter and redshift-dependence of the requested delta quantities, with an
    optional comparison against a second delta set.

    Args:
        plot_delta_path (str): Output directory for the figures.
        delta_path (str): Input delta directory.
        delta_plot_config: ``[delta plot]`` config section.
        delta_config: ``[delta convert]`` config section (for the cuts).
        main_config: ``[main]`` config section.
    """
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
    """Run the tomographic solver (``[tomography launching]`` stage).

    Wraps :class:`lelantos.task_manager.TomographyManager` to launch the solver
    (e.g. Dachshund) on the pixel file, wait for completion, copy the outputs
    back and clean the temporary run directory.

    Args:
        tomo_path (str): Run/output directory for the solver.
        pixel_path (str): Directory holding the pixel input and launch pickle.
        software_config: ``[tomography software]`` config section.
        tomography_config: ``[tomography launching]`` config section.
        main_config: ``[main]`` config section.
    """
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
    """Post-process solver output (``[tomography process]`` stage).

    Optionally merges the per-chunk solver outputs into a single reconstructed
    map and/or builds the voxel-to-nearest-line-of-sight distance map used for
    void cuts and plot masking.

    Args:
        tomo_abs_path (str): Absolute path to the tomography run directory.
        pixel_path (str): Directory holding the pixel input and property file.
        tomography_process_config: ``[tomography process]`` config section.
        delta_config: ``[delta convert]`` config section.
        software_config: ``[tomography software]`` config section.
        main_config: ``[main]`` config section.
    """
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
    """Run the void finder on the map (``[void finder]`` stage).

    Wraps :class:`lelantos.voidfinder.VoidFinder` to detect voids
    (SPHERICAL or WATERSHED) on the reconstructed map and write a void catalog.

    Args:
        void_path (str): Output directory for the void catalog.
        tomo_abs_path (str): Absolute path to the tomography directory (map).
        pixel_path (str): Directory holding the pixel input and property file.
        void_finder_config: ``[void finder]`` config section.
        tomography_process_config: ``[tomography process]`` config section.
        delta_config: ``[delta convert]`` config section.

    Returns:
        str: Path of the void catalog produced by the finder.
    """
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
    """Post-process the void catalog (``[void process]`` stage).

    Optionally computes additional per-void statistics, applies the configured
    cuts (CROSSING / RADIUS / BORDER / DIST) writing a cut catalog, and emits a
    QSO-like catalog for cross-correlation tools.

    Args:
        void_path (str): Void catalog directory.
        tomo_abs_path (str): Absolute path to the tomography directory.
        pixel_path (str): Directory holding the pixel input and property file.
        void_process_config: ``[void process]`` config section.
        delta_config: ``[delta convert]`` config section.
        tomography_process_config: ``[tomography process]`` config section.
        software_config: ``[tomography software]`` config section.
        void_catalog_name_default (str | None): Catalog produced by ``find_void``.

    Returns:
        str: Path of the (possibly cut) void catalog to use downstream.
    """
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
    """Stack the map at catalog positions (``[void stack]`` stage).

    Wraps :class:`lelantos.tomography.TomographyStack` to stack map cut-outs at
    each void/QSO/galaxy position and write the stack binary and its property
    file.

    Args:
        stack_void_path (str): Output directory for the stack products.
        tomo_abs_path (str): Absolute path to the tomography directory (map).
        void_path (str): Void catalog directory.
        pixel_path (str): Directory holding the property file.
        stack_void_config: ``[void stack]`` config section.
        tomography_process_config: ``[tomography process]`` config section.
        delta_config: ``[delta convert]`` config section.
        void_catalog_name_default (str | None): Catalog produced earlier.

    Returns:
        tuple[str, str]: ``(stack_binary_path, stack_property_file_path)``.
    """
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
    """Produce void-catalog diagnostic plots (``[void plot]`` stage).

    Wraps :class:`lelantos.voidfinder.PlotVoid` to plot histograms and
    redshift-dependence of void quantities (radius, redshift, central value,
    LOS distance), with an optional comparison catalog.

    Args:
        plot_void_path (str): Output directory for the figures.
        void_path (str): Void catalog directory.
        void_plot_config: ``[void plot]`` config section.
        main_config: ``[main]`` config section.
        void_catalog_name_default (str | None): Catalog produced earlier.
    """
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
    """Produce tomographic-map plots (``[tomography plot]`` stage).

    Wraps :class:`lelantos.tomography.TomographyPlot` to draw map slices with
    optional QSO/void/galaxy overlays, delta histograms (single and comparison),
    the redshift-integrated map, and catalog-centred maps, according to the
    enabled ``plot_*`` flags.

    Args:
        plot_tomography_path (str): Output directory for the figures.
        tomo_abs_path (str): Absolute path to the tomography directory (map).
        pixel_path (str): Directory holding the pixel input and property file.
        void_path (str): Void catalog directory.
        tomography_plot_config: ``[tomography plot]`` config section.
        tomography_process_config: ``[tomography process]`` config section.
        software_config: ``[tomography software]`` config section.
        delta_config: ``[delta convert]`` config section.
        main_config: ``[main]`` config section.
        void_catalog_name_default (str | None): Catalog produced earlier.
    """
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

    if tomography_plot_config.getstr("void_catalog") is not None:
        void_catalog_name_plot_tomo = os.path.join(
            void_path, tomography_plot_config.getstr("void_catalog")
        )
    else:
        void_catalog_name_plot_tomo = void_catalog_name_default

    if tomography_plot_config.getboolean("plot_map"):
        try:
            center_mpc = tomography_plot_config.getfloat("center_mpc")
        except:
            center_mpc = tomography_plot_config.getstr("center_mpc")
        if tomography_process_config.getstr("name_dist_map") is not None:
            name_dist_map = os.path.join(
                tomo_abs_path, tomography_process_config.getstr("name_dist_map")
            )
        else:
            name_dist_map = None
        Treat.plot(
            main_config.getstr("name"),
            tomography_plot_config.getstr("direction"),
            tomography_plot_config.getfloat("space"),
            center_mpc,
            qso=os.path.join(pixel_path, delta_config.getstr("return_qso_catalog")),
            void=void_catalog_name_plot_tomo,
            galaxy=tomography_plot_config.getstr("galaxy"),
            distance_mask=name_dist_map,
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
            cut_plot=tomography_plot_config.gettuplefloat("cut_plot"),
            void=void_catalog_name_plot_tomo,
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
    """Plot the void stack (``[stack void plot]`` stage).

    Wraps :meth:`lelantos.tomography.TomographyStack.plot_stack` to plot the
    stacked map, optionally with an ellipticity fit and QSO-distance overlay.

    Note:
        This stage reads the keys ``stack_name`` and ``property_stack_name``,
        which differ from the ``name_stack`` / ``property_file_stack`` keys
        written by the ``[void stack]`` stage (see interface_explanatory.ini).

    Args:
        plot_stack_void_path (str): Output directory for the figures.
        stack_void_path (str): Directory holding the stack products.
        stack_void_plot_config: ``[stack void plot]`` config section.
        stack_void_name_default (str | None): Stack binary produced earlier.
        property_stack_void_name_default (str | None): Stack property file
            produced earlier.
    """
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


def print_approximate_shape_size_from_interface_file(input_file):
    """Estimate and log the map shape/size implied by a config file.

    Reads the ``[delta convert]`` geometry (footprint, redshift window,
    Omega_m, chunking, overlap, sub-map shape) and calls
    :func:`lelantos.cosmology.compute_shape_size_parallel_from_interface` to
    report the approximate full-map pixel shape, physical size (Mpc.h^-1) and
    the Mpc-per-pixel / pixel-per-Mpc resolution. Useful to size a run before
    launching it.

    Args:
        input_file (str): Path to the interface ``.ini`` configuration file.
    """
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
    delta_config = config["delta convert"]

    ramin = delta_config.getfloat("ra_cut_min")
    ramax = delta_config.getfloat("ra_cut_max")
    decmin = delta_config.getfloat("dec_cut_min")
    decmax = delta_config.getfloat("dec_cut_max")
    zmin = delta_config.getfloat("z_cut_min")
    zmax = delta_config.getfloat("z_cut_max")
    coordinate_transform = delta_config.getstr("coordinate_transform")
    number_chunks = delta_config.gettupleint("number_chunks")
    overlaping = delta_config.getfloat("overlaping")
    shape_sub_map = delta_config.gettupleint("shape_sub_map")
    Omega_m = delta_config.getfloat("Omega_m")

    shape, size = cosmology.compute_shape_size_parallel_from_interface(
        ramin,
        ramax,
        decmin,
        decmax,
        zmin,
        zmax,
        Omega_m,
        coordinate_transform,
        number_chunks,
        overlaping,
        shape_sub_map,
        N_coord_edge=100,
    )

    log = utils.create_log()
    log.add(f"Approximate shape of the associated map: {shape}")
    log.add(f"Approximate size of the associated map: {size}")
    log.add(f"Mpc per pixels: {utils.mpc_per_pixel(size, shape)}")
    log.add(f"Pixels per mpc: {utils.pixel_per_mpc(size, shape)}")
    log.close()
