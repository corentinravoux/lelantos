#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Corentin Ravoux

Description : Classes used to treat input and output of the Dachshund software.
Can be used to plot the slices of a Tomographic map or to return a Paraview
file for 3D visualization. Can be also used to compare a simulation to a data
analysis.
"""


#############################################################################
#############################################################################
########################## MODULE IMPORTATION ###############################
#############################################################################
#############################################################################

import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from lelantos import tomographic_objects, utils


def create_merged_map(submap_directory, launching_file_name, map_name, property_file):
    """Merge per-chunk solver output maps into a single map and write it.

    Args:
        submap_directory (str): Directory holding the per-chunk maps.
        launching_file_name (str): Launch pickle describing the chunk layout.
        map_name (str): Output merged-map file name.
        property_file (str): Map property/pickle file.
    """
    map_merged = tomographic_objects.TomographicMap.init_by_merging(
        submap_directory, launching_file_name, map_name, property_file
    )
    map_merged.write()


def rebin_map(
    map_name, property_file, new_shape, new_name, new_prop_name, operation="mean"
):
    """Rebin a tomographic map to a new pixel shape and write it out.

    Args:
        map_name (str): Input map file.
        property_file (str): Input map property file.
        new_shape (tuple[int]): Target pixel shape.
        new_name (str): Output map file.
        new_prop_name (str): Output property file.
        operation (str, optional): Rebin reduction (``"mean"``, ``"sum"`` ...).
    """
    map_class = tomographic_objects.TomographicMap.init_from_property_files(
        property_file, name=map_name
    )
    map_class.read()
    map_class.rebin_map(new_shape, operation=operation)
    map_class.name = new_name
    map_class.write_property_file(new_prop_name)
    map_class.write()


def create_distance_map(
    map_name, pixel_file, property_file, nb_process=1, radius_local=50
):
    """Compute and write a voxel-to-nearest-line-of-sight distance map.

    Args:
        map_name (str): Output distance-map file.
        pixel_file (str): Pixel (line-of-sight) file.
        property_file (str): Map property file defining the grid.
        nb_process (int, optional): Number of parallel processes.
        radius_local (float, optional): Local search radius (Mpc.h^-1).
    """
    pixel = tomographic_objects.Pixel(name=pixel_file)
    pixel.read()
    tomographic_map = tomographic_objects.TomographicMap.init_from_property_files(
        property_file
    )
    distance_map = tomographic_objects.DistanceMap.init_by_computing(
        pixel,
        tomographic_map,
        map_name,
        nb_process=nb_process,
        radius_local=radius_local,
    )
    distance_map.write()


def convert_to_vtk(map_name, property_file, new_name):
    """Export a tomographic map to VTK for 3D visualisation (needs pyevtk).

    Args:
        map_name (str): Input map file.
        property_file (str): Map property file.
        new_name (str): Output VTK base name.
    """
    map_class = tomographic_objects.TomographicMap.init_from_property_files(
        property_file, name=map_name
    )
    map_class.read()
    map_class.name = new_name
    map_class.write_in_vtk()


def mask_map_to_3d(map_name, property_file, new_name, distance_map, distance):
    """Mask a map beyond a distance-to-LOS threshold and write it out.

    Args:
        map_name (str): Input map file.
        property_file (str): Map property file.
        new_name (str): Output map file.
        distance_map (str): Distance-map file used for masking.
        distance (float): Distance-to-LOS threshold (Mpc.h^-1).
    """
    map_class = tomographic_objects.TomographicMap.init_from_property_files(
        property_file, name=map_name
    )
    map_class.read()
    map_class.name = new_name
    map_class.mask_map(distance_map, distance)
    map_class.write()


def pixel_to_3d(pixel_name, new_name):
    """Export a pixel (line-of-sight) file to text for 3D visualisation.

    Args:
        pixel_name (str): Input pixel file.
        new_name (str): Output text file.
    """
    pixel = tomographic_objects.Pixel(name=pixel_name)
    pixel.read()
    pixel.writetxt(new_name)


#############################################################################
#############################################################################
############################### CLASSES #####################################
#############################################################################
#############################################################################


class TomographyPlot(object):
    """Plot slices, histograms and derived products of a tomographic map.

    Handles loading the map plus optional QSO/void/galaxy catalogs and distance
    mask, and drawing map slices (with overlays and an optional redshift axis),
    delta histograms, integrated maps, 2D shells, catalog-centred maps and 3D
    power spectra.
    """

    def __init__(
        self,
        pwd,
        map_name=None,
        map_shape=None,
        pixel_name=None,
        property_file=None,
        **kwargs,
    ):
        """Store paths and styling for the plots.

        Args:
            pwd (str): Output directory for the figures.
            map_name (str, optional): Map binary file.
            map_shape (tuple[int], optional): Map pixel shape (if no property file).
            pixel_name (str, optional): Pixel (line-of-sight) file to overlay.
            property_file (str, optional): Map property/pickle file.
            **kwargs: Default plot styling (``plot_args``); ``style`` selects a
                matplotlib style.
        """
        self.pwd = pwd
        self.map_name = map_name
        self.map_shape = map_shape
        self.pixel_name = pixel_name
        self.property_file = property_file
        self.kwargs = kwargs
        style = utils.return_key(kwargs, "style", None)
        if style is not None:
            plt.style.use(style)

    def load_tomographic_objects(
        self, qso=None, void=None, galaxy=None, distance_mask=None, cut_plot=None
    ):
        """Load the map and any requested overlay catalogs / distance mask.

        Args:
            qso (str, optional): QSO catalog file.
            void (str, optional): Void catalog file.
            galaxy (str, optional): Galaxy catalog file.
            distance_mask (str, optional): Distance-map file.
            cut_plot (sequence, optional): Fractional sub-box crop applied to
                every loaded object.

        Returns:
            tuple: ``(tomographic_map, pixel, quasar_catalog, void_catalog,
            galaxy_catalog, dist_map)`` (unused entries are None).
        """
        tomographic_map = tomographic_objects.TomographicMap.init_classic(
            name=self.map_name, shape=self.map_shape, property_file=self.property_file
        )
        tomographic_map.read()

        pixel, quasar_catalog, void_catalog, galaxy_catalog, dist_map = (
            None,
            None,
            None,
            None,
            None,
        )
        if self.pixel_name is not None:
            if self.property_file is not None:
                pixel = tomographic_objects.Pixel.init_from_property_files(
                    self.property_file, name=self.pixel_name
                )
            else:
                pixel = tomographic_objects.Pixel(name=self.pixel_name)
            pixel.read()
        if qso is not None:
            quasar_catalog = tomographic_objects.QSOCatalog.init_from_fits(qso)
        if void is not None:
            void_catalog = tomographic_objects.VoidCatalog.init_from_fits(void)
        if galaxy is not None:
            galaxy_catalog = tomographic_objects.GalaxyCatalog.init_from_fits(galaxy)
        if distance_mask is not None:
            dist_map = tomographic_objects.DistanceMap.init_from_tomographic_map(
                tomographic_map, name=distance_mask
            )
            dist_map.read()

        if cut_plot is not None:
            self.cut_objects(
                cut_plot,
                tomographic_map,
                pixel,
                quasar_catalog,
                void_catalog,
                galaxy_catalog,
                dist_map,
            )

        return (
            tomographic_map,
            pixel,
            quasar_catalog,
            void_catalog,
            galaxy_catalog,
            dist_map,
        )

    def cut_objects(
        self,
        cut_plot,
        tomographic_map,
        pixel,
        quasar_catalog,
        void_catalog,
        galaxy_catalog,
        dist_map,
    ):
        """Crop the map and all overlay objects to a fractional sub-box.

        Args:
            cut_plot (sequence): Per-axis fractions in ``[0, 1]`` defining the
                sub-box to keep.
            tomographic_map: The map (cropped in place).
            pixel: Pixel object (masked in place) or None.
            quasar_catalog: QSO catalog (masked in place) or None.
            void_catalog: Void catalog (masked in place) or None.
            galaxy_catalog: Galaxy catalog (masked in place) or None.
            dist_map: Distance map (cropped in place) or None.
        """
        size_map = np.array(cut_plot) * tomographic_map.size
        shape_map = np.round(np.array(cut_plot) * tomographic_map.shape, 0).astype(int)
        tomographic_map.shape = shape_map
        tomographic_map.size = size_map
        tomographic_map.map_array = tomographic_map.map_array[
            0 : shape_map[0], 0 : shape_map[1], 0 : shape_map[2]
        ]

        if dist_map is not None:
            dist_map.shape = shape_map
            dist_map.size = size_map
            dist_map.map_array = dist_map.map_array[
                0 : shape_map[0], 0 : shape_map[1], 0 : shape_map[2]
            ]

        if pixel is not None:
            mask = (
                (pixel.pixel_array[:, 0] < size_map[0])
                & (pixel.pixel_array[:, 1] < size_map[1])
                & (pixel.pixel_array[:, 2] < size_map[2])
            )
            pixel.pixel_array = pixel.pixel_array[mask]

        if quasar_catalog is not None:
            mask = (
                (quasar_catalog.coord[:, 0] < size_map[0])
                & (quasar_catalog.coord[:, 1] < size_map[1])
                & (quasar_catalog.coord[:, 2] < size_map[2])
            )
            quasar_catalog.apply_mask(mask)

        if void_catalog is not None:
            mask = (
                (void_catalog.coord[:, 0] < size_map[0])
                & (void_catalog.coord[:, 1] < size_map[1])
                & (void_catalog.coord[:, 2] < size_map[2])
            )
            void_catalog.apply_mask(mask)

        if galaxy_catalog is not None:
            mask = (
                (galaxy_catalog.coord[:, 0] < size_map[0])
                & (galaxy_catalog.coord[:, 1] < size_map[1])
                & (galaxy_catalog.coord[:, 2] < size_map[2])
            )
            galaxy_catalog.apply_mask(mask)

    def mask_tomographic_objects(
        self,
        void_catalog,
        dist_map,
        tomographic_map,
        criteria_distance_mask=None,
        minimal_void_crossing=None,
    ):
        """Mask the map by distance-to-LOS and flag well-crossed voids.

        Args:
            void_catalog: Void catalog or None.
            dist_map: Distance map or None.
            tomographic_map: The map (masked in place).
            criteria_distance_mask (float, optional): Distance-to-LOS threshold.
            minimal_void_crossing (float, optional): Minimum void crossing to
                keep a void in the "in" set.

        Returns:
            numpy.ndarray | None: Boolean mask of the kept voids, or None.
        """
        mask_void = None
        if void_catalog is not None:
            if minimal_void_crossing is not None:
                crossing_param = void_catalog.crossing_param
                mask_void = crossing_param > minimal_void_crossing
        if dist_map is not None:
            mask_distance = dist_map.get_mask_distance(criteria_distance_mask)
            tomographic_map.mask_map_from_mask(mask_distance)
        return mask_void

    def delta_histogram(self, listdeltas, nb_bins, norm=True, gauss_fit=True, alpha=1):
        """Histogram map delta values and optionally overlay a Gaussian fit.

        Args:
            listdeltas (array-like): Flattened map delta values.
            nb_bins (int): Number of histogram bins.
            norm (bool, optional): Normalise to a density.
            gauss_fit (bool, optional): Fit and draw a Gaussian.
            alpha (float, optional): Histogram transparency.

        Returns:
            tuple | None: ``(mu, sigma)`` of the Gaussian fit if ``gauss_fit``.
        """
        data, bins, patches = plt.hist(
            listdeltas, nb_bins, density=norm, alpha=alpha, range=(-0.4, 0.4)
        )
        if gauss_fit:
            bin_centers = np.array(
                [0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)]
            )
            fit_function = lambda x, A, mu, sigma: A * np.exp(
                -1.0 * (x - mu) ** 2 / (2 * sigma**2)
            )
            popt, pcov = curve_fit(
                fit_function, xdata=bin_centers, ydata=data, p0=[1, 0.0, 0.1]
            )
            x = np.linspace(min(bins), max(bins), 1000)
            y = fit_function(x, *popt)
            mu, sigma = popt[1], popt[2]
            plt.plot(x, y, "r--", linewidth=2)
            return (mu, sigma)

    def select_pixel_in(self, center_mpc, space_mpc, pixel, index_direction):
        """Mask pixels within half a slab thickness of the slice centre.

        Args:
            center_mpc (float): Slice centre along the slice-normal axis.
            space_mpc (float): Slab thickness (Mpc.h^-1).
            pixel: Pixel object.
            index_direction (int): Slice-normal axis index.

        Returns:
            numpy.ndarray: Boolean mask of the selected pixels.
        """
        mask = pixel.pixel_array[:, index_direction] < center_mpc + space_mpc / 2
        mask &= pixel.pixel_array[:, index_direction] >= center_mpc - space_mpc / 2
        return mask

    def select_qso_in(self, center_mpc, space_mpc, quasar_catalog, index_direction):
        """Mask QSOs within half a slab thickness of the slice centre.

        Args:
            center_mpc (float): Slice centre along the slice-normal axis.
            space_mpc (float): Slab thickness (Mpc.h^-1).
            quasar_catalog: QSO catalog.
            index_direction (int): Slice-normal axis index.

        Returns:
            numpy.ndarray: Boolean mask of the selected QSOs.
        """
        mask = quasar_catalog.coord[:, index_direction] < center_mpc + space_mpc / 2
        mask &= quasar_catalog.coord[:, index_direction] >= center_mpc - space_mpc / 2
        return mask

    def select_void_in(self, center_mpc, space_mpc, void, index_direction):
        """Mask voids whose sphere intersects the slice plane.

        Args:
            center_mpc (float): Slice centre along the slice-normal axis.
            space_mpc (float): Slab thickness (unused; radius is used instead).
            void: Void catalog.
            index_direction (int): Slice-normal axis index.

        Returns:
            numpy.ndarray: Boolean mask of the intersecting voids.
        """
        mask = void.coord[:, index_direction] - void.radius < center_mpc
        mask &= void.coord[:, index_direction] + void.radius >= center_mpc
        return mask

    def select_galaxy_in(self, center_mpc, space_mpc, galaxy, index_direction):
        """Mask galaxies within half a slab thickness of the slice centre.

        Args:
            center_mpc (float): Slice centre along the slice-normal axis.
            space_mpc (float): Slab thickness (Mpc.h^-1).
            galaxy: Galaxy catalog.
            index_direction (int): Slice-normal axis index.

        Returns:
            numpy.ndarray: Boolean mask of the selected galaxies.
        """
        mask = galaxy.coord[:, index_direction] < center_mpc + space_mpc / 2
        mask &= galaxy.coord[:, index_direction] >= center_mpc - space_mpc / 2
        return mask

    @staticmethod
    def get_direction_informations(direction, rotate, size_map):
        """Return plotting metadata for a slicing direction.

        Args:
            direction (str): Slicing axis (``x``/``y``/``z`` or ra/dec/redshift).
            rotate (bool): Whether the slice is rotated.
            size_map (array-like): Map physical size per axis.

        Returns:
            tuple: ``(x_index, y_index, index_direction, extentmap, xlab, ylab)``.
        """
        x_index, y_index, index_direction, index_dict = utils.get_direction_indexes(
            direction, rotate
        )

        if rotate:
            extentmap = [0, size_map[x_index], 0, size_map[y_index]]
        else:
            extentmap = [0, size_map[x_index], size_map[y_index], 0]

        xlab = (
            f"Relative comoving distance in the {index_dict[f'{index_dict[x_index]}_lab']} direction ["
            + r"$h^{-1}$"
            + r"$\cdot$"
            + "Mpc"
            + "]"
        )
        ylab = (
            f"Relative comoving distance in the {index_dict[f'{index_dict[y_index]}_lab']} direction ["
            + r"$h^{-1}$"
            + r"$\cdot$"
            + "Mpc"
            + "]"
        )

        return (x_index, y_index, index_direction, extentmap, xlab, ylab)

    def print_one_slice(
        self,
        name,
        tomographic_map,
        direction,
        space_mpc,
        center_mpc,
        pixel=None,
        quasar_catalog=None,
        void_catalog=None,
        mask_void=None,
        galaxy_catalog=None,
        rotate=False,
        redshift_axis=False,
    ):
        """Extract one map slice and draw it with the selected overlays.

        Selects the pixels/QSOs/voids/galaxies within the slab, computes each
        void's effective (in-plane) radius and calls :meth:`plot_slice`.

        Args:
            name (str): Base output name.
            tomographic_map: The map.
            direction (str): Slicing axis.
            space_mpc (float): Slab thickness (Mpc.h^-1).
            center_mpc (float): Slice centre along the slice-normal axis.
            pixel: Pixel object or None.
            quasar_catalog: QSO catalog or None.
            void_catalog: Void catalog or None.
            mask_void (numpy.ndarray, optional): Which voids are "primary".
            galaxy_catalog: Galaxy catalog or None.
            rotate (bool, optional): Rotate the slice.
            redshift_axis (bool, optional): Draw the secondary redshift axis.
        """
        size_map = tomographic_map.size
        pixel_per_mpc = tomographic_map.pixel_per_mpc

        (
            x_index,
            y_index,
            index_direction,
            extentmap,
            xlab,
            ylab,
        ) = TomographyPlot.get_direction_informations(direction, rotate, size_map)

        center_pix = int(round(center_mpc * pixel_per_mpc[index_direction], 0))
        if direction == "x":
            map_slice = tomographic_map.map_array[center_pix, :, :]
            if rotate:
                map_slice = np.transpose(
                    np.flip(tomographic_map.map_array[center_pix, :, :], axis=1)
                )
        elif direction == "y":
            map_slice = tomographic_map.map_array[:, center_pix, :]
            if rotate:
                map_slice = np.transpose(
                    np.flip(tomographic_map.map_array[:, center_pix, :], axis=1)
                )
        elif direction == "z":
            map_slice = tomographic_map.map_array[:, :, center_pix]
            if rotate:
                map_slice = np.transpose(
                    np.flip(tomographic_map.map_array[:, :, center_pix], axis=1)
                )

        pixel_in, pixel_bis_in, qso_in, qso_bis_in, void_in, void_bis_in, galaxy_in = (
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        if pixel is not None:
            mask_pixel_in = self.select_pixel_in(
                center_mpc, space_mpc, pixel, index_direction
            )
            pixel_in = pixel.pixel_array[mask_pixel_in]
            mask_pixel_bis_in = (
                self.select_pixel_in(center_mpc, 2 * space_mpc, pixel, index_direction)
            ) & (~mask_pixel_in)
            pixel_bis_in = pixel.pixel_array[mask_pixel_bis_in]
        if quasar_catalog is not None:
            mask_qso_in = self.select_qso_in(
                center_mpc, space_mpc, quasar_catalog, index_direction
            )
            qso_in = quasar_catalog.coord[mask_qso_in]
            mask_qso_bis_in = (
                self.select_qso_in(
                    center_mpc, 2 * space_mpc, quasar_catalog, index_direction
                )
            ) & (~mask_qso_in)
            qso_bis_in = quasar_catalog.coord[mask_qso_bis_in]
        if void_catalog is not None:
            mask_void_in = self.select_void_in(
                center_mpc, space_mpc, void_catalog, index_direction
            )
            void_coord_in = void_catalog.coord[mask_void_in]
            void_radius_in = void_catalog.radius[mask_void_in]
            effective_radius = np.sqrt(
                void_radius_in**2
                - (void_coord_in[:, index_direction] - center_mpc) ** 2
            )
            void_in = np.transpose(
                np.vstack(
                    [
                        void_coord_in[:, 0],
                        void_coord_in[:, 1],
                        void_coord_in[:, 2],
                        effective_radius,
                    ]
                )
            )
            if mask_void is not None:
                void_bis_in = void_in[~mask_void]
                void_in = void_in[mask_void]
        if galaxy_catalog is not None:
            mask_galaxy_in = self.select_galaxy_in(
                center_mpc, space_mpc, galaxy_catalog, index_direction
            )
            galaxy_in = np.transpose(
                np.vstack(
                    [
                        galaxy_catalog.coord[:, 0],
                        galaxy_catalog.coord[:, 1],
                        galaxy_catalog.coord[:, 2],
                        galaxy_catalog.standard_deviation,
                    ]
                )
            )[mask_galaxy_in]

        name_plot = f"{name}_direction_{direction}_mpc_{center_mpc}"
        TomographyPlot.plot_slice(
            self.pwd,
            map_slice,
            extentmap,
            xlab,
            ylab,
            name_plot,
            x_index,
            y_index,
            pixel_in=pixel_in,
            pixel_bis_in=pixel_bis_in,
            qso_in=qso_in,
            qso_bis_in=qso_bis_in,
            void_in=void_in,
            void_bis_in=void_bis_in,
            galaxy_in=galaxy_in,
            redshift_axis=redshift_axis,
            tomographic_map=tomographic_map,
            rotate=rotate,
            **self.kwargs,
        )

    @staticmethod
    def plot_slice(
        pwd,
        map_slice,
        extentmap,
        xlab,
        ylab,
        name,
        x_index,
        y_index,
        pixel_in=None,
        pixel_bis_in=None,
        qso_in=None,
        qso_bis_in=None,
        void_in=None,
        void_bis_in=None,
        galaxy_in=None,
        redshift_axis=False,
        tomographic_map=None,
        rotate=False,
        save_fig=True,
        **kwargs,
    ):
        """Render a map slice image with overlays, colour bar and axes.

        Args:
            pwd (str): Output directory.
            map_slice (numpy.ndarray): The 2D slice to display.
            extentmap (list): ``imshow`` extent.
            xlab, ylab (str): Axis labels.
            name (str): Output figure base name.
            x_index, y_index (int): In-plane axis indexes for overlays.
            pixel_in, pixel_bis_in, qso_in, qso_bis_in, void_in, void_bis_in,
                galaxy_in (numpy.ndarray, optional): Overlay coordinate arrays.
            redshift_axis (bool, optional): Add a secondary redshift axis.
            tomographic_map (optional): Needed for the redshift axis.
            rotate (bool, optional): Rotated layout.
            save_fig (bool, optional): Save the figure to disk.
            **kwargs: Styling options.
        """
        plt.figure()
        fig = plt.gcf()
        ax = plt.gca()
        size = fig.get_size_inches()
        fig.set_size_inches(1.75 * size)

        label_size = utils.return_key(kwargs, "label_size", 15)
        font_size = utils.return_key(kwargs, "font_size", 15)

        plt.xlabel(xlab, fontsize=font_size)
        plt.ylabel(ylab, fontsize=font_size)

        im = TomographyPlot.add_elements(
            map_slice,
            extentmap,
            x_index,
            y_index,
            pixel_in=pixel_in,
            pixel_bis_in=pixel_bis_in,
            qso_in=qso_in,
            qso_bis_in=qso_bis_in,
            void_in=void_in,
            void_bis_in=void_bis_in,
            galaxy_in=galaxy_in,
            **kwargs,
        )

        orientation_color_bar = utils.return_key(
            kwargs, "color_bar_orientation", "horizontal" if rotate else "vertical"
        )
        cbar = plt.colorbar(
            im,
            ax=ax,
            orientation=orientation_color_bar,
            fraction=utils.return_key(kwargs, "color_bar_fraction", 0.1),
        )
        cbar.set_label(
            utils.return_key(
                kwargs,
                "color_bar_label",
                "Reconstructed Ly" + r"$\alpha$" + " contrast " + r"$\delta_{Fmap}$",
            ),
            fontsize=font_size,
        )

        xlim_min = utils.return_key(kwargs, "map_xlim_min", extentmap[0])
        xlim_max = utils.return_key(kwargs, "map_xlim_max", extentmap[1])
        ylim_min = utils.return_key(kwargs, "map_ylim_min", extentmap[2])
        ylim_max = utils.return_key(kwargs, "map_ylim_max", extentmap[3])
        plt.xlim([xlim_min, xlim_max])
        plt.ylim([ylim_min, ylim_max])
        plt.gca().tick_params(axis="x", labelsize=label_size)
        plt.gca().tick_params(axis="y", labelsize=label_size)
        cbar.ax.tick_params(labelsize=label_size)
        if redshift_axis:
            TomographyPlot.add_reshift_axe(tomographic_map, rotate=rotate, **kwargs)

        if save_fig:
            format = utils.return_key(kwargs, "map_format", "pdf")
            plt.savefig(
                os.path.join(pwd, f"{name}.{format}"),
                format=format,
                dpi=utils.return_key(kwargs, "map_dpi", "figure"),
            )
            plt.close()

    @staticmethod
    def add_elements(
        map_slice,
        extentmap,
        x_index,
        y_index,
        pixel_in=None,
        pixel_bis_in=None,
        qso_in=None,
        qso_bis_in=None,
        void_in=None,
        void_bis_in=None,
        galaxy_in=None,
        **kwargs,
    ):
        """Draw the slice image and overlay pixels/QSOs/voids/galaxies.

        Args:
            map_slice (numpy.ndarray): The 2D slice.
            extentmap (list): ``imshow`` extent.
            x_index, y_index (int): In-plane axis indexes for overlays.
            pixel_in, pixel_bis_in, qso_in, qso_bis_in, void_in, void_bis_in,
                galaxy_in (numpy.ndarray, optional): Overlay coordinate arrays
                (``_bis_`` = secondary/faded series).
            **kwargs: Marker/colour styling options.

        Returns:
            matplotlib.image.AxesImage: The displayed slice image.
        """
        im = plt.imshow(
            map_slice,
            interpolation=utils.return_key(kwargs, "map_interpolation", "bilinear"),
            cmap=utils.return_key(kwargs, "map_color", "jet_r"),
            vmin=utils.return_key(kwargs, "map_delta_min", -1.0),
            vmax=utils.return_key(kwargs, "map_delta_max", 0.5),
            extent=extentmap,
        )
        if pixel_in is not None:
            if utils.return_key(kwargs, "pixel_on", True):
                plt.plot(
                    pixel_in[:, x_index],
                    pixel_in[:, y_index],
                    markersize=utils.return_key(kwargs, "pixel_marker_size", 2),
                    marker=utils.return_key(kwargs, "pixel_marker", "."),
                    markeredgewidth=utils.return_key(
                        kwargs, "pixel_marker_edge_size", 1
                    ),
                    color=utils.return_key(kwargs, "pixel_marker_color", "k"),
                    linestyle="None",
                )
        if pixel_bis_in is not None:
            if utils.return_key(kwargs, "pixel_bis_on", True):
                plt.plot(
                    pixel_bis_in[:, x_index],
                    pixel_bis_in[:, y_index],
                    markersize=utils.return_key(
                        kwargs,
                        "pixel_bis_marker_size",
                        utils.return_key(kwargs, "pixel_marker_size", 2),
                    ),
                    marker=utils.return_key(
                        kwargs,
                        "pixel_bis_marker",
                        utils.return_key(kwargs, "pixel_marker", "."),
                    ),
                    markeredgewidth=utils.return_key(
                        kwargs,
                        "pixel_bis_marker_edge_size",
                        utils.return_key(kwargs, "pixel_marker_edge_size", 1),
                    ),
                    color=utils.return_key(kwargs, "pixel_bis_grey", "0.5"),
                    alpha=utils.return_key(kwargs, "pixel_bis_transparency", 0.5),
                    linestyle="None",
                )
        if qso_in is not None:
            if utils.return_key(kwargs, "qso_on", True):
                plt.plot(
                    qso_in[:, x_index],
                    qso_in[:, y_index],
                    marker=utils.return_key(kwargs, "qso_marker", "*"),
                    markersize=utils.return_key(kwargs, "qso_marker_size", 8),
                    markeredgewidth=utils.return_key(kwargs, "qso_marker_edge_size", 1),
                    color=utils.return_key(kwargs, "qso_marker_color", "k"),
                    linestyle="None",
                )
        if qso_bis_in is not None:
            if utils.return_key(kwargs, "qso_bis_on", True):
                plt.plot(
                    qso_bis_in[:, x_index],
                    qso_bis_in[:, y_index],
                    marker=utils.return_key(
                        kwargs,
                        "qso_bis_marker",
                        utils.return_key(kwargs, "qso_marker", "*"),
                    ),
                    markersize=utils.return_key(
                        kwargs,
                        "qso_bis_marker_size",
                        utils.return_key(kwargs, "qso_marker_size", 8),
                    ),
                    markeredgewidth=utils.return_key(
                        kwargs,
                        "qso_bis_marker_edge_size",
                        utils.return_key(kwargs, "qso_marker_edge_size", 1),
                    ),
                    color=utils.return_key(
                        kwargs,
                        "qso_bis_marker_color",
                        utils.return_key(kwargs, "qso_marker_color", "k"),
                    ),
                    linestyle="None",
                    fillstyle="none",
                )
        if void_in is not None:
            if utils.return_key(kwargs, "void_on", True):
                for i in range(len(void_in)):
                    circle = plt.Circle(
                        (void_in[i, x_index], void_in[i, y_index]),
                        void_in[i, 3],
                        fill=False,
                        color=utils.return_key(kwargs, "void_marker_color", "r"),
                    )
                    plt.gcf().gca().add_artist(circle)
        if void_bis_in is not None:
            if utils.return_key(kwargs, "void_bis_on", True):
                for i in range(len(void_bis_in)):
                    circle = plt.Circle(
                        (void_bis_in[i, x_index], void_bis_in[i, y_index]),
                        void_bis_in[i, 3],
                        fill=False,
                        color=utils.return_key(kwargs, "void_bis_marker_color", "k"),
                    )
                    plt.gcf().gca().add_artist(circle)
        if galaxy_in is not None:
            if utils.return_key(kwargs, "galaxy_on", True):
                plt.plot(galaxy_in[:, x_index], galaxy_in[:, y_index], "rx")
                plt.errorbar(
                    galaxy_in[:, x_index],
                    galaxy_in[:, y_index],
                    xerr=2 * galaxy_in[:, 3],
                    capsize=0.01,
                    ecolor="red",
                    fmt="none",
                )
        return im

    @staticmethod
    def add_reshift_axe(tomographic_map, rotate=False, **kwargs):
        """Add a secondary redshift axis to the current slice figure.

        Only applies to the ``"middle"`` coordinate transform; otherwise a no-op.

        Args:
            tomographic_map: The map (provides Omega_m and cartesian bounds).
            rotate (bool, optional): Rotated layout (redshift on the y-axis).
            **kwargs: Axis position/label styling.
        """
        if tomographic_map.coordinate_transform != "middle":
            return ()

        label_size = utils.return_key(kwargs, "label_size", 15)
        font_size = utils.return_key(kwargs, "font_size", 15)

        fig = plt.gcf()
        ax1 = fig.axes[0]
        ax2 = fig.add_axes(ax1.get_position(), frameon=False)

        if rotate:
            bounds = ax1.get_ybound() + tomographic_map.boundary_cartesian_coord[0][2]
        else:
            bounds = ax1.get_xbound() + tomographic_map.boundary_cartesian_coord[0][2]

        z_array = np.linspace(bounds[0], bounds[1], 1000)
        (rcomov, distang, inv_rcomov, inv_distang) = utils.get_cosmo_function(
            tomographic_map.Omega_m
        )
        redshifts = utils.convert_z_cartesian_to_sky_middle(z_array, inv_rcomov)
        redshift_to_plot = np.unique(np.around(redshifts, decimals=1))
        tick_position = np.array(
            [
                np.argmin(np.abs(redshifts - redshift_to_plot[i]))
                for i in range(len(redshift_to_plot))
            ]
        ) / (len(redshifts) - 1)

        ax1.tick_params(labelbottom="on", labelleft="on", bottom="on", left="on")
        if rotate:
            ax2.set_yticks(tick_position)
            ax2.set_yticklabels(redshift_to_plot)
            ax2.tick_params(
                labelright="on",
                right="on",
                labelbottom=None,
                labelleft=None,
                bottom=None,
                left=None,
            )
            ax2.set_ylabel("Redshift $z$", fontsize=font_size)
        else:
            ax2.set_xticks(tick_position)
            ax2.set_xticklabels(redshift_to_plot)
            ax2.tick_params(
                labeltop="on",
                top="on",
                labelbottom=None,
                labelleft=None,
                bottom=None,
                left=None,
            )
            ax2.set_xlabel("Redshift $z$", fontsize=font_size)

        position_redshift_axe = utils.return_key(
            kwargs, "position_redshift_axe", "other"
        )
        outward_redshift_axe = utils.return_key(kwargs, "outward_redshift_axe", 50)
        if position_redshift_axe == "other":
            if rotate:
                ax2.yaxis.set_label_position("right")
                ax2.yaxis.set_ticks_position("right")
            else:
                ax2.xaxis.set_label_position("top")
                ax2.xaxis.set_ticks_position("top")
        elif position_redshift_axe == "same":
            if rotate:
                ax2.yaxis.set_label_position("left")
                ax2.yaxis.set_ticks_position("left")
                ax2.spines["left"].set_position(
                    ("outward", outward_redshift_axe)
                )  # put redshift axis at the bottom
            else:
                ax2.xaxis.set_label_position("bottom")
                ax2.xaxis.set_ticks_position("bottom")
                ax2.spines["bottom"].set_position(("outward", outward_redshift_axe))
        ax2.tick_params(axis="x", labelsize=label_size)
        ax2.tick_params(axis="y", labelsize=label_size)

    def plot_one_slice(
        self,
        name,
        direction,
        space,
        center_mpc,
        qso=None,
        void=None,
        galaxy=None,
        distance_mask=None,
        criteria_distance_mask=None,
        rotate=False,
        minimal_void_crossing=None,
        redshift_axis=False,
        cut_plot=None,
    ):
        """Load objects and plot a single map slice at a given position.

        Args:
            name (str): Base output name.
            direction (str): Slicing axis.
            space (float): Slab thickness (Mpc.h^-1).
            center_mpc (float): Slice centre along the slice-normal axis.
            qso, void, galaxy (str, optional): Overlay catalog files.
            distance_mask (str, optional): Distance-map file for masking.
            criteria_distance_mask (float, optional): Distance-to-LOS threshold.
            rotate (bool, optional): Rotate the slice.
            minimal_void_crossing (float, optional): Void crossing threshold.
            redshift_axis (bool, optional): Draw the secondary redshift axis.
            cut_plot (sequence, optional): Fractional sub-box crop.
        """
        (
            tomographic_map,
            pixel,
            quasar_catalog,
            void_catalog,
            galaxy_catalog,
            dist_map,
        ) = self.load_tomographic_objects(
            qso=qso,
            void=void,
            galaxy=galaxy,
            distance_mask=distance_mask,
            cut_plot=cut_plot,
        )
        mask_void = self.mask_tomographic_objects(
            void_catalog,
            dist_map,
            tomographic_map,
            criteria_distance_mask=criteria_distance_mask,
            minimal_void_crossing=minimal_void_crossing,
        )
        self.print_one_slice(
            name,
            tomographic_map,
            direction,
            space,
            center_mpc,
            pixel=pixel,
            quasar_catalog=quasar_catalog,
            void_catalog=void_catalog,
            mask_void=mask_void,
            galaxy_catalog=galaxy_catalog,
            rotate=rotate,
            redshift_axis=redshift_axis,
        )

    def plot_all_slice(
        self,
        name,
        direction,
        space,
        qso=None,
        void=None,
        galaxy=None,
        distance_mask=None,
        criteria_distance_mask=None,
        rotate=False,
        minimal_void_crossing=None,
        redshift_axis=False,
        cut_plot=None,
    ):
        """Plot successive slices spanning the whole box along a direction.

        Args:
            name (str): Base output name.
            direction (str): Slicing axis.
            space (float): Slab thickness and step between slices (Mpc.h^-1).
            qso, void, galaxy (str, optional): Overlay catalog files.
            distance_mask (str, optional): Distance-map file for masking.
            criteria_distance_mask (float, optional): Distance-to-LOS threshold.
            rotate (bool, optional): Rotate the slices.
            minimal_void_crossing (float, optional): Void crossing threshold.
            redshift_axis (bool, optional): Draw the secondary redshift axis.
            cut_plot (sequence, optional): Fractional sub-box crop.
        """
        (
            tomographic_map,
            pixel,
            quasar_catalog,
            void_catalog,
            galaxy_catalog,
            dist_map,
        ) = self.load_tomographic_objects(
            qso=qso,
            void=void,
            galaxy=galaxy,
            distance_mask=distance_mask,
            cut_plot=cut_plot,
        )
        mask_void = self.mask_tomographic_objects(
            void_catalog,
            dist_map,
            tomographic_map,
            criteria_distance_mask=criteria_distance_mask,
            minimal_void_crossing=minimal_void_crossing,
        )
        center_mpc = space / 2
        index_dict = {
            "x": 0,
            "y": 1,
            "z": 2,
            "ra": 0,
            "dec": 1,
            "redshift": 2,
            "x_lab": "X",
            "y_lab": "Y",
            "z_lab": "Z",
        }
        while center_mpc + space / 2 < tomographic_map.size[index_dict[direction]]:
            self.print_one_slice(
                name,
                tomographic_map,
                direction,
                space,
                center_mpc,
                pixel=pixel,
                quasar_catalog=quasar_catalog,
                void_catalog=void_catalog,
                mask_void=mask_void,
                galaxy_catalog=galaxy_catalog,
                rotate=rotate,
                redshift_axis=redshift_axis,
            )
            center_mpc += space

    def plot(
        self,
        name,
        direction,
        space,
        center_mpc,
        qso=None,
        void=None,
        galaxy=None,
        distance_mask=None,
        criteria_distance_mask=None,
        rotate=False,
        minimal_void_crossing=None,
        redshift_axis=False,
        cut_plot=None,
    ):
        """Plot one slice, or every slice, depending on ``center_mpc``.

        Args:
            name (str): Base output name.
            direction (str): Slicing axis.
            space (float): Slab thickness (Mpc.h^-1).
            center_mpc (float | str): A numeric slice position, or ``"all"`` to
                sweep the whole box.
            qso, void, galaxy (str, optional): Overlay catalog files.
            distance_mask (str, optional): Distance-map file for masking.
            criteria_distance_mask (float, optional): Distance-to-LOS threshold.
            rotate (bool, optional): Rotate the slices.
            minimal_void_crossing (float, optional): Void crossing threshold.
            redshift_axis (bool, optional): Draw the secondary redshift axis.
            cut_plot (sequence, optional): Fractional sub-box crop.

        Raises:
            ValueError: If ``center_mpc`` is neither numeric nor ``"all"``.
        """
        if (type(center_mpc) == float) | (type(center_mpc) == int):
            self.plot_one_slice(
                name,
                direction,
                space,
                center_mpc,
                qso=qso,
                void=void,
                galaxy=galaxy,
                distance_mask=distance_mask,
                criteria_distance_mask=criteria_distance_mask,
                rotate=rotate,
                minimal_void_crossing=minimal_void_crossing,
                redshift_axis=redshift_axis,
                cut_plot=cut_plot,
            )
        elif center_mpc.lower() == "all":
            self.plot_all_slice(
                name,
                direction,
                space,
                qso=qso,
                void=void,
                galaxy=galaxy,
                distance_mask=distance_mask,
                criteria_distance_mask=criteria_distance_mask,
                rotate=rotate,
                minimal_void_crossing=minimal_void_crossing,
                redshift_axis=redshift_axis,
                cut_plot=cut_plot,
            )
        else:
            raise ValueError(
                "Please give the distance of the slice you want to print or all"
            )

    def plot_integrate_image(self, zmin, zmax, name, void=None, cut_plot=None):
        """Plot the map averaged over a radial (z) slab, with void overlays.

        Args:
            zmin (float): Lower radial bound of the slab (Mpc.h^-1).
            zmax (float): Upper radial bound of the slab (Mpc.h^-1).
            name (str): Base output name.
            void (str, optional): Void catalog to overlay.
            cut_plot (sequence, optional): Fractional sub-box crop.
        """
        (
            tomographic_map,
            _,
            _,
            void_catalog,
            _,
            _,
        ) = self.load_tomographic_objects(
            void=void,
            cut_plot=cut_plot,
        )
        (
            x_index,
            y_index,
            index_direction,
            extentmap,
            xlab,
            ylab,
        ) = TomographyPlot.get_direction_informations("z", False, tomographic_map.size)
        map_data = tomographic_map.map_array
        i_pix_begin = int(
            round((zmin * tomographic_map.pixel_per_mpc[index_direction]), 0)
        )
        i_pix_end = int(
            round((zmax * tomographic_map.pixel_per_mpc[index_direction]), 0)
        )
        integrated_map = np.mean(map_data[:, :, i_pix_begin:i_pix_end], axis=2)

        void_in = None
        if void is not None:
            mask_void_in = void_catalog.coord[:, 2] < zmax
            mask_void_in &= void_catalog.coord[:, 2] >= zmin
            void_coord_in = void_catalog.coord[mask_void_in]
            void_radius_in = void_catalog.radius[mask_void_in]
            void_in = np.transpose(
                np.vstack(
                    [
                        void_coord_in[:, 0],
                        void_coord_in[:, 1],
                        void_coord_in[:, 2],
                        void_radius_in,
                    ]
                )
            )

        name = f"{name}_integrated_map"
        TomographyPlot.plot_slice(
            self.pwd,
            integrated_map,
            extentmap,
            xlab,
            ylab,
            name,
            x_index,
            y_index,
            rotate=False,
            void_in=void_in,
            **self.kwargs,
        )

    def plot_2d_shell(
        self,
        ra_array,
        dec_array,
        redshift,
        interpolation_method,
        void=None,
        rebin_shell=None,
        **kwargs,
    ):
        """Interpolate and plot a 2D (RA, Dec) shell of the map at fixed z.

        Builds an (RA, Dec, z) grid, converts it to map pixel coordinates,
        interpolates the map onto it, averages over the redshift thickness and
        displays the resulting shell (optionally rebinned), overlaying voids.

        Args:
            ra_array (array-like): RA grid (degrees).
            dec_array (array-like): Dec grid (degrees).
            redshift (array-like): Redshift slab sampling.
            interpolation_method (str): Map interpolation mode.
            void (str, optional): Void catalog to overlay.
            rebin_shell (tuple[int], optional): Rebin the shell to this shape.
            **kwargs: Colour-scale styling.

        Returns:
            tuple: ``(shell_mean, ra_dec_grid)``.
        """
        (
            tomographic_map,
            _,
            _,
            void_catalog,
            _,
            _,
        ) = self.load_tomographic_objects(
            void=void,
        )
        tomographic_map = self.load_tomographic_objects()[0]
        coords_ra_dec = np.moveaxis(
            np.array(np.meshgrid(ra_array, dec_array, redshift, indexing="ij")), 0, -1
        )
        coords_cartesian = np.zeros(coords_ra_dec.shape)
        suplementary_parameters = utils.return_suplementary_parameters(
            tomographic_map.coordinate_transform, property=tomographic_map
        )
        (rcomov, distang, _, _) = utils.get_cosmo_function(tomographic_map.Omega_m)

        (
            coords_cartesian[:, :, :, 0],
            coords_cartesian[:, :, :, 1],
            coords_cartesian[:, :, :, 2],
        ) = utils.convert_sky_to_cartesian(
            np.radians(coords_ra_dec[:, :, :, 0]),
            np.radians(coords_ra_dec[:, :, :, 1]),
            coords_ra_dec[:, :, :, 2],
            tomographic_map.coordinate_transform,
            rcomov=rcomov,
            distang=distang,
            suplementary_parameters=suplementary_parameters,
        )

        (
            coords_cartesian[:, :, :, 0],
            coords_cartesian[:, :, :, 1],
            coords_cartesian[:, :, :, 2],
        ) = (
            coords_cartesian[:, :, :, 0]
            - tomographic_map.boundary_cartesian_coord[0][0],
            coords_cartesian[:, :, :, 1]
            - tomographic_map.boundary_cartesian_coord[0][1],
            coords_cartesian[:, :, :, 2]
            - tomographic_map.boundary_cartesian_coord[0][2],
        )

        coords_pixel = np.zeros(coords_cartesian.shape)

        (
            coords_pixel[:, :, :, 0],
            coords_pixel[:, :, :, 1],
            coords_pixel[:, :, :, 2],
        ) = (
            coords_cartesian[:, :, :, 0] / tomographic_map.mpc_per_pixel[0],
            coords_cartesian[:, :, :, 1] / tomographic_map.mpc_per_pixel[1],
            coords_cartesian[:, :, :, 2] / tomographic_map.mpc_per_pixel[2],
        )

        shell = utils.interpolate_map(
            interpolation_method, tomographic_map.map_array, coords_pixel
        )
        shell_mean = np.mean(shell, axis=2)

        extent = [
            np.min(ra_array),
            np.max(ra_array),
            np.min(dec_array),
            np.max(dec_array),
        ]
        plt.figure(figsize=(20, 7))

        if rebin_shell:
            shell_mean = utils.bin_ndarray(shell_mean, rebin_shell, operation="mean")

        im = plt.imshow(
            np.flip(np.transpose(shell_mean), axis=0),
            extent=extent,
            cmap=utils.return_key(kwargs, "map_color", "jet_r"),
            vmin=utils.return_key(kwargs, "map_delta_min", -0.2),
            vmax=utils.return_key(kwargs, "map_delta_max", 0.2),
        )
        plt.colorbar(im, fraction=0.01, pad=0.04)
        plt.xlabel("RA")
        plt.ylabel("DEC")

        if void is not None:
            void_catalog.cut_catalog_void(
                [],
                coord_min=(np.min(ra_array), np.min(dec_array), np.min(redshift)),
                coord_max=(np.max(ra_array), np.max(dec_array), np.max(redshift)),
            )
            for i in range(len(void_catalog.coord)):
                plt.scatter(
                    void_catalog.coord[i, 0],
                    void_catalog.coord[i, 1],
                    color="r",
                    marker="x",
                )

        return shell_mean, coords_ra_dec[:, :, 0, :2]

    def plot_catalog_centered_maps(
        self,
        direction,
        name,
        space,
        void,
        nb_plot,
        radius_centered,
        qso=None,
        rotate=False,
    ):
        """Plot slices centred on the largest voids of a catalog.

        Args:
            direction (str): Slicing axis.
            name (str): Base output name.
            space (float): Slab thickness (Mpc.h^-1).
            void (str): Void catalog file.
            nb_plot (int): Number of (largest) voids to centre on.
            radius_centered (float): Half-size of each cut-out (Mpc.h^-1).
            qso (str, optional): QSO catalog to overlay.
            rotate (bool, optional): Rotate the slices.
        """
        load = self.load_tomographic_objects(void=void, qso=qso)
        tomographic_map, pixel, void_catalog = load[0], load[1], load[3]
        if qso is not None:
            qso = load[2]
        arg = np.array(void_catalog.radius).argsort()[-nb_plot:][::-1]
        coords = np.array(void_catalog.coord)[arg]
        (
            x_index,
            y_index,
            index_direction,
            extentmap,
            xlab,
            ylab,
        ) = TomographyPlot.get_direction_informations(
            direction, rotate, tomographic_map.size
        )
        for i in range(len(coords)):
            xlim = [
                coords[i][x_index] - radius_centered,
                coords[i][x_index] + radius_centered,
            ]
            ylim = [
                coords[i][y_index] - radius_centered,
                coords[i][y_index] + radius_centered,
            ]
            self.kwargs.update(
                {
                    "map_xlim_min": xlim[0],
                    "map_xlim_max": xlim[1],
                    "map_ylim_min": ylim[0],
                    "map_ylim_max": ylim[1],
                }
            )
            center_mpc = coords[i][index_direction]
            name_plot = "{}_number{}_radius{}".format(
                name, i, np.array(void_catalog.radius)[arg][i]
            )
            self.print_one_slice(
                name_plot,
                tomographic_map,
                direction,
                space,
                center_mpc,
                pixel=pixel,
                quasar_catalog=qso,
                void_catalog=void_catalog,
                rotate=rotate,
            )

    def plot_delta_histogram(
        self,
        name,
        nb_bins,
        gauss_fit=True,
        norm=True,
        distance_mask=None,
        criteria_distance_mask=None,
        log_scale=True,
        cut_plot=None,
    ):
        """Plot the histogram of the map delta values (optionally masked).

        Args:
            name (str): Output figure name.
            nb_bins (int): Number of histogram bins.
            gauss_fit (bool, optional): Overlay a Gaussian fit.
            norm (bool, optional): Normalise to a density.
            distance_mask (str, optional): Distance-map file for masking.
            criteria_distance_mask (float, optional): Distance-to-LOS threshold.
            log_scale (bool, optional): Log-scale the y-axis.
            cut_plot (sequence, optional): Fractional sub-box crop.
        """
        (
            tomographic_map,
            pixel,
            quasar_catalog,
            void_catalog,
            galaxy_catalog,
            dist_map,
        ) = self.load_tomographic_objects(
            distance_mask=distance_mask, cut_plot=cut_plot
        )
        self.mask_tomographic_objects(
            void_catalog,
            dist_map,
            tomographic_map,
            criteria_distance_mask=criteria_distance_mask,
        )
        listdeltas = tomographic_map.map_array.ravel()
        plt.figure()
        self.delta_histogram(listdeltas, nb_bins, norm=norm, gauss_fit=gauss_fit)
        if log_scale:
            plt.yscale("log")
        plt.grid()
        plt.savefig(os.path.join(self.pwd, "{}.pdf".format(name)), format="pdf")

    def plot_delta_histogram_comparison(
        self,
        name,
        name_second_map,
        nb_bins,
        legend,
        gauss_fit=True,
        norm=True,
        distance_mask=None,
        distance_mask2=None,
        criteria_distance_mask=None,
        log_scale=True,
        cut_plot=None,
    ):
        """Overlay the delta histograms of two maps for comparison.

        Args:
            name (str): Output figure name.
            name_second_map (str): Second map file.
            nb_bins (int): Number of histogram bins.
            legend (list[str]): Legend labels for the two maps.
            gauss_fit (bool, optional): Overlay Gaussian fits.
            norm (bool, optional): Normalise to a density.
            distance_mask (str, optional): Distance-map for the first map.
            distance_mask2 (str, optional): Distance-map for the second map.
            criteria_distance_mask (float, optional): Distance-to-LOS threshold.
            log_scale (bool, optional): Log-scale the y-axis.
            cut_plot (sequence, optional): Fractional sub-box crop.
        """
        (
            tomographic_map,
            pixel,
            quasar_catalog,
            void_catalog,
            galaxy_catalog,
            dist_map,
        ) = self.load_tomographic_objects(
            distance_mask=distance_mask, cut_plot=cut_plot
        )
        self.mask_tomographic_objects(
            void_catalog,
            dist_map,
            tomographic_map,
            criteria_distance_mask=criteria_distance_mask,
        )
        listdeltas = tomographic_map.map_array.ravel()

        tomo_plot = TomographyPlot(
            self.pwd,
            map_name=name_second_map,
            map_shape=self.map_shape,
            pixel_name=self.pixel_name,
            property_file=self.property_file,
        )
        (
            tomographic_map2,
            pixel,
            quasar_catalog,
            void_catalog,
            galaxy_catalog,
            dist_map2,
        ) = tomo_plot.load_tomographic_objects(
            distance_mask=distance_mask2, cut_plot=cut_plot
        )
        tomo_plot.mask_tomographic_objects(
            void_catalog,
            dist_map2,
            tomographic_map2,
            criteria_distance_mask=criteria_distance_mask,
        )
        listdeltas2 = tomographic_map2.map_array.ravel()

        plt.figure()
        self.delta_histogram(
            listdeltas, nb_bins, norm=norm, gauss_fit=gauss_fit, alpha=0.5
        )
        self.delta_histogram(
            listdeltas2, nb_bins, norm=norm, gauss_fit=gauss_fit, alpha=0.5
        )
        if log_scale:
            plt.yscale("log")
        plt.legend(legend, fontsize=15)
        plt.xlabel(r"$\delta_{F\mathrm{map}}$", fontsize=15)
        plt.ylabel("#", fontsize=15)
        plt.gca().tick_params(axis="x", labelsize=15)
        plt.gca().tick_params(axis="y", labelsize=15)
        plt.grid()
        plt.savefig(os.path.join(self.pwd, "{}.pdf".format(name)), format="pdf")

    def compare_two_map(
        self,
        name,
        name_second_map,
        distance_mask,
        distance_second_mask,
        dist_extremum,
        bin_dist,
        legend,
        shuffle_map=None,
        cut_plot=None,
    ):
        """Plot the correlation between two maps versus distance-to-LOS.

        For a range of distance-to-LOS thresholds, computes the Pearson
        correlation of the two maps over voxels closer than the threshold in
        both, optionally including a shuffled map as a null reference.

        Args:
            name (str): Output figure name.
            name_second_map (str): Second map file.
            distance_mask (str): Distance-map for the first map.
            distance_second_mask (str): Distance-map for the second map.
            dist_extremum (sequence): ``(min, max)`` distance range.
            bin_dist (int): Number of distance thresholds.
            legend (list[str]): Legend labels.
            shuffle_map (str, optional): Shuffled map for the null reference.
            cut_plot (sequence, optional): Fractional sub-box crop.

        Returns:
            tuple: ``(corr, dist_range)``.
        """
        (
            tomographic_map,
            pixel,
            quasar_catalog,
            void_catalog,
            galaxy_catalog,
            dist_map,
        ) = self.load_tomographic_objects(
            distance_mask=distance_mask, cut_plot=cut_plot
        )

        tomo_plot = TomographyPlot(
            self.pwd,
            map_name=name_second_map,
            map_shape=self.map_shape,
            pixel_name=self.pixel_name,
            property_file=self.property_file,
        )
        (
            tomographic_map2,
            pixel,
            quasar_catalog,
            void_catalog,
            galaxy_catalog,
            dist_map2,
        ) = tomo_plot.load_tomographic_objects(
            distance_mask=distance_second_mask, cut_plot=cut_plot
        )

        corr = []
        if shuffle_map is not None:
            tomo_plot_shuffle = TomographyPlot(
                self.pwd,
                map_name=shuffle_map,
                map_shape=self.map_shape,
                pixel_name=self.pixel_name,
                property_file=self.property_file,
            )
            tomographic_shuffle = tomo_plot_shuffle.load_tomographic_objects(
                cut_plot=cut_plot
            )[0]
            corr_shuffle = []

        dist_range = np.linspace(dist_extremum[0], dist_extremum[1], bin_dist)
        for i in range(len(dist_range)):
            mask = (dist_map < dist_range[i]) & (dist_map2 < dist_range[i])
            corr.append(
                np.corrcoef(
                    tomographic_map.map_array[mask], tomographic_map2.map_array[mask]
                )[0][1]
            )
            if shuffle_map is not None:
                corr_shuffle.append(
                    np.corrcoef(
                        tomographic_map.map_array[mask],
                        tomographic_shuffle.map_array[mask],
                    )[0][1]
                )
        plt.plot(dist_range, corr)
        plt.plot(dist_range, corr_shuffle)
        plt.legend(legend)
        plt.xlabel("Distance to the nearest los [" + r"$\mathrm{h^{-1}Mpc}$" + "]")
        plt.ylabel("Correlation coefficient")
        plt.grid()
        plt.savefig(os.path.join(self.pwd, "{}.pdf".format(name)), format="pdf")
        return (corr, dist_range)

    def plot_pk3D(
        self,
        name_map,
        name_prop,
        n_k,
        kmin,
        kmax,
        log=False,
        distance_map=None,
        criteria_distance_mask=None,
    ):
        """Compute and plot the 3D power spectrum of a map.

        Args:
            name_map (str): Map file.
            name_prop (str): Map property file.
            n_k (int): Number of k bins.
            kmin (float): Minimum wavenumber (h/Mpc).
            kmax (float): Maximum wavenumber (h/Mpc).
            log (bool, optional): Use log-spaced k / semilog plot.
            distance_map (str, optional): Distance-map for masking.
            criteria_distance_mask (float, optional): Distance-to-LOS threshold.

        Returns:
            tuple: ``(pk_3D, k_space)``.
        """
        tomographic_map = tomographic_objects.TomographicMap.init_from_property_files(
            name_prop, name=name_map
        )
        (k_space_final, pk_3D_final) = tomographic_map.compute_pk3d(
            kmin,
            kmax,
            n_k,
            distance_map=distance_map,
            criteria_distance_mask=criteria_distance_mask,
            log=log,
        )
        if log:
            plt.semilogx(k_space_final, pk_3D_final)
        else:
            plt.plot(k_space_final, pk_3D_final)
        plt.grid()
        plt.xlabel("Comoving wavevector in h.Mpc-1")
        plt.ylabel("Pk 3D")
        plt.savefig(os.path.join(self.pwd, "Pk_3D.pdf"), format="pdf")
        return (pk_3D_final, k_space_final)


class TomographyStack(object):
    """Stack a tomographic map at catalog positions and plot the stack.

    Stacks fixed-size cut-outs of the map around each catalog object
    (void/QSO/galaxy), optionally radius-normalised, and provides plotting and
    ellipticity-analysis helpers on the resulting mean stack.
    """

    def __init__(
        self,
        map_name,
        catalog_name,
        type_catalog,
        property_file_stack,
        size_stack,
        name_stack,
        shape_stack=None,
        map_shape=None,
        map_size=None,
        property_file=None,
        coordinate_convert=None,
        interpolation_method="NEAREST",
        normalized=False,
    ):
        """Load the map and catalog and set the stack geometry.

        Args:
            map_name (str): Map file.
            catalog_name (str): Catalog file to stack on.
            type_catalog (str): Catalog type (``"void"``/``"qso"``/``"galaxy"``).
            property_file_stack (str): Output stack property file.
            size_stack (float): Half-size of the stack cube (Mpc.h^-1).
            name_stack (str): Output stack binary file.
            shape_stack (tuple[int], optional): Explicit stack pixel shape.
            map_shape (tuple[int], optional): Map pixel shape (if no property file).
            map_size (tuple[float], optional): Map physical size.
            property_file (str, optional): Map property/pickle file.
            coordinate_convert (str, optional): Coordinate conversion before
                stacking.
            interpolation_method (str, optional): ``"NEAREST"`` or ``"LINEAR"``.
            normalized (bool, optional): Radius-normalise each cut-out.
        """
        self.size_stack = size_stack
        self.name_stack = name_stack
        self.property_file_stack = property_file_stack
        self.coordinate_convert = coordinate_convert
        self.interpolation_method = interpolation_method
        self.normalized = normalized

        self.tomographic_map = tomographic_objects.TomographicMap.init_classic(
            name=map_name, shape=map_shape, size=map_size, property_file=property_file
        )
        self.tomographic_map.read()
        self.catalog = tomographic_objects.Catalog.init_catalog_from_fits(
            catalog_name, type_catalog
        )

        if shape_stack is None:
            self.shape_stack = tuple(
                np.around(
                    utils.get_map_shape(size_stack, self.tomographic_map.mpc_per_pixel),
                    decimals=0,
                ).astype(int)
            )
        else:
            self.shape_stack = shape_stack

    def stack(self):
        """Compute the stack over all catalog positions and write it to disk."""
        stack = tomographic_objects.StackMap.init_by_tomographic_map(
            self.tomographic_map,
            self.catalog,
            self.size_stack,
            self.shape_stack,
            self.property_file_stack,
            interpolation_method=self.interpolation_method,
            name=self.name_stack,
            normalized=self.normalized,
            coordinate_convert=self.coordinate_convert,
        )
        stack.write()
        self.stack = stack

    def merge_stack(
        self,
        stack_name,
        property_stack_name,
        name,
        property_name,
        ellipticity_calculation=False,
    ):
        """Merge several stacks (optionally with jackknife ellipticity errors).

        Args:
            stack_name (list[str]): Stack binary files to merge.
            property_stack_name (list[str]): Their property files.
            name (str): Output merged stack file.
            property_name (str): Output merged property file.
            ellipticity_calculation (bool, optional): Compute the stack
                ellipticity and its jackknife errors.
        """
        merge_stack = tomographic_objects.StackMap.init_by_merging(
            stack_name, property_stack_name, name, property_name
        )
        if ellipticity_calculation:
            merge_stack.compute_stack_ellipticity()
            ellipticities = []
            for i in range(len(stack_name)):
                jack_knife_list = stack_name.copy()
                jack_knife_property_list = property_stack_name.copy()
                jack_knife_list.remove(stack_name[i])
                jack_knife_property_list.remove(property_stack_name[i])
                stack = tomographic_objects.StackMap.StackMap.init_by_merging(
                    jack_knife_list, jack_knife_property_list, None, None
                )
                stack.compute_stack_ellipticity()
                ellipticities.append(stack.ellipticity)
            merge_stack.add_ellipticity_errors(ellipticities)
        stack.write()
        self.stack = stack

    @staticmethod
    def plot_stack(
        pwd,
        stack_name,
        stack_property,
        name_plot,
        rotate=False,
        ellipticity=False,
        pixel_file_qso_distance=None,
        **kwargs,
    ):
        """Plot the central slices of a stack along x, y and z.

        Args:
            pwd (str): Output directory.
            stack_name (str): Stack binary file.
            stack_property (str): Stack property file.
            name_plot (str): Output figure base name.
            rotate (bool, optional): Rotate the slices.
            ellipticity (bool, optional): Overlay the ellipticity fit.
            pixel_file_qso_distance (str, optional): Pixel file for the mean
                QSO-distance overlay.
            **kwargs: Styling options.
        """
        stack = tomographic_objects.StackMap.init_classic(
            name=stack_name, property_file=stack_property
        )
        stack.read()
        for direction in ["x", "y", "z"]:
            TomographyStack.plot_stack_direction(
                pwd,
                stack,
                name_plot,
                direction,
                rotate=rotate,
                ellipticity=ellipticity,
                pixel_file_qso_distance=pixel_file_qso_distance,
                **kwargs,
            )

    @staticmethod
    def plot_stack_direction(
        pwd,
        stack,
        name_plot,
        direction,
        rotate=False,
        ellipticity=False,
        pixel_file_qso_distance=None,
        **kwargs,
    ):
        """Plot the central stack slice along one direction.

        Args:
            pwd (str): Output directory.
            stack: The loaded stack object.
            name_plot (str): Output figure base name.
            direction (str): Slicing axis (``x``/``y``/``z``).
            rotate (bool, optional): Rotate the slice.
            ellipticity (bool, optional): Overlay the ellipticity fit.
            pixel_file_qso_distance (str, optional): Pixel file for the mean
                QSO-distance overlay.
            **kwargs: Styling options.
        """
        (
            x_index,
            y_index,
            index_direction,
            extentmap,
            xlab,
            ylab,
        ) = TomographyPlot.get_direction_informations(direction, rotate, stack.size)
        if direction == "x":
            stack_slice = stack.map_array[stack.shape[0] // 2, :, :]
            if rotate:
                stack_slice = np.transpose(np.flip(stack_slice, axis=1))
        elif direction == "y":
            stack_slice = stack.map_array[:, stack.shape[1] // 2, :]
            if rotate:
                stack_slice = np.transpose(np.flip(stack_slice, axis=1))
        elif direction == "z":
            stack_slice = stack.map_array[:, :, stack.shape[2] // 2]
            if rotate:
                stack_slice = np.transpose(np.flip(stack_slice, axis=1))

        extentmap = [
            -stack.size[0] / 2,
            +stack.size[0] / 2,
            -stack.size[0] / 2,
            +stack.size[0] / 2,
        ]
        TomographyPlot.plot_slice(
            pwd,
            np.transpose(stack_slice),
            extentmap,
            xlab,
            ylab,
            name_plot,
            x_index,
            y_index,
            save_fig=False,
            **kwargs,
        )

        if pixel_file_qso_distance is not None:
            TomographyStack.plot_mean_los_distance(
                direction, stack, pixel_file_qso_distance=pixel_file_qso_distance
            )

        if ellipticity:
            if stack.ellipticity is None:
                stack.compute_stack_ellipticity()
            TomographyStack.plot_ellipticity(
                stack_slice,
                stack.ellipticity,
                stack.mpc_per_pixel[y_index],
                stack.mpc_per_pixel[x_index],
                direction,
                rotate,
                **kwargs,
            )
        plt.plot([0], [0], "kx")
        plt.savefig(
            f"{name_plot}_{direction}.pdf",
            format="pdf",
            dpi=utils.return_key(kwargs, "map_dpi", "figure"),
        )
        plt.close()

    @staticmethod
    def plot_mean_los_distance(direction, stack, pixel_file_qso_distance=None):
        """Draw the mean QSO line-of-sight distance as a marker on the stack.

        Args:
            direction (str): Slicing axis.
            stack: The loaded stack object.
            pixel_file_qso_distance (str, optional): Pixel file to compute the
                mean LOS distance from.
        """
        if stack.mean_los_distance is None:
            stack.compute_distance_to_los(pixel_file_qso_distance)
        if stack.mean_los_distance < stack.size:
            if (direction == "x") | (direction == "y"):
                plt.plot([-stack.size, -stack.mean_los_distance], [0, 0], "r-")

    @staticmethod
    def plot_ellipticity(
        stack_slice, ellipticity, mpx, mpy, direction, rotate, **kwargs
    ):
        """Overlay Gaussian-fit ellipticity contours on a stack slice.

        Args:
            stack_slice (numpy.ndarray): The 2D stack slice.
            ellipticity (dict): Per-direction ellipticity fit parameters.
            mpx (float): Mpc.h^-1 per pixel along the x display axis.
            mpy (float): Mpc.h^-1 per pixel along the y display axis.
            direction (str): Slicing axis.
            rotate (bool): Must be False (rotation is unsupported here).
            **kwargs: Contour styling.

        Raises:
            NotImplementedError: If ``rotate`` is True.
        """
        if rotate:
            raise NotImplementedError(
                "showing ellipticity with a rotate figure is not implemented, please put rotate to False"
            )
        x, y = np.indices((stack_slice.shape[0], stack_slice.shape[0]), dtype=np.float)
        xcenter = (x - x.shape[0] // 2) * mpx
        ycenter = -(y - y.shape[1] // 2) * mpy
        gauss = utils.gaussian_fitter_2d()
        gaussian = gauss.Gaussian2D(*ellipticity[direction + "_gauss"])
        data_fitted = gaussian(xcenter, ycenter)
        levels = utils.return_key(
            kwargs,
            "levels",
            [
                1 - np.exp(-((1) ** 2) / 2),
                1 - np.exp(-((2) ** 2) / 2),
                1 - np.exp(-((3) ** 2) / 2),
            ],
        )
        plt.contour(
            np.transpose(xcenter),
            np.transpose(ycenter),
            np.transpose(
                data_fitted.reshape(stack_slice.shape[0], stack_slice.shape[0])
            ),
            levels,
            linewidths=utils.return_key(kwargs, "linewidths", 2),
            colors=utils.return_key(kwargs, "colors", "w"),
        )
        ticks = utils.return_key(kwargs, "ticks", False)
        if ticks:
            elname = f"""{ellipticity[direction + "_order"]} = {str(np.round(ellipticity[direction],2))}"""
            plt.text(
                0.1,
                0.9,
                elname,
                ha="center",
                va="center",
                transform=plt.gca().transAxes,
            )
