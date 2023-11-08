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


import os, pickle, glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from random import sample
from matplotlib.lines import Line2D
from lelantos import utils
from lelantos import tomographic_objects


#############################################################################
#############################################################################
############################### FUNCTIONS ###################################
#############################################################################
#############################################################################


def get_delta_list(delta_path):
    if type(delta_path) == str:
        delta_list = np.sort(glob.glob(os.path.join(delta_path, "delta-*.fits*")))
    elif type(delta_path) == list:
        delta_list = []
        for i in range(len(delta_path)):
            delta_list = delta_list + list(
                np.sort(glob.glob(os.path.join(delta_path[i], "delta-*.fits*")))
            )
    if len(delta_list) == 0:
        raise KeyError("No delta file was found")
    return delta_list


# CR - need to rethink delta class 
def preselect_deltas(
    namefile,
    ramin=None,
    ramax=None,
    decmin=None,
    decmax=None,
    center_ra=True,
    pk1d_type=True,
):
    subset_namefile = []
    if (ramin is None) & (ramax is None) & (decmin is None) & (decmax is None):
        return namefile
    for i in range(len(namefile)):
        delta_tomo = tomographic_objects.Delta(name=namefile[i], pk1d_type=pk1d_type)
        delta_tomo.read()
        ra, dec, redshift, redshift_qso, id, sigma, delta = delta_tomo.return_params(
            center_ra=center_ra
        )
        mask = np.full(ra.shape, True)
        if ramin is not None:
            mask &= ra > ramin
        if ramax is not None:
            mask &= ra < ramax
        if decmin is not None:
            mask &= dec > decmin
        if ramin is not None:
            mask &= dec < decmax
        if len(mask[mask]) != 0:
            subset_namefile.append(namefile[i])
    print(len(subset_namefile))
    if len(subset_namefile) == 0:
        raise KeyError("Select window does not contain any delta files")
    return subset_namefile


def get_deltas(namefile, center_ra=True, pk1d_type=True):
    """Extract delta properties"""
    ras, decs, redshifts, redshift_qsos, ids, sigmas, deltas = (
        [],
        [],
        [],
        [],
        [],
        [],
        [],
    )
    for i in range(len(namefile)):
        delta_tomo = tomographic_objects.Delta(name=namefile[i], pk1d_type=pk1d_type)
        delta_tomo.read()
        ra, dec, redshift, redshift_qso, id, sigma, delta = delta_tomo.return_params(
            center_ra=center_ra
        )
        ras.append(ra)
        decs.append(dec)
        redshift_qsos.append(redshift_qso)
        ids.append(id)
        redshifts = redshifts + redshift
        sigmas = sigmas + sigma
        deltas = deltas + delta
    return (
        np.concatenate(ras),
        np.concatenate(decs),
        redshifts,
        np.concatenate(redshift_qsos),
        np.concatenate(ids),
        sigmas,
        deltas,
    )


def get_merged_multiple_exposure_deltas(namefile):
    """Merge deltas with repeated observation"""
    # Pack LOS by Id in the dict Deltas
    ra, dec, z, deltas, sigmas, zqso = [], [], [], [], [], []
    (Deltas, ids) = get_id_list(namefile)

    # For each pack of LOS
    for i in range(len(ids)):
        # Get the data
        zqso.append(tomographic_objects.Delta.z_qso(Deltas[ids[i]][0]))
        ra.append(tomographic_objects.Delta.ra(Deltas[ids[i]][0]))
        dec.append(tomographic_objects.Delta.dec(Deltas[ids[i]][0]))
        listsigmas, listz, listdelta = [], [], []
        for j in range(len(Deltas[ids[i]])):
            listsigmas.append(1 / np.sqrt(np.asarray(Deltas[ids[i]][j].ivar)))
            listdelta.append(Deltas[ids[i]][j].delta)
            listz.append(
                ((10 ** np.asarray(Deltas[ids[i]][j].log_lambda) / utils.lambdaLy) - 1)
            )

        (listz, listsigmas, listdelta) = delete_los_extrema(
            listz, listsigmas, listdelta
        )

        (listz, listsigmas, listdelta, lenlists) = delete_missing_pixels(
            listz, listsigmas, listdelta
        )

        # Weighted merging of the LOSs
        zmerged, sigmamerged, deltamerged = listz[0], [], []
        for k in range(len(listz[0])):
            sigma = 0
            delta = 0
            sumsigma = 0
            for m in range(len(listz)):
                sumsigma = sumsigma + 1 / (listsigmas[m][k] ** 2)
            for m in range(len(listz)):
                delta = delta + listdelta[m][k] / ((listsigmas[m][k] ** 2) * sumsigma)

                sigma = sigma + 1 / ((listsigmas[m][k] ** 2) * (sumsigma**2))
            sigma = np.sqrt(sigma / len(listz))
            sigmamerged.append(sigma)
            deltamerged.append(delta)
        deltas.append(deltamerged)
        z.append(zmerged)
        sigmas.append(sigmamerged)
    return (
        np.array(ra),
        np.array(dec),
        np.asarray(z),
        np.array(zqso),
        np.array(ids),
        np.asarray(sigmas),
        np.asarray(deltas),
    )


def get_id_list(namefile):
    ids = []
    for i in range(len(namefile)):
        delta_tomo = tomographic_objects.Delta(name=namefile[i], pk1d_type=True)
        delta_tomo.read()
        id = []
        for i in range(len(delta_tomo.delta_array)):
            id.append(tomographic_objects.Delta.primary_key(delta_tomo.delta_array[i]))
        ids.append(id)
    ids = np.concatenate(ids)
    ids = list(set(ids))
    Deltas = {}
    for i in range(len(ids)):
        Deltas[ids[i]] = []
        for j in range(len(namefile)):
            delta_tomo = tomographic_objects.Delta(name=namefile[j], pk1d_type=True)
            delta_tomo.read()
            for k in range(len(delta_tomo.delta_array)):
                if (
                    tomographic_objects.Delta.primary_key(delta_tomo.delta_array[k])
                    == ids[i]
                ):
                    Deltas[ids[i]].append(delta_tomo.delta_array[k])
    return (Deltas, ids)


def delete_los_extrema(listz, listsigmas, listdelta):
    # Get the list of common elements in the pack to have minimum and maximum redshifts
    zmin, zmax = 0, 10**10
    zcommon = listz[0]
    for j in range(1, len(listz)):
        zcommon = list(set(zcommon).intersection(listz[j]))
    zmin = np.min(zcommon)
    zmax = np.max(zcommon)

    # Deleting pixels at the beginning and at the end of LOSs
    for j in range(len(listz)):
        lineToDeleteFirst = 0
        k = 0
        while listz[j][k] != zmin:
            lineToDeleteFirst = lineToDeleteFirst + 1
            k = k + 1
        lineToDeleteLast = 0
        k = -1
        while listz[j][k] != zmax:
            lineToDeleteLast = lineToDeleteLast + 1
            k = k - 1
        listz[j] = listz[j][lineToDeleteFirst : len(listz[j]) - lineToDeleteLast]
        listsigmas[j] = listsigmas[j][
            lineToDeleteFirst : len(listsigmas[j]) - lineToDeleteLast
        ]
        listdelta[j] = listdelta[j][
            lineToDeleteFirst : len(listdelta[j]) - lineToDeleteLast
        ]
    return (listz, listsigmas, listdelta)


def delete_missing_pixels(listz, listsigmas, listdelta):
    # Ensuring that all LOS have the same lenght
    lenlists = []
    for j in range(len(listz)):
        lenlists.append(len(listz[j]))

    # Selection of the pixels to delete in case of a missing pixel along one LOS + Deletion
    while np.max(lenlists) != np.min(lenlists):
        eltTodelete = [[] for j in range(len(listz))]
        mi = 10**10
        mins = []
        for j in range(len(lenlists)):
            if lenlists[j] < mi:
                mi = lenlists[j]
                mins = [j]
            elif lenlists[j] == mi:
                mins.append(j)
        for j in range(len(mins)):
            for k in range(len(listz)):
                for m in range(len(listz[k])):
                    if (k != mins[j]) & (
                        np.isin([listz[k][m]], listz[mins[j]]) == False
                    ):
                        eltTodelete[k].append(m)
        newlistz, newlistsigma, newlistdelta = (
            [[] for n in range(len(listz))],
            [[] for n in range(len(listz))],
            [[] for n in range(len(listz))],
        )
        for j in range(len(eltTodelete)):
            eltTodeletej = list(set(eltTodelete[j]))
            for k in range(len(listz[j])):
                if np.isin([k], eltTodeletej) == False:
                    newlistz[j].append(listz[j][k])
                    newlistsigma[j].append(listsigmas[j][k])
                    newlistdelta[j].append(listdelta[j][k])
        listz = newlistz
        listsigmas = newlistsigma
        listdelta = newlistdelta
        lenlists = []
        for j in range(len(listz)):
            lenlists.append(len(listz[j]))
        return (listz, listsigmas, listdelta, lenlists)


def compute_shape_size_parallel(
    extremum_coord, number_chunks, overlaping, shape_sub_map
):
    if overlaping is None:
        overlaping = 0.0
    minx, maxx = extremum_coord[0], extremum_coord[1]
    miny, maxy = extremum_coord[2], extremum_coord[3]
    minz, maxz = extremum_coord[4], extremum_coord[5]
    intervalx = maxx - minx
    intervaly = maxy - miny
    intervalz = maxz - minz
    subIntervalx = intervalx / number_chunks[0]
    subIntervaly = intervaly / number_chunks[1]
    shape_x = number_chunks[0] * shape_sub_map[0]
    shape_y = number_chunks[1] * shape_sub_map[1]
    remove_shape_x, remove_shape_y = 0, 0
    for i in range(number_chunks[0]):
        for j in range(number_chunks[1]):
            if (i == number_chunks[0] - 1) & (i == 0):
                intervalxChunk = [i * subIntervalx, (i + 1) * subIntervalx]
            elif i == 0:
                intervalxChunk = [
                    i * subIntervalx,
                    (i + 1) * subIntervalx + overlaping,
                ]
            elif i == number_chunks[0] - 1:
                intervalxChunk = [i * subIntervalx - overlaping, intervalx]
            else:
                intervalxChunk = [
                    i * subIntervalx - overlaping,
                    (i + 1) * subIntervalx + overlaping,
                ]
            if (j == number_chunks[1] - 1) & (j == 0):
                intervalyChunk = [j * subIntervaly, (j + 1) * subIntervaly]
            elif j == 0:
                intervalyChunk = [
                    j * subIntervaly,
                    (j + 1) * subIntervaly + overlaping,
                ]
            elif j == number_chunks[1] - 1:
                intervalyChunk = [j * subIntervaly - overlaping, intervaly]
            else:
                intervalyChunk = [
                    j * subIntervaly - overlaping,
                    (j + 1) * subIntervaly + overlaping,
                ]
            size = (
                intervalxChunk[1] - intervalxChunk[0],
                intervalyChunk[1] - intervalyChunk[0],
                intervalz,
            )
            pixel_to_remove = np.around(
                utils.pixel_per_mpc(size, shape_sub_map) * overlaping, 0
            ).astype(int)
            if number_chunks[0] != 1:
                if (i == 0) | (i == number_chunks[0] - 1):
                    remove_shape_x = remove_shape_x + pixel_to_remove[0]
                else:
                    remove_shape_x = remove_shape_x + 2 * pixel_to_remove[0]
            if number_chunks[1] != 1:
                if (j == 0) | (j == number_chunks[1] - 1):
                    remove_shape_y = remove_shape_y + pixel_to_remove[1]
                else:
                    remove_shape_y = remove_shape_y + 2 * pixel_to_remove[1]
    shape_x = shape_x - remove_shape_x // number_chunks[1]
    shape_y = shape_y - remove_shape_y // number_chunks[0]

    size = (maxx - minx, maxy - miny, maxz - minz)
    shape = (shape_x, shape_y, shape_sub_map[2])
    return (shape, size)


def compute_shape_size_parallel_from_interface(
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
):
    ra = np.linspace(ramin, ramax, N_coord_edge)
    dec = np.linspace(decmin, decmax, N_coord_edge)
    z = np.linspace(zmin, zmax, N_coord_edge)

    cube_edge_coord = np.concatenate(
        [
            np.array([[ra[i], decmin, zmin] for i in range(N_coord_edge)]),
            np.array([[ramin, dec[i], zmin] for i in range(N_coord_edge)]),
            np.array([[ramin, decmin, z[i]] for i in range(N_coord_edge)]),
            np.array([[ramax, dec[i], zmin] for i in range(N_coord_edge)]),
            np.array([[ramax, decmin, z[i]] for i in range(N_coord_edge)]),
            np.array([[ra[i], decmax, zmin] for i in range(N_coord_edge)]),
            np.array([[ramin, decmax, z[i]] for i in range(N_coord_edge)]),
            np.array([[ra[i], decmin, zmax] for i in range(N_coord_edge)]),
            np.array([[ramin, dec[i], zmax] for i in range(N_coord_edge)]),
            np.array([[ra[i], decmax, zmax] for i in range(N_coord_edge)]),
            np.array([[ramax, dec[i], zmax] for i in range(N_coord_edge)]),
            np.array([[ramax, decmax, z[i]] for i in range(N_coord_edge)]),
        ]
    )

    (rcomov, distang, _, _) = utils.get_cosmo_function(Omega_m)
    suplementary_parameters = utils.return_suplementary_parameters(
        coordinate_transform, zmin=zmin, zmax=zmax
    )

    cartesian_cube_edge_coord = np.zeros_like(cube_edge_coord)

    (
        cartesian_cube_edge_coord[:, 0],
        cartesian_cube_edge_coord[:, 1],
        cartesian_cube_edge_coord[:, 2],
    ) = utils.convert_sky_to_cartesian(
        np.radians(cube_edge_coord[:, 0]),
        np.radians(cube_edge_coord[:, 1]),
        cube_edge_coord[:, 2],
        coordinate_transform,
        rcomov=rcomov,
        distang=distang,
        suplementary_parameters=suplementary_parameters,
    )

    Xmin = np.min(cartesian_cube_edge_coord[:, 0])
    Xmax = np.max(cartesian_cube_edge_coord[:, 0])
    Ymin = np.min(cartesian_cube_edge_coord[:, 1])
    Ymax = np.max(cartesian_cube_edge_coord[:, 1])
    Zmin = np.min(cartesian_cube_edge_coord[:, 2])
    Zmax = np.max(cartesian_cube_edge_coord[:, 2])
    extremum_coord = [Xmin, Xmax, Ymin, Ymax, Zmin, Zmax]

    shape, size = compute_shape_size_parallel(
        extremum_coord,
        number_chunks,
        overlaping,
        shape_sub_map,
    )
    return shape, size


#############################################################################
#############################################################################
############################### CLASSES #####################################
#############################################################################
#############################################################################


class DeltaModifier(object):
    def __init__(self, pwd, delta_path):
        self.pwd = pwd
        self.delta_path = delta_path

    # CR - need to rethink the way Delta class is done (+ add shuffle inside)

    def shuffle_deltas(self, other_delta_path=None, other_path_out=None, seed=None):
        namefile = get_delta_list(self.delta_path)
        namefile_other = None
        if other_delta_path is not None:
            namefile_other = get_delta_list(other_delta_path)

        (delta, ivar, delta_other, weight_other) = self.get_delta_sigma_array(
            namefile, namefile_other=namefile_other
        )
        if seed is None:
            seed = np.random.randint(10000000)
        np.random.seed(seed)
        ivar_rand = np.random.permutation(ivar)
        delta_rand = np.random.permutation(delta)
        if other_delta_path is not None:
            np.random.seed(seed)
            weight_other_rand = np.random.permutation(weight_other)
            delta_other_rand = np.random.permutation(delta_other)

        self.write_delta_sigma_array(
            delta_rand,
            ivar_rand,
            namefile,
            other_delta_path=other_delta_path,
            namefile_other=namefile_other,
            delta_other_rand=delta_other_rand,
            weight_other_rand=weight_other_rand,
            other_path_out=other_path_out,
        )

    def shuffle_deltas_cut_z(
        self, n_cut, zmin, zmax, other_delta_path=None, other_path_out=None, seed=None
    ):
        namefile = get_delta_list(self.delta_path)
        namefile_other = None
        if other_delta_path is not None:
            namefile_other = get_delta_list(other_delta_path)
        redshift_cut = np.linspace(zmin, zmax, n_cut + 1)
        (delta, ivar, delta_other, weight_other) = self.get_delta_sigma_array_cut_z(
            namefile, redshift_cut, namefile_other=namefile_other
        )

        if seed is None:
            seed = np.random.randint(10000000)
        ivar_rand, delta_rand = [], []
        np.random.seed(seed)
        for k in range(n_cut):
            ivar_rand.append(np.random.permutation(ivar[k]))
        np.random.seed(seed)
        for k in range(n_cut):
            delta_rand.append(np.random.permutation(delta[k]))
        if other_delta_path is not None:
            weight_other_rand, delta_other_rand = [], []
            np.random.seed(seed)
            for k in range(n_cut):
                weight_other_rand.append(np.random.permutation(weight_other[k]))
            np.random.seed(seed)
            for k in range(n_cut):
                delta_other_rand.append(np.random.permutation(delta_other[k]))

        self.write_delta_sigma_array_cut_z(
            delta_rand,
            ivar_rand,
            namefile,
            n_cut,
            redshift_cut,
            other_delta_path=other_delta_path,
            namefile_other=namefile_other,
            delta_other_rand=delta_other_rand,
            weight_other_rand=weight_other_rand,
            other_path_out=other_path_out,
        )

    def get_delta_sigma_array(self, namefile, namefile_other=None):
        weight_other, delta_other = None, None
        if namefile_other is not None:
            weight_other, delta_other = [], []
        ivar, delta = [], []
        for i in range(len(namefile)):
            delta_tomo = tomographic_objects.Delta(name=namefile[i], pk1d_type=True)
            delta_tomo.read()
            for j in range(len(delta_tomo.delta_array)):
                ivar.append(delta_tomo.delta_array[j].ivar)
                delta.append(delta_tomo.delta_array[j].delta)
            if namefile_other is not None:
                delta_tomo_other = tomographic_objects.Delta(
                    name=namefile_other[i], pk1d_type=False
                )
                delta_tomo_other.read()
                for j in range(len(delta_tomo_other.delta_array)):
                    weight_other.append(delta_tomo_other.delta_array[j].weights)
                    delta_other.append(delta_tomo_other.delta_array[j].delta)
        ivar = np.concatenate(ivar, axis=0)
        delta = np.concatenate(delta, axis=0)
        if namefile_other is not None:
            weight_other = np.concatenate(weight_other, axis=0)
            delta_other = np.concatenate(delta_other, axis=0)
        return (delta, ivar, delta_other, weight_other)

    def write_delta_sigma_array(
        self,
        delta_rand,
        ivar_rand,
        namefile,
        other_delta_path=None,
        namefile_other=None,
        delta_other_rand=None,
        weight_other_rand=None,
        other_path_out=None,
    ):
        ibegin = 0
        for i in range(len(namefile)):
            delta_tomo = tomographic_objects.Delta(name=namefile[i], pk1d_type=True)
            delta_tomo.read()
            delta_obj_list = []
            for j in range(len(delta_tomo.delta_array)):
                delta_tomo.delta_array[j].de = delta_rand[
                    ibegin : ibegin + len(delta_tomo.delta_array[j].de)
                ]
                delta_tomo.delta_array[j].iv = ivar_rand[
                    ibegin : ibegin + len(delta_tomo.delta_array[j].iv)
                ]
                ibegin = ibegin + len(delta_tomo.delta_array[j].de)
                delta_obj_list.append(delta_tomo.delta_array[j])
            name_delta = os.path.join(
                self.pwd, namefile[i].split(self.delta_path)[-1].split("/")[-1]
            )
            new_delta_tomo = tomographic_objects.Delta(name=name_delta, pk1d_type=True)
            new_delta_tomo.delta_array = delta_obj_list
            new_delta_tomo.write()

        if other_delta_path is not None:
            ibegin = 0
            for i in range(len(namefile_other)):
                delta_tomo_other = tomographic_objects.Delta(
                    name=namefile_other[i], pk1d_type=False
                )
                delta_tomo_other.read()
                delta_obj_list = []
                for j in range(len(delta_tomo_other.delta_array)):
                    delta_tomo_other.delta_array[j].de = delta_other_rand[
                        ibegin : ibegin + len(delta_tomo_other.delta_array[j].de)
                    ]
                    delta_tomo_other.delta_array[j].we = weight_other_rand[
                        ibegin : ibegin + len(delta_tomo_other.delta_array[j].we)
                    ]
                    ibegin = ibegin + len(delta_tomo_other.delta_array[j].de)
                    delta_obj_list.append(delta_tomo_other.delta_array[j])
                name_delta = os.path.join(
                    other_path_out,
                    namefile_other[i].split(other_delta_path)[-1].split("/")[-1],
                )
                new_delta_tomo = tomographic_objects.Delta(
                    name=name_delta, pk1d_type=False
                )
                new_delta_tomo.delta_array = delta_obj_list
                new_delta_tomo.write()

    def get_delta_sigma_array_cut_z(self, namefile, redshift_cut, namefile_other=None):
        n_cut = len(redshift_cut) - 1
        weight_other, delta_other = None, None
        if namefile_other is not None:
            weight_other, delta_other = [[] for i in range(n_cut)], [
                [] for i in range(n_cut)
            ]
        ivar, delta = [[] for i in range(n_cut)], [[] for i in range(n_cut)]
        for i in range(len(namefile)):
            delta_tomo = tomographic_objects.Delta(name=namefile[i], pk1d_type=True)
            delta_tomo.read()
            for j in range(len(delta_tomo.delta_array)):
                redshift = (
                    10 ** delta_tomo.delta_array[j].log_lambda / utils.lambdaLy
                ) - 1
                for k in range(n_cut):
                    mask = (redshift >= redshift_cut[k]) & (
                        redshift < redshift_cut[k + 1]
                    )
                    ivar[k].append(delta_tomo.delta_array[j].ivar[mask])
                    delta[k].append(delta_tomo.delta_array[j].delta[mask])
            if namefile_other is not None:
                delta_tomo_other = tomographic_objects.Delta(
                    name=namefile_other[i], pk1d_type=False
                )
                delta_tomo_other.read()
                for j in range(len(delta_tomo_other.delta_array)):
                    redshift = (
                        10 ** delta_tomo_other.delta_array[j].log_lambda
                        / utils.lambdaLy
                    ) - 1
                    for k in range(n_cut):
                        mask = (redshift >= redshift_cut[k]) & (
                            redshift < redshift_cut[k + 1]
                        )
                        weight_other[k].append(
                            delta_tomo_other.delta_array[j].weights[mask]
                        )
                        delta_other[k].append(
                            delta_tomo_other.delta_array[j].delta[mask]
                        )
        for k in range(n_cut):
            ivar[k] = np.concatenate(ivar[k], axis=0)
            delta[k] = np.concatenate(delta[k], axis=0)
        if namefile_other is not None:
            for k in range(n_cut):
                weight_other[k] = np.concatenate(weight_other[k], axis=0)
                delta_other[k] = np.concatenate(delta_other[k], axis=0)
        return (delta, ivar, delta_other, weight_other)

    def write_delta_sigma_array_cut_z(
        self,
        delta_rand,
        ivar_rand,
        namefile,
        n_cut,
        redshift_cut,
        other_delta_path=None,
        namefile_other=None,
        delta_other_rand=None,
        weight_other_rand=None,
        other_path_out=None,
    ):
        ibegin = [0 for i in range(n_cut)]
        for i in range(len(namefile)):
            delta_tomo = tomographic_objects.Delta(name=namefile[i], pk1d_type=True)
            delta_tomo.read()
            delta_obj_list = []
            for j in range(len(delta_tomo.delta_array)):
                redshift = (
                    10 ** delta_tomo.delta_array[j].log_lambda / utils.lambdaLy
                ) - 1
                for k in range(n_cut):
                    mask = (redshift >= redshift_cut[k]) & (
                        redshift < redshift_cut[k + 1]
                    )
                    delta_tomo.delta_array[j].de[mask] = delta_rand[k][
                        ibegin[k] : ibegin[k] + len(delta_tomo.delta_array[j].de[mask])
                    ]
                    delta_tomo.delta_array[j].iv[mask] = ivar_rand[k][
                        ibegin[k] : ibegin[k] + len(delta_tomo.delta_array[j].iv[mask])
                    ]
                    ibegin[k] = ibegin[k] + len(delta_tomo.delta_array[j].de[mask])
                delta_obj_list.append(delta_tomo.delta_array[j])
            name_delta = os.path.join(
                self.pwd, namefile[i].split(self.delta_path)[-1].split("/")[-1]
            )
            new_delta_tomo = tomographic_objects.Delta(name=name_delta, pk1d_type=True)
            new_delta_tomo.delta_array = delta_obj_list
            new_delta_tomo.write()

        if other_delta_path is not None:
            ibegin = [0 for i in range(n_cut)]
            for i in range(len(namefile_other)):
                delta_tomo_other = tomographic_objects.Delta(
                    name=namefile_other[i], pk1d_type=False
                )
                delta_tomo_other.read()
                delta_obj_list = []
                for j in range(len(delta_tomo_other.delta_array)):
                    redshift = (
                        10 ** delta_tomo_other.delta_array[j].log_lambda
                        / utils.lambdaLy
                    ) - 1
                    for k in range(n_cut):
                        mask = (redshift >= redshift_cut[k]) & (
                            redshift < redshift_cut[k + 1]
                        )
                        delta_tomo_other.delta_array[j].de[mask] = delta_other_rand[k][
                            ibegin[k] : ibegin[k]
                            + len(delta_tomo_other.delta_array[j].de[mask])
                        ]
                        delta_tomo_other.delta_array[j].we[mask] = weight_other_rand[k][
                            ibegin[k] : ibegin[k]
                            + len(delta_tomo_other.delta_array[j].we[mask])
                        ]
                        ibegin[k] = ibegin[k] + len(
                            delta_tomo_other.delta_array[j].de[mask]
                        )
                    delta_obj_list.append(delta_tomo_other.delta_array[j])
                name_delta = os.path.join(
                    other_path_out,
                    namefile_other[i].split(other_delta_path)[-1].split("/")[-1],
                )
                new_delta_tomo = tomographic_objects.Delta(
                    name=name_delta, pk1d_type=False
                )
                new_delta_tomo.delta_array = delta_obj_list
                new_delta_tomo.write()

    def get_new_healpix(
        self,
        number_cut,
        random_density_parameter=None,
        number_repeat=1,
        iterative_selection_parameters=None,
        ra_cut_min=None,
        ra_cut_max=None,
        dec_cut_min=None,
        dec_cut_max=None,
        z_cut_min=None,
        z_cut_max=None,
        center_ra=True,
    ):
        if random_density_parameter is None:
            deltas = self.create_healpix(
                number_cut,
                random=None,
                return_len_ra=False,
                ra_cut_min=ra_cut_min,
                ra_cut_max=ra_cut_max,
                dec_cut_min=dec_cut_min,
                dec_cut_max=dec_cut_max,
                center_ra=center_ra,
            )
        else:
            if number_repeat == 1:
                deltas = self.create_healpix(
                    number_cut,
                    random=random_density_parameter,
                    return_len_ra=False,
                    ra_cut_min=ra_cut_min,
                    ra_cut_max=ra_cut_max,
                    dec_cut_min=dec_cut_min,
                    dec_cut_max=dec_cut_max,
                    center_ra=center_ra,
                )
            else:
                if iterative_selection_parameters is None:
                    return KeyError(
                        "Please dictionary parameter for the iterative selection"
                    )
                deltas = self.iterate_healpix_creation(
                    number_cut,
                    random_density_parameter,
                    number_repeat,
                    iterative_selection_parameters,
                    ra_cut_min=ra_cut_min,
                    ra_cut_max=ra_cut_max,
                    dec_cut_min=dec_cut_min,
                    dec_cut_max=dec_cut_max,
                    z_cut_min=z_cut_min,
                    z_cut_max=z_cut_max,
                )
        return deltas

    def create_healpix(
        self,
        number_cut,
        random=None,
        return_len_ra=False,
        ra_cut_min=None,
        ra_cut_max=None,
        dec_cut_min=None,
        dec_cut_max=None,
        center_ra=True,
    ):
        namefile = get_delta_list(self.delta_path)
        deltas = {}
        ra_array = []
        dec_array = []
        for cut in range(number_cut):
            deltas[cut] = []
        for i in range(len(namefile)):
            delta_tomo = tomographic_objects.Delta(name=namefile[i], pk1d_type=True)
            delta_tomo.read()
            for j in range(len(delta_tomo.delta_array)):
                if center_ra:
                    if delta_tomo.delta_array[j].ra * 180 / np.pi > 180:
                        ra = (delta_tomo.delta_array[j].ra * 180 / np.pi) - 360
                    else:
                        ra = delta_tomo.delta_array[j].ra * 180 / np.pi
                else:
                    ra = delta_tomo.delta_array[j].ra * 180 / np.pi
                dec = delta_tomo.delta_array[j].dec * 180 / np.pi
                if (
                    (ra > ra_cut_min)
                    & (ra < ra_cut_max)
                    & (dec > dec_cut_min)
                    & (dec < dec_cut_max)
                ):
                    ra_array.append(ra)
                    dec_array.append(dec)
                    for cut in range(number_cut):
                        interval_ra = ((cut) / (number_cut)) * (
                            ra_cut_max - ra_cut_min
                        ) + ra_cut_min, ((cut + 1) / (number_cut)) * (
                            ra_cut_max - ra_cut_min
                        ) + ra_cut_min
                        if (ra > interval_ra[0]) & (ra <= interval_ra[1]):
                            if center_ra:
                                delta_tomo.delta_array[j].ra = (
                                    delta_tomo.delta_array[j].ra - 2 * np.pi
                                )
                            deltas[cut].append(delta_tomo.delta_array[j])
        if random is not None:
            deltas = self.randomize_choice_of_los(
                deltas,
                random,
                number_cut,
                len(ra_array),
                ra_cut_min=ra_cut_min,
                ra_cut_max=ra_cut_max,
                dec_cut_min=dec_cut_min,
                dec_cut_max=dec_cut_max,
            )
        if return_len_ra:
            return (deltas, len(ra_array))
        else:
            return deltas

    def randomize_choice_of_los(
        self,
        deltas_dict,
        random,
        number_cut,
        number_ra,
        ra_cut_min=None,
        ra_cut_max=None,
        dec_cut_min=None,
        dec_cut_max=None,
    ):
        deltas = deltas_dict.copy()
        density = number_ra / ((ra_cut_max - ra_cut_min) * (dec_cut_max - dec_cut_min))
        utils.Logger.add("density before random choice =" + str(density))
        random_cut = random / density
        ra_random = []
        dec_random = []
        for cut in range(number_cut):
            number_of_delta_to_select = int(round(len(deltas[cut]) * random_cut, 0))
            deltas[cut] = sample(deltas[cut], number_of_delta_to_select)
            for i in range(len(deltas[cut])):
                ra_random.append(deltas[cut][i].ra * 180 / np.pi)
                dec_random.append(deltas[cut][i].dec * 180 / np.pi)
        density = len(ra_random) / (
            (ra_cut_max - ra_cut_min) * (dec_cut_max - dec_cut_min)
        )
        utils.Logger.add("density after random choice =" + str(density))
        return deltas

    def iterate_healpix_creation(
        self,
        number_cut,
        random_density_parameter,
        number_repeat,
        iterative_selection_parameters,
        property_file_name,
        ra_cut_min=None,
        ra_cut_max=None,
        dec_cut_min=None,
        dec_cut_max=None,
        z_cut_min=None,
        z_cut_max=None,
        center_ra=True,
    ):
        density_names = iterative_selection_parameters["density_names"]
        dperp_names = iterative_selection_parameters["separation_names"]
        Om = iterative_selection_parameters["Om"]
        (rcomov, distang, inv_rcomov, inv_distang) = utils.get_cosmo_function(Om)
        coordinate_transform = iterative_selection_parameters["coordinate_transform"]
        suplementary_parameters = utils.return_suplementary_parameters(
            coordinate_transform, zmin=z_cut_min, zmax=z_cut_max
        )
        density_ref = {}
        dperp_ref = {}
        for cut in range(number_cut):
            density_ref[cut] = PixelAnalizer.read_density_file(density_names[cut])[1]
            dperp_ref[cut] = PixelAnalizer.read_dperp_file(dperp_names[cut])[1]
        min_density, min_dperp = np.inf, np.inf
        delta_to_keep = {}
        utils.Logger.add("Beginning of the random iteration selection")
        deltas_dict, number_ra = self.create_healpix(
            number_cut, random=False, return_len_ra=True, center_ra=center_ra
        )
        deltas_random = {}
        for i in range(number_repeat):
            utils.Logger.add("repeat " + str(i))
            utils.Logger.add("deltas " + str(i) + " computed")
            diff_dperp, diff_density = [], []
            deltas_random = self.randomize_choice_of_los(
                deltas_dict,
                random_density_parameter,
                number_cut,
                number_ra,
                ra_cut_min=ra_cut_min,
                ra_cut_max=ra_cut_max,
                dec_cut_min=dec_cut_min,
                dec_cut_max=dec_cut_max,
            )
            for cut in range(number_cut):
                delta_file = tomographic_objects.Delta(
                    delta_file=None,
                    delta_array=deltas_random[cut],
                    name="delta_" + str(cut) + "_to_test.pickle",
                )
                delta_file.write()
                namefile = "delta_" + str(cut) + "_to_test.pickle"
                (ra, dec, z, zqso, ids, sigmas, deltas) = get_deltas(namefile)
                sky_deltas = np.array(
                    [
                        [ra[i], dec[i], z[i][j], sigmas[i][j], deltas[i][j]]
                        for i in range(len(ra))
                        for j in range(len(z[i]))
                    ]
                )
                sky_deltas = sky_deltas[
                    utils.cut_sky_catalog(
                        sky_deltas[:, 0],
                        sky_deltas[:, 1],
                        sky_deltas[:, 2],
                        ramin=ra_cut_min,
                        ramax=ra_cut_max,
                        decmin=dec_cut_min,
                        decmax=dec_cut_max,
                        zmin=z_cut_min,
                        zmax=z_cut_max,
                    )
                ]
                cartesian_deltas = np.zeros(sky_deltas.shape)
                (
                    cartesian_deltas[:, 0],
                    cartesian_deltas[:, 1],
                    cartesian_deltas[:, 2],
                ) = utils.convert_sky_to_cartesian(
                    sky_deltas[:, 0],
                    sky_deltas[:, 1],
                    sky_deltas[:, 2],
                    coordinate_transform,
                    rcomov=rcomov,
                    distang=distang,
                    suplementary_parameters=suplementary_parameters,
                )
                pixel = tomographic_objects.Pixel.init_from_property_files(
                    property_file_name, pixel_array=cartesian_deltas, name=None
                )
                pixel_analyzer = PixelAnalizer(pixel=pixel)
                (
                    zpar,
                    dperpz,
                    densityz,
                ) = pixel_analyzer.compute_plot_mean_distance_density("", plot=False)
                diff_dperp.append(
                    np.mean(abs(np.array(dperpz) - np.array(dperp_ref[cut])))
                )
                diff_density.append(
                    np.mean(abs(np.array(densityz) - np.array(density_ref[cut])))
                )
            utils.Logger.add(
                "Mean difference in term LOS density : {}".format(np.mean(diff_density))
            )
            utils.Logger.add(
                "Mean difference in term of Mean LOS separation : {}".format(
                    np.mean(diff_dperp)
                )
            )
            if (np.mean(diff_density) < min_density) & (
                np.mean(diff_dperp) < min_dperp
            ):
                utils.Logger.add("Better at the repeat " + str(i))
                min_density = np.mean(diff_density)
                min_dperp = np.mean(diff_dperp)
                delta_to_keep = deltas_random
        for cut in range(number_cut):
            os.remove("delta_" + str(cut) + "_to_test.pickle")
        utils.Logger.add("End of the random iteration selection")
        return delta_to_keep

    def save_deltas(self, deltas, name_out, number_cut):
        for cut in range(number_cut):
            delta = tomographic_objects.Delta(
                name=os.path.join(self.pwd, f"{name_out}_{cut}.fits"),
                delta_array=deltas[cut],
            )
            delta.write()

    def subsample_deltas(
        self,
        name_out,
        number_cut,
        random_density_parameter=None,
        number_repeat=1,
        iterative_selection_parameters=None,
        ra_cut_min=None,
        ra_cut_max=None,
        dec_cut_min=None,
        dec_cut_max=None,
        z_cut_min=None,
        z_cut_max=None,
        center_ra=True,
    ):
        deltas = self.get_new_healpix(
            number_cut,
            random_density_parameter=random_density_parameter,
            number_repeat=number_repeat,
            iterative_selection_parameters=iterative_selection_parameters,
            ra_cut_min=ra_cut_min,
            ra_cut_max=ra_cut_max,
            dec_cut_min=dec_cut_min,
            dec_cut_max=dec_cut_max,
            z_cut_min=z_cut_min,
            z_cut_max=z_cut_max,
            center_ra=center_ra,
        )
        self.save_deltas(deltas, name_out, number_cut)


class DeltaConverter:
    def __init__(
        self,
        pwd,
        Omega_m,
        delta_path,
        coordinate_transform,
        plot_pixel_properties,
        software,
        return_qso_catalog=None,
        return_dla_catalog=None,
        dla_catalog=None,
        return_sky_catalogs=False,
        repeat=False,
        center_ra=True,
    ):
        self.pwd = pwd
        self.delta_path = delta_path
        self.Omega_m = Omega_m
        self.coordinate_transform = coordinate_transform
        self.plot_pixel_properties = plot_pixel_properties
        self.return_qso_catalog = return_qso_catalog
        self.return_dla_catalog = return_dla_catalog
        self.dla_catalog = dla_catalog
        self.return_sky_catalogs = return_sky_catalogs
        self.repeat = repeat
        self.software = software
        self.center_ra = center_ra

    def transform_delta_to_pixel_file(
        self,
        rebin=None,
        sigma_min=None,
        sigma_max=None,
        z_cut_min=None,
        z_cut_max=None,
        dec_cut_min=None,
        dec_cut_max=None,
        ra_cut_min=None,
        ra_cut_max=None,
    ):
        namefile = get_delta_list(self.delta_path)
        # namefile = preselect_deltas(namefile,
        #                             ramin=ra_cut_min,
        #                             ramax=ra_cut_max,
        #                             decmin=dec_cut_min,
        #                             decmax=dec_cut_max)
        properties_map_pixels = {}
        (rcomov, distang, inv_rcomov, inv_distang) = utils.get_cosmo_function(
            self.Omega_m
        )
        if self.repeat:
            (
                ra,
                dec,
                z,
                zqso,
                ids,
                sigmas,
                deltas,
            ) = get_merged_multiple_exposure_deltas(namefile)
        else:
            (ra, dec, z, zqso, ids, sigmas, deltas) = get_deltas(
                namefile, center_ra=self.center_ra
            )

        if rebin is not None:
            (z, deltas, sigmas) = self.rebin_data(z, deltas, sigmas, rebin)

        if self.return_dla_catalog is not None:
            zdlas, z_qso_dlas = [], []
            if self.dla_catalog is None:
                raise KeyError(
                    "Please give a DLA catalog name or turn off the return_dla_catalog option"
                )
            dla_catalog = tomographic_objects.DLACatalog.init_from_fits(
                self.dla_catalog
            )
            for i in range(len(ids)):
                mask = dla_catalog.primary_key == ids[i]
                if z_cut_min is not None:
                    mask &= dla_catalog.coord_z > z_cut_min
                if z_cut_max is not None:
                    mask &= dla_catalog.coord_z < z_cut_max
                zdlas.append(dla_catalog.coord_z[mask])
                z_qso_dlas.append(dla_catalog.z_qso[mask])

        sky_deltas = np.array(
            [
                [ra[i], dec[i], z[i][j], sigmas[i][j], deltas[i][j]]
                for i in range(len(ra))
                for j in range(len(z[i]))
            ]
        )
        sky_deltas = sky_deltas[
            utils.cut_sky_catalog(
                sky_deltas[:, 0],
                sky_deltas[:, 1],
                sky_deltas[:, 2],
                ramin=ra_cut_min,
                ramax=ra_cut_max,
                decmin=dec_cut_min,
                decmax=dec_cut_max,
                zmin=z_cut_min,
                zmax=z_cut_max,
            )
        ]
        suplementary_parameters = utils.return_suplementary_parameters(
            self.coordinate_transform,
            zmin=np.min(sky_deltas[:, 2]),
            zmax=np.max(sky_deltas[:, 2]),
        )
        cartesian_deltas = np.zeros(sky_deltas.shape)
        (
            cartesian_deltas[:, 0],
            cartesian_deltas[:, 1],
            cartesian_deltas[:, 2],
        ) = utils.convert_sky_to_cartesian(
            sky_deltas[:, 0],
            sky_deltas[:, 1],
            sky_deltas[:, 2],
            self.coordinate_transform,
            rcomov=rcomov,
            distang=distang,
            suplementary_parameters=suplementary_parameters,
        )
        cartesian_deltas[:, 3], cartesian_deltas[:, 4] = (
            sky_deltas[:, 3],
            sky_deltas[:, 4],
        )
        (
            properties_map_pixels["minx"],
            properties_map_pixels["miny"],
            properties_map_pixels["minz"],
        ) = (
            np.min(cartesian_deltas[:, 0]),
            np.min(cartesian_deltas[:, 1]),
            np.min(cartesian_deltas[:, 2]),
        )
        (
            properties_map_pixels["maxx"],
            properties_map_pixels["maxy"],
            properties_map_pixels["maxz"],
        ) = (
            np.max(cartesian_deltas[:, 0]),
            np.max(cartesian_deltas[:, 1]),
            np.max(cartesian_deltas[:, 2]),
        )
        (
            properties_map_pixels["minra"],
            properties_map_pixels["mindec"],
            properties_map_pixels["minredshift"],
        ) = (
            np.min(sky_deltas[:, 0]),
            np.min(sky_deltas[:, 1]),
            np.min(sky_deltas[:, 2]),
        )
        (
            properties_map_pixels["maxra"],
            properties_map_pixels["maxdec"],
            properties_map_pixels["maxredshift"],
        ) = (
            np.max(sky_deltas[:, 0]),
            np.max(sky_deltas[:, 1]),
            np.max(sky_deltas[:, 2]),
        )
        cartesian_deltas = cartesian_deltas - np.array(
            [
                properties_map_pixels["minx"],
                properties_map_pixels["miny"],
                properties_map_pixels["minz"],
                0,
                0,
            ]
        )

        if sigma_min is not None:
            cartesian_deltas[:, 3][cartesian_deltas[:, 3] < sigma_min] = sigma_min
        if sigma_max is not None:
            cartesian_deltas[:, 3][cartesian_deltas[:, 3] > sigma_max] = sigma_max

        if self.return_qso_catalog is not None:
            sky_qso_catalog = np.array(
                [[ra[i], dec[i], zqso[i], ids[i]] for i in range(len(ra))]
            )
            sky_qso_catalog = sky_qso_catalog[
                utils.cut_sky_catalog(
                    sky_qso_catalog[:, 0],
                    sky_qso_catalog[:, 1],
                    sky_qso_catalog[:, 2],
                    ramin=ra_cut_min,
                    ramax=ra_cut_max,
                    decmin=dec_cut_min,
                    decmax=dec_cut_max,
                    zmin=z_cut_min,
                    zmax=z_cut_max,
                ),
                :,
            ]
            cartesian_qso_catalog = np.zeros(sky_qso_catalog.shape)
            (
                cartesian_qso_catalog[:, 0],
                cartesian_qso_catalog[:, 1],
                cartesian_qso_catalog[:, 2],
            ) = utils.convert_sky_to_cartesian(
                sky_qso_catalog[:, 0],
                sky_qso_catalog[:, 1],
                sky_qso_catalog[:, 2],
                self.coordinate_transform,
                rcomov=rcomov,
                distang=distang,
                suplementary_parameters=suplementary_parameters,
            )
            cartesian_qso_catalog[:, 3] = sky_qso_catalog[:, 3]
            cartesian_qso_catalog = cartesian_qso_catalog - np.array(
                [
                    properties_map_pixels["minx"],
                    properties_map_pixels["miny"],
                    properties_map_pixels["minz"],
                    0,
                ]
            )

        else:
            sky_qso_catalog, cartesian_qso_catalog = None, None

        if self.return_dla_catalog is not None:
            sky_dla_catalog = np.array(
                [
                    [ra[i], dec[i], zdlas[i][j], z_qso_dlas[i][j]]
                    for i in range(len(ra))
                    for j in range(len(zdlas[i]))
                ]
            )
            sky_dla_catalog = sky_dla_catalog[
                utils.cut_sky_catalog(
                    sky_dla_catalog[:, 0],
                    sky_dla_catalog[:, 1],
                    sky_dla_catalog[:, 2],
                    ramin=ra_cut_min,
                    ramax=ra_cut_max,
                    decmin=dec_cut_min,
                    decmax=dec_cut_max,
                    zmin=z_cut_min,
                    zmax=z_cut_max,
                )
            ]
            cartesian_dla_catalog = np.zeros(sky_dla_catalog.shape)
            (
                cartesian_dla_catalog[:, 0],
                cartesian_dla_catalog[:, 1],
                cartesian_dla_catalog[:, 2],
            ) = utils.convert_sky_to_cartesian(
                sky_dla_catalog[:, 0],
                sky_dla_catalog[:, 1],
                sky_dla_catalog[:, 2],
                self.coordinate_transform,
                rcomov=rcomov,
                distang=distang,
                suplementary_parameters=suplementary_parameters,
            )
            cartesian_qso_catalog[:, 3] = sky_qso_catalog[:, 3]
            cartesian_dla_catalog = cartesian_dla_catalog - np.array(
                [
                    properties_map_pixels["minx"],
                    properties_map_pixels["miny"],
                    properties_map_pixels["minz"],
                    0,
                ]
            )
        else:
            sky_dla_catalog, cartesian_dla_catalog = None, None

        return (
            cartesian_deltas,
            cartesian_qso_catalog,
            cartesian_dla_catalog,
            sky_deltas,
            sky_qso_catalog,
            sky_dla_catalog,
            properties_map_pixels,
        )

    def rebin_data(self, z, deltas, sigmas, bin_pixel, method="gauss"):
        for i in range(len(z)):
            if len(z[i]) > 1:
                if len(z[i]) <= bin_pixel:
                    z[i] = [np.mean(z[i])]
                    deltas[i] = [np.mean(deltas[i])]
                    sigmas[i] = [np.mean(sigmas[i])]
                else:
                    new_shape = len(z[i]) // bin_pixel
                    first_coord = len(z[i]) - (len(z[i]) // bin_pixel) * bin_pixel
                    if first_coord == 0:
                        z[i] = utils.bin_ndarray(
                            np.array(z[i])[:], [new_shape], operation=method
                        )
                        deltas[i] = utils.bin_ndarray(
                            np.array(deltas[i])[:], [new_shape], operation=method
                        )
                        sigmas[i] = utils.bin_ndarray(
                            np.array(sigmas[i])[:], [new_shape], operation=method
                        )
                    else:
                        z[i] = np.concatenate(
                            [
                                [np.mean(np.array(z[i])[:first_coord])],
                                utils.bin_ndarray(
                                    np.array(z[i])[first_coord:],
                                    [new_shape],
                                    operation=method,
                                ),
                            ]
                        )
                        deltas[i] = np.concatenate(
                            [
                                [np.mean(np.array(deltas[i])[:first_coord])],
                                utils.bin_ndarray(
                                    np.array(deltas[i])[first_coord:],
                                    [new_shape],
                                    operation=method,
                                ),
                            ]
                        )
                        sigmas[i] = np.concatenate(
                            [
                                [np.mean(np.array(sigmas[i])[:first_coord])],
                                utils.bin_ndarray(
                                    np.array(sigmas[i])[first_coord:],
                                    [new_shape],
                                    operation=method,
                                ),
                            ]
                        )
        return (z, deltas, sigmas)

    def create_input_files(
        self, coordinates_to_write, properties, name_pixel, create_launcher=None
    ):
        if self.software.lower() == "dachshund":
            self.create_dachshund_input_files(
                coordinates_to_write,
                properties,
                name_pixel,
                create_launcher=create_launcher,
            )

    def create_dachshund_input_files(
        self, coordinates_to_write, properties, name_pixel, create_launcher=None
    ):
        pixel = tomographic_objects.Pixel(
            name=os.path.join(self.pwd, name_pixel), pixel_array=coordinates_to_write
        )
        pixel.write()
        if create_launcher is not None:
            self.create_dachshund_launcher(
                np.max(coordinates_to_write[:, 0]),
                np.max(coordinates_to_write[:, 1]),
                np.max(coordinates_to_write[:, 2]),
                len(coordinates_to_write),
                properties["shape"][0],
                properties["shape"][1],
                properties["shape"][2],
                properties["sigma_f"],
                properties["lperp"],
                properties["lpar"],
                properties["name_pixel"],
                properties["name_map"],
                create_launcher,
            )

    def create_dachshund_launcher(
        self,
        lx,
        ly,
        lz,
        npix,
        nx,
        ny,
        nz,
        sigmaf,
        lperp,
        lpar,
        namepixel,
        namemap,
        nameinput,
    ):
        f = open(os.path.join(self.pwd, f"{nameinput}.cfg"), "w")
        f.write("#lx, ly, lz: the domain size in each direction.\n")
        f.write("#num_pixels: the *total* number of pixels.\n")
        f.write(
            "#map_nx, map_ny, map_nz: the number of map points. The map points are arbitrary but for now these n's are used to setup a uniform grid across the domain given above.\n"
        )
        f.write("#corr_var_s: the signal cov prefactor sigma_f^2\n")
        f.write("#corr_l_perp: the signal cov perp scale.\n")
        f.write("#corr_l_para: the signal cov para scale.\n")
        f.write(
            "#pcg_max_iter: the PCG max number of iterations. 100 should be good.\n"
        )
        f.write(
            "#pcg_tol: the PCG stopping tolerance. I found 1.0e-3 is good enough. Set it very small if you want the most accurate map.\n"
        )
        f.write("lx = {}\n".format(lx))
        f.write("ly = {}\n".format(ly))
        f.write("lz = {}\n".format(lz))
        f.write("\n")
        f.write("# From output of GEN_DACH_INPUT.PRO\n")
        f.write("num_pixels = {}\n".format(npix))
        f.write("\n")
        f.write("map_nx = {}\n".format(nx))
        f.write("map_ny = {}\n".format(ny))
        f.write("map_nz = {}\n".format(nz))
        f.write("\n")
        f.write("corr_var_s = {}\n".format(sigmaf))
        f.write("corr_l_perp = {}\n".format(lperp))
        f.write("corr_l_para = {}\n".format(lpar))
        f.write("\n")
        f.write("pcg_max_iter = 500\n")
        f.write("pcg_tol = 1.0e-3\n")
        f.write("#pcg_step_r = 1\n")
        f.write("\n")
        f.write("option_map_covar = 0\n")
        f.write("option_noise_covar = 0\n")
        f.write("pixel_data_path = {}\n".format(namepixel))
        f.write("map_path = {}\n".format(namemap))
        f.close()

    def create_dachshund_map_pixel_property_file(
        self,
        name_out,
        cartesian_coordinates,
        sky_coordinates,
        shape,
        properties_map_pixels,
    ):
        size = (
            np.max(cartesian_coordinates[:, 0]),
            np.max(cartesian_coordinates[:, 1]),
            np.max(cartesian_coordinates[:, 2]),
        )
        coordinate_transform = self.coordinate_transform
        boundary_cartesian_coord = (
            (
                properties_map_pixels["minx"],
                properties_map_pixels["miny"],
                properties_map_pixels["minz"],
            ),
            (
                properties_map_pixels["maxx"],
                properties_map_pixels["maxy"],
                properties_map_pixels["maxz"],
            ),
        )
        boundary_sky_coord = (
            (
                properties_map_pixels["minra"],
                properties_map_pixels["mindec"],
                properties_map_pixels["minredshift"],
            ),
            (
                properties_map_pixels["maxra"],
                properties_map_pixels["maxdec"],
                properties_map_pixels["maxredshift"],
            ),
        )
        property_file = tomographic_objects.MapPixelProperty(
            name=os.path.join(self.pwd, name_out),
            size=size,
            shape=shape,
            boundary_cartesian_coord=boundary_cartesian_coord,
            boundary_sky_coord=boundary_sky_coord,
            coordinate_transform=coordinate_transform,
            Omega_m=self.Omega_m,
        )
        return property_file

    def create_serial_input(self, nameout, properties, cartesian_deltas, sky_deltas):
        self.create_input_files(
            cartesian_deltas,
            properties,
            properties["name_pixel"],
            create_launcher=nameout,
        )
        if self.return_sky_catalogs:
            self.create_input_files(
                sky_deltas,
                properties,
                "{}_sky_coordinates".format(properties["name_pixel"]),
                create_launcher=None,
            )
        return properties["shape"]

    def cut_in_chunks(self, cartesian_deltas, number_chunks, overlaping, shape_sub_map):
        if overlaping is None:
            overlaping = 0.0
        minx, maxx = np.min(cartesian_deltas[:, 0]), np.max(cartesian_deltas[:, 0])
        miny, maxy = np.min(cartesian_deltas[:, 1]), np.max(cartesian_deltas[:, 1])
        minz, maxz = np.min(cartesian_deltas[:, 2]), np.max(cartesian_deltas[:, 2])
        intervalx = maxx - minx
        intervaly = maxy - miny
        intervalz = maxz - minz
        subIntervalx = intervalx / number_chunks[0]
        subIntervaly = intervaly / number_chunks[1]
        Chunks = {}
        shape_x = number_chunks[0] * shape_sub_map[0]
        shape_y = number_chunks[1] * shape_sub_map[1]
        remove_shape_x, remove_shape_y = 0, 0
        for i in range(number_chunks[0]):
            for j in range(number_chunks[1]):
                filename = f"{i:03d}" + f"{j:03d}"
                Chunks[filename] = {}
                if (i == number_chunks[0] - 1) & (i == 0):
                    intervalxChunk = [i * subIntervalx, (i + 1) * subIntervalx]
                elif i == 0:
                    intervalxChunk = [
                        i * subIntervalx,
                        (i + 1) * subIntervalx + overlaping,
                    ]
                elif i == number_chunks[0] - 1:
                    intervalxChunk = [i * subIntervalx - overlaping, intervalx]
                else:
                    intervalxChunk = [
                        i * subIntervalx - overlaping,
                        (i + 1) * subIntervalx + overlaping,
                    ]
                if (j == number_chunks[1] - 1) & (j == 0):
                    intervalyChunk = [j * subIntervaly, (j + 1) * subIntervaly]
                elif j == 0:
                    intervalyChunk = [
                        j * subIntervaly,
                        (j + 1) * subIntervaly + overlaping,
                    ]
                elif j == number_chunks[1] - 1:
                    intervalyChunk = [j * subIntervaly - overlaping, intervaly]
                else:
                    intervalyChunk = [
                        j * subIntervaly - overlaping,
                        (j + 1) * subIntervaly + overlaping,
                    ]
                mask = (cartesian_deltas[:, 0] < intervalxChunk[1]) & (
                    cartesian_deltas[:, 0] >= intervalxChunk[0]
                )
                mask &= (cartesian_deltas[:, 1] < intervalyChunk[1]) & (
                    cartesian_deltas[:, 1] >= intervalyChunk[0]
                )
                chunks_deltas = []
                chunks_deltas.append(cartesian_deltas[:, 0][mask] - intervalxChunk[0])
                chunks_deltas.append(cartesian_deltas[:, 1][mask] - intervalyChunk[0])
                chunks_deltas.append(cartesian_deltas[:, 2][mask])
                chunks_deltas.append(cartesian_deltas[:, 3][mask])
                chunks_deltas.append(cartesian_deltas[:, 4][mask])
                chunks_deltas = np.transpose(np.stack(chunks_deltas))
                Chunks[filename]["coord"] = chunks_deltas
                Chunks[filename]["limits"] = [
                    intervalxChunk[0],
                    intervalxChunk[1],
                    intervalyChunk[0],
                    intervalyChunk[1],
                    np.min(cartesian_deltas[:, 2]),
                    np.max(cartesian_deltas[:, 2]),
                ]
                size = (
                    intervalxChunk[1] - intervalxChunk[0],
                    intervalyChunk[1] - intervalyChunk[0],
                    intervalz,
                )
                pixel_to_remove = np.around(
                    utils.pixel_per_mpc(size, shape_sub_map) * overlaping, 0
                ).astype(int)
                if number_chunks[0] != 1:
                    if (i == 0) | (i == number_chunks[0] - 1):
                        remove_shape_x = remove_shape_x + pixel_to_remove[0]
                    else:
                        remove_shape_x = remove_shape_x + 2 * pixel_to_remove[0]
                if number_chunks[1] != 1:
                    if (j == 0) | (j == number_chunks[1] - 1):
                        remove_shape_y = remove_shape_y + pixel_to_remove[1]
                    else:
                        remove_shape_y = remove_shape_y + 2 * pixel_to_remove[1]
        shape_x = shape_x - remove_shape_x // number_chunks[1]
        shape_y = shape_y - remove_shape_y // number_chunks[0]
        Chunks["overlaping"] = overlaping
        return (Chunks, (shape_x, shape_y))

    def create_parallel_input(
        self, properties, cartesian_deltas, number_chunks, overlaping, shape_sub_map
    ):
        chunks, shape = self.cut_in_chunks(
            cartesian_deltas, number_chunks, overlaping, shape_sub_map
        )
        shape = (shape[0], shape[1], shape_sub_map[2])
        filename = []
        parallel_launcher_params = []
        for i in range(len(list(chunks.keys()))):
            key = list(chunks.keys())[i]
            if key != "overlaping":
                filename.append(key)
                parallel_launcher_params.append({})
                parallel_launcher_params[i]["maxx"] = chunks[key]["limits"][1]
                parallel_launcher_params[i]["minx"] = chunks[key]["limits"][0]
                parallel_launcher_params[i]["maxy"] = chunks[key]["limits"][3]
                parallel_launcher_params[i]["miny"] = chunks[key]["limits"][2]
                parallel_launcher_params[i]["maxz"] = chunks[key]["limits"][5]
                parallel_launcher_params[i]["minz"] = chunks[key]["limits"][4]
                parallel_launcher_params[i]["lx"] = (
                    chunks[key]["limits"][1] - chunks[key]["limits"][0]
                )
                parallel_launcher_params[i]["ly"] = (
                    chunks[key]["limits"][3] - chunks[key]["limits"][2]
                )
                parallel_launcher_params[i]["lz"] = (
                    chunks[key]["limits"][5] - chunks[key]["limits"][4]
                )
                parallel_launcher_params[i]["npix"] = len(chunks[key]["coord"])
                parallel_launcher_params[i]["nx"] = shape_sub_map[0]
                parallel_launcher_params[i]["ny"] = shape_sub_map[1]
                parallel_launcher_params[i]["nz"] = shape_sub_map[2]
                parallel_launcher_params[i]["sigmaf"] = properties["sigma_f"]
                parallel_launcher_params[i]["lperp"] = properties["lperp"]
                parallel_launcher_params[i]["lpar"] = properties["lpar"]
                parallel_launcher_params[i]["namepixel"] = "{}_{}".format(
                    properties["name_pixel"], key
                )
                parallel_launcher_params[i]["namemap"] = "map_{}_{}".format(
                    properties["name_pixel"], key
                )
                parallel_launcher_params[i]["nameinput"] = "input_{}.cfg".format(key)
        return (parallel_launcher_params, filename, chunks, shape)

    def write_parallel_input(
        self,
        cartesian_deltas,
        parallel_launcher_params,
        filename,
        chunks,
        properties,
        nameout,
        number_chunks,
        overlaping,
    ):
        self.create_input_files(
            cartesian_deltas, properties, properties["name_pixel"], create_launcher=None
        )
        for i in range(len(list(chunks.keys()))):
            key = list(chunks.keys())[i]
            if key != "overlaping":
                self.create_input_files(
                    chunks[key]["coord"],
                    properties,
                    f"{properties['name_pixel']}_{key}",
                    create_launcher=None,
                )
        pickle.dump(
            [filename, parallel_launcher_params, number_chunks, overlaping],
            open(os.path.join(self.pwd, f"{nameout}_launch_data.pickle"), "wb"),
        )

    def create_additional_catalogs(
        self,
        cartesian_qso_catalog,
        cartesian_dla_catalog,
        sky_qso_catalog,
        sky_dla_catalog,
        properties_map_pixels,
    ):
        boundary_cartesian_coord = (
            (
                properties_map_pixels["minx"],
                properties_map_pixels["miny"],
                properties_map_pixels["minz"],
            ),
            (
                properties_map_pixels["maxx"],
                properties_map_pixels["maxy"],
                properties_map_pixels["maxz"],
            ),
        )
        boundary_sky_coord = (
            (
                properties_map_pixels["minra"],
                properties_map_pixels["mindec"],
                properties_map_pixels["minredshift"],
            ),
            (
                properties_map_pixels["maxra"],
                properties_map_pixels["maxdec"],
                properties_map_pixels["maxredshift"],
            ),
        )
        if self.return_dla_catalog is not None:
            dla_catalog_cartesian = (
                tomographic_objects.DLACatalog.init_from_pixel_catalog(
                    cartesian_dla_catalog,
                    name=os.path.join(self.pwd, self.return_dla_catalog),
                    coordinate_transform=self.coordinate_transform,
                    Omega_m=self.Omega_m,
                    boundary_cartesian_coord=boundary_cartesian_coord,
                    boundary_sky_coord=boundary_sky_coord,
                )
            )
            dla_catalog_cartesian.write()
            if self.return_sky_catalogs:
                dla_catalog_sky = (
                    tomographic_objects.DLACatalog.init_from_pixel_catalog(
                        sky_dla_catalog,
                        name=os.path.join(
                            self.pwd, f"{self.return_dla_catalog}_sky_coordinates"
                        ),
                        coordinate_transform=self.coordinate_transform,
                        Omega_m=self.Omega_m,
                        boundary_cartesian_coord=boundary_cartesian_coord,
                        boundary_sky_coord=boundary_sky_coord,
                    )
                )
                dla_catalog_sky.write()
        if self.return_qso_catalog is not None:
            quasar_catalog_cartesian = (
                tomographic_objects.QSOCatalog.init_from_pixel_catalog(
                    cartesian_qso_catalog,
                    name=os.path.join(self.pwd, self.return_qso_catalog),
                    coordinate_transform=self.coordinate_transform,
                    Omega_m=self.Omega_m,
                    boundary_cartesian_coord=boundary_cartesian_coord,
                    boundary_sky_coord=boundary_sky_coord,
                )
            )
            quasar_catalog_cartesian.write()
            if self.return_sky_catalogs:
                quasar_catalog_sky = (
                    tomographic_objects.QSOCatalog.init_from_pixel_catalog(
                        sky_qso_catalog,
                        name=os.path.join(
                            self.pwd, f"{self.return_qso_catalog}_sky_coordinates"
                        ),
                        coordinate_transform=self.coordinate_transform,
                        Omega_m=self.Omega_m,
                        catalog_type="sky",
                        boundary_cartesian_coord=boundary_cartesian_coord,
                        boundary_sky_coord=boundary_sky_coord,
                    )
                )
                quasar_catalog_sky.write()

    def write_additional_catalogs(
        self,
        dla_catalog_sky,
        dla_catalog_cartesian,
        quasar_catalog_sky,
        quasar_catalog_cartesian,
    ):
        if dla_catalog_sky is not None:
            dla_catalog_sky.write()
        if dla_catalog_cartesian is not None:
            dla_catalog_cartesian.write()
        if quasar_catalog_sky is not None:
            quasar_catalog_sky.write()
        if quasar_catalog_cartesian is not None:
            quasar_catalog_cartesian.write()

    # CR - Modification : dissociate map properties and launcher properties

    def transform_delta(
        self,
        mode,
        nameout,
        properties,
        property_file_name,
        rebin=False,
        sigma_min=None,
        sigma_max=None,
        z_cut_min=None,
        z_cut_max=None,
        dec_cut_min=None,
        dec_cut_max=None,
        ra_cut_min=None,
        ra_cut_max=None,
        number_chunks=None,
        overlaping=None,
        shape_sub_map=None,
    ):
        (
            cartesian_deltas,
            cartesian_qso_catalog,
            cartesian_dla_catalog,
            sky_deltas,
            sky_qso_catalog,
            sky_dla_catalog,
            properties_map_pixels,
        ) = self.transform_delta_to_pixel_file(
            rebin=rebin,
            sigma_min=sigma_min,
            sigma_max=sigma_max,
            z_cut_min=z_cut_min,
            z_cut_max=z_cut_max,
            dec_cut_min=dec_cut_min,
            dec_cut_max=dec_cut_max,
            ra_cut_min=ra_cut_min,
            ra_cut_max=ra_cut_max,
        )
        if mode.lower() == "serial":
            shape = self.create_serial_input(
                nameout, properties, cartesian_deltas, sky_deltas
            )
        elif mode.lower() == "parallel":
            (
                parallel_launcher_params,
                filename,
                chunks,
                shape,
            ) = self.create_parallel_input(
                properties, cartesian_deltas, number_chunks, overlaping, shape_sub_map
            )
            self.write_parallel_input(
                cartesian_deltas,
                parallel_launcher_params,
                filename,
                chunks,
                properties,
                nameout,
                number_chunks,
                overlaping,
            )
        else:
            raise KeyError("Please choose a mode between serial and parallel")
        property_file = self.create_dachshund_map_pixel_property_file(
            property_file_name,
            cartesian_deltas,
            sky_deltas,
            shape,
            properties_map_pixels,
        )
        property_file.write()
        self.create_additional_catalogs(
            cartesian_qso_catalog,
            cartesian_dla_catalog,
            sky_qso_catalog,
            sky_dla_catalog,
            properties_map_pixels,
        )
        if self.plot_pixel_properties:
            pixel_analyzer = PixelAnalizer(
                pixel=os.path.join(self.pwd, properties["name_pixel"]),
                property_file=os.path.join(self.pwd, property_file_name),
            )
            pixel_analyzer.analyze_pixels(
                False,
                True,
                name_dperp=os.path.join(self.pwd, nameout),
                coupled_plot=True,
            )


class PixelAnalizer(object):
    def __init__(self, pixel=None, property_file=None):
        if type(pixel) == str:
            pixel_class = tomographic_objects.Pixel.init_from_property_files(
                property_file, name=pixel
            )
            pixel_class.read()
        else:
            pixel_class = pixel
        self.pixel = pixel_class

    @staticmethod
    def write_density_file(z, dperp, name):
        pickle.dump([z, dperp], open(name, "wb"))

    @staticmethod
    def write_dperp_file(z, dperp, name):
        pickle.dump([z, dperp], open(name, "wb"))

    @staticmethod
    def read_density_file(name):
        a = pickle.load(open(name, "rb"))
        return (a[0], a[1])

    @staticmethod
    def read_dperp_file(name):
        a = pickle.load(open(name, "rb"))
        return (a[0], a[1])

    @staticmethod
    def plot_histogram_mean_distance(zpar, dmin, name_histo, nb_bins=50):
        plt.figure()
        plt.hist(dmin, nb_bins)
        plt.xlabel("minimal distance histogram at Z={}".format(zpar))
        plt.savefig(f"{name_histo}_at_Z{zpar}.pdf", format="pdf")

    @staticmethod
    def plot_mean_distance_density(
        zpar,
        dperpz,
        densityz,
        nameout,
        coupled_plot=False,
        comparison=False,
        dperp_comparison=None,
        density_comparison=None,
        zpar_comparison=None,
        legend=None,
        dperp_other=None,
        density_other=None,
        **kwargs,
    ):
        if coupled_plot:
            PixelAnalizer.plot_mean_distance_density_coupled(
                zpar,
                dperpz,
                densityz,
                nameout,
                comparison=comparison,
                dperp_comparison=dperp_comparison,
                density_comparison=density_comparison,
                zpar_comparison=zpar_comparison,
                legend=legend,
                dperp_other=dperp_other,
                density_other=density_other,
                **kwargs,
            )
        else:
            PixelAnalizer.plot_mean_distance_density_not_coupled(
                zpar,
                dperpz,
                densityz,
                nameout,
                comparison=comparison,
                dperp_comparison=dperp_comparison,
                density_comparison=density_comparison,
                zpar_comparison=zpar_comparison,
                legend=legend,
                dperp_other=dperp_other,
                density_other=density_other,
                **kwargs,
            )

    @staticmethod
    def plot_mean_distance_density_not_coupled(
        zpar,
        dperpz,
        densityz,
        nameout,
        comparison=False,
        dperp_comparison=None,
        density_comparison=None,
        zpar_comparison=None,
        legend=None,
        dperp_other=None,
        density_other=None,
        **kwargs,
    ):
        plt.figure()
        plt.xlabel("Redshift")
        plt.ylabel("Mean los separation [" + r"$\mathrm{Mpc\cdot h^{-1}}$" + "]")
        plt.grid()
        plt.plot(zpar, dperpz)
        if comparison:
            plt.plot(zpar_comparison, dperp_comparison)
            plt.legend(legend)
            if dperp_other is not None:
                for i in range(len(dperp_other)):
                    plt.plot(dperp_other[i][0], dperp_other[i][1])
        plt.savefig(f"{nameout}_separation.pdf", format="pdf")

        plt.figure()
        plt.xlabel("Redshift")
        plt.ylabel("Density [" + r"$\mathrm{deg^{-2}}$" + "]")
        plt.grid()
        plt.plot(zpar, densityz)
        if comparison:
            plt.plot(zpar_comparison, density_comparison)
            plt.legend(legend)
            if density_other is not None:
                for i in range(len(density_other)):
                    plt.plot(density_other[i][0], density_other[i][1])
        plt.savefig(f"{nameout}_density.pdf", format="pdf")

    @staticmethod
    def plot_mean_distance_density_coupled(
        zpar,
        dperpz,
        densityz,
        nameout,
        comparison=False,
        dperp_comparison=None,
        density_comparison=None,
        zpar_comparison=None,
        legend=None,
        dperp_other=None,
        density_other=None,
        **kwargs,
    ):
        figsize = utils.return_key(kwargs, "figsize", (8, 6))
        grid = utils.return_key(kwargs, "grid", True)
        fontsize = utils.return_key(kwargs, "fontsize", 13)
        fontsize_scale = utils.return_key(kwargs, "fontscalesize", 13)
        ylabel1 = utils.return_key(
            kwargs,
            "ylabel1",
            "Mean separation between nearest los ["
            + r"$h^{-1}$"
            + r"$\cdot$"
            + "Mpc"
            + "]",
        )
        ylabel2 = utils.return_key(
            kwargs, "ylabel2", "Density [" + r"$\mathrm{deg}^{-2}$" + "]"
        )
        xlabel = utils.return_key(kwargs, "xlabel", r"Redshift $z$")

        fig, ax1 = plt.subplots(1, 1, figsize=figsize)
        if grid:
            ax1.grid()
        line = ["dotted", "dashdot", "densely dashdotdotted"]

        color = "C0"
        ax1.set_xlabel(xlabel, fontsize=fontsize)
        ax1.set_ylabel(ylabel1, color=color, fontsize=fontsize)
        ax1.tick_params(axis="x", labelsize=fontsize_scale)
        ax1.tick_params(axis="y", labelsize=fontsize_scale)
        ax1.plot(zpar, dperpz, color=color)
        if comparison:
            ax1.plot(zpar_comparison, dperp_comparison, color=color, linestyle="--")
            if dperp_other is not None:
                for i in range(len(dperp_other)):
                    ax1.plot(
                        dperp_other[i][0],
                        dperp_other[i][1],
                        color=color,
                        linestyle=line[i],
                    )
        ax1.tick_params(axis="y", labelcolor=color)

        ax2 = ax1.twinx()
        if grid:
            ax2.grid()
        color = "C1"
        ax2.set_ylabel(ylabel2, color=color, fontsize=fontsize)
        ax2.tick_params(axis="y", labelsize=fontsize_scale)
        ax2.plot(zpar, densityz, color=color)
        if comparison:
            ax2.plot(zpar_comparison, density_comparison, color=color, linestyle="--")
            if density_other is not None:
                for i in range(len(density_other)):
                    ax2.plot(
                        density_other[i][0],
                        density_other[i][1],
                        color=color,
                        linestyle=line[i],
                    )
        ax2.tick_params(axis="y", labelcolor=color)

        fig.tight_layout()
        if comparison:
            legend_elements = [
                Line2D([0], [0], color="k", lw=1, label=legend[0]),
                Line2D([0], [0], color="k", linestyle="--", lw=1, label=legend[1]),
            ]
            if dperp_other is not None:
                for i in range(len(dperp_other)):
                    legend_elements.append(
                        Line2D(
                            [0],
                            [0],
                            color="k",
                            linestyle=line[i],
                            lw=1,
                            label=legend[i + 2],
                        )
                    )

            ax1.legend(handles=legend_elements, loc="upper center")
        plt.savefig(f"{nameout}_separation_density.pdf", format="pdf")

    @staticmethod
    def plot_histo_mean_distance_comparison(
        dpername1,
        dpername2,
        densityname1,
        densityname2,
        nameout,
        legend,
        coupled_plot=False,
        **kwargs,
    ):
        if type(dpername1) is str:
            zpar, dperp = PixelAnalizer.read_dperp_file(dpername1)
            zpar_comparison, dperp_comparison = PixelAnalizer.read_dperp_file(dpername2)
            zpar, density = PixelAnalizer.read_density_file(densityname1)
            zpar_comparison, density_comparison = PixelAnalizer.read_density_file(
                densityname2
            )
        else:
            zpar, dperp = np.mean(
                [
                    PixelAnalizer.read_dperp_file(dpername1[i])
                    for i in range(len(dpername1))
                ],
                axis=0,
            )
            zpar_comparison, dperp_comparison = np.mean(
                [
                    PixelAnalizer.read_dperp_file(dpername2[i])
                    for i in range(len(dpername2))
                ],
                axis=0,
            )
            zpar, density = np.mean(
                [
                    PixelAnalizer.read_density_file(densityname1[i])
                    for i in range(len(densityname1))
                ],
                axis=0,
            )
            zpar_comparison, density_comparison = np.mean(
                [
                    PixelAnalizer.read_density_file(densityname2[i])
                    for i in range(len(densityname2))
                ],
                axis=0,
            )
        PixelAnalizer.plot_mean_distance_density(
            zpar,
            dperp,
            density,
            nameout,
            coupled_plot=coupled_plot,
            comparison=True,
            dperp_comparison=dperp_comparison,
            density_comparison=density_comparison,
            zpar_comparison=zpar_comparison,
            legend=legend,
            **kwargs,
        )

    def compute_plot_histo_mean_distance(self, zpar, name_histo):
        dmin = self.pixel.compute_mean_distance_histogram(zpar)
        PixelAnalizer.plot_histogram_mean_distance(zpar, dmin, name_histo, nb_bins=50)

    def compute_plot_mean_distance_density(self, nameout, coupled=False, plot=True):
        (zpar, dperpz, densityz) = self.pixel.compute_mean_distance_density()
        if plot:
            PixelAnalizer.write_dperp_file(zpar, dperpz, f"{nameout}_dperp_file.pickle")
            PixelAnalizer.write_density_file(
                zpar, densityz, f"{nameout}_density_file.pickle"
            )
            PixelAnalizer.plot_mean_distance_density(
                zpar, dperpz, densityz, nameout, coupled_plot=coupled
            )
        return (zpar, dperpz, densityz)

    def analyze_pixels(
        self,
        compute_histo,
        compute_mean_distance_density,
        histo_zpar=None,
        name_histo="histogram_dperp",
        name_dperp="density_mean_distance",
        coupled_plot=False,
    ):
        if compute_histo:
            self.compute_plot_histo_mean_distance(histo_zpar, name_histo)
        if compute_mean_distance_density:
            self.compute_plot_mean_distance_density(name_dperp, coupled=coupled_plot)


class DeltaAnalyzer(object):
    def __init__(
        self,
        pwd,
        delta_path,
        center_ra=True,
        z_cut_min=None,
        z_cut_max=None,
        dec_cut_min=None,
        dec_cut_max=None,
        ra_cut_min=None,
        ra_cut_max=None,
        degree=True,
        pk1d_type=True,
    ):
        self.pwd = pwd
        self.delta_path = delta_path
        self.center_ra = center_ra
        self.z_cut_min = z_cut_min
        self.z_cut_max = z_cut_max
        self.dec_cut_min = dec_cut_min
        self.dec_cut_max = dec_cut_max
        self.ra_cut_min = ra_cut_min
        self.ra_cut_max = ra_cut_max
        self.degree = degree
        self.pk1d_type = pk1d_type

    def get_ra_dec(self, delta_path):
        """Obtain arrays of RA and DEC coordinates from a list or a name of a delta file in pickle, fits or ascii format"""
        namefile = get_delta_list(delta_path)
        # namefile = preselect_deltas(namefile,
        #                             ramin=self.ra_cut_min,
        #                             ramax=self.ra_cut_max,
        #                             decmin=self.dec_cut_min,
        #                             decmax=self.dec_cut_max)
        (ra, dec, z, zqso, ids, sigmas, deltas) = get_deltas(
            namefile, center_ra=self.center_ra, pk1d_type=self.pk1d_type
        )
        pixel_coord = np.array(
            [
                [ra[i], dec[i], z[i][j], sigmas[i][j], deltas[i][j], zqso[i], ids[i]]
                for i in range(len(ra))
                for j in range(len(z[i]))
            ]
        )
        pixel_coord = pixel_coord[
            utils.cut_sky_catalog(
                pixel_coord[:, 0],
                pixel_coord[:, 1],
                pixel_coord[:, 2],
                ramin=self.ra_cut_min,
                ramax=self.ra_cut_max,
                decmin=self.dec_cut_min,
                decmax=self.dec_cut_max,
                zmin=self.z_cut_min,
                zmax=self.z_cut_max,
            )
        ]
        (redshift, redshift_qso, id, sigma, delta) = (
            pixel_coord[:, 2],
            pixel_coord[:, 5],
            pixel_coord[:, 6],
            pixel_coord[:, 3],
            pixel_coord[:, 4],
        )
        unique_coord = np.unique(pixel_coord[:, 0:2], axis=0)
        ra = unique_coord[:, 0]
        dec = unique_coord[:, 1]
        if self.degree:
            ra, dec = np.degrees(ra), np.degrees(dec)
        snr = np.abs((delta + 1) / sigma)
        return (ra, dec, redshift, redshift_qso, id, sigma, delta, snr)

    def get_ra_dec_comparison(self, comparison):
        if comparison is None:
            return (None, None, None, None, None, None, None, None)
        (
            ra_comp,
            dec_comp,
            redshift_comp,
            redshift_qso_comp,
            id_comp,
            sigma_comp,
            delta_comp,
            snr_comp,
        ) = ([], [], [], [], [], [], [], [])
        for i in range(len(comparison)):
            (ra, dec, redshift, redshift_qso, id, sigma, delta, snr) = self.get_ra_dec(
                comparison[i]
            )
            ra_comp.append(ra)
            dec_comp.append(dec)
            redshift_comp.append(redshift)
            redshift_qso_comp.append(redshift_qso)
            id_comp.append(id)
            sigma_comp.append(sigma)
            delta_comp.append(delta)
            snr_comp.append(snr)
        return (
            ra_comp,
            dec_comp,
            redshift_comp,
            redshift_qso_comp,
            id_comp,
            sigma_comp,
            delta_comp,
            snr_comp,
        )

    def plot(
        self,
        value_names,
        name,
        comparison=None,
        comparison_legend=None,
        histo=True,
        mean_z_dependence=True,
        z_dependence=True,
        ra_dec_plots=True,
        print_stats=False,
        **kwargs,
    ):
        style = utils.return_key(kwargs, "style", None)
        if style is not None:
            plt.style.use(style)

        (ra, dec, redshift, redshift_qso, id, sigma, delta, snr) = self.get_ra_dec(
            self.delta_path
        )
        (
            ra_comp,
            dec_comp,
            redshift_comp,
            redshift_qso_comp,
            id_comp,
            sigma_comp,
            delta_comp,
            snr_comp,
        ) = self.get_ra_dec_comparison(comparison)

        for value_name in value_names:
            value = locals()[value_name]
            comparison_value = locals()[value_name + "_comp"]
            lambda_rest = utils.return_key(kwargs, f"{value_name}_lambda_rest", False)
            if lambda_rest:
                kwargs[f"{value_name}_redshift_qso"] = redshift_qso
            if histo:
                utils.save_histo(
                    self.pwd,
                    value,
                    value_name,
                    name,
                    comparison=comparison_value,
                    comparison_legend=comparison_legend,
                    **kwargs,
                )
            if (mean_z_dependence) & (value_name not in ["redshift", "ra", "dec"]):
                utils.save_mean_redshift_dependence(
                    self.pwd,
                    value,
                    redshift,
                    value_name,
                    name,
                    comparison=comparison_value,
                    comparison_redshift=redshift_comp,
                    comparison_legend=None,
                    **kwargs,
                )

            if (z_dependence) & (value_name not in ["redshift", "ra", "dec"]):
                utils.save_redshift_dependence(
                    self.pwd,
                    value,
                    redshift,
                    value_name,
                    name,
                    comparison=comparison_value,
                    comparison_redshift=redshift_comp,
                    comparison_legend=None,
                    **kwargs,
                )
        if ra_dec_plots:
            self.plot_ra_dec(
                ra, dec, name, comparison=None, comparison_legend=None, **kwargs
            )
        if (comparison is not None) & (print_stats):
            self.stat_comparison(
                redshift, redshift_comp, sigma, sigma_comp, delta, delta_comp
            )

    def plot_ra_dec(
        self,
        ra,
        dec,
        name,
        comparison_ra=None,
        comparison_dec=None,
        comparison_legend=None,
        **kwargs,
    ):
        utils.save_ra_dec(
            self.pwd,
            ra,
            dec,
            name,
            comparison_ra=comparison_ra,
            comparison_dec=comparison_dec,
            comparison_legend=comparison_legend,
            **kwargs,
        )

        DeltaAnalyzer.plot_los_density(self.pwd, ra, dec, name, **kwargs)

    def stat_comparison(
        self, redshift, redshift_comp, sigma, sigma_comp, delta, delta_comp
    ):
        print("Redshift interval =", np.max(redshift), np.min(redshift))
        for i in range(len(redshift_comp)):
            print(
                f"Redshift interval for comp {i} =",
                np.max(redshift_comp[i]),
                np.min(redshift_comp[i]),
            )
        print("Maximal sigma =", np.max(sigma))
        for i in range(len(sigma_comp)):
            print(f"Maximal sigma for comp {i} =", np.max(sigma_comp[i]))
        print("Mean sigma =", np.mean(sigma))
        for i in range(len(sigma_comp)):
            print(f"Mean sigma for comp {i} =", np.mean(sigma_comp[i]))
        print("Median sigma =", np.median(sigma))
        for i in range(len(sigma_comp)):
            print(f"Median sigma for comp {i} =", np.median(sigma_comp[i]))
        print("Mean delta =", np.mean(delta))
        for i in range(len(delta_comp)):
            print(f"Mean delta for comp {i} =", np.mean(delta_comp[i]))

    @staticmethod
    def plot_los_density(pwd, ra, dec, plot_name, **kwargs):
        nb_interval = utils.return_key(kwargs, "nb_interval", 20)
        different_sign_region = utils.return_key(kwargs, "different_sign_region", False)

        ra_interval = np.linspace(np.min(ra), np.max(ra), nb_interval)
        ra_size = abs((np.max(ra) - np.min(ra)) / nb_interval)
        ra_array = []
        for i in range(len(ra_interval) - 1):
            ra_array.append((ra_interval[i] + ra_interval[i + 1]) / 2)
        ra_array = np.asarray(ra_array)
        density_array = np.zeros(ra_array.shape)
        density_array_plus = np.zeros(ra_array.shape)
        density_array_minus = np.zeros(ra_array.shape)
        for i in range(len(ra_interval) - 1):
            mask = (ra > ra_interval[i]) & (ra < ra_interval[i + 1])
            density_array[i] = len(ra[mask])
            density_array_plus[i] = len(ra[mask & (dec >= 0)])
            density_array_minus[i] = len(ra[mask & (dec < 0)])
        maxdec = np.max(dec)
        mindec = np.min(dec)
        plt.figure()
        plt.plot(ra_array, density_array / abs(ra_size * (maxdec - mindec)))
        plt.title("LOS density in function of RA")
        plt.grid()
        plt.savefig(os.path.join(pwd, f"{plot_name}_los_density.pdf"), format="pdf")
        if different_sign_region:
            plt.figure()
            plt.plot(ra_array, density_array_plus / abs(ra_size * maxdec))
            plt.title("LOS density in function of RA for DEC >= 0")
            plt.grid()
            plt.savefig(
                os.path.join(pwd, f"{plot_name}_los_density_dec_positive.pdf"),
                format="pdf",
            )
            plt.figure()
            plt.plot(ra_array, density_array_minus / abs(ra_size * mindec))
            plt.title("LOS density in function of RA for DEC < 0")
            plt.grid()
            plt.savefig(
                os.path.join(pwd, f"{plot_name}_los_density_dec_negative.pdf"),
                format="pdf",
            )

    @staticmethod
    def plot_delta_gaussian_fit(pwd, delta, plot_name, **kwargs):
        (name, data, bins, patches) = utils.plot_histo(delta, "delta", "", **kwargs)
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
        plt.plot(x, y, "r--", linewidth=2)
        mu, sigma = popt[1], popt[2]
        plt.text(0.8, 2 * np.max(data) / 6, "mu = " + str(round(mu, 8)))
        plt.text(0.8, 1.5 * np.max(data) / 6, "sigma = " + str(round(sigma, 8)))

        plt.savefig(
            os.path.join(pwd, f"{plot_name}_histo_delta_gaussian_fit.pdf"), format="pdf"
        )
