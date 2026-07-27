"""Unit tests for :mod:`lelantos.utils` (geometry and coordinate helpers)."""

import numpy as np
import pytest

from lelantos import utils


# --- Pixel / size geometry -------------------------------------------------


def test_mpc_pixel_roundtrip():
    size = (100.0, 200.0, 300.0)
    shape = (11, 21, 31)
    mpp = utils.mpc_per_pixel(size, shape)
    ppm = utils.pixel_per_mpc(size, shape)
    np.testing.assert_allclose(mpp * ppm, 1.0)
    np.testing.assert_allclose(mpp, np.array(size) / (np.array(shape) - 1))


def test_get_map_shape_size_inverse():
    size = (100.0, 200.0, 300.0)
    mpp = np.array([2.0, 4.0, 6.0])
    shape = utils.get_map_shape(size, mpp)
    np.testing.assert_allclose(utils.get_map_size(shape, mpp), size)


# --- Coordinate transforms -------------------------------------------------


@pytest.fixture(scope="module")
def cosmo_functions():
    return utils.get_cosmo_function(0.3147)


def test_get_cosmo_function_monotonic(cosmo_functions):
    rcomov, distang, inv_rcomov, inv_distang = cosmo_functions
    z = np.linspace(0.1, 3.0, 20)
    r = rcomov(z)
    # comoving distance is strictly increasing with redshift
    assert np.all(np.diff(r) > 0)
    # inverse maps distance back to redshift
    np.testing.assert_allclose(inv_rcomov(r), z, rtol=1e-4)


@pytest.mark.parametrize("method", ["full_angle", "full", "middle"])
def test_sky_cartesian_roundtrip(method, cosmo_functions):
    rcomov, distang, inv_rcomov, inv_distang = cosmo_functions
    ra = np.array([0.01, -0.02, 0.03])
    dec = np.array([0.02, 0.01, -0.015])
    z = np.array([2.2, 2.5, 2.8])
    supp = utils.return_suplementary_parameters(method, zmin=2.0, zmax=3.0)
    X, Y, Z = utils.convert_sky_to_cartesian(
        ra, dec, z, method, rcomov=rcomov, distang=distang,
        suplementary_parameters=supp,
    )
    ra2, dec2, z2 = utils.convert_cartesian_to_sky(
        X, Y, Z, method, inv_rcomov=inv_rcomov, distang=distang,
        suplementary_parameters=supp,
    )
    np.testing.assert_allclose(ra2, ra, atol=1e-6)
    np.testing.assert_allclose(dec2, dec, atol=1e-6)
    # the "full" (small-angle) transform carries a ~0.1% redshift approximation
    # because the transverse extent leaks into the radial distance.
    np.testing.assert_allclose(z2, z, rtol=5e-3)


def test_return_suplementary_parameters_middle():
    supp = utils.return_suplementary_parameters("middle", zmin=2.0, zmax=3.0)
    assert supp == [2.5]
    assert utils.return_suplementary_parameters("full") is None


def test_cut_sky_catalog_degrees():
    # bounds are in degrees, coordinates in radians
    ra = np.radians(np.array([10.0, 20.0, 30.0]))
    dec = np.radians(np.array([-5.0, 0.0, 5.0]))
    z = np.array([2.0, 2.5, 3.0])
    mask = utils.cut_sky_catalog(
        ra, dec, z, ramin=15.0, ramax=35.0, decmin=-1.0, decmax=6.0,
        zmin=2.2, zmax=3.5,
    )
    np.testing.assert_array_equal(mask, [False, True, True])


# --- Array helpers ---------------------------------------------------------


def test_bin_ndarray_sum_matches_docstring():
    m = np.arange(0, 100, 1).reshape((10, 10))
    n = utils.bin_ndarray(m, new_shape=(5, 5), operation="sum")
    expected = np.array(
        [
            [22, 30, 38, 46, 54],
            [102, 110, 118, 126, 134],
            [182, 190, 198, 206, 214],
            [262, 270, 278, 286, 294],
            [342, 350, 358, 366, 374],
        ]
    )
    np.testing.assert_array_equal(n, expected)


def test_bin_ndarray_mean_and_errors():
    m = np.ones((4, 4))
    np.testing.assert_allclose(utils.bin_ndarray(m, (2, 2), operation="mean"), 1.0)
    with pytest.raises(ValueError):
        utils.bin_ndarray(m, (2, 2), operation="bogus")
    with pytest.raises(ValueError):
        utils.bin_ndarray(m, (2, 2, 2), operation="mean")


def test_gaussian_smoothing_preserves_shape_and_mean():
    rng = np.random.default_rng(0)
    data = rng.normal(size=(8, 8, 8))
    smoothed = utils.gaussian_smoothing(data, sigma=1.0)
    assert smoothed.shape == data.shape
    # a Gaussian filter conserves the array mean (reflect boundary)
    np.testing.assert_allclose(smoothed.mean(), data.mean(), rtol=1e-6)


def test_interpolate_map_nearest():
    arr = np.arange(27, dtype=float).reshape((3, 3, 3))
    coord = np.array([[[[0, 0, 0], [1, 1, 1]]]], dtype=float)
    out = utils.interpolate_map("NEAREST", arr, coord)
    np.testing.assert_array_equal(out.ravel(), [arr[0, 0, 0], arr[1, 1, 1]])
    with pytest.raises(ValueError):
        utils.interpolate_map("BOGUS", arr, coord)


# --- Misc helpers ----------------------------------------------------------


def test_return_key():
    d = {"a": 1}
    assert utils.return_key(d, "a", 99) == 1
    assert utils.return_key(d, "b", 99) == 99


def test_latex_float():
    assert utils.latex_float(100.0) == "1 \\times 10^{2}"
    assert "times" not in utils.latex_float(1.5)


def test_get_direction_indexes():
    assert utils.get_direction_indexes("x", False)[:3] == (2, 1, 0)
    assert utils.get_direction_indexes("z", False)[:3] == (1, 0, 2)
    # rotate swaps the in-plane axes
    x_index, y_index, _, _ = utils.get_direction_indexes("x", True)
    assert (x_index, y_index) == (1, 2)
