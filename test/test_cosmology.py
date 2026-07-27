"""Tests for :mod:`lelantos.cosmology` delta reading and conversion.

These exercise the delta->pixel conversion against the committed reference
delta file (``test/data/delta/delta-19.fits.gz``), checking both the low-level
readers and the invariants of :class:`DeltaConverter`.
"""

import numpy as np
import pytest

from lelantos import cosmology, tomographic_objects

from conftest import DELTA_CUTS, DELTA_DIR, OMEGA_M


# --- Low-level delta readers ----------------------------------------------


def test_get_delta_list_finds_reference_file():
    files = cosmology.get_delta_list(str(DELTA_DIR))
    assert len(files) == 1
    assert files[0].endswith(".fits.gz")


def test_get_deltas_shapes_and_ranges():
    files = cosmology.get_delta_list(str(DELTA_DIR))
    ra, dec, z, zqso, ids, sigmas, deltas = cosmology.get_deltas(files, center_ra=True)
    n = len(ra)
    assert n == 154
    # per-line-of-sight arrays are all the same length as the LOS list
    for arr in (dec, z, zqso, ids, sigmas, deltas):
        assert len(arr) == n
    # each LOS has matching redshift / delta / sigma pixel counts
    for i in range(n):
        assert len(z[i]) == len(deltas[i]) == len(sigmas[i])
    # ids are unique quasar identifiers
    assert len(np.unique(ids)) == n


# --- DeltaConverter --------------------------------------------------------


def test_transform_delta_invariants(converted_deltas):
    cart, cart_qso, cart_dla, sky, sky_qso, sky_dla, props = converted_deltas
    # cartesian pixels are (N, 5): x, y, z, sigma, delta
    assert cart.ndim == 2 and cart.shape[1] == 5
    assert cart.shape[0] > 0
    # box is shifted so the minimum cartesian corner sits at the origin
    np.testing.assert_allclose(cart[:, 0:3].min(axis=0), [0.0, 0.0, 0.0], atol=1e-9)
    # sigma has been clipped from below at sigma_min
    assert cart[:, 3].min() >= DELTA_CUTS["sigma_min"] - 1e-12
    # the property dict exposes the cartesian and sky bounding boxes
    for key in ("minx", "maxx", "minra", "maxredshift"):
        assert key in props
    # redshift window respected
    assert props["minredshift"] >= DELTA_CUTS["z_cut_min"]
    assert props["maxredshift"] <= DELTA_CUTS["z_cut_max"]


def test_transform_qso_catalog_consistency(converted_deltas):
    cart, cart_qso, cart_dla, sky, sky_qso, sky_dla, props = converted_deltas
    assert cart_qso is not None and sky_qso is not None
    # same number of quasars in sky and cartesian views, sharing primary keys
    assert cart_qso.shape == sky_qso.shape
    np.testing.assert_array_equal(cart_qso[:, 3], sky_qso[:, 3])
    # cartesian quasars share the box origin shift with the pixels
    assert cart_qso[:, 0:3].min() >= -1e-9


def test_transform_is_deterministic(delta_converter):
    a = delta_converter.transform_delta_to_pixel_file(**DELTA_CUTS)[0]
    b = delta_converter.transform_delta_to_pixel_file(**DELTA_CUTS)[0]
    np.testing.assert_array_equal(a, b)


def test_empty_footprint_raises(delta_converter):
    # a footprint that excludes every line of sight leaves nothing to reduce
    with pytest.raises(ValueError):
        delta_converter.transform_delta_to_pixel_file(
            sigma_min=0.1, z_cut_min=2.1, z_cut_max=3.2,
            dec_cut_min=51.0, dec_cut_max=55.0,
            ra_cut_min=10.0, ra_cut_max=20.0,  # data sits near -148 deg
        )


# --- Serial solver input ---------------------------------------------------


def test_create_serial_input_writes_files(workdir, delta_converter, converted_deltas):
    cart, cart_qso, cart_dla, sky, sky_qso, sky_dla, props = converted_deltas
    properties = {
        "sigma_f": 0.1,
        "lperp": 13,
        "lpar": 13,
        "shape": (10, 10, 50),
        "name_pixel": "pixel.bin",
        "name_map": "map.bin",
    }
    shape = delta_converter.create_serial_input("launch", properties, cart, sky)
    assert shape == (10, 10, 50)
    assert (workdir / "launch.cfg").is_file()
    assert (workdir / "pixel.bin").is_file()
    # the written pixel binary round-trips to the in-memory cartesian pixels
    pixel = tomographic_objects.Pixel(name=str(workdir / "pixel.bin"))
    pixel.read()
    np.testing.assert_allclose(pixel.pixel_array, cart)


def test_create_map_pixel_property_file(workdir, delta_converter, converted_deltas):
    cart, cart_qso, cart_dla, sky, sky_qso, sky_dla, props = converted_deltas
    prop = delta_converter.create_dachshund_map_pixel_property_file(
        "property_file.pickle", cart, sky, (10, 10, 50), props
    )
    prop.write()
    read = tomographic_objects.MapPixelProperty(name=str(workdir / "property_file.pickle"))
    read.read()
    assert read.shape == (10, 10, 50)
    assert read.Omega_m == OMEGA_M
    assert read.coordinate_transform == "middle"
