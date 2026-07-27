"""Unit tests for :mod:`lelantos.tomographic_objects` (I/O round-trips)."""

import numpy as np
import pytest

from lelantos import tomographic_objects

from conftest import OMEGA_M, PIXEL_DIR, make_synthetic_property


# --- MapPixelProperty ------------------------------------------------------


def test_map_pixel_property_roundtrip(workdir):
    name = workdir / "prop.pickle"
    prop = tomographic_objects.MapPixelProperty(
        name=str(name),
        size=(80.0, 80.0, 120.0),
        shape=(40, 40, 60),
        boundary_cartesian_coord=((0.0, 0.0, 0.0), (80.0, 80.0, 120.0)),
        boundary_sky_coord=((-0.01, -0.01, 2.1), (0.01, 0.01, 2.3)),
        coordinate_transform="middle",
        Omega_m=OMEGA_M,
    )
    prop.write()
    read = tomographic_objects.MapPixelProperty(name=str(name))
    read.read()
    assert read.shape == prop.shape
    assert read.size == prop.size
    assert read.coordinate_transform == "middle"
    assert read.Omega_m == OMEGA_M
    np.testing.assert_allclose(read.mpc_per_pixel, prop.mpc_per_pixel)


def test_init_false_prop_middle(workdir):
    prop = make_synthetic_property(workdir / "prop.pickle")
    assert prop.coordinate_transform == "middle"
    # cartesian box origin is the requested corner (0, 0, minz)
    (minx, miny, minz), (maxx, maxy, maxz) = prop.boundary_cartesian_coord
    assert (minx, miny) == (0.0, 0.0)
    np.testing.assert_allclose((maxx - minx, maxy - miny, maxz - minz), prop.size)
    # sky redshift bounds bracket the input min redshift 2.2
    assert prop.boundary_sky_coord[0][2] == 2.2
    assert prop.boundary_sky_coord[1][2] > 2.2


def test_init_false_prop_rejects_non_middle(workdir):
    with pytest.raises(NotImplementedError):
        tomographic_objects.MapPixelProperty.init_false_prop(
            (10, 10, 10), (10.0, 10.0, 10.0), OMEGA_M, 0.0, 0.0, 2.2, "full"
        )


# --- TomographicMap --------------------------------------------------------


def test_tomographic_map_roundtrip(workdir):
    shape = (5, 6, 7)
    rng = np.random.default_rng(1)
    arr = rng.normal(size=shape)
    name = workdir / "map.bin"
    tomo = tomographic_objects.TomographicMap(
        map_array=arr, name=str(name), shape=shape, size=(10.0, 12.0, 14.0)
    )
    tomo.write()
    read = tomographic_objects.TomographicMap(name=str(name), shape=shape)
    read.read()
    np.testing.assert_allclose(read.map_array, arr)


def test_tomographic_map_read_requires_shape(workdir):
    tomo = tomographic_objects.TomographicMap(name=str(workdir / "x.bin"))
    with pytest.raises(ValueError):
        tomo.read()


def test_tomographic_map_write_property_file(workdir):
    shape = (4, 4, 4)
    tomo = tomographic_objects.TomographicMap(
        map_array=np.zeros(shape),
        name=str(workdir / "map.bin"),
        shape=shape,
        size=(8.0, 8.0, 8.0),
        coordinate_transform="middle",
        Omega_m=OMEGA_M,
    )
    prop_name = str(workdir / "prop.pickle")
    tomo.write_property_file(prop_name)
    prop = tomographic_objects.MapPixelProperty(name=prop_name)
    prop.read()
    assert prop.shape == shape
    assert prop.Omega_m == OMEGA_M


def test_tomographic_map_rebin_and_smooth(workdir):
    shape = (8, 8, 8)
    tomo = tomographic_objects.TomographicMap(
        map_array=np.ones(shape), name=str(workdir / "m.bin"),
        shape=shape, size=(8.0, 8.0, 8.0),
    )
    tomo.rebin_map((4, 4, 4), operation="mean")
    assert tomo.shape == (4, 4, 4)
    np.testing.assert_allclose(tomo.map_array, 1.0)
    tomo.smooth(2.0)  # must run and keep the shape
    assert tomo.map_array.shape == (4, 4, 4)


# --- Pixel -----------------------------------------------------------------


def test_pixel_roundtrip(workdir):
    rng = np.random.default_rng(2)
    arr = rng.normal(size=(30, 5))
    name = workdir / "pixel.bin"
    pixel = tomographic_objects.Pixel(pixel_array=arr, name=str(name))
    pixel.write()
    read = tomographic_objects.Pixel(name=str(name))
    read.read()
    np.testing.assert_allclose(read.pixel_array, arr)


def test_pixel_init_from_property_files(workdir):
    prop = make_synthetic_property(workdir / "prop.pickle")
    pixel = tomographic_objects.Pixel.init_from_property_files(
        str(workdir / "prop.pickle"), pixel_array=np.zeros((3, 5)), name="p.bin"
    )
    assert pixel.coordinate_transform == "middle"
    assert pixel.boundary_cartesian_coord == prop.boundary_cartesian_coord


def test_pixel_write_requires_array(workdir):
    pixel = tomographic_objects.Pixel(name=str(workdir / "p.bin"))
    with pytest.raises(ValueError):
        pixel.write()


def test_pixel_repack_by_los():
    # two lines of sight at (0,0) and (1,1), two pixels each
    arr = np.array(
        [
            [0.0, 0.0, 10.0, 0.1, 0.0],
            [0.0, 0.0, 20.0, 0.1, 0.0],
            [1.0, 1.0, 10.0, 0.1, 0.0],
            [1.0, 1.0, 20.0, 0.1, 0.0],
        ]
    )
    pixel = tomographic_objects.Pixel(pixel_array=arr)
    x, y, z = pixel.repack_by_los()
    assert len(x) == 2
    assert all(len(zi) == 2 for zi in z)


# --- Reference-data reads --------------------------------------------------


def test_read_reference_pixel_binary():
    """The committed reference pixel binary reads as an (N, 5) array."""
    prop = str(PIXEL_DIR / "property_file.pickle")
    pixel = tomographic_objects.Pixel.init_from_property_files(
        prop, name=str(PIXEL_DIR / "pixel_stripe82_DR16.bin")
    )
    pixel.read()
    assert pixel.pixel_array.ndim == 2
    assert pixel.pixel_array.shape[1] == 5
    assert pixel.pixel_array.shape[0] > 0


def test_read_reference_qso_catalog():
    """The committed reference QSO catalog loads with cartesian coordinates."""
    qso = tomographic_objects.QSOCatalog.init_from_fits(
        str(PIXEL_DIR / "qso_catalog.fits")
    )
    assert qso.coord.shape[1] == 3
    assert len(qso.primary_key) == qso.coord.shape[0]


# --- VoidCatalog round-trip ------------------------------------------------


def test_void_catalog_roundtrip(workdir):
    radius = np.array([5.0, 8.0, 12.0])
    coord = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
    mean_value = np.array([0.2, 0.3, 0.4])
    central_value = np.array([0.25, 0.35, 0.45])
    name = str(workdir / "voids.fits")
    void = tomographic_objects.VoidCatalog.init_from_dictionary(
        name,
        radius,
        coord,
        "cartesian",
        "middle",
        OMEGA_M,
        ((0.0, 0.0, 0.0), (10.0, 10.0, 10.0)),
        ((-0.01, -0.01, 2.1), (0.01, 0.01, 2.3)),
        other_array=[mean_value, central_value],
        other_array_name=["MEAN", "VALUE"],
    )
    void.write()
    read = tomographic_objects.VoidCatalog.init_from_fits(name)
    np.testing.assert_allclose(np.sort(read.radius), np.sort(radius))
    np.testing.assert_allclose(read.coord, coord)
    np.testing.assert_allclose(read.mean_value, mean_value)
    np.testing.assert_allclose(read.central_value, central_value)
    assert read.catalog_type == "cartesian"
    assert read.Omega_m == OMEGA_M


def test_void_catalog_merge(workdir):
    def build(name, radius):
        coord = np.zeros((len(radius), 3))
        v = tomographic_objects.VoidCatalog.init_from_dictionary(
            str(name), np.array(radius), coord, "cartesian", "middle", OMEGA_M,
            ((0.0, 0.0, 0.0), (10.0, 10.0, 10.0)),
            ((-0.01, -0.01, 2.1), (0.01, 0.01, 2.3)),
        )
        v.write()
        return str(name)

    a = build(workdir / "a.fits", [5.0, 6.0])
    b = build(workdir / "b.fits", [7.0])
    merged = tomographic_objects.VoidCatalog.init_by_merging(
        [a, b], name=str(workdir / "merged.fits")
    )
    assert len(merged.radius) == 3
