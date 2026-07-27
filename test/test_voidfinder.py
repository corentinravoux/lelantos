"""Unit and integration tests for :mod:`lelantos.voidfinder`."""

import numpy as np
import pytest

from lelantos import tomographic_objects, voidfinder

from conftest import (
    SPHERICAL_VOID_PARAMS,
    make_synthetic_void_map,
)


@pytest.fixture
def finder(workdir):
    """A VoidFinder whose map is never read (for testing pure helpers)."""
    return voidfinder.VoidFinder(
        str(workdir),
        "map.bin",
        dict(SPHERICAL_VOID_PARAMS),
        map_shape=(10, 10, 10),
        map_size=(10.0, 10.0, 10.0),
        number_core=1,
        delete_option="CLUSTERS",
    )


# --- Pure geometry helpers -------------------------------------------------


def test_convert_to_mpc(finder):
    coord = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    out = finder.convert_to_Mpc((2.0, 2.0, 2.0), coord)
    np.testing.assert_allclose(out, coord * 2.0)
    # empty input passes through untouched
    empty = np.empty((0, 3))
    assert finder.convert_to_Mpc((2.0, 2.0, 2.0), empty).shape == (0, 3)


def test_find_overlapers(finder):
    mpp = np.array([1.0, 1.0, 1.0])
    coord = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    radius = np.array([2.0, 2.0, 2.0])
    mask = finder.find_overlapers(mpp, coord[0], radius[0], radius, coord, mpc=True)
    # the first two overlap (dist 1 < 4), the far one does not
    np.testing.assert_array_equal(mask, [True, True, False])


def test_delete_voids_clusters_merges_overlap(finder):
    mpp = np.array([1.0, 1.0, 1.0])
    # two overlapping voids + one isolated -> two survivors
    coord = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [20.0, 0.0, 0.0]])
    radius = np.array([3.0, 2.0, 4.0])
    new_coord, new_radius = finder.delete_voids(mpp, radius, coord, mpc=True)
    assert len(new_radius) == 2
    # the larger radius of the overlapping pair is kept
    assert np.isclose(np.max(new_radius), 4.0)
    assert 3.0 in np.round(new_radius, 5)


def test_delete_voids_none_keeps_all(workdir):
    vf = voidfinder.VoidFinder(
        str(workdir), "map.bin", dict(SPHERICAL_VOID_PARAMS),
        map_shape=(10, 10, 10), map_size=(10.0, 10.0, 10.0),
        number_core=1, delete_option="NONE",
    )
    coord = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    radius = np.array([3.0, 2.0])
    new_coord, new_radius = vf.delete_voids(
        np.array([1.0, 1.0, 1.0]), radius, coord, mpc=True
    )
    np.testing.assert_array_equal(new_radius, radius)


def test_delete_voids_bad_option(workdir):
    vf = voidfinder.VoidFinder(
        str(workdir), "map.bin", dict(SPHERICAL_VOID_PARAMS),
        map_shape=(10, 10, 10), map_size=(10.0, 10.0, 10.0),
        number_core=1, delete_option="BOGUS",
    )
    with pytest.raises(ValueError):
        vf.delete_voids(np.array([1.0, 1.0, 1.0]), np.array([1.0]),
                        np.array([[0.0, 0.0, 0.0]]), mpc=True)


def test_get_name_catalog(finder):
    name = finder.get_name_catalog()
    assert name.startswith("Voids_SPHERICAL")
    assert "CLUSTERS_deletion" in name


# --- Integration: spherical finder on a synthetic map ----------------------


def test_find_voids_sphere_detects_injected_void(workdir):
    tomo, centre = make_synthetic_void_map(
        workdir / "property_file.pickle", workdir / "map.bin"
    )
    vf = voidfinder.VoidFinder(
        str(workdir),
        str(workdir / "map.bin"),
        dict(SPHERICAL_VOID_PARAMS),
        map_property_file=str(workdir / "property_file.pickle"),
        number_core=1,
        delete_option="CLUSTERS",
    )
    name = vf.find_voids()
    catalog = tomographic_objects.VoidCatalog.init_from_fits(name)
    # exactly one void, of a non-trivial radius, near the injected centre
    assert len(catalog.radius) == 1
    assert catalog.radius[0] > SPHERICAL_VOID_PARAMS["minimal_radius"]
    mpp = np.array(tomo.size) / (np.array(tomo.shape) - 1)
    expected_centre_mpc = np.array(centre) * mpp
    np.testing.assert_allclose(catalog.coord[0], expected_centre_mpc, atol=3 * mpp.max())


def test_find_voids_empty_map_raises(workdir):
    """Characterise the current behaviour on a map with no over-threshold pixel.

    The spherical finder does not currently short-circuit when no seed pixel is
    found: it indexes an empty ``(0,)`` coordinate array as if it were 2D and
    raises ``IndexError``.  This test pins that behaviour so a future fix (e.g.
    returning an empty catalog) is noticed.
    """
    from conftest import make_synthetic_property

    prop = make_synthetic_property(workdir / "property_file.pickle")
    tomo = tomographic_objects.TomographicMap.init_from_property_files(
        str(workdir / "property_file.pickle"), name=str(workdir / "map.bin")
    )
    tomo.map_array = np.zeros(prop.shape)  # all below threshold
    tomo.write()
    vf = voidfinder.VoidFinder(
        str(workdir), str(workdir / "map.bin"), dict(SPHERICAL_VOID_PARAMS),
        map_property_file=str(workdir / "property_file.pickle"),
        number_core=1, delete_option="CLUSTERS",
    )
    with pytest.raises(IndexError):
        vf.find_voids()
