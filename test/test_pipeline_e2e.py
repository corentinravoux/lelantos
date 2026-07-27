"""End-to-end pipeline tests.

The lelantos pipeline is:

    picca deltas -> pixel binary (+ launcher, property file)
                 -> [dachshund Wiener-filter solver] -> tomographic map
                 -> void catalog

The Wiener-filter solver (``dachshund.exe``) is not shipped with the source, so
the two ends of the pipeline are exercised separately:

* ``test_delta_to_solver_input`` covers delta -> solver input (real data);
* ``test_map_to_void_catalog`` covers map -> void catalog (synthetic map);
* ``test_full_chain_with_synthetic_map`` stitches them, substituting a
  fabricated map for the missing solver step.
"""

import numpy as np

from lelantos import cosmology, tomographic_objects, voidfinder

from conftest import (
    DELTA_CUTS,
    DELTA_DIR,
    OMEGA_M,
    SPHERICAL_VOID_PARAMS,
    make_synthetic_void_map,
)


def _converter():
    return cosmology.DeltaConverter(
        "./", OMEGA_M, str(DELTA_DIR), "middle", False, "dachshund",
        return_qso_catalog="qso_catalog.fits", return_sky_catalogs=True, repeat=False,
    )


def test_delta_to_solver_input(workdir):
    """Real deltas produce a pixel binary, a launcher and a property file."""
    converter = _converter()
    out = converter.transform_delta_to_pixel_file(**DELTA_CUTS)
    cart, cart_qso, _, sky, sky_qso, _, props = out

    properties = {
        "sigma_f": 0.1, "lperp": 13, "lpar": 13, "shape": (12, 12, 40),
        "name_pixel": "pixel.bin", "name_map": "map.bin",
    }
    shape = converter.create_serial_input("launch", properties, cart, sky)
    prop = converter.create_dachshund_map_pixel_property_file(
        "property_file.pickle", cart, sky, shape, props
    )
    prop.write()

    # all solver inputs are on disk
    for fname in ("launch.cfg", "pixel.bin", "property_file.pickle"):
        assert (workdir / fname).is_file()

    # the launcher declares the requested map grid
    cfg = (workdir / "launch.cfg").read_text()
    assert "map_nx = 12" in cfg
    assert "map_nz = 40" in cfg

    # the property file and pixel binary are mutually consistent
    pixel = tomographic_objects.Pixel.init_from_property_files(
        str(workdir / "property_file.pickle"), name=str(workdir / "pixel.bin")
    )
    pixel.read()
    np.testing.assert_allclose(pixel.pixel_array, cart)
    assert pixel.Omega_m == OMEGA_M
    assert pixel.coordinate_transform == "middle"


def test_map_to_void_catalog(workdir):
    """A tomographic map is turned into a void catalog by the finder."""
    make_synthetic_void_map(workdir / "property_file.pickle", workdir / "map.bin")
    finder = voidfinder.VoidFinder(
        str(workdir), str(workdir / "map.bin"), dict(SPHERICAL_VOID_PARAMS),
        map_property_file=str(workdir / "property_file.pickle"),
        number_core=1, delete_option="CLUSTERS",
    )
    name = finder.find_voids()
    catalog = tomographic_objects.VoidCatalog.init_from_fits(name)
    assert len(catalog.radius) == 1
    assert catalog.Omega_m == OMEGA_M


def test_full_chain_with_synthetic_map(workdir):
    """delta -> solver input -> (synthetic) map -> void catalog, end to end."""
    converter = _converter()
    out = converter.transform_delta_to_pixel_file(**DELTA_CUTS)
    cart, cart_qso, _, sky, sky_qso, _, props = out
    properties = {
        "sigma_f": 0.1, "lperp": 13, "lpar": 13, "shape": (12, 12, 40),
        "name_pixel": "pixel.bin", "name_map": "map.bin",
    }
    shape = converter.create_serial_input("launch", properties, cart, sky)
    prop = converter.create_dachshund_map_pixel_property_file(
        "property_file.pickle", cart, sky, shape, props
    )
    prop.write()

    # --- substitute the missing dachshund solver with a fabricated map -----
    # (same shape/geometry as the real solver would output) holding one void.
    prop.read()
    arr = np.zeros(shape)
    mpp = np.array(prop.size) / (np.array(shape) - 1)
    centre = (shape[0] // 2, shape[1] // 2, shape[2] // 2)
    ii, jj, kk = np.indices(shape)
    dist = np.sqrt(
        ((ii - centre[0]) * mpp[0]) ** 2
        + ((jj - centre[1]) * mpp[1]) ** 2
        + ((kk - centre[2]) * mpp[2]) ** 2
    )
    void_radius = 3 * mpp.min()
    arr[dist < void_radius] = 0.5
    tomo = tomographic_objects.TomographicMap.init_from_property_files(
        str(workdir / "property_file.pickle"), name=str(workdir / "map.bin")
    )
    tomo.map_array = arr
    tomo.write()

    # --- run the finder on the reconstructed-map stand-in ------------------
    # Scale the finder radii to the (coarse, anisotropic) pixel size of this
    # real-data box so that radii expressed in Mpc.h^-1 map onto a few pixels.
    params = dict(SPHERICAL_VOID_PARAMS)
    params["minimal_radius"] = mpp.min()
    params["maximal_radius"] = 6 * mpp.min()
    params["radius_step"] = mpp.min() / 2
    finder = voidfinder.VoidFinder(
        str(workdir), str(workdir / "map.bin"), params,
        map_property_file=str(workdir / "property_file.pickle"),
        number_core=1, delete_option="CLUSTERS",
    )
    name = finder.find_voids()
    catalog = tomographic_objects.VoidCatalog.init_from_fits(name)
    assert len(catalog.radius) >= 1
    assert catalog.catalog_type == "cartesian"
    assert catalog.Omega_m == OMEGA_M
