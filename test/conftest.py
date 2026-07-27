"""Shared pytest fixtures and helpers for the lelantos test suite.

The suite is built from three complementary sources:

* the committed reference data in ``test/data`` (one picca delta file and the
  matching pixel / property / QSO reference files);
* short pipeline runs driven directly through the public API
  (:class:`lelantos.cosmology.DeltaConverter`,
  :class:`lelantos.voidfinder.VoidFinder`);
* smoke imports / ``--help`` runs of every ``scripts/`` and ``launchers/``
  entry point.

The heavy Wiener-filter solver (``dachshund.exe``) is not shipped with the
source tree, so end-to-end tests that would need a reconstructed map fabricate
a controlled synthetic map instead of invoking the solver.
"""

import os
import sys
from pathlib import Path

import numpy as np
import pytest

# --- Locations -------------------------------------------------------------

TEST_DIR = Path(__file__).resolve().parent
# Make ``lelantos`` importable even when pytest is not launched from the
# package root and the package is not pip-installed.
PACKAGE_ROOT = TEST_DIR.parent
for _path in (str(PACKAGE_ROOT), str(TEST_DIR)):
    if _path not in sys.path:
        sys.path.insert(0, _path)
DATA_DIR = TEST_DIR / "data"
DELTA_DIR = DATA_DIR / "delta"
PIXEL_DIR = DATA_DIR / "pixel"

# Fiducial cosmology used to generate the committed reference data.
OMEGA_M = 0.3147

# Sky / redshift footprint of the committed ``delta-19.fits.gz`` line of sight
# bundle.  RA/Dec cut bounds are expressed in DEGREES (compared internally
# against the radian coordinates via ``np.radians``); the delta file sits
# around RA ~ -148 deg, Dec ~ 53 deg.
DELTA_CUTS = dict(
    sigma_min=0.1,
    sigma_max=None,
    z_cut_min=2.1,
    z_cut_max=3.2,
    dec_cut_min=51.0,
    dec_cut_max=55.0,
    ra_cut_min=-150.0,
    ra_cut_max=-143.0,
)


# --- Fixtures --------------------------------------------------------------


@pytest.fixture
def workdir(tmp_path, monkeypatch):
    """Change into an isolated temp directory for the duration of a test.

    Many pipeline helpers write their outputs relative to ``pwd`` / the current
    working directory; running each test in its own directory keeps them
    hermetic and leaves the repository clean.
    """
    monkeypatch.chdir(tmp_path)
    return tmp_path


@pytest.fixture
def delta_converter():
    """A :class:`DeltaConverter` reading the committed reference delta file."""
    from lelantos import cosmology

    return cosmology.DeltaConverter(
        "./",
        OMEGA_M,
        str(DELTA_DIR),
        "middle",
        False,
        "dachshund",
        return_qso_catalog="qso_catalog.fits",
        return_dla_catalog=None,
        dla_catalog=None,
        return_sky_catalogs=True,
        repeat=False,
    )


@pytest.fixture
def converted_deltas(delta_converter):
    """Run the delta->pixel conversion once and return all of its outputs."""
    return delta_converter.transform_delta_to_pixel_file(**DELTA_CUTS)


def make_synthetic_property(name, shape=(40, 40, 60), size=(80.0, 80.0, 120.0)):
    """Build and write a controlled :class:`MapPixelProperty` (middle transform).

    The box is intentionally small so that ``mpc_per_pixel`` is O(2 Mpc.h^-1)
    and void radii expressed in Mpc map onto several pixels.
    """
    from lelantos import tomographic_objects

    prop = tomographic_objects.MapPixelProperty.init_false_prop(
        shape, size, OMEGA_M, 0.0, 0.0, 2.2, "middle", name=str(name)
    )
    prop.write()
    return prop


def make_synthetic_void_map(property_name, map_name, radius_mpc=8.0):
    """Fabricate a tomographic map holding a single spherical high-flux void.

    Returns the written :class:`TomographicMap` and the (pixel-space) centre of
    the injected void.
    """
    from lelantos import tomographic_objects

    prop = make_synthetic_property(property_name)
    shape = prop.shape
    size = prop.size
    mpp = np.array(size) / (np.array(shape) - 1)
    centre = (shape[0] // 2, shape[1] // 2, shape[2] // 2)
    ii, jj, kk = np.indices(shape)
    dist = np.sqrt(
        ((ii - centre[0]) * mpp[0]) ** 2
        + ((jj - centre[1]) * mpp[1]) ** 2
        + ((kk - centre[2]) * mpp[2]) ** 2
    )
    arr = np.zeros(shape)
    arr[dist < radius_mpc] = 0.5
    tomo = tomographic_objects.TomographicMap.init_from_property_files(
        str(property_name), name=str(map_name)
    )
    tomo.map_array = arr
    tomo.write()
    return tomo, centre


SPHERICAL_VOID_PARAMS = dict(
    method="SPHERICAL",
    threshold=0.1,
    average=0.1,
    minimal_radius=2.0,
    maximal_radius=30.0,
    radius_step=1.0,
    dist_clusters=2.0,
)
