# LELANTOS — Large-scalE LymAN-alpha TOmography Suite

LELANTOS — Greek titan associated with hunting and air (name found after
discussion with [E. Chaussidon](https://github.com/echaussidon)).

Package to create and manipulate large-scale Lyman-α tomographic maps.

It takes the flux-contrast files (**delta**) produced by
[picca](https://github.com/igmhub/picca), converts them into the input of a
tomographic solver, runs the reconstruction, and post-processes the resulting
3D map into line-of-sight (**pixel**) products, **tomographic maps**, **void
catalogs** and diagnostic plots. Every stage can be run on its own (via a
launcher, a script or the API) or all at once through a single configuration
file.

📖 **Full documentation (API reference + user guide): the MkDocs site built
from this repository** — see [Documentation](#documentation) below.

---

## Pipeline overview

```
picca deltas ──DeltaConverter──▶ pixel binary ──TomographyManager──▶ tomographic map
                                                     (Dachshund)              │
                                                                              ├─▶ VoidFinder ──▶ void catalog ──▶ stacks / QSO-like x-corr catalog
                                                                              └─▶ map & delta plots
```

---

## Install

This package can be installed with pip:

```
pip install . --user
```

The minimal required python packages are: `fitsio`, `numpy`, `scipy`,
`matplotlib`, `picca` (installed automatically with pip).

To convert a tomographic box into a format for 3D visualization, you need
`pyevtk` as an additional package.

This package uses code parts from the public packages
[picca](https://github.com/igmhub/picca) and
[SaclayMocks](https://github.com/igmhub/SaclayMocks) — it is advised to install
them independently (a vendored SaclayMocks fallback ships under
`lelantos/saclaymocks/`).

To run tomographic mapping, place the executable of the tomographic algorithm
you want to use in the `lelantos/lelantos/exec` directory. For now the only
algorithm available is [dachshund](https://github.com/caseywstark/dachshund).

---

## Running the whole pipeline

Most of the launchers below can be run in one go by editing a configuration
file and running:

```
lelantos_interface.py interface_example.ini
```

- `scripts/interface_example.ini` — a ready-to-edit example.
- `scripts/interface_explanatory.ini` — **a fully commented reference** that
  documents *every* option and whether it is mandatory or optional.

The `[main]` section is a switchboard of booleans, one per stage; each stage
reads its own config section. See the documentation
[Pipeline & config](#documentation) page for the full description.

---

## Launchers

The main launchers available in the `launchers/` folder are:

- **lelantos_subsample_delta.py** — Subsampling of flux-contrast files (delta), output of the picca delta generation routine.
- **lelantos_shuffle_delta.py** — Shuffling of flux-contrast files (delta).
- **lelantos_convert_delta.py** — Conversion of flux-contrast files (delta) to a pixel file, input of the tomographic mapping procedure.
- **lelantos_tomography_manager.py** — Running the tomography algorithm.
- **lelantos_tomography_process.py** — Processing the tomographic map output.
- **lelantos_void_finder.py** — Launching the void-finding routine on the tomographic map.
- **lelantos_void_process.py** — Post-processing of the void catalog.
- **lelantos_plot_delta.py** — Plotting routines for flux-contrast files (delta).
- **lelantos_plot_pixel.py** — Plotting routines for the pixel file.
- **lelantos_plot_tomography.py** — Plotting routines for the tomographic map.
- **lelantos_plot_void.py** — Plotting routines for the void catalog.

Additional launchers for specific studies:

- **lelantos_stack_void.py** — Stacking of void positions on a tomographic map.
- **lelantos_stack_qso.py** — Stacking of quasar positions on a tomographic map.
- **lelantos_plot_stack_void.py** — Plotting routines for the void stack.
- **lelantos_plot_stack_qso.py** — Plotting routines for the QSO stack.
- **lelantos_tomography_3d.py** — Conversion of tomographic objects into a format for 3D visualization.
- **lelantos_boxdm.py** — Extraction of the underlying dark-matter box of mocks associated with a given box geometry.

Each launcher is configured by editing the variables at the top of the file,
then running it directly (`python launchers/<name>.py`).

---

## Scripts

The `scripts/` folder contains small utilities using `argparse` or a config
parser to quickly modify or show information about tomographic objects
(print a property file, cut a QSO/DLA catalog, show void statistics, estimate a
map shape/size, …), plus `lelantos_interface.py` for the full pipeline.

---

## Package modules

| Module | Role |
|---|---|
| `lelantos.cosmology` | Delta ingestion, shuffling/subsampling, solver-input generation, delta/pixel diagnostics |
| `lelantos.tomographic_objects` | Core data objects: maps, distance maps, stacks, pixels, deltas, property files, catalogs |
| `lelantos.task_manager` | Run the tomographic solver locally (bash) or on a cluster (NERSC / Irene) |
| `lelantos.tomography` | Map post-processing, slice/histogram plots, stacking |
| `lelantos.voidfinder` | Spherical / watershed void detection and post-processing |
| `lelantos.boxdm` | Extract the underlying SaclayMocks dark-matter field onto a map geometry |
| `lelantos.utils` | Coordinate transforms, interpolation, plotting, logging helpers |
| `lelantos.interface` | Config-file driven end-to-end pipeline |

---

## Conventions

- RA / Dec handled in **radians** internally (degrees at the plotting
  boundary); comoving distances and box sizes in **Mpc·h⁻¹**; wavenumbers in
  h/Mpc.
- The sky ↔ cartesian transform is chosen per run via `coordinate_transform`
  (`middle`, `full_angle`/`full`, `cartesian`); most downstream code assumes the
  `middle` tangent-plane transform.
- The fiducial cosmology enters **only** at delta conversion, through `Omega_m`.

---

## Documentation

Documentation is written with [MkDocs Material](https://squidfunk.github.io/mkdocs-material/)
and [mkdocstrings](https://mkdocstrings.github.io/) (the API reference is
generated from the in-code docstrings) and is published to **GitHub Pages** by
the `.github/workflows/docs.yml` workflow on every push to `main`
(enable *Settings → Pages → Source: GitHub Actions* once).

Build it locally with:

```
pip install -r requirements-docs.txt
mkdocs serve      # live preview at http://127.0.0.1:8000
mkdocs build      # static site into site/
```
