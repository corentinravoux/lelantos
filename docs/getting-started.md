# Getting started

## Installation

```bash
pip install .            # from the repository root
# or, for development:
pip install -e .
```

Minimal required packages (installed automatically): `fitsio`, `numpy`,
`scipy`, `matplotlib`, `picca`.

Optional:

- `pyevtk` — to export a tomographic box to VTK for 3D visualisation.
- [SaclayMocks](https://github.com/igmhub/SaclayMocks) — used by
  `lelantos.boxdm`; a vendored fallback copy ships under
  `lelantos/saclaymocks/`, but installing SaclayMocks directly is recommended.

### The tomographic solver

To run the reconstruction you must place the executable of the tomographic
algorithm in `lelantos/lelantos/exec/`. The only algorithm currently supported
is [Dachshund](https://github.com/caseywstark/dachshund) (a Wiener-filter
reconstruction), expected as `dachshund.exe`.

## Two ways to run

### 1. The all-in-one interface

Copy and edit `scripts/interface_example.ini` (a fully commented reference is
provided as `scripts/interface_explanatory.ini`), then run:

```bash
lelantos_interface.py interface_example.ini
```

Each stage is switched on/off by a boolean in the `[main]` section. See the
[Pipeline & config](pipeline.md) page for the meaning of every option.

### 2. Individual launchers

The `launchers/` folder holds standalone scripts that expose one stage each.
Edit the variables at the top of the file and run it directly, e.g.:

```bash
python launchers/lelantos_convert_delta.py
python launchers/lelantos_void_finder.py
```

The `scripts/` folder additionally provides small `argparse` / config-parser
utilities to inspect or modify tomographic objects (print a property file, cut
a QSO/DLA catalog, show void statistics, …).

## Coordinate & unit conventions

- RA / Dec are handled in **radians** internally (degrees at the plotting
  boundary), comoving distances and box sizes in **Mpc·h⁻¹**, velocities in
  km/s and wavenumbers in h/Mpc.
- The sky ↔ cartesian transform is selected per run by `coordinate_transform`
  (`middle`, `full_angle`/`full`, or `cartesian`). Most downstream code assumes
  the `middle` tangent-plane transform.
- The fiducial cosmology enters **only** at the delta-conversion step, through
  `Omega_m`.
