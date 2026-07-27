# Pipeline & configuration

The end-to-end pipeline is driven by `lelantos.interface.main`, which reads an
`.ini` configuration file. A fully commented reference version — documenting
**every option** and whether it is mandatory — ships as
[`scripts/interface_explanatory.ini`](https://github.com/corentinravoux/lelantos/blob/main/scripts/interface_explanatory.ini).

## Stages

The `[main]` section is a switchboard: each boolean turns one stage on or off.
Stages run in this order, each feeding the next:

| Order | `[main]` gate | Section(s) read | What it does | Implemented by |
|---|---|---|---|---|
| 1 | `transform_delta` | `[delta transform]` | Shuffle / redshift-cut the deltas | [`DeltaModifier`][lelantos.cosmology.DeltaModifier] |
| 2 | `plot_delta` | `[delta plot]` | Diagnostic plots of the deltas | [`DeltaAnalyzer`][lelantos.cosmology.DeltaAnalyzer] |
| 3 | `convert_delta` | `[delta convert]`, `[tomography software]` | Deltas → solver pixel input + property file | [`DeltaConverter`][lelantos.cosmology.DeltaConverter] |
| 4 | `launch_tomography` | `[tomography launching]`, `[tomography software]` | Run the solver | [`TomographyManager`][lelantos.task_manager.TomographyManager] |
| 5 | `process_tomography` | `[tomography process]` | Merge sub-maps, build distance map | [`create_merged_map`][lelantos.tomography.create_merged_map], [`create_distance_map`][lelantos.tomography.create_distance_map] |
| 6 | `find_void` | `[void finder]` | Detect voids | [`VoidFinder`][lelantos.voidfinder.VoidFinder] |
| 7 | `process_void` | `[void process]` | Stats, cuts, cross-corr catalog | [`cut_catalog`][lelantos.voidfinder.cut_catalog] |
| 8 | `stack_void` | `[void stack]` | Stack the map at void positions | [`TomographyStack`][lelantos.tomography.TomographyStack] |
| 9 | `plot_void` | `[void plot]` | Void-catalog plots | [`PlotVoid`][lelantos.voidfinder.PlotVoid] |
| 10 | `plot_tomography` | `[tomography plot]` | Map slices / histograms | [`TomographyPlot`][lelantos.tomography.TomographyPlot] |
| 11 | `plot_stack_void` | `[stack void plot]` | Void-stack plots | [`TomographyStack.plot_stack`][lelantos.tomography.TomographyStack.plot_stack] |

## Reading the config

Every value is parsed through a typed getter:

- `getstr` / `getint` / `getfloat` — scalars;
- `gettupleint` / `gettuplefloat` — comma-separated lists (`2,2,300`);
- `gettuplestr` — comma-separated strings;
- `getdict` — a Python-literal dict (`{"n": 8}`);
- `getboolean` — `True` / `False`.

The literal string `None` is parsed to Python `None` for every scalar getter —
use it to explicitly disable an option.

!!! note "Mandatory vs optional"
    In `interface_explanatory.ini` an option marked **MANDATORY** must be
    present *and* hold a real (non-`None`) value whenever its stage is enabled;
    an **OPTIONAL** option must still be present but may be set to `None` to
    disable the associated behaviour. A missing key raises
    `configparser.NoOptionError` at runtime.

!!! warning "Known key-name mismatch"
    The `[stack void plot]` stage reads `stack_name` and `property_stack_name`,
    which differ from the `name_stack` / `property_file_stack` keys written by
    the `[void stack]` stage. Use the names documented in
    `interface_explanatory.ini`.

## Sizing a run before launching

Use
[`print_approximate_shape_size_from_interface_file`][lelantos.interface.print_approximate_shape_size_from_interface_file]
(exposed by `scripts/lelantos_print_map_prop.py`) to estimate the resulting map
pixel shape and physical size from the footprint and chunking configuration,
before committing to a full reconstruction.
