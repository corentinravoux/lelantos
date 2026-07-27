# LELANTOS

**L**arge-scal**E** **L**ym**AN**-alpha **TO**mography **S**uite.

LELANTOS creates and manipulates large-scale Lyman-α tomographic maps. It takes
the flux-contrast files (**delta**) produced by
[picca](https://github.com/igmhub/picca), converts them into the input of a
tomographic solver, runs the reconstruction, and post-processes the resulting
3D map into **line-of-sight (pixel)** products, **tomographic maps**, **void
catalogs** and diagnostic plots.

The name is that of the Greek titan associated with hunting and air (found after
discussion with [E. Chaussidon](https://github.com/echaussidon)).

## The pipeline at a glance

```mermaid
flowchart LR
    delta["picca deltas"] -->|DeltaConverter| pixel["pixel binary"]
    pixel -->|TomographyManager<br/>(Dachshund)| map["tomographic map"]
    map -->|VoidFinder| voids["void catalog"]
    map --> plots["map / delta plots"]
    voids --> stack["void stacks"]
    voids --> xcorr["QSO-like catalog<br/>(cross-correlation)"]
```

Every stage can be run individually through a **launcher**, a **script**, or
the lelantos API — or all at once through a single configuration file with
`lelantos_interface.py`.

## Where to go next

- **[Getting started](getting-started.md)** — install and run your first map.
- **[Pipeline & config](pipeline.md)** — the end-to-end `interface.ini` workflow
  and every configuration option.
- **[API reference](api/index.md)** — the documented modules, classes and
  functions.

## Package layout

| Component | Role |
|---|---|
| `lelantos.cosmology` | Delta ingestion, shuffling/subsampling, solver-input generation, delta/pixel diagnostics |
| `lelantos.tomographic_objects` | Core data objects: maps, pixels, deltas, catalogs, property files |
| `lelantos.task_manager` | Launch the tomographic solver locally or on a cluster |
| `lelantos.tomography` | Map post-processing, slicing/plots, stacking |
| `lelantos.voidfinder` | Spherical / watershed void detection and post-processing |
| `lelantos.boxdm` | Extract the underlying SaclayMocks dark-matter field |
| `lelantos.utils` | Coordinate transforms, interpolation, plotting, logging helpers |
| `lelantos.interface` | Config-file driven end-to-end pipeline |
