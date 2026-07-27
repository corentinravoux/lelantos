# API reference

Every public module, class and function of lelantos is documented from its
in-code docstrings. Use the navigation to browse a module, or jump straight to
a topic below.

- **[interface](interface.md)** — config-file driven end-to-end pipeline.
- **[cosmology](cosmology.md)** — delta ingestion, shuffling/subsampling,
  solver-input generation and delta/pixel diagnostics.
- **[tomographic_objects](tomographic_objects.md)** — the core data objects:
  maps, distance maps, stacks, pixels, deltas, property files and catalogs
  (QSO / DLA / galaxy / void).
- **[tomography](tomography.md)** — map post-processing, slice/histogram plots
  and stacking.
- **[voidfinder](voidfinder.md)** — spherical / watershed void detection and
  void-catalog post-processing and plotting.
- **[task_manager](task_manager.md)** — running the tomographic solver locally
  (bash) or on a cluster (NERSC / Irene).
- **[boxdm](boxdm.md)** — extracting the underlying SaclayMocks dark-matter
  field onto a map geometry.
- **[utils](utils.md)** — coordinate transforms, interpolation, plotting and
  logging helpers.
- **[saclaymocks](saclaymocks.md)** — vendored SaclayMocks geometry/constants.
