# lelantos test suite

Run from the package root:

```bash
pip install -e ".[dev]"   # installs pytest alongside the runtime deps
pytest -q                 # or: pytest test/
```

The suite mixes fast unit tests with data-driven integration and end-to-end
tests built from the reference data in [`data/`](data) (one picca delta file
plus the matching pixel / property / QSO reference files).

| File | Scope | What it covers |
|---|---|---|
| `test_utils.py` | unit | Pixel/size geometry, sky↔cartesian transforms (round-trip for `full_angle`/`full`/`middle`), footprint cuts, `bin_ndarray`, smoothing, interpolation, misc helpers |
| `test_tomographic_objects.py` | unit | Read/write round-trips for `TomographicMap`, `Pixel`, `MapPixelProperty`, `QSOCatalog`, `VoidCatalog`; `init_false_prop`; rebin/smooth; reference-data reads |
| `test_interface.py` | unit | The `interface` config-parsing helpers (`parse_int_tuple`, `parse_dict`, …) |
| `test_cosmology.py` | unit + integration | Delta readers (`get_delta_list`, `get_deltas`), `DeltaConverter` invariants and determinism, serial solver-input writing |
| `test_voidfinder.py` | unit + integration | Overlap geometry (`find_overlapers`, `delete_voids`), spherical finder on a synthetic map, empty-map edge case |
| `test_task_manager.py` | unit | Log parsing, Dachshund `.cfg` writing, manager directory helpers, backend selection (no solver execution) |
| `test_pipeline_e2e.py` | end-to-end | delta → solver input; map → void catalog; and the full chain with a fabricated map standing in for the solver |
| `test_scripts_smoke.py` | smoke | Byte-compiles every `scripts/` and `launchers/` entry point; runs `--help` on the argparse CLIs |

## Notes

* **No solver.** The Wiener-filter solver (`dachshund.exe`) is not shipped with
  the source tree, so end-to-end tests substitute a fabricated tomographic map
  (matching the real box geometry) for the reconstruction step.
* **Hermetic.** Tests that write files run inside a temporary working directory
  (the `workdir` fixture in [`conftest.py`](conftest.py)); the repository is
  never modified.
* **Reference footprint.** RA/Dec cut bounds are given in **degrees**; the
  committed delta file sits around RA ≈ −148°, Dec ≈ 53°, z ∈ [1.9, 3.5].
