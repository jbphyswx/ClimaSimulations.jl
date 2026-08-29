# MOSAiC_AYiL

Single-column ClimaAtmos runs of the MOSAiC "A Year in LES" days: 190 DALES large-eddy
simulations over the drifting MOSAiC ice camp, 2019-10-16 to 2020-09-11, forced by ERA5
Lagrangian-trajectory `scm_in` files.

Schnierstein et al. (2024), JAMES, [10.1029/2024MS004296](https://doi.org/10.1029/2024MS004296);
data at [10.5281/zenodo.10491362](https://doi.org/10.5281/zenodo.10491362).

`docs/design.md` is the data contract — file roles, the species and scalar mapping, the three
moisture conventions, the DALES constants and how far they sit from CliMA's, the nudging ramp,
the surface flux blending, and every open discrepancy, each marked for whether it was verified
against the data or taken from a document. Read it before changing this package.

## Running a day

```julia
using MOSAiC_AYiL: MOSAiC_AYiL as MA

c = MA.case("20200503")
MA.run_case(c; output_dir = "output/20200503")
```

`mosaic_simulation(FT, case; output_dir, kwargs...)` returns the `AtmosSimulation` without
solving it. Every default is a keyword, so the grid, the parameters, `t_end`, `dt`, the
diagnostics and the model components can each be replaced at the call site.

`examples/run_days.jl` drives several days in one session.

## Data

The 190 days are a lazy `Pkg` artifact fetched from Zenodo on first use (~911 MB), built by
`gen/build_data_artifact.jl`. Set `MOSAIC_AYIL_DATA_ROOT` to point at a local tree instead.

```julia
MA.available_dates()      # the days whose directory is present
MA.case("20200503")       # one day
```

## Grid

The LES grid is 286 cells to a top face of 11 949.728 m. A column can use it as it stands, or
truncate to the day's usable top and coarsen:

```julia
faces = MA.truncate_faces_to_top(MA.native_faces(c), MA.best_simulation_top(c))
grid = MA.mosaic_grid(Float64, c; faces = MA.coarsen_faces_to_dz_min(faces, 200))
```

## Tests

```bash
julia --project=test test/runtests.jl     # no forward run
julia --project=test test/integration.jl  # builds a column and steps it; gated
```

`runtests.jl` is the default suite. `integration.jl` is a separate gated job: it constructs a
real column and solves a short window, so it costs minutes and several GB.

## Calibration

`calibration/` holds the EKI layer and does not participate in forward runs.

```bash
julia --project=calibration/test calibration/test/runtests.jl
julia --project=calibration/examples calibration/examples/run_calibration.jl
```
