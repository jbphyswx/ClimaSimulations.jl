# SOCRATES

Single-column simulations of the SOCRATES field campaign, set up to reproduce the
large-eddy simulations of [Atlas et al. (2020)](https://doi.org/10.1029/2020MS002205),
and to calibrate microphysics against them.

Eleven cases: research flights 1, 9, 10, 11, 12 and 13, each forced either from
the observed soundings (`Obs`, a 12 h run) or from ERA5 (`ERA5`, 14 h). Flight 11
has no `Obs` artifact. Forcing and reference data come from
[SOCRATESSingleColumnForcings.jl](https://github.com/jbphyswx/SOCRATESSingleColumnForcings.jl)
as lazily-downloaded artifacts.

## Running a case

Nothing has to be generated or configured first.

```julia
using SOCRATES: SOCRATES

case = SOCRATES.case("RF09_Obs")
dir = SOCRATES.run_case(case; output_dir = "output/RF09_Obs")

SOCRATES.print_comparison(SOCRATES.compare_to_les(case, dir))
```

`run_case` builds the case on its own Atlas vertical grid and runs it for the
Atlas run length. `SOCRATES.all_cases()` lists all eleven.

## What the case is

  - **Grid**: the Atlas LES levels for the flight — 320 levels to 6.1 km for
    flight 9, 320 to 4.8 km for 1, 10 and 11, 192 to 2.9 km for 12 and 13.
  - **Forcing**: horizontal advection of temperature and total water, relaxation
    of temperature, total water and horizontal wind toward the case profiles, and
    large-scale subsidence. SSCF carries no vertical eddy tendency; vertical
    transport comes from subsidence acting on the model's own profiles.
  - **Surface**: Monin-Obukhov fluxes over the prescribed sea surface temperature.
  - **Radiation**: interactive RRTMGP at a zenith angle held fixed at the case
    reference time, following Atlas et al. (2020) section 4.2.
  - **Physics**: prognostic EDMFX with one updraft, non-equilibrium 1-moment
    microphysics with a quadrature cloud.

## Layout

| path | holds |
|---|---|
| `src/casedata.jl` | readers for `cases/` and `configs/` |
| `src/cases.jl` | the case type, its grid, its parameter set |
| `src/forcing.jl` | `SOCRATESForcing` and its cache and tendency |
| `src/model.jl` | scaled sedimentation, prescribed surface temperature and insolation, the initial condition, and the simulation |
| `src/scoring.jl` | Atlas LES reference, scored region, normalized misfit |
| `cases/socrates.toml` | per-case data, each value with its source |
| `configs/socrates.toml` | ClimaParams overrides |
| `calibration/` | the EKI layer, a separate package |

Case values and model components are all keyword arguments. `socrates_model` and
`socrates_simulation` forward `kwargs...` to `AtmosModel` and `AtmosSimulation`,
so any field of either can be replaced without editing this package:

```julia
SOCRATES.run_case(
    case;
    output_dir = "output/coarse",
    grid = SOCRATES.socrates_grid(Float64, case; dz_min = 100.0),
    t_end = 3600,
    model_kwargs = (; cloud_model = SOCRATES.CA.GridScaleCloud()),
)
```

## Tests

```julia
using Pkg
Pkg.test("SOCRATES")
```

The suite covers the case registry, the grids, the parameter composition, the
forcing on every case, the Atlas reference, and a short simulation that has to
reach its end time.
