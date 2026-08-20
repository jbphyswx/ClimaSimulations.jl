# SOCRATES Single-Column Model (SCM) & Calibration Suite

This package provides a standalone, production-ready implementation of the Southern Ocean Clouds, Radiation, Aerosol Transport Experimental Study (**SOCRATES**) Single-Column Model (SCM) simulation, Atlas Large Eddy Simulation (LES) profile scoring, and Ensemble Kalman Inversion (EKI) calibration suite within `ClimaSimulations.jl`.

---

## Architecture

The codebase is organized into two strictly decoupled layers:

### Layer 1: Core Forward Model & LES Scoring (`src/`)
- **Independent Package (`SOCRATES.jl`)**: Contains all forward simulation logic, in-memory SSCF forcing, diagnostic registration, grid management, and Atlas LES profile scoring.
- **Zero Calibration Dependencies**: Does not import `ClimaCalibrate`, `EnsembleKalmanProcesses`, or `JLD2`.

Key files:
- [`src/cases.jl`](src/cases.jl): `SocratesCase` registry, 11 flight configurations (5 Obs, 6 ERA5), vertical grids, parameter merging.
- [`src/model.jl`](src/model.jl): In-memory `ColumnMemoryTimeVaryingInput`, `ScaledTerminalVelocity`, `socrates_simulation`, and execution runners (`run_case`, `run_cases`).
- [`src/scoring.jl`](src/scoring.jl): Atlas LES reference loading, `ScoreTransform`, vertical domain bounding, and normalized profile scoring (`compare_to_les`).

### Layer 2: Outer Calibration Layer (`calibration/`)
- **Decoupled Package (`SOCRATESCalibration.jl`)**: Implements `ClimaCalibrate.AbstractModelInterface`, SVD+diagonal noise covariance construction, `GEnsembleBuilder` observation maps, and EKI driver loops with staged `T_stops`.

Key files:
- [`calibration/src/calibration.jl`](calibration/src/calibration.jl): `SocratesInterface`, observation mapping, priors, and flat-parallel `calibrate` driver.

---

## Quickstart

### 1. Forward Simulation & LES Scoring (Layer 1)

```julia
using SOCRATES: SOCRATES

# Select case (e.g. Flight 09 with observational forcing)
case = SOCRATES.case("RF09_Obs")

# Run forward SCM simulation
output_dir = SOCRATES.run_case(case; FT = Float64, output_dir = "output/RF09_Obs")

# Compare against Atlas LES reference
comparison = SOCRATES.compare_to_les(case, output_dir)
SOCRATES.print_comparison(comparison)
```

### 2. EKI Parameter Calibration (Layer 2)

```julia
using ClimaCalibrate: ClimaCalibrate
using Distributed: Distributed
using SOCRATES: SOCRATES
using SOCRATESCalibration: SOCRATESCalibration

# Select cases for joint calibration
cases = [SOCRATES.case("RF09_Obs"), SOCRATES.case("RF10_Obs")]

# Build interface, prior, and EKP
interface = SOCRATESCalibration.SocratesInterface(; cases, output_dir = "output/calibration_run")
prior = SOCRATESCalibration.default_prior()
ekp = SOCRATESCalibration.build_ekp(interface, prior; ensemble_size = 4)

# Run calibration
backend = ClimaCalibrate.WorkerBackend(Distributed.default_worker_pool())
SOCRATESCalibration.calibrate(backend, ekp, interface; prior, n_iterations = 5)
```

---

## Running Tests

Run the unit test suite:

```julia
using Pkg
Pkg.test("SOCRATES")
```
