# Single-column cases

One package per case. Each is self-contained.

| package | state |
|---|---|
| `SOCRATES` | Atlas LES of SOCRATES research flights (Atlas et al. 2020), 11 cases |
| `MOSAiC_AYiL` | MOSAiC "A Year in LES" (Schnierstein et al. 2024), 190 DALES days |
| `SwirlLMCloudBench` | Swirl-LM CloudBench (Chammas et al. 2026), 500 sites × 4 months × 5 scenarios |
| `BOMEX` | trade-wind cumulus, one case; profiles from `AtmosphericProfilesLibrary` |
| `DYCOMS`, `TRMM_LBA`, `VARANAL`, `ShenGCMForcedLES` | not started |

## Layout
```
<PKG>/
  Project.toml  README.md
  src/          <PKG>.jl  cases.jl  forcing.jl  model.jl  io.jl  parallel.jl
                [+ case-specific]
  cases/        per-case data as declarative files, each value carrying its source
  configs/      parameter overrides
  test/         Project.toml  runtests.jl  scheduler.jl  [+ integration.jl]
  examples/     Project.toml  run_single_case.jl
  calibration/
    Project.toml
    src/        <PKG>Calibration.jl  calibration.jl  scoring.jl  reference.jl
                postprocess.jl
    test/       Project.toml  runtests.jl  postprocess.jl
    examples/   run_calibration.jl
  plots/
    Project.toml
    src/        <PKG>Plots.jl  run_plots.jl  calibration_plots.jl
    test/       Project.toml  runtests.jl
```

`src/` is the model and the forward run: the case types, the forcing, the grid, the
`AtmosModel`, the `AtmosSimulation`, reading a run's own output (`io.jl`'s `run_outputvars`), and `parallel.jl`, which runs a list of cases over a worker pool.

`calibration/` is the EKI layer: which variables are compared, over what window, with what
weighting, the prior, the observations, and the `ClimaCalibrate` interface.

`plots/` is the figures