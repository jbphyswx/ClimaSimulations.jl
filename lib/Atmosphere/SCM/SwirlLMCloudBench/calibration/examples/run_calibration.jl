"""
    run_calibration.jl

Ensemble Kalman Inversion of CloudBench microphysics parameters against a case's published
output.

```julia
include("calibration/examples/run_calibration.jl")

run_calibration()                                       
run_calibration(; run_kwargs = (; t_end = 3600.0))       # short members, for a smoke test
run_calibration(; ensemble_size = 4, n_iterations = 1)
run_calibration(; cases = [S.CloudBenchInstance(7, 4, :amip_p4k)])
```

Every argument is a keyword with a default; none of it is baked in. `run_kwargs` is forwarded
to `SwirlLMCloudBenchSim.run_case`
"""

using ClimaCalibrate: ClimaCalibrate
using SwirlLMCloudBench: Simulation as S
using SwirlLMCloudBenchCalibration: SwirlLMCloudBenchCalibration as CBC

function run_calibration(;
    cases = [S.CloudBenchInstance(0, 1, :amip)],
    output_dir = joinpath(@__DIR__, "output", "calibration_demo"),
    ensemble_size::Int = 10,
    n_iterations::Int = 10,
    backend = ClimaCalibrate.JuliaBackend(),
    prior = CBC.default_prior(),
    vars = collect(String, CBC.DEFAULT_SCORED_VARS),
    run_kwargs = (;),
    interface_kwargs = (;),
    ekp_kwargs = (;),
)
    mkpath(output_dir)
    @info "calibrating" cases = length(cases) ensemble_size n_iterations output_dir

    interface = CBC.CloudBenchInterface(;
        cases, output_dir, vars, run_kwargs, interface_kwargs...,
    )
    ekp = CBC.build_ekp(interface, prior; ensemble_size, ekp_kwargs...)

    return ClimaCalibrate.calibrate(
        backend, ekp, interface, n_iterations, prior, output_dir,
    )
end