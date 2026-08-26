"""
    run_calibration.jl

Run an AYiL EKI calibration.

```
> include("calibration/run_calibration.jl")
```

Everything the run chooses is here: which days, which variables, the grid, the
ensemble size, and the EKP process. The interface and the scoring hold no run
policy of their own.
"""

using ClimaCalibrate: ClimaCalibrate
using EnsembleKalmanProcesses: EnsembleKalmanProcesses as EKP
using MOSAiC_AYiL: MOSAiC_AYiL
using Random: Random

calibration_vars = ("cli", "clw", "husra", "hussn", "iwp", "lwp", "rwp", "swp")
# calibration_vars = ("ql_all", "qi_all", "qi_all_wp", "ql_all_wp")

include(joinpath(@__DIR__, "src", "MOSAiCAYiLCalibration.jl"))
include(joinpath(@__DIR__, "configs", "prior.jl"))
const MC = MOSAiCAYiLCalibration

function parse_args(args = ARGS)
    out = Dict{String, Any}(
        "output-dir" => joinpath(pwd(), "ayil_calibration"),
        "ensemble-size" => 20,
        "iterations" => 5,
        "minibatch-size" => 8,
        "dz-min" => 50,
        "dates" => nothing,
    )
    i = 1
    while i <= length(args)
        key = lstrip(args[i], '-')
        haskey(out, key) || error("Unknown argument `$(args[i])`")
        i + 1 <= length(args) || error("`$(args[i])` needs a value")
        raw = args[i + 1]
        out[key] = key in ("output-dir",) ? raw :
                   key == "dates" ? split(raw, ',') : parse(Int, raw)
        i += 2
    end
    return out
end

function calibrate(;
    cases = MC.default_calibration_cases(),
    outdir::AbstractSring,
    calibration_vars = ("cli", "clw", "husra", "hussn", "iwp", "lwp", "rwp", "swp"),
    ensemble_size = 20, 
    iterations = 10, 
    T_stops = [1., 10., 100., 1000.],
    minibatch_size = length(cases), 
    grid = MC.default_mosaic_ayil_grid(; dz_min = 100),
    dates = MC.default_calibration_dates(),
    prior = default_ayil_prior(),
    seed = 1234,
)
    interface = MC.MOSAiCInterface(;
        output_dir = outdir,
        cases,
        vars = collect(String, calibration_vars),
        grid,
    )
    series = MC.observation_series(interface; minibatch_size)


    rng = Random.MersenneTwister(seed)
    # `Inversion` takes the initial ensemble as its first argument; only the
    # unscented processes are built from the observation series alone
    ekp = EKP.EnsembleKalmanProcess(
        EKP.construct_initial_ensemble(rng, prior, ensemble_size),
        series,
        EKP.Inversion();
        rng,
        scheduler = EKP.DataMisfitController(; terminate_at = first(T_stops)),
        localization_method = EKP.SECNice(),
        accelerator = EKP.NesterovAccelerator(),
        # a member whose runs failed comes back as a NaN column from
        # `observation_map`, and this is what draws from the successful ones instead
        # of propagating it
        failure_handler_method = EKP.SampleSuccGauss(),
        verbose = true,
    )

    @info "AYiL calibration" days = length(cases) vars = interface.vars ensemble_size iterations minibatch_size T_stops output =
        interface.output_dir

    ekp = MC.calibrate(
        ClimaCalibrate.WorkerBackend(),
        ekp,
        interface,
        prior,
        interface.output_dir;
        T_stops,
        max_iterations = iterations,
    )
    return (; interface, ekp)
end