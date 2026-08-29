"""
    run_calibration.jl

Ensemble Kalman Inversion of the BOMEX microphysics and entrainment parameters.

```julia
include("calibration/examples/run_calibration.jl")

run_calibration(; reference = "my_reference.nc")
run_calibration(; reference, nworkers = 16)         # more members in flight
run_calibration(; reference, nworkers = 0)          # use the existing pool as it is
run_calibration(; reference, startup_waves = 4, startup_wave_pause = 60)
```

`reference` is required and has no default: BOMEX ships no observational reference, so what the
ensemble is scored against is yours to supply. It may be a NetCDF path or a `Dict` of profiles —
see `BOMEXCalibration.read_reference`.

Each member is a full ClimaAtmos column, so expect several GB per worker.
"""

using BOMEX: BOMEX
using BOMEXCalibration: BOMEXCalibration as BC
using ClimaCalibrate: ClimaCalibrate
using Distributed: Distributed

# Recorded, not just broadcast, so workers that connect later also load the model code.
ClimaCalibrate.@worker_setup begin
    using BOMEX: BOMEX
    using BOMEXCalibration: BOMEXCalibration
end

"""
    worker_pool(; nworkers, workers, device, cluster, kwargs...)

The pool the calibration runs on.

Pass `workers` to use a specific set of existing workers, `nworkers = 0` to take the pool as it
already is, or `nworkers > 0` to add that many.
"""
function worker_pool(;
    nworkers::Int = 4,
    workers = nothing,
    device::Symbol = :cpu,
    cluster::Symbol = :auto,
    kwargs...,
)
    if !isnothing(workers)
        ids = collect(Int, workers)
        live = Distributed.workers()
        missing_ids = setdiff(ids, live)
        isempty(missing_ids) ||
            error("workers $missing_ids are not running; live workers are $live")
        return Distributed.WorkerPool(ids)
    end
    nworkers > 0 && wait(ClimaCalibrate.add_workers(nworkers; device, cluster, kwargs...))
    pool = ClimaCalibrate.default_worker_pool()
    isempty(pool.workers) && error(
        "no workers available; pass `nworkers` to add some or `workers` to name existing ones",
    )
    return pool
end

function run_calibration(;
    reference,
    output_dir = joinpath(@__DIR__, "output", "calibration_demo"),
    ensemble_size::Int = 4,
    n_iterations::Int = 2,
    startup_waves::Integer = 0,
    startup_wave_pause::Real = 0,
    interface_kwargs = (;),
    ekp_kwargs = (;),
    pool_kwargs...,
)
    pool = worker_pool(; pool_kwargs...)
    @info "calibrating BOMEX" ensemble_size n_iterations workers =
        collect(pool.workers) output_dir

    interface = BC.BOMEXInterface(;
        reference,
        output_dir,
        startup_waves,
        startup_wave_pause,
        interface_kwargs...,
    )
    prior = BC.default_prior()
    ekp = BC.build_ekp(interface, prior; ensemble_size, ekp_kwargs...)
    backend = ClimaCalibrate.WorkerBackend(; worker_pool = pool)

    return ClimaCalibrate.calibrate(
        backend, ekp, interface, n_iterations, prior, interface.output_dir,
    )
end
