"""
    run_calibration.jl

Ensemble Kalman Inversion of the SOCRATES microphysics parameters against the
Atlas LES.

```julia
include("examples/run_calibration.jl")

run_calibration()                              # 4 workers on whatever this is
run_calibration(; nworkers = 16)               # more members in flight
run_calibration(; nworkers = 1)                # one at a time, for a memory limit
run_calibration(; workers = [2, 3])            # reuse two workers already running
run_calibration(; nworkers = 0)                # use the existing pool as it is
run_calibration(; cluster = :slurm, device = :cpu, time = 120)
```

Each ensemble member is a full SOCRATES case, so a worker holds a whole
ClimaAtmos simulation: expect several GB per worker.
"""

using ClimaCalibrate: ClimaCalibrate
using Distributed: Distributed
using SOCRATES: SOCRATES
using SOCRATESCalibration: SOCRATESCalibration as SC

# Recorded, not just broadcast, so workers that connect later also load the model
# code. `@everywhere` would skip them.
ClimaCalibrate.@worker_setup begin
    using SOCRATES: SOCRATES
    using SOCRATESCalibration: SOCRATESCalibration
end

"""
    worker_pool(; nworkers, workers, device, cluster, kwargs...)

The pool the calibration runs on.

Pass `workers` to use a specific set of existing workers, `nworkers = 0` to take
the pool as it already is, or `nworkers > 0` to add that many. `device` and
`cluster` are forwarded to `ClimaCalibrate.add_workers`, which detects SLURM, PBS
or a local session.
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
    nworkers > 0 && wait(
        ClimaCalibrate.add_workers(nworkers; device, cluster, kwargs...),
    )
    pool = ClimaCalibrate.default_worker_pool()
    isempty(pool.workers) && error(
        "no workers available; pass `nworkers` to add some or `workers` to name \
         existing ones",
    )
    return pool
end

function run_calibration(;
    cases = [SOCRATES.case("RF09_Obs"), SOCRATES.case("RF10_Obs")],
    output_dir = joinpath(@__DIR__, "output", "calibration_demo"),
    ensemble_size::Int = 4,
    n_iterations::Int = 2,
    interface_kwargs = (;),
    ekp_kwargs = (;),
    pool_kwargs...,
)
    pool = worker_pool(; pool_kwargs...)
    @info "calibrating" cases = length(cases) ensemble_size n_iterations workers =
        collect(pool.workers) output_dir

    interface = SC.SocratesInterface(; cases, output_dir, interface_kwargs...)
    prior = SC.default_prior()
    ekp = SC.build_ekp(interface, prior; ensemble_size, ekp_kwargs...)
    backend = ClimaCalibrate.WorkerBackend(; worker_pool = pool)

    return SC.calibrate(backend, ekp, interface; prior, n_iterations)
end
