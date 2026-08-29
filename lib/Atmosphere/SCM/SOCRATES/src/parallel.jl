"""
Executing independent case runs.

`run_tasks` is the only task-scheduling mechanism here: a case sweep and a calibration
iteration are the same call with a different task list.
"""

using Distributed: Distributed

"""How a list of independent tasks is executed. Concrete: [`SerialExecutor`](@ref), [`WorkerPoolExecutor`](@ref)."""
abstract type AbstractExecutor end

"""Run tasks one at a time in the current process."""
struct SerialExecutor <: AbstractExecutor end

"""Seconds to wait for a worker, with nothing in flight, before erroring."""
const DEFAULT_EMPTY_POOL_TIMEOUT = 7200

"""
    WorkerPoolExecutor(pool; empty_pool_timeout, startup_waves, startup_wave_pause)

Run tasks across `pool`, handing each worker a new task as it frees up.

`empty_pool_timeout` bounds the wait for a worker while nothing is in flight; workers may join
after the run starts, so an empty pool with work in flight is not a stall.

`startup_waves` groups the workers by order of first dispatch and delays each group's first
task by `(wave - 1) * startup_wave_pause` seconds, spreading the memory peak of a first run
across time. Tasks after a worker's first are never delayed.
"""
struct WorkerPoolExecutor{P <: Distributed.AbstractWorkerPool} <: AbstractExecutor
    pool::P
    empty_pool_timeout::Int
    startup_waves::Int
    startup_wave_pause::Float64
end

WorkerPoolExecutor(
    pool::Distributed.AbstractWorkerPool;
    empty_pool_timeout::Integer = DEFAULT_EMPTY_POOL_TIMEOUT,
    startup_waves::Integer = 0,
    startup_wave_pause::Real = 0,
) = WorkerPoolExecutor(
    pool,
    Int(empty_pool_timeout),
    Int(startup_waves),
    Float64(startup_wave_pause),
)

"""
    topology_key(grid)

The affinity key of a run on `grid`: in this case, its number of vertical levels
"""
topology_key(grid) = length(socrates_z(grid))

function _take_warm_first!(pending::Vector{Int}, affinity_keys, warm)
    if !isempty(warm)
        at = findfirst(i -> affinity_keys[i] in warm, pending)
        isnothing(at) || return popat!(pending, at)
    end
    return popfirst!(pending)
end

"""Seconds the `ordinal`-th worker to be dispatched waits before its first task."""
_startup_delay(ordinal::Integer, startup_waves::Integer, startup_wave_pause::Real) =
    startup_waves > 0 ? ((ordinal - 1) % startup_waves) * Float64(startup_wave_pause) : 0.0

_check_keys(affinity_keys, tasks) =
    length(affinity_keys) == length(tasks) || error(
        "$(length(affinity_keys)) affinity keys for $(length(tasks)) tasks; there must be one \
         per task.",
    )

"""
    run_tasks(f, tasks, [affinity_keys], executor)

Apply `f` to each of `tasks`, returning results in `tasks` order. `f` runs on a worker, so it
may only close over serializable data. A task that throws is logged and its result is
`nothing`; the remaining tasks still run.

`affinity_keys[i]` is the affinity key of `tasks[i]`: a worker that has completed a key is
preferred for the remaining tasks carrying it, so compiled specializations are reused rather
than rediscovered. A worker with no warm task pending takes the next task regardless, so the
keys affect neither which tasks run nor the order of the results.
"""
function run_tasks(f, tasks, affinity_keys, ::SerialExecutor)
    _check_keys(affinity_keys, tasks)
    results = Vector{Any}(nothing, length(tasks))
    for (i, task) in enumerate(tasks)
        try
            results[i] = f(task)
        catch e
            @error "Task failed" task exception = (e, catch_backtrace())
        end
    end
    return results
end

function run_tasks(f, tasks, affinity_keys, executor::WorkerPoolExecutor)
    _check_keys(affinity_keys, tasks)
    (; pool, empty_pool_timeout, startup_waves, startup_wave_pause) = executor
    results = Vector{Any}(nothing, length(tasks))
    pending = collect(eachindex(tasks))
    warm = Dict{Int, Set{Any}}()
    started = Set{Int}()
    guard = ReentrantLock()
    inflight = Threads.Atomic{Int}(0)
    t_last_available = time()
    @sync while !isempty(pending)
        if isempty(pool.workers)
            inflight[] > 0 && (t_last_available = time())
            waited = time() - t_last_available
            inflight[] == 0 &&
                waited > empty_pool_timeout &&
                error(
                    "No workers available for $(round(Int, waited)) s (timeout \
                     $(empty_pool_timeout) s) with nothing in flight.",
                )
            sleep(1)
            continue
        end
        t_last_available = time()
        worker = take!(pool)
        i, delay = lock(guard) do
            index = _take_warm_first!(
                pending,
                affinity_keys,
                get!(() -> Set{Any}(), warm, worker),
            )
            worker in started && return (index, 0.0)
            push!(started, worker)
            return (
                index,
                _startup_delay(length(started), startup_waves, startup_wave_pause),
            )
        end
        Threads.atomic_add!(inflight, 1)
        @async try
            delay > 0 && Distributed.remotecall_fetch(sleep, worker, delay)
            results[i] = Distributed.remotecall_fetch(f, worker, tasks[i])
            lock(guard) do
                push!(warm[worker], affinity_keys[i])
            end
        catch e
            @error "Task failed on worker $worker" task = tasks[i] exception = e
        finally
            Threads.atomic_sub!(inflight, 1)
            put!(pool, worker)
        end
    end
    return results
end

run_tasks(f, tasks, executor::AbstractExecutor = SerialExecutor()) =
    run_tasks(f, tasks, fill(nothing, length(tasks)), executor)

"""
    run_cases(cases; FT, output_dir, executor, grids, kwargs...)

Run each case into its own subdirectory of `output_dir`, returning each run's directory — or
`nothing` for a case that failed — in `cases` order.
"""
function run_cases(
    cases::AbstractVector{<:SOCRATESCase};
    FT::Type{<:AbstractFloat} = Float64,
    output_dir::AbstractString,
    executor::AbstractExecutor = SerialExecutor(),
    grids::AbstractVector = [socrates_grid(FT, c) for c in cases],
    kwargs...,
)
    run_one =
        ((c, grid),) -> run_case(
            c;
            FT,
            grid,
            output_dir = joinpath(output_dir, case_name(c)),
            kwargs...,
        )
    return run_tasks(run_one, collect(zip(cases, grids)), map(topology_key, grids), executor)
end