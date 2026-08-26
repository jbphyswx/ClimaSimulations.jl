"""
    run_days.jl

Run AYiL days as simulations — no calibration, no ensemble.

From the REPL, which is how it is meant to be used, since building one simulation
compiles ClimaAtmos and RRTMGP:

```julia
include("examples/run_days.jl")
run_days()                                    # one day at the reference resolution
run_days(; dates = ("20200503", "20200210"))  # several
run_days(; dz_min = 100, t_end = 3600)        # coarser and shorter
compare_day("20200503", out)                  # the misfit table against the LES
```

Or as a script:

```
julia --project=. examples/run_days.jl 20200503 20200210
```

Each day is truncated at its own [`best_simulation_top`], above which the reference's
ice is not reproducible, and written to its own directory under `output_dir`.
"""

using ClimaAtmos: ClimaAtmos as CA
using MOSAiC_AYiL: MOSAiC_AYiL as MA
using Printf: Printf

"""
    run_days(; dates, output_dir, FT, dz_min, t_end, dt, kwargs...)

Run each of `dates` and return `date => output directory`.

`dz_min` coarsens the reference's own grid, which is what a scientific comparison
does not want and a quick look does. `nothing` keeps every level. A day that fails is
reported and the rest still run, since one bad day should not lose the others.
"""
function run_simulation(;
    cases = MA.default_cases()
    output_dir::AbstractString = joinpath(pwd(), "ayil_runs"),
    FT::Type{<:AbstractFloat} = Float64,
    t_end = nothing,
    grid = MA.mosaic_grid(
        FT, c; faces = MA.coarsen_faces_to_dz_min(faces, 10),
    ),
    diagnostics = MA.default_mosaic_diagnostics(; n_levels),
    kwargs...,
)
    out = Dict{String, String}()
    failed = Dict{String, Any}()
    for case in cases
        n_levels = length(MA.mosaic_z(grid))
        dir = joinpath(abspath(output_dir), MA.case_name(c))
        @info "running $(MA.case_name(c))" levels = n_levels top =
            MA.best_simulation_top(c) output = dir
        try
            out[date] = MA.run_case(
                c;
                FT,
                output_dir = dir,
                grid,
                diagnosics,
                (isnothing(t_end) ? (;) : (; t_end))...,
                kwargs...,
            )
        catch e
            @error "$(MA.case_name(c)) failed" exception = (e, catch_backtrace())
            failed[date] = e
        end
    end
    isempty(failed) ||
        @warn "$(length(failed)) of $(length(dates)) days failed" dates = keys(failed)
    return out
end

"""
    compare_day(date, output_dir; kwargs...)

The misfit table for one finished run, against the DALES reference.

Loads the calibration side, because what is compared and how it is weighted are
calibration choices rather than properties of the case.
"""
function compare_day(date, output_dir::AbstractString; kwargs...)
    cal = joinpath(@__DIR__, "..", "calibration", "src")
    isdefined(Main, :compare_to_les) || for f in ("ayil_info.jl", "scoring.jl")
        Base.include(Main, joinpath(cal, f))
    end
    c = MA.case(date)
    comparison = Main.compare_to_les(c, output_dir; kwargs...)
    Main.print_comparison(comparison)
    return comparison
end

if abspath(PROGRAM_FILE) == @__FILE__
    dates = isempty(ARGS) ? ("20200503",) : Tuple(ARGS)
    run_days(; dates)
end