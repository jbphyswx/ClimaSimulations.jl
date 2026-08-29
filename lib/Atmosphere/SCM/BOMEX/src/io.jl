"""
    run_outputvars(output_dir, vars; period, reduction)

The named diagnostics of a finished run, as `ClimaAnalysis.OutputVar`s.

`period` and `reduction` default to `nothing`, letting `ClimaAnalysis` discover them; it
errors if a run wrote more than one
"""
function run_outputvars(
    output_dir::AbstractString,
    vars = default_diagnostic_vars;
    period::Union{String, Nothing} = nothing,
    reduction::Union{String, Nothing} = nothing,
)
    simdir = ClimaAnalysis.SimDir(active_output_dir(output_dir))
    return Dict{String, ClimaAnalysis.OutputVar}(
        name => ClimaAnalysis.get(simdir; short_name = name, reduction, period)
        for name in vars
    )
end

"""A run writes into `output_active` while it is going; afterwards `output_dir` is the tree."""
active_output_dir(dir::AbstractString) =
    isdir(joinpath(dir, "output_active")) ? joinpath(dir, "output_active") : dir