"""
    run_profiles(dir, vars; path, period, ncols)

One time-averaged profile panel per variable from the run in `dir`.

Drawn with `ClimaAnalysis.Visualize.line_plot1D!`, so each panel takes its labels and units from
the diagnostic itself.
"""
function run_profiles(
    dir::AbstractString,
    vars = BOMEX.default_diagnostic_vars;
    path::AbstractString,
    period::Union{String, Nothing} = nothing,
    ncols::Int = 3,
)
    outputs = BOMEX.run_outputvars(dir, vars; period)
    names = collect(String, vars)
    fig = CairoMakie.Figure(; size = (380 * ncols, 330 * cld(length(names), ncols)))
    for (i, name) in enumerate(names)
        row, col = fldmod1(i, ncols)
        ClimaAnalysis.Visualize.line_plot1D!(
            fig,
            ClimaAnalysis.average_time(outputs[name]);
            p_loc = (row, col),
            more_kwargs = Dict(:axis => Dict(:dim_on_y => true), :plot => Dict()),
        )
    end
    return save_figure(fig, path)
end

"""
    run_evolution(dir, vars; path, period, ncols)

One height-time heatmap per variable from the run in `dir`, so the whole record is visible
rather than a single average.
"""
function run_evolution(
    dir::AbstractString,
    vars = BOMEX.default_diagnostic_vars;
    path::AbstractString,
    period::Union{String, Nothing} = nothing,
    ncols::Int = 2,
)
    outputs = BOMEX.run_outputvars(dir, vars; period)
    names = [n for n in collect(String, vars) if length(outputs[n].dims) == 2]
    isempty(names) && error(
        "None of $(collect(vars)) is a height-time field; `run_evolution` has nothing to draw.",
    )
    fig = CairoMakie.Figure(; size = (460 * ncols, 340 * cld(length(names), ncols)))
    for (i, name) in enumerate(names)
        row, col = fldmod1(i, ncols)
        ClimaAnalysis.Visualize.heatmap2D!(fig, outputs[name]; p_loc = (row, col))
    end
    return save_figure(fig, path)
end

"""
    run_vs_reference(dir, reference; path, vars, levels, period, ncols)

The run in `dir` against a reference, one time-averaged profile panel per variable, both put on
`levels` by `BOMEXCalibration.reference_on_levels`.

`reference` is whatever `BOMEXCalibration.read_reference` accepts — this case ships no
observational reference, so it has to be supplied.
"""
function run_vs_reference(
    dir::AbstractString,
    reference;
    path::AbstractString,
    vars = BC.default_calibration_vars,
    levels,
    period::Union{String, Nothing} = nothing,
    ncols::Int = 2,
)
    names = collect(String, vars)
    model = BC.model_values(dir, names, levels; period)
    truth = BC.reference_values(reference, names, levels)
    fig = CairoMakie.Figure(; size = (440 * ncols, 340 * cld(length(names), ncols)))
    for (i, name) in enumerate(names)
        row, col = fldmod1(i, ncols)
        ax = CairoMakie.Axis(
            fig[row, col]; title = name, ylabel = col == 1 ? "z [m]" : "",
        )
        CairoMakie.lines!(
            ax, model[name], levels; color = :steelblue, linewidth = 3, label = "model",
        )
        CairoMakie.lines!(
            ax, truth[name], levels; color = :black, linewidth = 3.5, label = "reference",
        )
        i == 1 && CairoMakie.axislegend(ax; position = :rt, framevisible = false)
    end
    return save_figure(fig, path)
end
