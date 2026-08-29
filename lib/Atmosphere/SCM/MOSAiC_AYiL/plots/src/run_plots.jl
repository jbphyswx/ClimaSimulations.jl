"""
    run_profiles(dir, vars; path, period, ncols)

One time-averaged profile panel per variable from the run in `dir`.

Drawn with `ClimaAnalysis.Visualize.line_plot1D!`, so each panel takes its labels and units from
the diagnostic itself.
"""
function run_profiles(
    dir::AbstractString,
    vars = MOSAiC_AYiL.DEFAULT_DIAGNOSTIC_VARS;
    path::AbstractString,
    period::Union{String, Nothing} = nothing,
    ncols::Int = 3,
)
    outputs = MOSAiC_AYiL.run_outputvars(dir, vars; period)
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
    vars = MOSAiC_AYiL.DEFAULT_DIAGNOSTIC_VARS;
    path::AbstractString,
    period::Union{String, Nothing} = nothing,
    ncols::Int = 2,
)
    outputs = MOSAiC_AYiL.run_outputvars(dir, vars; period)
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
    run_vs_reference(dir, case, vars; path, period, ncols)

The run in `dir` against that day's DALES reference, one panel per variable, both put on the
levels and window the day is scored over by `MOSAiCAYiLCalibration.calibrated_values`.

A scored water path is a single number, so it is drawn as a horizontal line across the panel.
"""
function run_vs_reference(
    dir::AbstractString,
    c::MOSAiC_AYiL.MOSAiCAYiLCase,
    vars = MC.default_calibration_vars;
    path::AbstractString,
    period::Union{String, Nothing} = nothing,
    ncols::Int = 2,
)
    window = MC.score_window(c)
    levels = MC.scored_levels(
        MOSAiC_AYiL.mosaic_z(MOSAiC_AYiL.mosaic_grid(Float64, c)), MC.z_bounds(c),
    )
    reference =
        MC.calibrated_values(c, levels, MC.reference_fetch(c, levels); vars, window)
    model = MC.calibrated_values(
        c, levels, MC.model_fetch(dir, levels; period); vars, window,
    )
    names = collect(String, vars)
    fig = CairoMakie.Figure(; size = (440 * ncols, 340 * cld(length(names), ncols)))
    for (i, name) in enumerate(names)
        row, col = fldmod1(i, ncols)
        ax = CairoMakie.Axis(
            fig[row, col]; title = name, ylabel = col == 1 ? "z [m]" : "",
        )
        if length(model[name]) == length(levels)
            CairoMakie.lines!(
                ax, model[name], levels; color = :steelblue, linewidth = 3, label = "model",
            )
            CairoMakie.lines!(
                ax, reference[name], levels;
                color = :black, linewidth = 3.5, label = "DALES",
            )
        else
            CairoMakie.hlines!(
                ax, only(model[name]); color = :steelblue, linewidth = 3, label = "model",
            )
            CairoMakie.hlines!(
                ax, only(reference[name]);
                color = :black, linewidth = 3.5, label = "DALES",
            )
        end
        i == 1 && CairoMakie.axislegend(ax; position = :rt, framevisible = false)
    end
    CairoMakie.Label(
        fig[0, :],
        "$(MOSAiC_AYiL.case_name(c)) — averaged over $(round.(Int, window)) s";
        fontsize = 16,
    )
    return save_figure(fig, path)
end
