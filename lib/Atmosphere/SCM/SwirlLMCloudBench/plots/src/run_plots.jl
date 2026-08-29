"""
    run_profiles(dir, vars; path, period, ncols)

One time-averaged profile panel per variable from the run in `dir`.

Drawn with `ClimaAnalysis.Visualize.line_plot1D!`, so each panel takes its labels and units from
the diagnostic itself.
"""
function run_profiles(
    dir::AbstractString,
    vars = CB.default_diagnostic_vars;
    path::AbstractString,
    period::Union{String, Nothing} = nothing,
    ncols::Int = 3,
)
    outputs = CB.run_outputvars(dir, vars; period)
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
    vars = CB.default_diagnostic_vars;
    path::AbstractString,
    period::Union{String, Nothing} = nothing,
    ncols::Int = 2,
)
    outputs = CB.run_outputvars(dir, vars; period)
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
    run_vs_reference(dir, case, vars; path, grid, period, ncols)

The run in `dir` against that case's published output, one panel per scored quantity, both sides
put on `grid`'s levels and reduced by `SwirlLMCloudBenchCalibration.scored_rows` — the same
reduction the calibration scores with, so a panel shows what a residual is computed from.

A height-resolved quantity is drawn as a profile and a scalar as a horizontal line. Reads the
published store, so the case's `data.zarr` has to be reachable.
"""
function run_vs_reference(
    dir::AbstractString,
    case,
    vars = CBC.DEFAULT_SCORED_VARS;
    path::AbstractString,
    grid,
    period::Union{String, Nothing} = nothing,
    ncols::Int = 3,
)
    names = collect(String, vars)
    levels = collect(Float64, CB.cloudbench_z(grid))
    needed = CBC.scored_model_vars(
        Dict(n => CBC.REFERENCE_PAIRINGS[n] for n in names),
    )
    run_vars = CB.run_outputvars(dir, needed; period)
    values = Dict{String, Any}(n => CBC.altitude_first(v) for (n, v) in run_vars)
    z_model = CBC.model_altitude(run_vars)
    ref = CBC.reference_comparable(case; names)

    fig = CairoMakie.Figure(; size = (400 * ncols, 340 * cld(length(names), ncols)))
    for (i, name) in enumerate(names)
        row, col = fldmod1(i, ncols)
        model = CBC.scored_rows(
            CBC.model_comparable(values, name), name; z_src = z_model, levels,
        )
        reference = CBC.scored_rows(ref.data[name], name; z_src = ref.z, levels)
        ax = CairoMakie.Axis(
            fig[row, col];
            title = name,
            xlabel = CBC.REFERENCE_PAIRINGS[name].units,
            ylabel = col == 1 ? "z [m]" : "",
        )
        if size(model, 1) == length(levels)
            CairoMakie.lines!(
                ax, vec(Statistics.mean(model; dims = 2)), levels;
                color = :steelblue, linewidth = 3, label = "model",
            )
            CairoMakie.lines!(
                ax, vec(Statistics.mean(reference; dims = 2)), levels;
                color = :black, linewidth = 3.5, label = "CloudBench LES",
            )
        else
            CairoMakie.hlines!(
                ax, only(Statistics.mean(model; dims = 2));
                color = :steelblue, linewidth = 3, label = "model",
            )
            CairoMakie.hlines!(
                ax, only(Statistics.mean(reference; dims = 2));
                color = :black, linewidth = 3.5, label = "CloudBench LES",
            )
        end
        i == 1 && CairoMakie.axislegend(ax; position = :rt, framevisible = false)
    end
    CairoMakie.Label(fig[0, :], CBC.case_name(case); fontsize = 16)
    return save_figure(fig, path)
end
