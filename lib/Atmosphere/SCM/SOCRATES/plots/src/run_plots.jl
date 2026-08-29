"""
    run_profiles(dir, vars; path, period, ncols)

One time-averaged profile panel per variable from the run in `dir`.

Drawn with `ClimaAnalysis.Visualize.line_plot1D!`, so each panel takes its labels and units from
the diagnostic itself.
"""
function run_profiles(
    dir::AbstractString,
    vars = SOCRATES.default_diagnostic_vars;
    path::AbstractString,
    period::Union{String, Nothing} = nothing,
    ncols::Int = 3,
)
    outputs = SOCRATES.run_outputvars(dir, vars; period)
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
    vars = SOCRATES.default_diagnostic_vars;
    path::AbstractString,
    period::Union{String, Nothing} = nothing,
    ncols::Int = 2,
)
    outputs = SOCRATES.run_outputvars(dir, vars; period)
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

The run in `dir` against that case's Atlas LES, one time-averaged profile panel per variable,
both restricted to the levels the case is scored on.
"""
function run_vs_reference(
    dir::AbstractString,
    c::SOCRATES.SOCRATESCase,
    vars = ("clw", "cli", "husra", "hussn");
    path::AbstractString,
    period::Union{String, Nothing} = nothing,
    ncols::Int = 2,
)
    window = SC.score_window(c)
    outputs = SOCRATES.run_outputvars(dir, vars; period)
    reference = SC.les_outputvars(c)
    names = collect(String, vars)
    levels = SC.scored_levels(
        outputs[first(names)].dims[ClimaAnalysis.altitude_name(outputs[first(names)])],
        SC.default_z_bounds(c),
    )
    fig = CairoMakie.Figure(; size = (440 * ncols, 340 * cld(length(names), ncols)))
    for (i, name) in enumerate(names)
        row, col = fldmod1(i, ncols)
        model = SC.windowed_time_mean(SC.restrict_to_levels(outputs[name], levels), window)
        les = SC.windowed_time_mean(SC.reference_on_levels(reference[name], levels), window)
        ax = CairoMakie.Axis(
            fig[row, col];
            title = name,
            xlabel = ClimaAnalysis.units(outputs[name]),
            ylabel = col == 1 ? "z [m]" : "",
        )
        CairoMakie.lines!(
            ax, vec(Array{Float64}(model.data)), levels;
            color = :steelblue, linewidth = 3, label = "model",
        )
        CairoMakie.lines!(
            ax, vec(Array{Float64}(les.data)), levels;
            color = :black, linewidth = 3.5, label = "Atlas LES",
        )
        i == 1 && CairoMakie.axislegend(ax; position = :rt, framevisible = false)
    end
    CairoMakie.Label(
        fig[0, :],
        "$(SOCRATES.case_name(c)) — averaged over $(round.(Int, window)) s";
        fontsize = 16,
    )
    return save_figure(fig, path)
end
