"""
    parameter_evolution(ekp, prior; path, logscale)

Every member's parameter value against iteration in physical units, one panel per parameter,
with the ensemble mean as a thick line.
"""
function parameter_evolution(ekp, prior; path::AbstractString, logscale::Bool = false)
    names = EKP.ParameterDistributions.get_name(prior)
    phi = EKP.get_ϕ(prior, ekp)
    n_iter = length(phi)
    n_par = length(names)
    fig = CairoMakie.Figure(; size = (450 * min(n_par, 2), 320 * cld(n_par, 2)))
    for (i, name) in enumerate(names)
        row, col = fldmod1(i, 2)
        ax = CairoMakie.Axis(
            fig[row, col];
            title = name,
            xlabel = "iteration",
            ylabel = "value",
            yscale = logscale ? log10 : identity,
        )
        values = [phi[k][i, :] for k in 1:n_iter]
        for member in eachindex(first(values))
            CairoMakie.lines!(
                ax, 1:n_iter, [v[member] for v in values];
                color = (:grey, 0.5), linewidth = 1,
            )
        end
        CairoMakie.lines!(
            ax, 1:n_iter, [Statistics.mean(v) for v in values];
            color = :firebrick, linewidth = 3, label = "ensemble mean",
        )
        i == 1 && CairoMakie.axislegend(ax; position = :rt, framevisible = false)
    end
    return save_figure(fig, path)
end

"""
    error_evolution(ekp; path)

The EKP misfit and the accumulated algorithmic time `T = ΣΔt` against iteration. `nothing`
before the first ensemble update, which is what records the misfit.
"""
function error_evolution(ekp; path::AbstractString)
    haskey(EKP.get_error_metrics(ekp), "loss") || begin
        @warn "No misfit recorded yet (no completed ensemble update); skipping" path
        return nothing
    end
    err = EKP.get_error(ekp)
    dt = EKP.get_Δt(ekp)
    isempty(err) && (@warn "Misfit history is empty; skipping" path; return nothing)
    fig = CairoMakie.Figure(; size = (900, 350))
    ax1 = CairoMakie.Axis(
        fig[1, 1];
        title = "misfit", xlabel = "iteration", ylabel = "error", yscale = log10,
    )
    CairoMakie.scatterlines!(ax1, 1:length(err), err; color = :firebrick, linewidth = 2)
    ax2 = CairoMakie.Axis(
        fig[1, 2];
        title = "accumulated T = ΣΔt", xlabel = "iteration", ylabel = "T", yscale = log10,
    )
    CairoMakie.scatterlines!(
        ax2, 1:length(dt), cumsum(dt); color = :steelblue, linewidth = 2,
    )
    return save_figure(fig, path)
end

"""
    profile_grid(z, reference_by_var, groups_by_var, vars; path, xlabels, title, reference_label, ncols)

One panel per variable: the reference as a thick black line, and for each `label => members`
group in `groups_by_var[var]` every member thin with the group mean thick.

`z`, the reference and every member must already be restricted to the scored levels and
averaged over the scoring window, so what is drawn is what is scored.
"""
function profile_grid(
    z,
    reference_by_var::AbstractDict,
    groups_by_var::AbstractDict,
    vars;
    path::AbstractString,
    xlabels::AbstractDict = Dict{String, String}(),
    title::AbstractString = "",
    reference_label::AbstractString = "Atlas LES",
    ncols::Int = 2,
)
    n = length(vars)
    fig = CairoMakie.Figure(; size = (460 * ncols, 340 * cld(n, ncols)))
    for (i, var) in enumerate(vars)
        row, col = fldmod1(i, ncols)
        ax = CairoMakie.Axis(
            fig[row, col];
            title = var, xlabel = get(xlabels, var, ""), ylabel = "z [m]",
        )
        for (g, (label, members)) in enumerate(groups_by_var[var])
            color = GROUP_COLORS[mod1(g, length(GROUP_COLORS))]
            for m in members
                CairoMakie.lines!(ax, m, z[var]; color = (color, 0.35), linewidth = 1)
            end
            isempty(members) && continue
            mean_profile =
                [Statistics.mean(m[k] for m in members) for k in eachindex(z[var])]
            CairoMakie.lines!(
                ax, mean_profile, z[var]; color, linewidth = 3, label = String(label),
            )
        end
        CairoMakie.lines!(
            ax, reference_by_var[var], z[var];
            color = :black, linewidth = 3.5, label = reference_label,
        )
        i == 1 && CairoMakie.axislegend(ax; position = :rt, framevisible = false)
    end
    isempty(title) || CairoMakie.Label(fig[0, :], title; fontsize = 18)
    return save_figure(fig, path)
end

"""
    scalar_comparison(labels, reference, groups; path, title, ylabel, reference_label)

Reference against each group's members for the scalar diagnostics named by `labels`.
"""
function scalar_comparison(
    labels::AbstractVector,
    reference::AbstractVector,
    groups::AbstractVector;
    path::AbstractString,
    title::AbstractString = "",
    ylabel::AbstractString = "",
    reference_label::AbstractString = "Atlas LES",
)
    fig = CairoMakie.Figure(; size = (620, 380))
    ax = CairoMakie.Axis(
        fig[1, 1];
        title, ylabel, yscale = log10,
        xticks = (1:length(labels), collect(String.(labels))),
    )
    for (g, (label, members)) in enumerate(groups)
        color = GROUP_COLORS[mod1(g, length(GROUP_COLORS))]
        for m in members
            CairoMakie.scatter!(
                ax, 1:length(labels), max.(m, eps());
                color = (color, 0.4), markersize = 7,
            )
        end
        isempty(members) && continue
        mean_values =
            [max(Statistics.mean(m[i] for m in members), eps()) for i in eachindex(labels)]
        CairoMakie.scatter!(
            ax, 1:length(labels), mean_values;
            color, markersize = 16, marker = :diamond, label = String(label),
        )
    end
    CairoMakie.scatter!(
        ax, 1:length(labels), max.(reference, eps());
        color = :black, markersize = 16, marker = :star5, label = reference_label,
    )
    CairoMakie.axislegend(ax; position = :rt, framevisible = false)
    return save_figure(fig, path)
end

"""
    observation_blocks(ekp)

One `(short_name, range, z, scale)` per scored block. `range` indexes the observation vector,
`z` is empty for a column integral, and `scale` is `sqrt(pool_var)` — observations are stored
divided by it, so multiplying returns physical units.
"""
function observation_blocks(ekp)
    md = EKP.get_metadata(first(EKP.get_observations(EKP.get_observation_series(ekp))))
    blocks = Tuple{String, UnitRange{Int}, Vector{Float64}, Float64}[]
    offset = 0
    for m in md
        dims = collect(values(m.dims))
        n = isempty(dims) ? 1 : prod(length(v) for v in dims)
        z = isempty(dims) ? Float64[] : collect(Float64, first(dims))
        scale = sqrt(parse(Float64, string(get(m.attributes, "pool_var", "1.0"))))
        push!(
            blocks,
            (get(m.attributes, "short_name", "?"), (offset + 1):(offset + n), z, scale),
        )
        offset += n
    end
    return blocks
end

"""
    prior_posterior_profiles(ekp, prior_g, posterior_g; path, n_sigma, ncols, reference_label, xlim_factor)

Reference against the prior and posterior ensembles, one panel per scored profile, with a band
at each multiple of the ensemble standard deviation in `n_sigma`.

`prior_g` and `posterior_g` are `(n_obs, n_members)` `G_ensemble` matrices. Each block is
scaled by its `pool_var`, so panels read in physical units rather than score units.
"""
function prior_posterior_profiles(
    ekp,
    prior_g::AbstractMatrix,
    posterior_g::AbstractMatrix;
    path::AbstractString,
    n_sigma = (1, 2),
    ncols::Int = 4,
    reference_label::AbstractString = "Atlas LES",
    xlim_factor::Union{Nothing, Real} = nothing,
)
    y = EKP.get_obs(ekp)
    blocks = filter(b -> !isempty(b[3]), observation_blocks(ekp))
    isempty(blocks) &&
        error("No profile variables to plot; every scored block is a scalar.")
    nrows = cld(length(blocks), ncols)
    fig = CairoMakie.Figure(; size = (420 * ncols, 300 * nrows))
    for (i, (name, range, z, scale)) in enumerate(blocks)
        row, col = fldmod1(i, ncols)
        ax = CairoMakie.Axis(
            fig[row, col];
            title = name, xlabel = "kg/kg", ylabel = col == 1 ? "z [m]" : "",
        )
        for (g, color, label) in
            ((prior_g, :firebrick, "prior"), (posterior_g, :steelblue, "posterior"))
            members = view(g, range, :)
            keep = [j for j in axes(members, 2) if !any(isnan, view(members, :, j))]
            isempty(keep) && continue
            μ = Statistics.mean(view(members, :, keep); dims = 2)[:] .* scale
            σ = Statistics.std(view(members, :, keep); dims = 2)[:] .* scale
            # Widest band first, so the narrower ones stay visible over it.
            for k in sort(collect(n_sigma); rev = true)
                CairoMakie.band!(
                    ax,
                    CairoMakie.Point2f.(μ .- k .* σ, z),
                    CairoMakie.Point2f.(μ .+ k .* σ, z);
                    color = (color, k == minimum(n_sigma) ? 0.30 : 0.15),
                )
            end
            CairoMakie.lines!(ax, μ, z; color, linewidth = 3, label)
        end
        reference = y[range] .* scale
        CairoMakie.lines!(
            ax, reference, z; color = :black, linewidth = 4, label = reference_label,
        )
        if !isnothing(xlim_factor)
            span = maximum(abs, reference)
            span > 0 &&
                CairoMakie.xlims!(ax, -0.1 * xlim_factor * span, xlim_factor * span)
        end
        i == 1 && CairoMakie.axislegend(ax; position = :rt, framevisible = false)
    end
    return save_figure(fig, path)
end

"""
    member_values(interface, iteration, case, levels; vars)

Every member's scored values for one case at one iteration, as
`MOSAiCAYiLCalibration.calibrated_values` returns them. Stops at the first member whose output
is absent.
"""
function member_values(
    interface::MC.MOSAiCInterface,
    iteration::Integer,
    c::MOSAiC_AYiL.MOSAiCAYiLCase,
    levels::AbstractVector;
    vars = interface.vars,
)
    window = MC.score_window(c)
    values = Dict{String, Vector{Float64}}[]
    member = 1
    while isdir(MC.case_output_dir(interface, iteration, member, c))
        dir = MC.case_output_dir(interface, iteration, member, c)
        push!(
            values,
            MC.calibrated_values(c, levels, MC.model_fetch(dir, levels); vars, window),
        )
        member += 1
    end
    return values
end

"""
    calibration_figures(interface, prior; output_dir, iterations, logscale)

Write every figure a finished calibration supports, returning the paths written.

The parameter and misfit figures need only the EKP files. The per-case figures read each
member's run output and are skipped with a warning when it is not there. A scored variable
whose values are as long as `levels` is drawn as a profile; the rest go into one scalar panel
per case.
"""
function calibration_figures(
    interface::MC.MOSAiCInterface,
    prior;
    output_dir::AbstractString = joinpath(interface.output_dir, "figures"),
    iterations = nothing,
    logscale::Bool = false,
)
    last_iteration = ClimaCalibrate.last_completed_iteration(interface.output_dir)
    last_iteration >= 1 ||
        error("No completed iterations in $(interface.output_dir); nothing to plot.")
    ekp = ClimaCalibrate.load_ekp_struct(interface.output_dir, last_iteration)
    wanted = isnothing(iterations) ? unique((1, last_iteration)) : collect(iterations)
    @info "plotting calibration" last_iteration n_ens = EKP.get_N_ens(ekp) output_dir

    written = String[]
    push!(
        written,
        parameter_evolution(
            ekp, prior; path = joinpath(output_dir, "parameters.png"), logscale,
        ),
    )
    err = error_evolution(ekp; path = joinpath(output_dir, "error.png"))
    isnothing(err) || push!(written, err)

    for c in interface.cases
        name = MOSAiC_AYiL.case_name(c)
        window = MC.score_window(c)
        levels = MC.case_levels(interface, c)
        have = [it for it in wanted if isdir(MC.case_output_dir(interface, it, 1, c))]
        if isempty(have)
            @warn "No member output for $name; skipping its per-case figures."
            continue
        end
        reference = MC.calibrated_values(
            c, levels, MC.reference_fetch(c, levels); vars = interface.vars, window,
        )
        by_iteration = Dict(it => member_values(interface, it, c, levels) for it in have)

        profile_vars = [v for v in interface.vars if length(reference[v]) == length(levels)]
        scalar_vars = [v for v in interface.vars if length(reference[v]) != length(levels)]

        if !isempty(profile_vars)
            z = Dict(v => levels for v in profile_vars)
            reference_by_var = Dict(v => reference[v] for v in profile_vars)
            groups_by_var = Dict(
                v => ["iteration $it" => [m[v] for m in by_iteration[it]] for it in have]
                for v in profile_vars
            )
            push!(
                written,
                profile_grid(
                    z, reference_by_var, groups_by_var, profile_vars;
                    path = joinpath(output_dir, "$(name)_profiles.png"),
                    reference_label = "DALES",
                    title = "$name — scored levels, averaged over $(round.(Int, window)) s",
                ),
            )
        end

        if !isempty(scalar_vars)
            push!(
                written,
                scalar_comparison(
                    scalar_vars,
                    [only(reference[v]) for v in scalar_vars],
                    [
                        "iteration $it" =>
                            [[only(m[v]) for v in scalar_vars] for m in by_iteration[it]]
                        for it in have
                    ];
                    path = joinpath(output_dir, "$(name)_scalars.png"),
                    title = "$name scored scalars",
                    reference_label = "DALES",
                ),
            )
        end
    end
    return written
end

