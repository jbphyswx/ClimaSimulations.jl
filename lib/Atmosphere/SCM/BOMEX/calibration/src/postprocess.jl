"""
    postprocess.jl

Rerunning selected ensemble members with the full diagnostics.


The selection metric is the one the EKP scheduler reacts to, `Φⱼ = ½ (gⱼ - y)ᵀ Γ⁻¹ (gⱼ - y)`.
"""

"""
    member_misfits(ekp, g_ensemble)

`Φⱼ` for each column of `g_ensemble`, `NaN` for a member that failed.
"""
function member_misfits(ekp, g_ensemble::AbstractMatrix)
    y = EKP.get_obs(ekp)
    Γ = Matrix(EKP.get_obs_noise_cov(ekp))
    factorized = LinearAlgebra.cholesky(LinearAlgebra.Symmetric(Γ))
    return map(axes(g_ensemble, 2)) do j
        d = view(g_ensemble, :, j) .- y
        any(isnan, d) ? NaN : 0.5 * LinearAlgebra.dot(d, factorized \ d)
    end
end

"""Where the observation map wrote `G_ensemble` for `iteration`."""
g_ensemble_path(output_dir::AbstractString, iteration::Integer) =
    joinpath(output_dir, "iteration_" * lpad(iteration, 3, '0'), "G_ensemble.jld2")

"""
    best_members(interface; last_iteration)

The `(iteration, member, misfit)` of the lowest-misfit member over all completed iterations, and
of the lowest within the last one. Iterations with no `G_ensemble` are skipped and failed members
ignored.
"""
function best_members(
    interface::BOMEXInterface;
    last_iteration::Integer = ClimaCalibrate.last_completed_iteration(
        interface.output_dir,
    ),
)
    last_iteration >= 1 ||
        error("No completed iterations in $(interface.output_dir); nothing to postprocess.")
    ekp = ClimaCalibrate.load_ekp_struct(interface.output_dir, last_iteration)
    candidates = NamedTuple{(:iteration, :member, :misfit), Tuple{Int, Int, Float64}}[]
    for iteration in 1:last_iteration
        path = g_ensemble_path(interface.output_dir, iteration)
        isfile(path) || continue
        for (member, φ) in enumerate(member_misfits(ekp, JLD2.load_object(path)))
            isnan(φ) || push!(candidates, (; iteration, member, misfit = φ))
        end
    end
    isempty(candidates) && error("Every ensemble member failed; there is nothing to rerun.")
    final = filter(c -> c.iteration == last_iteration, candidates)
    isempty(final) && error(
        "Every member of iteration $last_iteration failed, so there is no best final member.",
    )
    return (;
        best = argmin(c -> c.misfit, candidates),
        best_final = argmin(c -> c.misfit, final),
    )
end

"""
    rerun_member(interface, iteration, member; output_dir, vars, period_seconds)

Rerun one member with the full diagnostics, returning the directory written.

Uses the member's own `parameters.toml` and the grid the calibration used, so this reproduces
what the calibration ran.
"""
function rerun_member(
    interface::BOMEXInterface,
    iteration::Integer,
    member::Integer;
    output_dir::AbstractString,
    vars = (BOMEX.default_diagnostic_vars..., BOMEX.STATE_VARS...),
    period_seconds::Real = 600,
)
    parameter_file =
        ClimaCalibrate.parameter_path(interface.output_dir, iteration, member)
    isfile(parameter_file) ||
        error("No parameters.toml for iteration $iteration member $member: $parameter_file")
    return BOMEX.run_case(
        interface.case;
        FT = interface.float_type,
        grid = interface.grid,
        params = member_params(parameter_file; FT = interface.float_type),
        output_dir,
        verbose = false,
        diagnostics = BOMEX.bomex_diagnostics(
            vars;
            period_seconds,
            n_levels = length(BOMEX.bomex_z(interface.grid)),
        ),
        interface.run_kwargs...,
    )
end

"""
    postprocess_best_members(interface; output_dir, last_iteration, period_seconds)

Rerun the best and best-final members with full diagnostics, returning a `Dict` from label to
`(; iteration, member, misfit, dir)`.

Need to still add proper functioanlity, median/mean/best prior, meadian/mean/bset posterior, etc
"""
function postprocess_best_members(
    interface::BOMEXInterface;
    output_dir::AbstractString = joinpath(interface.output_dir, "postprocess"),
    last_iteration::Integer = ClimaCalibrate.last_completed_iteration(
        interface.output_dir,
    ),
    period_seconds::Real = 600,
)
    selected = best_members(interface; last_iteration)
    results = Dict{String, Any}()
    dirs = Dict{Tuple{Int, Int}, String}()
    for label in ("best", "best_final")
        pick = getproperty(selected, Symbol(label))
        key = (pick.iteration, pick.member)
        dir = get(dirs, key, nothing)
        if isnothing(dir)
            @info "Rerunning $label member with full diagnostics" pick.iteration pick.member pick.misfit
            dir = rerun_member(
                interface,
                pick.iteration,
                pick.member;
                output_dir = joinpath(output_dir, label),
                period_seconds,
            )
            dirs[key] = dir
        else
            @info "$label is the same member; reusing its rerun" pick.iteration pick.member dir
        end
        results[label] = (; pick.iteration, pick.member, pick.misfit, dir)
    end
    return results
end
