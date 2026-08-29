"""
    postprocess.jl

Rerunning selected ensemble members with the full diagnostics.

A calibration writes only what is scored, so a finished one cannot say what drove the result.
This picks members out by their misfit and reruns them asking for
[`MOSAiC_AYiL.DEFAULT_DIAGNOSTIC_VARS`](@ref).

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
    interface::MOSAiCInterface;
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
    rerun_member(interface, iteration, member; output_dir, vars, executor)

Rerun every case for one member with the full diagnostics, into `output_dir/<case_name>`,
returning the directories in `interface.cases` order.

Uses the member's own `parameters.toml` and the grid the calibration used, so this reproduces
what the calibration ran.
"""
function rerun_member(
    interface::MOSAiCInterface,
    iteration::Integer,
    member::Integer;
    output_dir::AbstractString,
    vars = MOSAiC_AYiL.DEFAULT_DIAGNOSTIC_VARS,
    executor::MOSAiC_AYiL.AbstractExecutor = MOSAiC_AYiL.SerialExecutor(),
)
    params = ClimaCalibrate.parameter_path(interface.output_dir, iteration, member)
    isfile(params) ||
        error("No parameters.toml for iteration $iteration member $member: $params")
    run_one = c -> begin
        grid = case_grid(interface, c)
        MOSAiC_AYiL.run_case(
            c;
            FT = interface.float_type,
            params,
            grid,
            output_dir = joinpath(output_dir, MOSAiC_AYiL.case_name(c)),
            verbose = false,
            diagnostics = MOSAiC_AYiL.mosaic_diagnostics(
                vars; n_levels = length(MOSAiC_AYiL.mosaic_z(grid)),
            ),
            interface.run_kwargs...,
        )
    end
    affinity_keys =
        [MOSAiC_AYiL.topology_key(case_grid(interface, c)) for c in interface.cases]
    return MOSAiC_AYiL.run_tasks(run_one, interface.cases, affinity_keys, executor)
end

"""
    postprocess_best_members(interface; output_dir, last_iteration, executor)

Rerun the best and best-final members with full diagnostics, returning a `Dict` from label to
`(; iteration, member, misfit, dirs)`. Both are rerun even when they are the same member, so a
caller can always index either.
"""
function postprocess_best_members(
    interface::MOSAiCInterface;
    output_dir::AbstractString = joinpath(interface.output_dir, "postprocess"),
    last_iteration::Integer = ClimaCalibrate.last_completed_iteration(
        interface.output_dir,
    ),
    executor::MOSAiC_AYiL.AbstractExecutor = MOSAiC_AYiL.SerialExecutor(),
)
    selected = best_members(interface; last_iteration)
    results = Dict{String, Any}()
    for label in ("best", "best_final")
        pick = getproperty(selected, Symbol(label))
        @info "Rerunning $label member with full diagnostics" pick.iteration pick.member pick.misfit
        dirs = rerun_member(
            interface,
            pick.iteration,
            pick.member;
            output_dir = joinpath(output_dir, label),
            executor,
        )
        results[label] = (; pick.iteration, pick.member, pick.misfit, dirs)
    end
    return results
end
