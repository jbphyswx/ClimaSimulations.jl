using EnsembleKalmanProcesses: EnsembleKalmanProcesses as EKP
using LinearAlgebra: LinearAlgebra
using Random: Random
using SOCRATESPlots: SOCRATESPlots as PL
using Test: Test

# Every figure here is drawn from data built in this file, so the assertions exercise the
# drawing itself rather than a calibration's output. A written PNG is the only evidence that a
# figure function works end to end; `isdefined` says nothing about whether it can render.

"""A profile with a single interior maximum, so a panel has something to show."""
profile(levels, amplitude) = amplitude .* exp.(-((levels .- 400.0) ./ 250.0) .^ 2)

# The do-block form removes the directory on the way out, rather than at process exit: these
# run on a long-lived REPL server that does not exit between runs.
mktempdir() do out
    Test.@testset "SOCRATESPlots" begin

        Test.@testset "profile_grid writes a figure" begin
            levels = collect(0.0:100.0:900.0)
            vars = ("clw", "cli")
            z = Dict(v => levels for v in vars)
            reference = Dict(v => profile(levels, 1.0e-4) for v in vars)
            groups = Dict(
                v => ["iteration 1" => [profile(levels, 1.0e-4 * (1 + 0.1k)) for k in 1:3]]
                for v in vars
            )
            path = PL.profile_grid(
                z, reference, groups, vars;
                path = joinpath(out, "profiles.png"),
                xlabels = Dict(v => "kg/kg" for v in vars),
                title = "synthetic",
            )
            Test.@test isfile(path)
            Test.@test filesize(path) > 1_000
        end

        Test.@testset "an empty member group still draws the reference" begin
            levels = collect(0.0:200.0:800.0)
            path = PL.profile_grid(
                Dict("clw" => levels),
                Dict("clw" => profile(levels, 1.0e-4)),
                Dict("clw" => ["iteration 1" => Vector{Float64}[]]),
                ("clw",);
                path = joinpath(out, "profiles_empty.png"),
            )
            Test.@test isfile(path)
        end

        Test.@testset "scalar_comparison writes a figure" begin
            path = PL.scalar_comparison(
                ["lwp", "iwp"],
                [1.0e-2, 2.0e-3],
                ["iteration 1" => [[1.1e-2, 2.2e-3], [0.9e-2, 1.8e-3]]];
                path = joinpath(out, "scalars.png"),
                title = "synthetic",
                ylabel = "kg/m^2",
            )
            Test.@test isfile(path)
            Test.@test filesize(path) > 1_000
        end

        Test.@testset "the EKI figures over a real EnsembleKalmanProcess" begin
            rng = Random.MersenneTwister(1234)
            constraint = EKP.ParameterDistributions.bounded(1.0e-2, 1.0e2)
            prior = EKP.ParameterDistributions.combine_distributions([
                EKP.ParameterDistributions.ParameterDistribution(
                    EKP.ParameterDistributions.Parameterized(
                        EKP.ParameterDistributions.Normal(0.0, 1.0),
                    ),
                    constraint,
                    name,
                ) for name in ("a", "b")
            ])
            n_ens = 6
            y = [1.0, 2.0, 3.0]
            observation = EKP.Observation(
                Dict(
                    "samples" => y,
                    "covariances" => LinearAlgebra.Diagonal(ones(length(y))),
                    "names" => "synthetic",
                ),
            )
            ekp = EKP.EnsembleKalmanProcess(
                EKP.construct_initial_ensemble(rng, prior, n_ens),
                observation,
                EKP.Inversion();
                rng,
                verbose = false,
            )

            # Documented behaviour: the misfit is recorded by the ensemble update, so there is
            # nothing to draw before one has happened.
            Test.@test isnothing(
                PL.error_evolution(ekp; path = joinpath(out, "error_before.png")),
            )
            Test.@test !isfile(joinpath(out, "error_before.png"))

            g = reduce(hcat, [y .+ 0.1 * j for j in 1:n_ens])
            EKP.update_ensemble!(ekp, g)

            params =
                PL.parameter_evolution(ekp, prior; path = joinpath(out, "parameters.png"))
            Test.@test isfile(params)
            Test.@test filesize(params) > 1_000

            err = PL.error_evolution(ekp; path = joinpath(out, "error.png"))
            Test.@test isfile(err)
            Test.@test filesize(err) > 1_000
        end

    end
end
