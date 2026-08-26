using EnsembleKalmanProcesses: EnsembleKalmanProcesses as EKP
using LinearAlgebra: LinearAlgebra
using SOCRATES: SOCRATES
using SOCRATESCalibration: SOCRATESCalibration as SC
using Test: Test

# What is checkable without running a forward model: the prior, the interface's validation, the
# misfit metric, and the T_stops ratchet.

Test.@testset "prior" begin
    prior = SC.default_prior()
    names = SC.prior_names()
    Test.@test length(names) == length(SC.DEFAULT_PRIOR_SPEC)
    Test.@test EKP.get_name(prior) == names
    # the four terminal-velocity factors the model reads must be among them, spelled identically
    for param in keys(SOCRATES.TERMINAL_VELOCITY_SCALING_PARAMS)
        Test.@test String(param) in names
    end

    # a mean outside its own bounds is a silent disaster in the unconstrained space, so it errors
    Test.@test_throws ErrorException SC.default_prior((; bad = (10.0, 0.0, 1.0, 1.0)))
    Test.@test_throws ErrorException SC.default_prior((; bad = (0.5, 0.0, 1.0, 0.0)))
    Test.@test_throws ErrorException SC.default_prior(NamedTuple())

    # every default mean lies inside its bounds
    for (name, (mean, lower, upper, σ)) in pairs(SC.DEFAULT_PRIOR_SPEC)
        Test.@test lower <= mean <= upper
        Test.@test σ > 0
    end
end

Test.@testset "interface validation" begin
    mktempdir() do out
        Test.@test_throws ErrorException SC.SocratesInterface(;
            output_dir = out, cases = SOCRATES.SocratesCase[],
        )
        Test.@test_throws ErrorException SC.SocratesInterface(;
            output_dir = out, cases = [SOCRATES.case("RF09_Obs")], vars = String[],
        )
        # the grid is the interface's business, so passing one through `run_kwargs` is refused rather
        # than silently overriding the observation's grid
        for key in (:grid, :dz_min, :faces)
            Test.@test_throws ErrorException SC.SocratesInterface(;
                output_dir = out,
                cases = [SOCRATES.case("RF09_Obs")],
                run_kwargs = (; key => 1),
            )
        end
    end
end

struct StubEKP
    y::Vector{Float64}
    Γ::Matrix{Float64}
end
EKP.get_obs(e::StubEKP) = e.y
EKP.get_obs_noise_cov(e::StubEKP) = e.Γ

Test.@testset "member misfit" begin
    y = [1.0, 2.0, 3.0]
    ekp = StubEKP(y, Matrix{Float64}(LinearAlgebra.I, 3, 3))
    g = [1.0 2.0 NaN; 2.0 2.0 3.0; 3.0 2.0 3.0]
    φ = SC.member_misfits(ekp, g)
    Test.@test φ[1] == 0                      # exact match
    Test.@test φ[2] ≈ 1.0                     # d = [1, 0, -1]
    Test.@test isnan(φ[3])                    # a failed member never scores best

    # Γ is inverted, not ignored: with Γ₁₁ = 4 a deviation of 2 costs ½·4/4
    ekp2 = StubEKP(y, [4.0 0 0; 0 1.0 0; 0 0 1.0])
    Test.@test only(SC.member_misfits(ekp2, reshape([3.0, 2.0, 3.0], 3, 1))) ≈ 0.5
end

Test.@testset "set_field" begin
    x = (a = 1, b = 2)
    # a NamedTuple is reconstructed positionally, so the wrong name would silently move a value
    Test.@test_throws ErrorException SC.set_field(x, :nope, 3)
end

struct StubScheduler
    terminate_at::Float64
end
struct StubTStopEKP
    scheduler::StubScheduler
    Δt::Vector{Float64}
end
EKP.get_Δt(e::StubTStopEKP) = e.Δt

Test.@testset "T_stops" begin
    stops = [1.0, 10.0, 100.0]
    at(Δt, terminate_at) = StubTStopEKP(StubScheduler(terminate_at), Δt)

    # accumulated T = 0.5 is below the first stop, so the bound stays there
    Test.@test SC.apply_terminate_at(at([0.5], 1.0), stops).scheduler.terminate_at == 1.0
    # past the first stop, the bound advances to the next
    Test.@test SC.apply_terminate_at(at([2.0], 1.0), stops).scheduler.terminate_at == 10.0
    Test.@test SC.apply_terminate_at(at([20.0], 10.0), stops).scheduler.terminate_at == 100.0
    # past the last stop there is nothing to advance to
    Test.@test SC.apply_terminate_at(at([200.0], 100.0), stops).scheduler.terminate_at == 100.0
    # and the bound never moves backwards
    Test.@test SC.apply_terminate_at(at([2.0], 50.0), stops).scheduler.terminate_at == 50.0

    # no schedule leaves the object alone
    unchanged = at([2.0], 1.0)
    Test.@test SC.apply_terminate_at(unchanged, nothing) === unchanged
    Test.@test SC.apply_terminate_at(unchanged, Float64[]) === unchanged

    Test.@test SC.accumulated_T(at([1.0, 2.0, 3.0], 1.0)) == 6.0
end

Test.@testset "G_ensemble location" begin
    Test.@test SC.g_ensemble_path("/out", 7) ==
               joinpath("/out", "iteration_007", "G_ensemble.jld2")
end
