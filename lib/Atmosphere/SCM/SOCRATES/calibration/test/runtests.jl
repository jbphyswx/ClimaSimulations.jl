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
        Test.@test_throws ErrorException SC.SOCRATESInterface(;
            output_dir = out, cases = SOCRATES.SOCRATESCase[],
        )
        Test.@test_throws ErrorException SC.SOCRATESInterface(;
            output_dir = out, cases = [SOCRATES.case("RF09_Obs")], vars = String[],
        )
        # the grid is the interface's business, so passing one through `run_kwargs` is refused rather
        # than silently overriding the observation's grid
        for key in (:grid, :dz_min, :faces)
            Test.@test_throws ErrorException SC.SOCRATESInterface(;
                output_dir = out,
                cases = [SOCRATES.case("RF09_Obs")],
                run_kwargs = (; key => 1),
            )
        end
    end
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

include("postprocess.jl")
include("reference.jl")
