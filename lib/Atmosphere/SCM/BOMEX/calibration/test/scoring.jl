using BOMEX: BOMEX
using BOMEXCalibration: BOMEXCalibration as BC
using Test: Test

Test.@testset "scoring" begin
    FT = Float64
    c = BOMEX.case(FT)

    Test.@testset "the window is the last hour of the run" begin
        w = BC.score_window(c)
        Test.@test last(w) == BOMEX.t_end(c)
        Test.@test last(w) - first(w) == 3600.0
    end

    Test.@testset "levels are the model's own, inside the compared region" begin
        z = collect(0.0:100.0:3000.0)
        levels = BC.scored_levels(z, (500.0, 1500.0))
        Test.@test first(levels) >= 500.0
        Test.@test last(levels) <= 1500.0
        Test.@test issubset(levels, z)
        # a region containing no model level is an error, not an empty comparison
        Test.@test_throws ErrorException BC.scored_levels(z, (10.0, 20.0))
    end

    Test.@testset "a reference must say what it is" begin
        # no heights
        Test.@test_throws ErrorException BC.read_reference(
            Dict("clw" => [1.0, 2.0]),
        )
        # a variable that does not match its own height axis
        Test.@test_throws ErrorException BC.read_reference(
            Dict("z" => [0.0, 100.0, 200.0], "clw" => [1.0, 2.0]),
        )
        Test.@test_throws ErrorException BC.read_reference("not_a_file.nc")

        ref = BC.read_reference(
            Dict("z" => [0.0, 100.0, 200.0], "clw" => [1.0, 2.0, 3.0]),
        )
        Test.@test ref.z == [0.0, 100.0, 200.0]
        Test.@test ref.data["clw"] == [1.0, 2.0, 3.0]
    end

    Test.@testset "the reference is resampled, never extrapolated" begin
        z = [0.0, 100.0, 200.0]
        values = [0.0, 10.0, 20.0]
        # linear in height, so a midpoint is the midpoint
        Test.@test BC.reference_on_levels(z, values, [50.0]) ≈ [5.0]
        Test.@test BC.reference_on_levels(z, values, [0.0, 200.0]) ≈ [0.0, 20.0]
        # reaching past the reference column would invent values
        Test.@test_throws ErrorException BC.reference_on_levels(z, values, [500.0])
    end

    Test.@testset "a variable the reference lacks is named, not skipped" begin
        source = Dict("z" => [0.0, 100.0], "clw" => [1.0, 2.0])
        Test.@test_throws ErrorException BC.reference_values(
            source, ("cl",), [0.0, 100.0],
        )
        got = BC.reference_values(source, ("clw",), [0.0, 100.0])
        Test.@test got["clw"] ≈ [1.0, 2.0]
    end

    Test.@testset "the prior is over parameters configs/BOMEX.toml actually sets" begin
        # calibrating a name the case never configures would sample a value nothing reads
        configured = keys(BOMEX._override_dict(BOMEX.reference_parameter_file()))
        for name in BC.prior_names()
            Test.@test name in configured
        end
        prior = BC.default_prior()
        Test.@test BC.EKP.ParameterDistributions.get_name(prior) == BC.prior_names()
        Test.@test BC.EKP.ParameterDistributions.ndims(prior) == length(BC.prior_names())
        # bounds are checked, so a mean outside them is refused rather than silently clipped
        Test.@test_throws ErrorException BC.default_prior((; bad = (5.0, 1.0, 2.0, 1.0)))
        Test.@test_throws ErrorException BC.default_prior((; bad = (1.5, 1.0, 2.0, -1.0)))
        Test.@test_throws ErrorException BC.default_prior((;))
    end
end
