using ClimaAtmos: ClimaAtmos as CA
using SwirlLMCloudBenchSim: SwirlLMCloudBenchSim as CB
using SwirlLMCloudBenchCalibration: SwirlLMCloudBenchCalibration as CBC
using Test: Test

# Reading a case's published output needs the network; these assertions do not.

Test.@testset "SwirlLMCloudBenchCalibration" begin

    Test.@testset "every pairing names real diagnostics with the units it declares" begin
        Test.@test CBC.assert_pairings_are_diagnostics() === nothing
        # a pairing that names something ClimaAtmos does not register is refused
        Test.@test_throws ErrorException CBC.assert_pairings_are_diagnostics(;
            pairings = Dict(
                "olr" => (
                    model = ("not_a_diagnostic",), combine = only,
                    to_model = identity, units = "W m^-2",
                ),
            ),
        )
        # and so is one whose declared units disagree with ClimaAtmos's
        Test.@test_throws ErrorException CBC.assert_pairings_are_diagnostics(;
            pairings = Dict(
                "olr" => (
                    model = ("rlut",), combine = only, to_model = identity,
                    units = "kg m^-2",
                ),
            ),
        )
    end

    Test.@testset "the model side is the arithmetic the pairing states" begin
        # `cl` is a (z, time) profile; the rest are one value per time
        cl = [1.0 2.0; 8.5 3.0; 0.0 40.0]
        run_vars = Dict(
            "rsdt" => [400.0], "rsut" => [130.0], "rlut" => [255.0],
            "rlutcs" => [265.0], "rsutcs" => [70.0], "lwp" => [0.08],
            "hfls" => [90.0], "hfss" => [12.0], "cl" => cl,
        )
        Test.@test CBC.model_comparable(run_vars, "lwp") == [0.08]
        Test.@test CBC.model_comparable(run_vars, "olr") == [255.0]
        # absorbed shortwave is incident less outgoing
        Test.@test CBC.model_comparable(run_vars, "asr") == [270.0]
        # longwave cloud effect warms, shortwave cools
        Test.@test CBC.model_comparable(run_vars, "cre_lw") == [10.0]
        Test.@test CBC.model_comparable(run_vars, "cre_sw") == [-60.0]
        # maximum overlap: the column is as cloudy as its cloudiest layer, per time
        Test.@test CBC.model_comparable(run_vars, "cloud_cover") == [8.5, 40.0]
        Test.@test CBC.maximum_overlap_cover(cl) == [8.5, 40.0]
        # the profile is reduced over height, not over time
        Test.@test length(CBC.model_comparable(run_vars, "cloud_cover")) == size(cl, 2)
        # a run missing a diagnostic the pairing needs fails by name
        Test.@test_throws ErrorException CBC.model_comparable(
            Dict("rsdt" => [400.0]), "asr",
        )
        # a store variable nobody paired is refused rather than guessed at
        Test.@test_throws ErrorException CBC.model_comparable(
            run_vars, "sfc_flux_rad_lw",
        )
    end

    Test.@testset "the scored diagnostics are what a run has to write" begin
        vars = CBC.scored_model_vars()
        Test.@test issorted(vars)
        Test.@test allunique(vars)
        for name in vars
            Test.@test haskey(CA.Diagnostics.ALL_DIAGNOSTICS, name)
        end
        # the clear-sky fluxes the cloud radiative effects need are in the set
        Test.@test "rlutcs" in vars && "rsutcs" in vars
        # a run configured from this package writes every one of them
        Test.@test CB.cloudbench_diagnostics(vars; n_levels = 10) !== nothing
    end


    include("postprocess.jl")

end
