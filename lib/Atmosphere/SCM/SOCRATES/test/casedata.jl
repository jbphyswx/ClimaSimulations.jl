using SOCRATES: SOCRATES
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF
using Test: Test

# `cases/socrates.toml` is hand-edited, so what is worth checking is its schema against
# SSCF's own case list and the few values a typo would silently change. That a timescale
# is positive or a droplet number is within some invented window would pass for any
# plausible edit and catch none.
Test.@testset "case data" begin

    data = SOCRATES.CASE_DATA

    Test.@testset "covers exactly the valid cases" begin
        Test.@test Set(e["number"] for e in data["flight"]) == Set(SSCF.flight_numbers)
        Test.@test length(data["flight"]) == length(SSCF.flight_numbers)

        tabulated = Set((e["flight"], e["forcing"]) for e in data["case"])
        expected = Set(
            (fl, String(SSCF.symbol(ft))) for ft in SSCF.forcing_types,
            fl in SSCF.flight_numbers if SSCF.is_valid_flight_number(ft, fl)
        )
        Test.@test tabulated == expected
        Test.@test length(data["case"]) == length(expected)  # no duplicates
    end

    Test.@testset "the scored region fits inside each case's column" begin
        for ft in SSCF.forcing_types, fl in SSCF.flight_numbers
            SSCF.is_valid_flight_number(ft, fl) || continue
            top = SOCRATES.scored_z_top(fl, SSCF.symbol(ft))
            Test.@test top <= SOCRATES.z_max_default(SOCRATES.SocratesCase(fl, ft))
        end
    end

    Test.@testset "values a typo would change" begin
        Test.@test SOCRATES.droplet_number(9) == 1.9e8
        Test.@test SOCRATES.run_duration(:Obs) == 12 * 3600.0
        Test.@test SOCRATES.run_duration(:ERA5) == 14 * 3600.0
        # the scored window has to close within the run it scores
        Test.@test last(SOCRATES.obs_score_window()) <= SOCRATES.run_duration(:Obs)
    end

    Test.@testset "an absent case errors rather than returning a default" begin
        Test.@test_throws ErrorException SOCRATES.droplet_number(7)
        Test.@test_throws ErrorException SOCRATES.scored_z_top(11, :Obs)  # no Obs artifact
        Test.@test_throws ErrorException SOCRATES.scored_z_top(7, :ERA5)
    end

end
