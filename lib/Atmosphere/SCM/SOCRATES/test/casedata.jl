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

    end

    Test.@testset "values a typo would change" begin
        Test.@test SOCRATES.droplet_number(9) == 1.9e8
        Test.@test SOCRATES.run_duration(:Obs) == 12 * 3600.0
        Test.@test SOCRATES.run_duration(:ERA5) == 14 * 3600.0
    end

    Test.@testset "an absent case errors rather than returning a default" begin
        Test.@test_throws ErrorException SOCRATES.droplet_number(7)
    end

end
