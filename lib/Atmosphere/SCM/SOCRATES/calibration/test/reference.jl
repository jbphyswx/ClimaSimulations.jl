using ClimaAnalysis: ClimaAnalysis
using ClimaAtmos: ClimaAtmos as CA
using SOCRATES: SOCRATES
using SOCRATESCalibration: SOCRATESCalibration as SC
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF
using Test: Test

Test.@testset "scoring" begin

    Test.@testset "Atlas LES reference loads for every case" begin
        for c in SOCRATES.all_cases()
            ref = SC.les_outputvars(c)
            for name in SOCRATES.SCORED_VARS
                Test.@test haskey(ref, name)
            end
            clw = ref["clw"]
            Test.@test ClimaAnalysis.has_altitude(clw)
            Test.@test ClimaAnalysis.has_time(clw)
            # profiles are (z, time); water paths are (time,)
            Test.@test ndims(clw.data) == 2
            Test.@test ndims(ref["lwp"].data) == 1
            # the reference time axis starts at zero after `_elapsed`
            Test.@test first(clw.dims["time"]) == 0.0
        end
    end

    Test.@testset "the reference is labelled in the units the model writes" begin
        # the conversion itself is the registry's and is checked in `atlas_registry.jl`;
        # what this file owns is that the `OutputVar` is labelled with ClimaAtmos's
        # spelling, since the comparison pairs them by name
        c = SOCRATES.case("RF09_Obs")
        ref = SC.les_outputvars(c)
        for (name, pair) in SC.MODEL_TO_ATLAS
            haskey(ref, name) || continue
            Test.@test ref[name].attributes["units"] == pair.model_units
        end
        # a short name with no Atlas counterpart is refused rather than skipped
        Test.@test_throws ErrorException SC.les_outputvars(c; vars = ("ta_not_real",))
    end

    Test.@testset "scored region" begin
        for c in SOCRATES.all_cases()
            lo, hi = SC.default_z_bounds(c)
            Test.@test hi <= SOCRATES.z_max(c)
            # `scored_levels` selects, so what matters is that it selects exactly the
            # levels inside the bounds and nothing else
            levels = SC.scored_levels(SOCRATES.native_z(c), (lo, hi))
            Test.@test levels == filter(l -> lo <= l <= hi, SOCRATES.native_z(c))
        end
        # a region with no levels in it is an error, not an empty comparison
        Test.@test_throws ErrorException SC.scored_levels([5000.0], (0.0, 100.0))
    end

    Test.@testset "reference resampling onto model levels" begin
        c = SOCRATES.case("RF09_Obs")
        ref = SC.les_outputvars(c; vars = ("clw",))["clw"]
        levels = SC.scored_levels(SOCRATES.native_z(c), SC.default_z_bounds(c))
        on_levels = SC.reference_on_levels(ref, levels)
        Test.@test on_levels.dims["z"] ≈ levels
        Test.@test size(on_levels.data, 1) == length(levels)
        # asking for a level above the reference column is an error, not a clamp
        Test.@test_throws ErrorException SC.reference_on_levels(ref, [1.0e6])
    end

    Test.@testset "ScoreTransform" begin
        t = SC.ScoreTransform()
        # a nonzero reference normalizes by its own magnitude
        pv = SC.pool_var(t, "clw", [1.0e-4, 2.0e-4, 1.5e-4])
        Test.@test pv > 0
        Test.@test isfinite(pv)
        # an all-zero reference falls back to the characteristic magnitude, so the
        # normalization can never divide by zero
        pv0 = SC.pool_var(t, "clw", zeros(4))
        Test.@test pv0 ≈ SC.DEFAULT_OBS_VAR_SCALING["clw"] *
                         SC.DEFAULT_CHARACTERISTIC["clw"]^2
        Test.@test pv0 > 0

        diag = SC.uncertainty_diagonal(t, "clw", fill(1.0e-4, 3, 5))
        Test.@test length(diag) == 3
        Test.@test all(>(0), diag)

        Test.@test_throws ErrorException SC.pool_var(t, "not_a_variable", [1.0])
    end

end

const PARAMS = SOCRATES.socrates_params(Float64, SOCRATES.case("RF09_Obs"))

Test.@testset "the two vocabularies cannot drift apart" begin
    diags = CA.Diagnostics.ALL_DIAGNOSTICS
    for (name, pair) in SC.MODEL_TO_ATLAS
        Test.@test haskey(SOCRATES.atlas_specs(PARAMS), pair.atlas)
        Test.@test SOCRATES.atlas_specs(PARAMS)[pair.atlas].units == pair.atlas_units
        Test.@test haskey(diags, name)
        Test.@test diags[name].units == pair.model_units
    end
    for name in SOCRATES.SCORED_VARS
        Test.@test haskey(SC.MODEL_TO_ATLAS, name)
    end
    # the one pairing whose two vocabularies differ numerically rather than in spelling
    Test.@test SC.MODEL_TO_ATLAS["cl"].to_model(0.25) == 25
end

