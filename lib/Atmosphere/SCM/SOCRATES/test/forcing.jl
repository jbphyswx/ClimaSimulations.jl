using ClimaAtmos: ClimaAtmos as CA
using SOCRATES: SOCRATES
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF
using Test: Test

const FT = Float64

# A `Y` on the case's own grid is all `external_forcing_cache` needs: it reads
# `axes(Y.c)`, the centre heights, and allocates with `similar(Y.c, FT)`.
function case_state(case)
    grid = SOCRATES.socrates_grid(FT, case)
    params = CA.ClimaAtmosParameters(FT)
    model = CA.AtmosModel()
    spaces = CA.get_spaces(grid)
    setup = CA.Setups.DecayingProfile(; perturb = false, params)
    Y = CA.Setups.initial_state(
        setup, params, model, spaces.center_space, spaces.face_space,
    )
    return (; grid, params, Y, z = SOCRATES.socrates_z(grid))
end

forcing_for(source) = SOCRATES.SOCRATESForcing(
    FT, source; scalar_nudge_timescale = 1200.0, wind_nudge_timescale = 1200.0,
)

Test.@testset "SOCRATESForcing" begin

    Test.@testset "construction validates" begin
        case = SOCRATES.case("RF09_Obs")
        f = forcing_for(case)
        Test.@test f.source == case
        Test.@test f.source.forcing_type isa SSCF.ObsForcing
        Test.@test f.inv_τ_scalar ≈ 1 / 1200
        Test.@test f.inv_τ_wind ≈ 1 / 1200
        Test.@test isbits(f)

        # flight 11 has no Obs artifact; built directly, since `case` rejects it first
        Test.@test_throws ErrorException forcing_for(
            SOCRATES.SocratesCase(11, SSCF.ObsForcing()),
        )
        # a path that does not exist is refused rather than deferred to the cache
        Test.@test_throws ErrorException forcing_for(joinpath(mktempdir(), "absent.nc"))
        # a non-positive timescale is a division by zero waiting to happen
        Test.@test_throws ErrorException SOCRATES.SOCRATESForcing(
            FT, case;
            scalar_nudge_timescale = 0.0, wind_nudge_timescale = 1200.0,
        )
        Test.@test_throws ErrorException SOCRATES.SOCRATESForcing(
            FT, case;
            scalar_nudge_timescale = 1200.0, wind_nudge_timescale = -60.0,
        )
    end

    Test.@testset "cache is sampled on the model levels" begin
        case = SOCRATES.case("RF09_Obs")
        (; params, Y, z) = case_state(case)
        cache = CA.external_forcing_cache(
            Y, forcing_for(case), params, SOCRATES.les_start_datetime(case),
        )
        for name in SOCRATES.SSCF_FORCING_VARS
            Test.@test haskey(cache.interpolants, name)
            Test.@test length(cache.interpolants[name]) == length(z)
        end
        for name in
            (:ᶜdTdt_hadv, :ᶜdqtdt_hadv, :ᶜT_nudge, :ᶜqt_nudge, :ᶜu_nudge, :ᶜv_nudge, :ᶜls_subsidence)
            Test.@test CA.CC.Spaces.nlevels(axes(getproperty(cache, name))) == length(z)
        end
    end

    Test.@testset "fill reproduces the interpolants exactly" begin
        case = SOCRATES.case("RF09_Obs")
        (; params, Y, z) = case_state(case)
        cache = CA.external_forcing_cache(
            Y, forcing_for(case), params, SOCRATES.les_start_datetime(case),
        )
        field = cache.ᶜT_nudge
        itps = cache.interpolants.T_nudge
        # on a node, between nodes, and at the end of the run
        for t in (0.0, 1800.0, 3600.0, 21_600.0, 43_200.0)
            SOCRATES._fill_column!(field, itps, t)
            got = vec(Array(parent(field)))
            want = [itps[k](t) for k in eachindex(itps)]
            Test.@test got == want
            Test.@test all(isfinite, got)
        end
        # the fill must not allocate: it runs every stage of every step
        SOCRATES._fill_column!(field, itps, 0.0)
        Test.@test (@allocated SOCRATES._fill_column!(field, itps, 600.0)) == 0
    end

    Test.@testset "fill rejects a geometry it cannot fill" begin
        case = SOCRATES.case("RF09_Obs")
        (; params, Y) = case_state(case)
        cache = CA.external_forcing_cache(
            Y, forcing_for(case), params, SOCRATES.les_start_datetime(case),
        )
        # wrong number of interpolants for the field
        Test.@test_throws ErrorException SOCRATES._fill_column!(
            cache.ᶜT_nudge, cache.interpolants.T_nudge[1:3], 0.0,
        )
    end

    Test.@testset "every case samples cleanly on its own grid" begin
        for c in SOCRATES.all_cases()
            (; params, Y, z) = case_state(c)
            cache = CA.external_forcing_cache(
                Y, forcing_for(c), params, SOCRATES.les_start_datetime(c),
            )
            t_end = SOCRATES.t_end(c)
            for name in SOCRATES.SSCF_FORCING_VARS
                itps = cache.interpolants[name]
                Test.@test length(itps) == length(z)
                # the case must be evaluable across its whole run: SSCF
                # interpolants error outside their node span rather than
                # extrapolating, so this is a real check, not a formality
                for t in (0.0, t_end / 2, t_end)
                    vals = [itps[k](t) for k in eachindex(itps)]
                    Test.@test all(isfinite, vals)
                end
            end
        end
    end

    Test.@testset "a written file reproduces the in-memory forcing" begin
        case = SOCRATES.case("RF09_Obs")
        (; params, Y, z) = case_state(case)
        thermo = CA.Parameters.thermodynamics_params(params)
        start = SOCRATES.les_start_datetime(case)

        path = joinpath(mktempdir(), "RF09_Obs_forcing.nc")
        SOCRATES.write_forcing_file(path, case, z; thermo_params = thermo)

        from_case = CA.external_forcing_cache(Y, forcing_for(case), params, start)
        from_file = CA.external_forcing_cache(Y, forcing_for(path), params, start)

        # `external_forcing_tendency!` reads the state only through these fields, so
        # equality here is equality of the tendency the two sources produce
        for name in SOCRATES.SSCF_FORCING_VARS, t in (0.0, 1234.5, 21_600.0, 43_200.0)
            a, b = similar(from_case.ᶜT_nudge), similar(from_file.ᶜT_nudge)
            SOCRATES._fill_column!(a, from_case.interpolants[name], t)
            SOCRATES._fill_column!(b, from_file.interpolants[name], t)
            Test.@test vec(Array(parent(a))) == vec(Array(parent(b)))
        end
        for name in SOCRATES.SSCF_SURFACE_VARS, t in (0.0, 1234.5, 43_200.0)
            Test.@test from_case.surface[name](t) == from_file.surface[name](t)
        end

        # the file belongs to its grid and to the thermodynamics it was derived with
        Test.@test_throws ErrorException SOCRATES.read_forcing_file(
            path, z[1:(end - 1)], thermo,
        )
        other = CA.Parameters.thermodynamics_params(
            SOCRATES.socrates_params(SOCRATES.socrates_toml_dict(FT, case)),
        )
        Test.@test_throws ErrorException SOCRATES.read_forcing_file(path, z, other)
    end

end
