using BOMEX: BOMEX
using ClimaAtmos: ClimaAtmos as CA
using Test: Test

# These assert the wiring between cases/bomex.toml and the generics ClimaAtmos dispatches on. A
# value read from the wrong key, or a generic left unimplemented so ClimaAtmos's `nothing`
# default applies, both surface here rather than as a silently unforced column.

Test.@testset "setup generics" begin
    FT = Float64
    c = BOMEX.case(FT)
    params = CA.ClimaAtmosParameters(FT)

    Test.@testset "subsidence is prescribed, not defaulted away" begin
        prof = CA.Setups.subsidence_forcing(c, FT)
        Test.@test !isnothing(prof)
        Test.@test isfinite(prof(0.0))
        # subsidence is downward through the trade inversion
        Test.@test prof(1000.0) < 0
        # ClimaAtmos wraps the profile; `nothing` here would mean an unforced column
        Test.@test CA.get_subsidence_model(FT; setup_type = c) isa CA.LargeScaleSubsidence
    end

    Test.@testset "large-scale advection supplies both tendencies" begin
        data = CA.Setups.large_scale_advection_forcing(c, FT)
        Test.@test !isnothing(data)
        Test.@test isfinite(data.prof_dTdt(1.0, 500.0))
        Test.@test isfinite(data.prof_dqtdt(500.0))
        Test.@test CA.get_large_scale_advection_model(FT; setup_type = c) isa
                   CA.LargeScaleAdvection
    end

    Test.@testset "Coriolis carries the case's own parameter" begin
        f = CA.Setups.coriolis_forcing(c, FT)
        Test.@test f.coriolis_param == FT(BOMEX.coriolis_parameter())
        Test.@test all(f.prof_vg(z) == 0 for z in (0.0, 500.0, 3000.0))
        Test.@test isfinite(f.prof_ug(500.0))
        Test.@test !isnothing(CA.get_scm_coriolis(FT; setup_type = c))
    end

    Test.@testset "the surface threads the case's prescribed fluxes through" begin
        s = CA.Setups.surface_condition(c, params)
        # `MoninObukhov` takes `z0` but stores a momentum and a scalar roughness, and holds the
        # two prescribed fluxes in a `θAndQFluxes` rather than as its own fields
        Test.@test s.flux_scheme.z0m == FT(BOMEX.surface_roughness())
        Test.@test s.flux_scheme.z0b == FT(BOMEX.surface_roughness())
        Test.@test s.flux_scheme.fluxes.θ_flux == FT(BOMEX.surface_theta_flux())
        Test.@test s.flux_scheme.fluxes.q_flux == FT(BOMEX.surface_q_flux())
        Test.@test s.flux_scheme.ustar == FT(BOMEX.surface_ustar())
        Test.@test s.overrides.p == FT(BOMEX.surface_override_pressure())
        Test.@test s.overrides.q_vap == FT(BOMEX.surface_override_q_vap())
    end

    Test.@testset "diagnostics are checked before a model is built" begin
        Test.@test_throws ErrorException BOMEX.bomex_diagnostics(
            ("not_a_diagnostic",); n_levels = 10,
        )
        Test.@test_throws ErrorException BOMEX.bomex_diagnostics((); n_levels = 10)
        for name in (BOMEX.default_diagnostic_vars..., BOMEX.STATE_VARS...)
            Test.@test haskey(CA.Diagnostics.ALL_DIAGNOSTICS, name)
        end
        # the case carries no radiation, so a radiative flux must not be among its defaults
        Test.@test !any(
            startswith(n, "rl") || startswith(n, "rs")
            for n in BOMEX.default_diagnostic_vars
        )
    end
end
