using ClimaAtmos: ClimaAtmos as CA
using SOCRATES: SOCRATES
using Test: Test

# The only thing this mode computes is the conversion of a prescribed heating rate into an
# energy tendency. Sampling is SSCF's, and dispatch is a property of the source file, so
# neither is worth asserting; the conversion is where a wrong heat capacity, a dropped
# density or a sign error would live.

Test.@testset "prescribed radiation applies rho*cv_m*dTdt" begin
    case = SOCRATES.case("RF09_Obs")
    params = SOCRATES.socrates_params(Float64, case)
    thermo = CA.Parameters.thermodynamics_params(params)
    grid = SOCRATES.socrates_grid(Float64, case)
    z = SOCRATES.socrates_z(grid)
    spaces = CA.get_spaces(grid)
    Y = CA.Setups.initial_state(
        CA.Setups.DecayingProfile(; perturb = false, params),
        params,
        CA.AtmosModel(),
        spaces.center_space,
        spaces.face_space,
    )

    rad = SOCRATES.SOCRATESPrescribedRadiation(Float64, case, z; thermo_params = thermo)
    cache = CA.radiation_model_cache(Y, rad)

    # a state with condensate, so `cv_m` is not just the dry value
    q_tot = similar(Y.c, Float64)
    q_liq = similar(Y.c, Float64)
    q_ice = similar(Y.c, Float64)
    parent(q_tot) .= 5.0e-3
    parent(q_liq) .= 3.0e-4
    parent(q_ice) .= 1.0e-5
    p = (;
        params,
        radiation = cache,
        precomputed = (;
            ᶜq_tot_nonneg = q_tot, ᶜq_liq = q_liq, ᶜq_ice = q_ice,
        ),
    )

    t = 21_600.0
    Yₜ = similar(Y)
    parent(Yₜ.c.ρe_tot) .= 0
    CA.radiation_tendency!(Yₜ, Y, p, t, rad)

    # the heating that should have been applied, level by level, from the archive itself
    dTdt = [rad.dTdt_rad[k](t) for k in eachindex(rad.dTdt_rad)]
    cv_m = CA.TD.cv_m(thermo, 5.0e-3, 3.0e-4, 1.0e-5)
    want = vec(Array(parent(Y.c.ρ))) .* cv_m .* dTdt
    got = vec(Array(parent(Yₜ.c.ρe_tot)))
    Test.@test got ≈ want

    # the density and the heat capacity both matter: dropping either changes every level
    Test.@test !(got ≈ cv_m .* dTdt)
    Test.@test !(got ≈ vec(Array(parent(Y.c.ρ))) .* dTdt)

    # a tendency is accumulated onto whatever is already there, not assigned over it
    CA.radiation_tendency!(Yₜ, Y, p, t, rad)
    Test.@test vec(Array(parent(Yₜ.c.ρe_tot))) ≈ 2 .* want
end
