using Test: Test
using MOSAiC_AYiL: MOSAiC_AYiL as MA
using ClimaAtmos: ClimaAtmos as CA

# 20200503 carries cloud liquid and no ice; 20200219 carries ice and, being far
# colder, no liquid. Both branches of the condensate split are therefore
# exercised, which one day alone cannot do.

function initial_column(date; dz_min = 50)
    FT = Float64
    c = MA.case(date)
    fd = MA.read_scm_in(c.date)
    grid = MA.mosaic_grid(
        FT, c;
        faces = MA.coarsen_faces_to_dz_min(MA.native_faces(c), dz_min),
    )
    z = MA.mosaic_z(grid)
    params = MA.mosaic_params(FT, c)
    setup = MA.MOSAiCSetup(FT, c; forcing = fd)
    lg = CA.CC.Fields.local_geometry_field(CA.get_spaces(grid).center_space)
    states = [CA.Setups.center_initial_condition(setup, lg[k], params) for k in 1:length(z)]
    return (; z, fd, params, states)
end

Test.@testset "condensate split" for date in ("20200503", "20200219")
    (; z, fd, params, states) = initial_column(date)
    thermo = CA.Parameters.thermodynamics_params(params)
    q_tot = [s.q_tot for s in states]
    q_liq = [s.q_liq for s in states]
    q_ice = [s.q_ice for s in states]
    T = [s.T for s in states]
    ρ = [s.ρ for s in states]

    Test.@test all(isfinite, q_tot) && all(isfinite, T) && all(isfinite, ρ)
    Test.@test all(>(0), ρ)
    Test.@test all(>=(0), q_liq)
    Test.@test all(>=(0), q_ice)
    Test.@test all(q_liq .+ q_ice .<= q_tot .+ 1.0e-14)

    # the liquid is the saturation excess of the vapour-plus-liquid pool alone —
    # ice excluded, or it would be counted twice
    for k in eachindex(z)
        q_sat_liq = CA.TD.q_vap_saturation(thermo, T[k], ρ[k], CA.TD.Liquid())
        Test.@test q_liq[k] ≈ max(0.0, (q_tot[k] - q_ice[k]) - q_sat_liq)
    end

    # the ice is the file's own profile, carried over rather than re-derived
    Test.@test maximum(q_ice) ≈ maximum(fd.qi) rtol = 0.05
end

Test.@testset "which branch each day exercises" begin
    liquid_day = initial_column("20200503")
    ice_day = initial_column("20200219")
    Test.@test count(>(1.0e-9), [s.q_liq for s in liquid_day.states]) > 0
    Test.@test all(iszero, [s.q_ice for s in liquid_day.states])
    Test.@test count(>(1.0e-9), [s.q_ice for s in ice_day.states]) > 0
end

Test.@testset "density is the reference run's own" begin
    c = MA.case("20200503")
    z_les, ρ_les = MA.les_density(c.date)
    (; z, states) = initial_column("20200503")
    ρ = [s.ρ for s in states]
    # sampled from `rhof`, so it spans the same range
    Test.@test minimum(ρ) >= minimum(ρ_les) - 1.0e-6
    Test.@test maximum(ρ) <= maximum(ρ_les) + 1.0e-6
end

Test.@testset "TKE seed" begin
    # the DALES cold-start profile, decaying away from the surface
    Test.@test MA.dales_tke_seed(0.0) ≈ (MA.DALES_CONSTANTS.e12_min + 1)^2
    Test.@test MA.dales_tke_seed(50.0) < MA.dales_tke_seed(0.0)
    Test.@test MA.dales_tke_seed(5000.0) ≈ MA.DALES_CONSTANTS.e12_min^2 rtol = 1.0e-6
    Test.@test MA.dales_tke_seed(100.0; decay_length = 25) <
               MA.dales_tke_seed(100.0; decay_length = 100)
end

Test.@testset "insolation" begin
    FT = Float64
    # early May at 82 N is daylight
    may = MA.MOSAiCInsolation(FT, MA.case("20200503"))
    Test.@test may.cos_zenith > 0
    Test.@test may.toa_flux > 1000

    # midwinter at 86 N is polar night: no flux, and a positive zenith cosine so
    # RRTMGP still has something to divide by
    december = MA.MOSAiCInsolation(FT, MA.case("20191219"))
    Test.@test december.toa_flux == 0
    Test.@test december.cos_zenith == eps(FT)

    # the reference time is 11:00 UTC on the case date
    Test.@test MA.reference_datetime(MA.case("20200503")) ==
               MA.Dates.DateTime(2020, 5, 3, 11)
end
