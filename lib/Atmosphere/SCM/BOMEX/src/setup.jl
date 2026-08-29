"""
    ClimaAtmos.Setups.center_initial_condition(c::BOMEXCase, local_geometry, params)

The initial centre-level state: temperature from the case's hydrostatic pressure against its
`θ_li` and `q_tot` profiles, the prescribed wind, and TKE either carried by the model from zero
or read from the prescribed profile.
"""
function CA.Setups.center_initial_condition(c::BOMEXCase, local_geometry, params)
    FT = eltype(params)
    thermo_params = CA.Parameters.thermodynamics_params(params)
    (; z) = local_geometry.coordinates
    (; prognostic_tke, profiles) = c

    q_tot = profiles.q_tot(z)
    p = profiles.p(z)
    T = CA.TD.air_temperature(thermo_params, CA.TD.pθ_li(), p, profiles.θ(z), q_tot)

    return CA.Setups.physical_state(;
        T,
        p,
        q_tot,
        u = profiles.u(z),
        tke = prognostic_tke ? FT(0) : profiles.tke(z),
        q_gas_A = FT(1),
    )
end

"""
    ClimaAtmos.Setups.surface_condition(c::BOMEXCase, params)

The surface: prescribed `θ` and `q` fluxes and a prescribed friction velocity, so the fluxes
are the case's rather than a bulk formula's, with the boundary state pinned to the case's
surface pressure and saturation humidity.
"""
function CA.Setups.surface_condition(::BOMEXCase, params)
    FT = eltype(params)
    flux_scheme = CA.SurfaceConditions.MoninObukhov(;
        z0 = FT(surface_roughness()),
        θ_flux = FT(surface_theta_flux()),
        q_flux = FT(surface_q_flux()),
        ustar = FT(surface_ustar()),
    )
    return (;
        flux_scheme,
        temperature = CA.SurfaceConditions.AnalyticTemperature(
            Returns(FT(surface_temperature())),
        ),
        overrides = CA.SurfaceConditions.SurfaceBoundaryOverrides(;
            p = FT(surface_override_pressure()),
            q_vap = FT(surface_override_q_vap()),
        ),
    )
end

"""Prescribed large-scale subsidence [m s^-1] against height."""
CA.Setups.subsidence_forcing(::BOMEXCase, ::Type{FT}) where {FT} =
    CA.APL.Bomex_subsidence(FT)

"""Prescribed large-scale drying and radiative cooling against height."""
CA.Setups.large_scale_advection_forcing(::BOMEXCase, ::Type{FT}) where {FT} =
    (; prof_dTdt = CA.APL.Bomex_dTdt(FT), prof_dqtdt = CA.APL.Bomex_dqtdt(FT))

"""Geostrophic wind and the case's Coriolis parameter; geostrophic `v` is identically zero."""
CA.Setups.coriolis_forcing(::BOMEXCase, ::Type{FT}) where {FT} = (;
    prof_ug = CA.APL.Bomex_geostrophic_u(FT),
    prof_vg = Returns(FT(0)),
    coriolis_param = FT(coriolis_parameter()),
)
