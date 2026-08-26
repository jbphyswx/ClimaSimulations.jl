"""
    radiation.jl

Radiative heating read from the forcing rather than computed, as an alternative to the
interactive RRTMGP `radiation_mode` [`socrates_model`](@ref) defaults to.

SSCF carries `:dTdt_rad` from the Atlas LES itself, so a column can be driven with the
heating the reference produced instead of recomputing it. Which of the two to use is the
caller's choice: pass `radiation_mode = SOCRATESPrescribedRadiation(...)` to
[`socrates_model`](@ref).
"""

"""
    SOCRATESPrescribedRadiation(FT, case, z; thermo_params, forcing_type)

Radiative heating prescribed from the Atlas LES, as a per-level vector of time
interpolants on levels `z`.

Grid-dependent, because the sampling happens here: `ClimaAtmos.radiation_model_cache`
receives only `Y`, so it cannot reach the thermodynamic parameters SSCF derives the profile
with. Build it on the same levels the simulation runs on.

With this mode ClimaAtmos installs no RRTMGP callback, so the `insolation` model is never
consulted and the solar geometry plays no part in the heating.

Unlike the nudging targets, which error outside their node span, `:dTdt_rad` is built on
the LES output axis with an extrapolating boundary condition. An Obs case runs 300 s past
that axis's last node at 42900 s, so its final step extrapolates: measured at most
1.02 K/day away from the last node, 0.0097 K/day at the median.
"""
struct SOCRATESPrescribedRadiation{V <: AbstractVector}
    dTdt_rad::V
end

function SOCRATESPrescribedRadiation(
    ::Type{FT},
    flight_number::Integer,
    forcing_type::SSCF.AbstractForcingType,
    z::AbstractVector;
    thermo_params,
) where {FT <: AbstractFloat}
    sampled = SSCF.get_column_forcing(
        Int(flight_number),
        forcing_type,
        (:dTdt_rad,);
        new_z = collect(Float64, z),
        thermodynamics_backend = thermo_params,
    )
    interpolants = sampled.dTdt_rad
    length(interpolants) == length(z) || error(
        "SSCF returned $(length(interpolants)) level interpolants for `dTdt_rad` on a \
         $(length(z))-level grid.",
    )
    return SOCRATESPrescribedRadiation(interpolants)
end

SOCRATESPrescribedRadiation(
    ::Type{FT},
    c::SOCRATESCase,
    z::AbstractVector;
    kwargs...,
) where {FT <: AbstractFloat} =
    SOCRATESPrescribedRadiation(FT, c.flight_number, c.forcing_type, z; kwargs...)

"""
    ClimaAtmos.radiation_model_cache(Y, mode::SOCRATESPrescribedRadiation)

Allocate the heating-rate field and the flux accumulators the diagnostics read.

The two accumulators stay zero: a prescribed heating rate carries no flux, which is also
how `ClimaAtmos.RadiationTRMM_LBA` behaves.
"""
function CA.radiation_model_cache(Y, mode::SOCRATESPrescribedRadiation)
    FT = CC.Spaces.undertype(axes(Y.c))
    return (;
        ᶜdTdt_rad = similar(Y.c, FT),
        net_energy_flux_toa = [CC.Geometry.WVector(FT(0))],
        net_energy_flux_sfc = [CC.Geometry.WVector(FT(0))],
    )
end

"""
    ClimaAtmos.radiation_tendency!(Yₜ, Y, p, t, mode::SOCRATESPrescribedRadiation)

Add the prescribed heating to `Yₜ.c.ρe_tot`, converted with the moist isochoric heat
capacity, `ρ cv_m dT/dt`. Applied at every stage, so there is no refresh interval.
"""
function CA.radiation_tendency!(Yₜ, Y, p, t, mode::SOCRATESPrescribedRadiation)
    thermo = CA.Parameters.thermodynamics_params(p.params)
    ᶜdTdt_rad = p.radiation.ᶜdTdt_rad
    _fill_column!(ᶜdTdt_rad, mode.dTdt_rad, _t_seconds(t))
    (; ᶜq_tot_nonneg, ᶜq_liq, ᶜq_ice) = p.precomputed
    @. Yₜ.c.ρe_tot +=
        Y.c.ρ * TD.cv_m(thermo, ᶜq_tot_nonneg, ᶜq_liq, ᶜq_ice) * ᶜdTdt_rad
    return nothing
end
