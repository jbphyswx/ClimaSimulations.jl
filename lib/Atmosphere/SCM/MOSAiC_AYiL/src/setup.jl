"""
    setup.jl

The MOSAiC AYiL initial condition: the state the reference DALES run started
from, on the model levels, with the condensate split as DALES splits it
(`docs/design.md` §6).
"""

"""
    dales_tke_seed(z; e12_min = DALES_CONSTANTS.e12_min, decay_length = 50)

DALES's cold-start turbulence seed as a specific kinetic energy [m²/s²]:
`(e12_min + exp(-z / decay_length))²`, from `modstartup.f90:449`.

An LES seeds a *subgrid* kinetic energy and grows the rest from resolved
perturbations a single column cannot have, while a single-column model's TKE
stands for the whole of the turbulence. So this is a spin-up choice with no
reference value: it is the default because it is the seed the reference used, and
it is a keyword because it is a choice.
"""
dales_tke_seed(z; e12_min = DALES_CONSTANTS.e12_min, decay_length = 50) =
    (e12_min + exp(-z / decay_length))^2

"""
    MOSAiCSetup(FT, case; kwargs...)

Initial state of one AYiL day, as a profile per prognostic quantity.

Every entry is a callable `z -> value`, which is what
[`Setups.center_initial_condition`](@ref) needs: it is handed one level's local
geometry at a time, never an index into a vector. Each is built once over the grid
its data arrives on, so a value is interpolated from the source and not from a
resampling of it.

Temperature, wind and moisture come from `scm_in`, which is what DALES initialized
from. The density is the reference run's own `rhof` ([`les_density`](@ref)) — its
real pressure written as a density, `presf/(R_d T_v)` — so, since the model
diagnoses pressure from `ρ` and `T`, the column starts at the pressure the run had.
`docs/design.md` §8 records what the alternatives cost. `q_ice` is prescribed rather
than diagnosed because the condensate split needs it (below), and `tke` is a
cold-start seed with no counterpart in the archive.

# Keyword Arguments

  - `root`, `time_index`: passed to [`read_scm_in`](@ref).
  - `forcing`: the `scm_in` state; pass it to share one read of the file.
  - `density`: `(z, ρ)` in [m] and [kg/m³]; [`les_density`](@ref) by default.
  - `tke`: callable `z -> TKE` [m²/s²]; [`dales_tke_seed`](@ref) by default.
"""
struct MOSAiCSetup{P <: NamedTuple}
    profiles::P
end

function MOSAiCSetup(
    ::Type{FT},
    c::MOSAiCAYiLCase;
    root = data_root(),
    time_index::Int = 1,
    forcing = read_scm_in(c.date; root, time_index),
    density = les_density(c.date; root),
    tke = dales_tke_seed,
) where {FT <: AbstractFloat}
    z_scm = collect(FT, forcing.z)
    from_scm(values) = _column_profile(z_scm, collect(FT, values))
    z_ρ, ρ = density
    return MOSAiCSetup((;
        T = from_scm(forcing.ta),
        u = from_scm(forcing.ua),
        v = from_scm(forcing.va),
        q_tot = from_scm(forcing.hus),
        q_ice = from_scm(forcing.qi),
        ρ = _column_profile(collect(FT, z_ρ), collect(FT, ρ)),
        tke,
    ))
end

"""
    ClimaAtmos.Setups.center_initial_condition(setup::MOSAiCSetup, local_geometry, params)

The initial state at one level, with cloud ice prescribed and the liquid
saturating the vapour-plus-liquid pool alone.

Ice is excluded from the pool the liquid adjustment sees, matching
`initclouds3`/`initiceclouds3`: adjusting the whole of `q_tot` toward liquid
saturation while separately pinning `q_ice` would count the ice twice.
"""
function CA.Setups.center_initial_condition(
    setup::MOSAiCSetup,
    local_geometry,
    params,
)
    (; z) = local_geometry.coordinates
    FT = typeof(z)
    p = setup.profiles
    T = FT(p.T(z))
    ρ = FT(p.ρ(z))
    q_tot = FT(p.q_tot(z))
    q_ice = FT(p.q_ice(z))
    thermo = CA.Parameters.thermodynamics_params(params)
    q_sat_liq = FT(TD.q_vap_saturation(thermo, T, ρ, TD.Liquid()))
    return CA.Setups.physical_state(;
        T,
        ρ,
        q_tot,
        q_liq = max(zero(FT), (q_tot - q_ice) - q_sat_liq),
        q_ice,
        u = FT(p.u(z)),
        v = FT(p.v(z)),
        tke = FT(p.tke(z)),
    )
end



# ============================================================================================= #