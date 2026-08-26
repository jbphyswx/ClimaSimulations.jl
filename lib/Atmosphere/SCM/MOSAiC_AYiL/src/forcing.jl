"""
    forcing.jl

The MOSAiC AYiL large-scale forcing: a forcing type owned by this package, with
`external_forcing_cache` / `external_forcing_tendency!` methods composing
ClimaAtmos's shared forcing kernels.

The closure is specified in `docs/design.md` §7, which cites the DALES source for
every term. Two parts of it cannot be expressed by ClimaAtmos's stock forcing
terms and are why this type exists: the relaxation rate is anchored on an
inversion height re-diagnosed from the model state every step, and the horizontal
advection of momentum has no canonical column variable.
"""

"""Lower bound of the inversion search [m] (DALES `tb_minzinv`)."""
const INVERSION_SEARCH_MIN = 100.0

"""Upper bound of the inversion search [m] (DALES `tb_maxzinv`)."""
const INVERSION_SEARCH_MAX = 5000.0

"""
Total water below which moisture is relaxed within one step [kg/kg] (DALES
`qtthres`, `modtestbed.f90:1551`).
"""
const DRY_AIR_NUDGE_THRESHOLD = 1.0e-6

"""
    MOSAiCForcing(FT, case; kwargs...)

Large-scale forcing for one AYiL day: relaxation of temperature, total water and
horizontal wind toward the ERA5 testbed profiles above a diagnosed inversion,
horizontal advection of heat, moisture and momentum, and large-scale subsidence.

The profiles are held as plain vectors on the `scm_in` levels and sampled onto the
model grid once in `ClimaAtmos.external_forcing_cache`. They carry no time axis:
each `scm_in` file is a single 05:00–11:00 UTC composite average, written twice
with bitwise-identical records, so `time_index` selects a record rather than
bracketing an interval (`docs/design.md` §1).

The geostrophic wind is not part of this term. DALES applies it as a pressure
gradient balanced against Coriolis, netting `-f × (u - u_g)`, which is what
ClimaAtmos's `scm_coriolis` computes; [`mosaic_scm_coriolis`](@ref) builds it.

# Keyword Arguments

  - `root`, `time_index`, `thermo_params`: passed to [`read_scm_in`](@ref).
  - `nudging`: `(; timescale, ramp_depth, z_min)`, by default the case's own
    namelist values.
  - `z_inv_min`, `z_inv_max`: bounds of the inversion search [m].
  - `q_tot_threshold`: total water below which moisture relaxes within one step
    [kg/kg]; pass `0` to disable that branch.
"""
struct MOSAiCForcing{FT, V <: AbstractVector{FT}}
    z::V
    ta::V
    hus::V
    ua::V
    va::V
    wa::V
    dTdt_hadv::V
    dqtdt_hadv::V
    dudt_hadv::V
    dvdt_hadv::V
    inv_τ::FT
    ramp_depth::FT
    z_inv_min::FT
    z_inv_max::FT
    q_tot_threshold::FT
end

function MOSAiCForcing(
    ::Type{FT},
    c::MOSAiCAYiLCase;
    root = data_root(),
    time_index::Int = 1,
    forcing = read_scm_in(c.date; root, time_index),
    nudging = nudging_parameters(c; root),
    z_inv_min::Real = INVERSION_SEARCH_MIN,
    z_inv_max::Real = INVERSION_SEARCH_MAX,
    q_tot_threshold::Real = DRY_AIR_NUDGE_THRESHOLD,
) where {FT <: AbstractFloat}
    # A positive onset height would be DALES's fixed-ramp mode, which this term
    # does not implement; it would be silently ignored rather than honoured.
    nudging.z_min < 0 || error(
        "`nudging.z_min = $(nudging.z_min)` asks for a fixed nudging onset \
         height. Only the diagnosed-inversion mode is implemented.",
    )

    f = forcing
    v(x) = collect(FT, x)
    return MOSAiCForcing{FT, Vector{FT}}(
        v(f.z),
        v(f.ta),
        v(f.hus),
        v(f.ua),
        v(f.va),
        v(f.wa),
        v(f.tntha),
        v(f.tnhusha),
        v(f.tnua),
        v(f.tnva),
        FT(1 / nudging.timescale),
        FT(nudging.ramp_depth),
        FT(z_inv_min),
        FT(z_inv_max),
        FT(q_tot_threshold),
    )
end

# A ClimaCore broadcast is pointwise over the space and never exposes the level
# index, so a sampled profile is written level by level. The guard requires one
# level per value and one column: for a wider geometry `length(field)` exceeds the
# level count, and writing one column's profile into a flat index would be wrong.
function _fill_column!(field, values)
    n = length(values)
    (CC.Spaces.nlevels(axes(field)) == n && length(field) == n) || error(
        "The MOSAiC forcing fills a single column of $n levels; got a field with \
         $(CC.Spaces.nlevels(axes(field))) levels and $(length(field)) points.",
    )
    FT = eltype(field)
    @inbounds for k in 1:n
        field[k] = FT(values[k])
    end
    return nothing
end

"""
    ClimaAtmos.external_forcing_cache(Y, forcing::MOSAiCForcing, params, start_date)

Sample the testbed profiles onto the model levels and allocate the working fields.

The profiles are constant in time, so they are sampled once here and never
refreshed; the per-step work is the inversion diagnosis and the two rate fields.
"""
function CA.external_forcing_cache(Y, forcing::MOSAiCForcing, params, start_date)
    FT = CC.Spaces.undertype(axes(Y.c))
    z = vec(Array(parent(CC.Fields.coordinate_field(Y.c).z)))
    sampled(values) = begin
        field = similar(Y.c, FT)
        _fill_column!(field, CA.interp_vertical_prof(z, forcing.z, values))
        field
    end
    thermo = CA.Parameters.thermodynamics_params(params)
    return (;
        z = collect(FT, z),
        Lv_over_cp = FT(
            TD.Parameters.LH_v0(thermo) / TD.Parameters.cp_d(thermo),
        ),
        Rd_over_cp = FT(TD.Parameters.R_d(thermo) / TD.Parameters.cp_d(thermo)),
        p_ref = FT(TD.Parameters.p_ref_theta(thermo)),
        ᶜT_nudge = sampled(forcing.ta),
        ᶜqt_nudge = sampled(forcing.hus),
        ᶜu_nudge = sampled(forcing.ua),
        ᶜv_nudge = sampled(forcing.va),
        ᶜls_subsidence = sampled(forcing.wa),
        ᶜdTdt_hadv = sampled(forcing.dTdt_hadv),
        ᶜdqtdt_hadv = sampled(forcing.dqtdt_hadv),
        ᶜdudt_hadv = sampled(forcing.dudt_hadv),
        ᶜdvdt_hadv = sampled(forcing.dvdt_hadv),
        ᶜθ_l = similar(Y.c, FT),
        ᶜinv_τ = similar(Y.c, FT),
        ᶜinv_τ_q = similar(Y.c, FT),
    )
end

"""
    inversion_height(ᶜθ_l, z, z_min, z_max)

The height [m] of the largest `∂θ_l/∂z` on the levels with `z_min < z < z_max`,
by the centred difference `(θ_l[k+1] - θ_l[k-1]) / (z[k+1] - z[k-1])`.

`0` when no level lies in the window. The centred difference places the maximum
one level below a sharp jump in `θ_l`, which is a property of the DALES algorithm
this reproduces (`modtestbed.f90:1521-1543`).
"""
function inversion_height(ᶜθ_l, z, z_min, z_max)
    n = length(z)
    FT = eltype(ᶜθ_l)
    best = typemin(FT)
    k_inv = 0
    @inbounds for k in 2:(n - 1)
        z[k] > z_min || continue
        z[k] < z_max || break
        gradient = (ᶜθ_l[k + 1] - ᶜθ_l[k - 1]) / FT(z[k + 1] - z[k - 1])
        if gradient > best
            best = gradient
            k_inv = k
        end
    end
    return k_inv == 0 ? zero(FT) : FT(z[k_inv])
end

"""
    nudge_ramp(z, z_inv, z_mid)

The DALES relaxation shape: `0` at and below `z_inv`, rising linearly to `1` at
`z_mid`, and `1` above. A zero-depth ramp is a step at `z_inv`.
"""
nudge_ramp(z, z_inv, z_mid) =
    z <= z_inv ? zero(z) :
    (z_mid > z_inv ? min(one(z), (z - z_inv) / (z_mid - z_inv)) : one(z))

"""
    ClimaAtmos.external_forcing_tendency!(Yₜ, Y, p, t, forcing::MOSAiCForcing)

Apply the AYiL forcing: relaxation above the inversion diagnosed from the current
state, horizontal advection of heat, moisture and momentum, and subsidence.
"""
function CA.external_forcing_tendency!(Yₜ, Y, p, t, forcing::MOSAiCForcing)
    (;
        z,
        Lv_over_cp,
        Rd_over_cp,
        p_ref,
        ᶜT_nudge,
        ᶜqt_nudge,
        ᶜu_nudge,
        ᶜv_nudge,
        ᶜls_subsidence,
        ᶜdTdt_hadv,
        ᶜdqtdt_hadv,
        ᶜdudt_hadv,
        ᶜdvdt_hadv,
        ᶜθ_l,
        ᶜinv_τ,
        ᶜinv_τ_q,
    ) = p.external_forcing
    (; ᶜT, ᶜp, ᶜq_liq) = p.precomputed
    FT = eltype(ᶜθ_l)
    thermo = CA.Parameters.thermodynamics_params(p.params)
    cp_d = TD.Parameters.cp_d(thermo)
    L_v = TD.Parameters.LH_v0(thermo)
    R_d = TD.Parameters.R_d(thermo)
    p_ref = TD.Parameters.p_ref_theta(thermo)

    @. ᶜθ_l = (ᶜT - FT(L_v / cp_d) * ᶜq_liq) * (FT(p_ref) / ᶜp)^FT(R_d / cp_d)
    z_inv = inversion_height(ᶜθ_l, z, forcing.z_inv_min, forcing.z_inv_max)
    z_mid = z_inv + forcing.ramp_depth

    # DALES never relaxes faster than one timestep, and relaxes moisture *at* the
    # timestep where the column is drier than the threshold.
    inv_τ = min(forcing.inv_τ, inv(p.dt))
    inv_τ_dry = inv(p.dt)
    ᶜz = CC.Fields.coordinate_field(axes(Y.c)).z
    ᶜq_tot = @. lazy(Y.c.ρq_tot / Y.c.ρ)
    @. ᶜinv_τ = nudge_ramp(ᶜz, z_inv, z_mid) * inv_τ
    @. ᶜinv_τ_q =
        nudge_ramp(ᶜz, z_inv, z_mid) *
        ifelse(ᶜq_tot < forcing.q_tot_threshold, inv_τ_dry, inv_τ)

    CA.nudge_uv!(Yₜ, Y, p, ᶜu_nudge, ᶜv_nudge, ᶜinv_τ)

    ᶜlg = CC.Fields.local_geometry_field(Y.c)
    ᶜuₕ_hadv = p.scratch.ᶜtemp_C12
    @. ᶜuₕ_hadv = CA.C12(CC.Geometry.UVVector(ᶜdudt_hadv, ᶜdvdt_hadv), ᶜlg)
    @. Yₜ.c.uₕ += ᶜuₕ_hadv

    # `nudge_Tq!` carries one rate for both scalars; the dry-air branch gives
    # moisture its own.
    ᶜdTdt = p.scratch.ᶜtemp_scalar
    ᶜdqtdt = p.scratch.ᶜtemp_scalar_2
    @. ᶜdTdt = -(ᶜT - ᶜT_nudge) * ᶜinv_τ + ᶜdTdt_hadv
    @. ᶜdqtdt = -(ᶜq_tot - ᶜqt_nudge) * ᶜinv_τ_q + ᶜdqtdt_hadv
    CA.apply_Tq_forcing!(Yₜ, Y, p, ᶜdTdt, ᶜdqtdt)

    CA.apply_subsidence_forcing!(Yₜ, Y, p, ᶜls_subsidence)
    return nothing
end

# Linear in height, flat outside the data — the interpolation `ColumnProfiles` and
# `interp_vertical_prof` use, as a scalar callable of `z`.
_column_profile(z, values) = Intp.extrapolate(
    Intp.interpolate((z,), values, Intp.Gridded(Intp.Linear())),
    Intp.Flat(),
)

"""
    mosaic_scm_coriolis(FT, case; params, root, time_index, forcing, latitude, omega)

ClimaAtmos's `scm_coriolis` for one AYiL day: the geostrophic wind profiles from
`scm_in` and the Coriolis parameter at the day's drift position.

Together with Coriolis this is `-f × (u - u_g)`, the net of DALES's geostrophic
pressure gradient and Coriolis terms (`docs/design.md` §7).
"""
function mosaic_scm_coriolis(
    ::Type{FT},
    c::MOSAiCAYiLCase;
    params,
    root = data_root(),
    time_index::Int = 1,
    forcing = read_scm_in(c.date; root, time_index),
    latitude::Real = forcing.surface.trajectory_latitude,
    omega = CA.Parameters.Omega(params),
) where {FT <: AbstractFloat}
    z = collect(FT, forcing.z)
    return (;
        prof_ug = _column_profile(z, collect(FT, forcing.ug)),
        prof_vg = _column_profile(z, collect(FT, forcing.vg)),
        coriolis_param = FT(2 * omega * sind(latitude)),
    )
end
