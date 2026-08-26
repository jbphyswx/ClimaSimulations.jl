"""
    atlas_archive.jl

The unit conventions of the Atlas LES archive, and the conversions onto SI and specific
humidity. [`atlas_registry.jl`](@ref) classifies every variable in terms of these.

SAM writes water as **mixing ratios**, per kilogram of dry air, on a time axis in days.
With `q_tot = w_all / (1 + w_all)` the specific total water, `1 - q_tot` is the dry-air
fraction, so a per-dry-air quantity becomes a per-total-mass one by multiplying by
`(1 - q_tot)`, and each water factor in a quantity contributes one such factor. Because
`q_tot` is held fixed, a *tendency* converts by the same factor as the quantity — which
is why `dwcdt_to_dqcdt` is `wc_to_qc`. Total water itself takes two factors, since
`dq_t/dw_t = (1 - q_t)^2`.

Number quantities carry no water factor: their per-kilogram basis is not established by
the Atlas metadata, so only their unambiguous time and volume units are converted.
"""

const SECONDS_PER_DAY = 86400

"""Sentinel the raw Atlas files declare for missing data."""
const ATLAS_MISSING_VALUE = -9999.0

"""
`SHFOBS`/`LHFOBS` carry this where they have no observation. It is not the declared
`missing_value`, so nothing masks it and it would otherwise be written as a flux.
"""
const ATLAS_UNDECLARED_SENTINEL = -99999.0

"""
    atlas_values(a)

`a` with both spellings of an absent Atlas sample — `missing` and the
[`ATLAS_MISSING_VALUE`](@ref) sentinel — mapped to `NaN`, in the archive's element type.
"""
function atlas_values(a::AbstractArray)
    FT = nonmissingtype(eltype(a))
    out = Array{FT}(undef, size(a))
    @inbounds for i in eachindex(a, out)
        v = a[i]
        out[i] = (v === missing || v == ATLAS_MISSING_VALUE) ? FT(NaN) : FT(v)
    end
    return out
end

"""
    ComposeFirst(f, g)

Callable applying `g` to the first argument only, then `f` to that result and any
remaining arguments: `ComposeFirst(f, g)(x, args...) == f(g(x), args...)`.

`∘` cannot express this, because the extra arguments must bypass `g`. It lets a unit
scaling compose with a conversion that additionally needs the local state, as in
`ComposeFirst(wc_to_qc, g_to_kg)(w_gkg, q_tot)`.
"""
struct ComposeFirst{F, G}
    f::F
    g::G
end
@inline (h::ComposeFirst)(x, args...) = h.f(h.g(x), args...)

"""A time axis rebased on its own first sample, so `t = 0` is the start of the run."""
_elapsed(t) = t .- first(t)

# --- unit scalings ---------------------------------------------------------- #

g_to_kg(x) = x / oftype(x, 1000)
perday_to_persec(x) = x / oftype(x, SECONDS_PER_DAY)
days_to_seconds(x) = x * oftype(x, SECONDS_PER_DAY)
hPa_to_Pa(x) = x * oftype(x, 100)
percm3_to_perm3(x) = x * oftype(x, 1.0e6)
km_to_m(x) = x * oftype(x, 1000)
km2_to_m2(x) = x * oftype(x, 1.0e6)
micron_to_m(x) = x * oftype(x, 1.0e-6)
permicron_to_perm(x) = x * oftype(x, 1.0e6)
percent_to_fraction(x) = x / oftype(x, 100)
g2kg2_to_kg2kg2(x) = x / oftype(x, 1_000_000)

"""Mask the undeclared `-99999` sentinel of the observation series."""
mask_undeclared_sentinel(x) =
    x == oftype(x, ATLAS_UNDECLARED_SENTINEL) ? oftype(x, NaN) : x

"""SAM stores static energies divided by `c_pd`, so `K` becomes J/kg."""
K_to_J_kg(x, cp_d) = x * oftype(x, cp_d)

"""Latent heat flux [W m^-2] to the liquid mass flux [kg m^-2 s^-1] carrying it."""
wm2_to_kgm2s_liq(x, L_v) = x / oftype(x, L_v)

"""Latent heat flux [W m^-2] to the ice mass flux [kg m^-2 s^-1] carrying it."""
wm2_to_kgm2s_ice(x, L_s) = x / oftype(x, L_s)

# --- mixing ratio to specific humidity -------------------------------------- #

"""Total-water mixing ratio to total-water specific humidity."""
wt_to_qt(w) = w / (one(w) + w)

"""Total-water specific humidity to total-water mixing ratio."""
qt_to_wt(q) = q / (one(q) - q)

"""A condensate mixing ratio to a specific content, against the column's total water."""
wc_to_qc(w, q_tot) = w * (one(q_tot) - q_tot)

"""A condensate specific content to a mixing ratio."""
qc_to_wc(q, q_tot) = q / (one(q_tot) - q_tot)

"""
A condensate mixing-ratio tendency to a specific-content tendency.

`q_tot` is fixed with respect to the derivative, so a tendency converts by the same
factor as the quantity does.
"""
const dwcdt_to_dqcdt = wc_to_qc

"""A total-water mixing-ratio tendency to a specific tendency; `dq_t/dw_t = (1 - q_t)^2`."""
dwtdt_to_dqtdt(dwdt, q_tot) = dwdt * (one(q_tot) - q_tot)^2

"""A condensate `(g/kg)^2` moment to `(kg/kg)^2`: two water factors."""
wc2_g2kg2_to_kg2kg2(var, q_tot) = g2kg2_to_kg2kg2(var) * (one(q_tot) - q_tot)^2

"""A total-water `(g/kg)^2` moment to `(kg/kg)^2`: four water factors."""
wt2_g2kg2_to_kg2kg2(var, q_tot) = g2kg2_to_kg2kg2(var) * (one(q_tot) - q_tot)^4

# Every conversion above is scalar; the registry applies them to whole arrays.
for f in (
    :g_to_kg, :perday_to_persec, :days_to_seconds, :hPa_to_Pa, :percm3_to_perm3,
    :km_to_m, :km2_to_m2, :micron_to_m, :permicron_to_perm, :percent_to_fraction,
    :g2kg2_to_kg2kg2, :mask_undeclared_sentinel, :wt_to_qt, :qt_to_wt,
)
    @eval $f(x::AbstractArray) = $f.(x)
end
for f in (:K_to_J_kg, :wm2_to_kgm2s_liq, :wm2_to_kgm2s_ice)
    @eval $f(x::AbstractArray, c) = $f.(x, c)
end
for f in (:wc_to_qc, :qc_to_wc, :dwtdt_to_dqtdt, :wc2_g2kg2_to_kg2kg2, :wt2_g2kg2_to_kg2kg2)
    @eval $f(x::AbstractArray, q_tot::AbstractArray) = $f.(x, q_tot)
end

# --- vertical operator ------------------------------------------------------ #

"""
    atlas_forward_flux_divergence_tendency(F, z, ρ)

Vertical tendency [kg kg^-1 s^-1] from a mass flux `F` [kg m^-2 s^-1]:
`-(F[k+1] - F[k]) / ((z[k+1] - z[k]) ρ[k])`, the top level copying the cell below.

The forward stencil is the one that reproduces Atlas's own `*SED` from its `*SDFLX`; no
centred or backward form does. Each difference is taken in `Float64` because neighbouring
fluxes are nearly equal, so differencing them at the archive's own precision would lose
most of the result's significant digits.
"""
function atlas_forward_flux_divergence_tendency(
    F::AbstractMatrix,
    z::AbstractVector,
    ρ::AbstractMatrix,
)
    nz, nt = size(F)
    length(z) == nz || throw(ArgumentError("length(z) must match size(F, 1)"))
    size(ρ) == size(F) || throw(ArgumentError("ρ must match F"))
    FT = eltype(F)
    out = similar(F)
    for t in 1:nt
        for k in 1:(nz - 1)
            dz = Float64(z[k + 1]) - Float64(z[k])
            out[k, t] =
                FT(-(Float64(F[k + 1, t]) - Float64(F[k, t])) / (dz * Float64(ρ[k, t])))
        end
        out[nz, t] = out[nz - 1, t]
    end
    return out
end