"""
    transport_diagnostics.jl

Vertical transport terms of the microphysics budgets.

ClimaAtmos registers the process rates but not the sedimentation, vertical advection and vertical
diffusion of the same tracers, so a budget assembled from its diagnostics alone does not close. These
register the missing terms through `add_diagnostic_variable!`, computing each with the model's own
operators on the model's own state, so each is the term the model applied.
"""

"""Tracer suffix → (`ρq` field, precomputed fall speed, whether `K_h` reaches it through `q_tot`)."""
const TRANSPORT_SPECIES = (
    ("q_lcl", :ρq_lcl, :ᶜwₗ, true),
    ("q_icl", :ρq_icl, :ᶜwᵢ, true),
    ("q_rai", :ρq_rai, :ᶜwᵣ, false),
    ("q_sno", :ρq_sno, :ᶜwₛ, false),
)

"""
    sedimentation_tendency(state, cache, ρq_name, w_name)

`∂q/∂t` from sedimentation [kg kg^-1 s^-1].

`implicit_vertical_advection_tendency!` subtracts this flux divergence from `ρq`, so the sign is
negative and the result is divided by `ρ`.
"""
function sedimentation_tendency(state, cache, ρq_name::Symbol, w_name::Symbol)
    ᶜρ = state.c.ρ
    ᶜJ = CC.Fields.local_geometry_field(axes(state.c)).J
    ᶠJ = CC.Fields.local_geometry_field(axes(state.f)).J
    ᶜw = getproperty(cache.precomputed, w_name)
    ᶜρq = getproperty(state.c, ρq_name)
    return @. CA.lazy(
        -CA.ᶜprecipdivᵥ(
            CA.ᶠinterp(ᶜρ * ᶜJ) / ᶠJ *
            CA.ᶠright_bias(CC.Geometry.WVector(-(ᶜw)) * CA.specific(ᶜρq, ᶜρ)),
        ) / ᶜρ,
    )
end

"""
    vertical_advection_tendency(state, cache, ρq_name)

`∂q/∂t` from resolved vertical advection [kg kg^-1 s^-1], using the `vertical_transport` and
`tracer_upwinding` the explicit tendency uses.
"""
function vertical_advection_tendency(state, cache, ρq_name::Symbol)
    ᶜρ = state.c.ρ
    ᶜρq = getproperty(state.c, ρq_name)
    ᶜχ = @. CA.lazy(CA.specific(ᶜρq, ᶜρ))
    vtt = CA.vertical_transport(
        ᶜρ,
        cache.precomputed.ᶠu³,
        ᶜχ,
        cache.dt,
        cache.atmos.numerics.tracer_upwinding,
    )
    return @. CA.lazy(vtt / ᶜρ)
end

"""
    vertical_diffusion_tendency(state, cache, ρq_name, gets_K_h)

`∂q/∂t` from SGS vertical diffusion [kg kg^-1 s^-1].

`edmfx_sgs_diffusive_flux_tendency!` diffuses a microphysics tracer in two pieces: the unified tracer
loop applies `K_entr` to the species' own gradient, and for cloud species only, the `K_h` diffusion of
`q_tot_eff` is distributed by the clipped mass ratio `q/q_tot_eff`. Rain and snow receive no `K_h`
transport, which `gets_K_h` selects.

Errors when the configuration routes tracer diffusion through a path this does not reproduce, so a
term this cannot compute never reads as an absent one.
"""
function vertical_diffusion_tendency(state, cache, ρq_name::Symbol, gets_K_h::Bool)
    cache.atmos.edmfx_model.sgs_diffusive_flux isa Val{true} || error(
        "sgs_diffusive_flux is off, so `edmfx_sgs_diffusive_flux_tendency!` applies no tracer \
         diffusion and this diagnostic would not be the term the model used.",
    )
    isnothing(cache.atmos.vertical_diffusion) || error(
        "vertical_diffusion = $(cache.atmos.vertical_diffusion) also diffuses tracers through \
         `vertical_diffusion_boundary_layer_tendency!`, which this diagnostic does not include.",
    )
    ᶜρ = state.c.ρ
    FT = eltype(ᶜρ)
    (; ᶠK_h, ᶠK_entr) = cache.precomputed
    ᶜρq = getproperty(state.c, ρq_name)
    ᶜχ = @. CA.lazy(CA.specific(ᶜρq, ᶜρ))
    ᶠρK_e = @. CA.lazy(CA.ᶠinterp(ᶜρ) * ᶠK_entr)
    ᶜflux = CA.ᶜdiffusive_flux_divergenceᵥ(ᶠρK_e, ᶜχ)
    gets_K_h || return @. CA.lazy(-(ᶜflux) / ᶜρ)

    ᶠρK_h = @. CA.lazy(CA.ᶠinterp(ᶜρ) * ᶠK_h)
    ᶜq_tot_eff = CA.ᶜdiffusing_water(state, cache)
    ᶜratio = @. CA.lazy(
        max(zero(FT), min(one(FT), CA.specific(ᶜρq, ᶜρ) / max(ᶜq_tot_eff, eps(FT)))),
    )
    ᶜflux_h = CA.ᶜdiffusive_flux_divergenceᵥ(ᶠρK_h, ᶜq_tot_eff)
    return @. CA.lazy(-(ᶜflux + ᶜratio * ᶜflux_h) / ᶜρ)
end

"""
    register_transport_diagnostics!()

Register `sed_`, `adv_` and `dif_` for every microphysics tracer.

Any previous entry is dropped first, so re-including this file replaces the closure registered
earlier in the session.
"""
function register_transport_diagnostics!()
    for (species, ρq_name, w_name, gets_K_h) in TRANSPORT_SPECIES
        for (prefix, description, compute) in (
            (
                "sed",
                "sedimentation",
                (state, cache, _) -> sedimentation_tendency(state, cache, ρq_name, w_name),
            ),
            (
                "adv",
                "resolved vertical advection",
                (state, cache, _) -> vertical_advection_tendency(state, cache, ρq_name),
            ),
            (
                "dif",
                "SGS vertical diffusion",
                (state, cache, _) ->
                    vertical_diffusion_tendency(state, cache, ρq_name, gets_K_h),
            ),
        )
            name = "$(prefix)_$(species)"
            delete!(CA.Diagnostics.ALL_DIAGNOSTICS, name)
            CA.Diagnostics.add_diagnostic_variable!(;
                short_name = name,
                long_name = "Tendency of $species from $description",
                units = "kg kg^-1 s^-1",
                comments = "Computed with the model's own vertical operators on the model state.",
                compute,
            )
        end
    end
    return nothing
end
