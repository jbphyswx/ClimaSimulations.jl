"""
    les_rates.jl

Atlas LES process and transport rates under the model's own tendency names. Units and the
mixing-ratio conversion come from [`atlas_var_specs`](@ref); this file only says which
Atlas variables make up each model rate.
"""

"""
Model rate name => the Atlas variables summing to it, with signs.

The process rates are the **unsigned** magnitudes ClimaAtmos registers;
[`MP1M_BUDGETS`](@ref) supplies each one's sign within a given budget. Atlas's own family
signs are uniform — growth and collection positive, evaporation, sublimation and melting
negative — so only melting is negated here to recover a magnitude. The transport terms are
already `∂q/∂t` and keep the archive's sign.

`sed_q_lcl` and `adv_q_lcl` read the derived `QCSED` and `QCADV`: `QTO` is vapour plus
cloud liquid and vapour does not sediment, so `QTOSED` is cloud liquid's sedimentation,
and Atlas reports cloud liquid's vertical advection as a flux rather than a tendency.

Absent: `dif_q_lcl`, for which the archive has no cloud-liquid diffusion; the freezing
terms, which need a nucleation the 1-moment scheme lacks; and `S_accr_icl_rai`,
`S_accr_lcl_sno_warm`, `S_accr_rai_sno_warm`, `S_accr_freeze_icl_rai`, `S_melt_icl_lcl`
and `S_accr_melt_*`, which have no single Atlas counterpart.
"""
const LES_RATE_SOURCES = Dict{String, Vector{Tuple{Symbol, Int}}}(
    "sed_q_lcl" => [(:QCSED, +1)],
    "sed_q_icl" => [(:QISED, +1)],
    "sed_q_rai" => [(:QRSED, +1)],
    "sed_q_sno" => [(:QSSED, +1)],
    "adv_q_lcl" => [(:QCADV, +1)],
    "adv_q_icl" => [(:QIADV, +1)],
    "adv_q_rai" => [(:QRADV, +1)],
    "adv_q_sno" => [(:QSADV, +1)],
    "dif_q_icl" => [(:QIDIFF, +1)],
    "dif_q_rai" => [(:QRDIFF, +1)],
    "dif_q_sno" => [(:QSDIFF, +1)],
    "S_phase_change_vap_lcl" => [(:PCC, +1)],
    "S_phase_change_vap_icl" => [(:PRDT, +1)],
    "S_phase_change_vap_sno" => [(:PRDST, +1)],
    "S_phase_change_vap_rai" => [(:PRE, +1)],
    "S_acnv_lcl_rai" => [(:PRC, +1)],
    "S_acnv_icl_sno" => [(:PRCI, +1), (:PITOSN, +1)],
    "S_accr_lcl_rai" => [(:PRA, +1)],
    "S_accr_lcl_sno_cold" => [(:PSACWS, +1)],
    "S_accr_icl_sno" => [(:PRAI, +1)],
    "S_accr_rai_sno_cold" => [(:PRACS, +1)],
    "S_melt_sno_rai" => [(:PSMLT, -1)], # Atlas signs melting negative; the budget wants the magnitude
)

"""Every Atlas variable [`LES_RATE_SOURCES`](@ref) names."""
les_rate_variables() =
    unique(v for terms in values(LES_RATE_SOURCES) for (v, _) in terms)

"""
    les_rate_profiles(case; params)

`(; z, time, data)` of every Atlas variable [`LES_RATE_SOURCES`](@ref) names, as
specific-content tendencies [kg kg^-1 s^-1] with `time` in elapsed seconds.
"""
les_rate_profiles(c::SocratesCase; params) =
    read_atlas(c, les_rate_variables(); params)

"""
    case_les_rates(case, z; params, window)

Atlas rates for `case` [kg kg^-1 s^-1], averaged over `window` and resampled onto `z`,
keyed by the model rate name.
"""
function case_les_rates(
    c::SocratesCase,
    z::AbstractVector;
    params,
    window = score_window(c),
)
    raw = les_rate_profiles(c; params)
    keep = findall(t -> first(window) <= t <= last(window), raw.time)
    isempty(keep) && error(
        "No LES times inside $(window) s for $(case_name(c)); the record spans \
         $(extrema(raw.time)) s.",
    )
    FT = eltype(raw.z)
    # linear in height with flat ends: the Atlas grid is finer than any model grid here
    onto(profile) = map(z) do zi
        k = searchsortedfirst(raw.z, zi)
        k <= 1 && return profile[1]
        k > length(raw.z) && return profile[end]
        z0, z1 = raw.z[k - 1], raw.z[k]
        w = (zi - z0) / (z1 - z0)
        (1 - w) * profile[k - 1] + w * profile[k]
    end
    rates = Dict{String, Vector{FT}}()
    for (rate, terms) in LES_RATE_SOURCES
        total = zeros(FT, length(raw.z))
        for (name, sign) in terms
            a = raw.data[name]
            total .+= sign .* [FT(nanmean(view(a, k, keep))) for k in axes(a, 1)]
        end
        rates[rate] = collect(FT, onto(total))
    end
    return rates
end
