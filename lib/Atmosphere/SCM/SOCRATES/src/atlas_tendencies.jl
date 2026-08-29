"""
    les_tendencies.jl

The Atlas LES budget under the model's own tendency names: the transport terms named here,
the process rates from [`ATLAS_PROCESS_RATES`](@ref). Units and the mixing-ratio conversion
come from [`atlas_var_specs`](@ref).
"""



"""
    atlas_process_rates.jl

Every microphysics process rate of the Atlas LES archive — 41 mass rates and 36 number
rates — with the archive's own description and the ClimaAtmos tendency it feeds.

Atlas reports all 77 as `?/kg/day`, so `long_name` is the only authority for what each one
is; it is recorded here verbatim (minus the leading `NAME , `) and `test/atlas_registry.jl`
asserts it against the file. `kind` gives the numerator the `?` stands for, which fixes the units
and whether the mixing-ratio conversion applies.
"""

"""
Why an Atlas rate has no ClimaAtmos counterpart. One of these stands in `model` wherever a
`(term, sign)` pair cannot.

  - `:number`               a number rate; `CA.Diagnostics` registers process-rate
                            diagnostics for `MP1M_SOURCE_TERMS` only, and
                            `CloudMicrophysics.BulkMicrophysicsTendencies` offers
                            `InstantaneousVerbose` for `Microphysics1Moment` alone — its
                            2-moment methods return aggregated `dn_*_dt` and no source terms
  - `:graupel`              changes graupel, a category the 1-moment scheme does not carry
  - `:freezing`             droplet or rain freezing; the 1-moment scheme's only path to ice
                            is vapour deposition
  - `:secondary_ice`        ice multiplication
  - `:activation`           droplet activation
  - `:liquid_ice_collision` cloud liquid collected by cloud ice
  - `:melt_evaporation`     evaporation of meltwater on snow
"""
const NO_MODEL_COUNTERPART = (
    :number, :graupel, :freezing, :secondary_ice, :activation, :liquid_ice_collision,
    :melt_evaporation,
)

"""
Atlas process rate => `(; kind, long_name, model)`.

`kind` is `:mass` or `:number`. `model` is either `(term, sign)` naming the ClimaAtmos
tendency the rate feeds and the sign bringing it onto that term's convention, or a symbol
from [`NO_MODEL_COUNTERPART`](@ref).

Signs were measured over all 11 cases, whole field, full precision. Atlas's family signs
are *not* uniform, so none of them is assumed: deposition (`PRD`, `PRDS`) is 100% positive
against sublimation (`EPRD`, `EPRDS`) 100% negative, `PRE` is 100% negative, and the two
melting rates disagree — `PSMLT` 100% negative, `QMELTI` 100% positive. `MP1M_SOURCE_TERMS`
are non-negative transfer magnitudes except the four `S_phase_change_*`, which are signed
positive towards the condensate
(`CloudMicrophysics.BulkMicrophysicsTendencies._aggregate_tendencies`), so only `PSMLT`
needs negating.

"""
const ATLAS_PROCESS_RATES = Dict{
    Symbol,
    @NamedTuple{
        kind::Symbol,
        long_name::String,
        model::Union{Tuple{Symbol, Int}, Tuple{Vararg{Tuple{Symbol, Int}}}, Symbol},
    },
}(
    # ---- mass rates with a 1-moment counterpart ----------------------------------------
    :PCC => (kind = :mass, long_name = "COND/EVAP DROPLETS",
        model = (:S_phase_change_vap_lcl, +1)),
    :PRD => (kind = :mass, long_name = "DEP CLOUD ICE",
        model = (:S_phase_change_vap_icl, +1)),
    :EPRD => (kind = :mass, long_name = "SUBLIMATION CLOUD ICE",
        model = (:S_phase_change_vap_icl, +1)),
    :PRDS => (kind = :mass, long_name = "DEP SNOW",
        model = (:S_phase_change_vap_sno, +1)),
    :EPRDS => (kind = :mass, long_name = "SUBLIMATION SNOW",
        model = (:S_phase_change_vap_sno, +1)),
    :PRE => (kind = :mass, long_name = "EVAP OF RAIN",
        model = (:S_phase_change_vap_rai, +1)),
    :PRC => (kind = :mass, long_name = "AUTOCONVERSION DROPLETS",
        model = (:S_acnv_lcl_rai, +1)),
    :PRCI => (kind = :mass, long_name = "CHANGE Q AUTOCONVERSION CLOUD ICE BY SNOW",
        model = (:S_acnv_icl_sno, +1)),
    :PITOSN => (kind = :mass, long_name = "CHANGE Q ICE TO SNOW DUE TO THRESHOLD",
        model = (:S_acnv_icl_sno, +1)),
    :PRA => (kind = :mass, long_name = "ACCRETION DROPLETS BY RAIN",
        model = (:S_accr_lcl_rai, +1)),
    :PSACWS => (kind = :mass, long_name = "CHANGE Q DROPLET ACCRETION BY SNOW",
        model = ((:S_accr_lcl_sno_cold, +1), (:S_accr_lcl_sno_warm, +1))),
    :PRACI => (kind = :mass, long_name = "CHANGE QI, ICE-RAIN COLLECTION",
        model = (:S_accr_icl_rai, +1)),
    :PRACIS => (kind = :mass, long_name = "CHANGE QI, ICE RAIN COLLISION, ADDED TO SNOW",
        model = (:S_accr_icl_rai, +1)),
    :PIACR => (kind = :mass, long_name = "CHANGE QR, ICE-RAIN COLLECTION",
        model = (:S_accr_freeze_icl_rai, +1)),
    :PIACRS => (kind = :mass, long_name = "CHANGE QR, ICE RAIN COLLISION, ADDED TO SNOW",
        model = (:S_accr_freeze_icl_rai, +1)),
    :PRAI => (kind = :mass, long_name = "CHANGE Q ACCRETION CLOUD ICE",
        model = (:S_accr_icl_sno, +1)),
    :PRACS => (kind = :mass, long_name = "CHANGE Q RAIN-SNOW COLLECTION",
        model = ((:S_accr_rai_sno_cold, +1), (:S_accr_rai_sno_warm, +1))),
    :QMELTI => (kind = :mass,
        long_name = "CHANGE Q DUE TO MELTING OF CLOUD ICE TO FORM RAIN",
        model = (:S_melt_icl_lcl, +1)),
    :PSMLT => (kind = :mass, long_name = "CHANGE Q MELTING SNOW TO RAIN",
        model = (:S_melt_sno_rai, -1)),

    # ---- mass rates of the graupel category --------------------------------------------
    :PRDG => (kind = :mass, long_name = "DEP OF GRAUPEL", model = :graupel),
    :EPRDG => (kind = :mass, long_name = "SUB OF GRAUPEL", model = :graupel),
    :PGMLT => (kind = :mass, long_name = "CHANGE Q MELTING OF GRAUPEL", model = :graupel),
    :EVPMG => (kind = :mass,
        long_name = "CHANGE Q MELTING OF GRAUPEL AND EVAPORATION", model = :graupel),
    :PRACG => (kind = :mass, long_name = "CHANGE IN Q COLLECTION RAIN BY GRAUPEL",
        model = :graupel),
    :PSACWG => (kind = :mass, long_name = "CHANGE IN Q COLLECTION DROPLETS BY GRAUPEL",
        model = :graupel),
    :PGRACS => (kind = :mass,
        long_name = "CONVERSION Q TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW",
        model = :graupel),
    :PGSACW => (kind = :mass,
        long_name = "CONVERSION Q TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW",
        model = :graupel),
    :PSACR => (kind = :mass, long_name = "CONVERSION DUE TO COLL OF SNOW BY RAIN",
        model = :graupel),

    # ---- mass rates of processes the 1-moment scheme does not have ----------------------
    :MNUCCC => (kind = :mass, long_name = "CHANGE Q DUE TO CONTACT FREEZ DROPLETS",
        model = :freezing),
    :MNUCCD => (kind = :mass,
        long_name = "CHANGE Q FREEZING AEROSOL (PRIM ICE NUCLEATION)", model = :freezing),
    :MNUCCI => (kind = :mass, long_name = "CHANGE Q DUE TO IMMERSION FREEZ DROPLETS",
        model = :freezing),
    :MNUCCR => (kind = :mass, long_name = "CHANGE Q DUE TO CONTACT FREEZ RAIN",
        model = :freezing),
    :QHOMOC => (kind = :mass,
        long_name = "CHANGE Q DUE TO HOMOGENEOUS FREEZING OF CLOUD WATER TO FORM CLOUD ICE",
        model = :freezing),
    :QHOMOR => (kind = :mass,
        long_name = "CHANGE Q DUE TO HOMOGENEOUS FREEZING OF RAIN TO FORM GRAUPEL",
        model = :freezing),
    :QMULTS => (kind = :mass, long_name = "CHANGE Q DUE TO ICE MULT DROPLETS/SNOW",
        model = :secondary_ice),
    :QMULTG => (kind = :mass,
        long_name = "CHANGE Q DUE TO ICE MULT DROPLETS/GRAUPEL, SINK OF QC, SOURCE OF QI",
        model = :secondary_ice),
    :QMULTR => (kind = :mass, long_name = "CHANGE Q DUE TO ICE RAIN/SNOW",
        model = :secondary_ice),
    :QMULTRG => (kind = :mass,
        long_name = "CHANGE Q DUE TO ICE MULT RAIN/GRAUPEL, SINK OF QR, SOURCE OF QI",
        model = :secondary_ice),
    :PCCN => (kind = :mass, long_name = "CHANGE Q DROPLET ACTIVATION", model = :activation),
    :PSACWI => (kind = :mass, long_name = "CHANGE Q DROPLET ACCRETION BY CLOUD ICE",
        model = :liquid_ice_collision),
    :EVPMS => (kind = :mass, long_name = "CHNAGE Q MELTING SNOW EVAPORATING",
        model = :melt_evaporation),

    # ---- number rates ------------------------------------------------------------------
    :NPRC => (kind = :number, long_name = "CHANGE NC AUTOCONVERSION DROPLETS",
        model = :number),
    :NPRC1 => (kind = :number, long_name = "CHANGE NR AUTOCONVERSION DROPLETS",
        model = :number),
    :NPRCI => (kind = :number, long_name = "CHANGE N AUTOCONVERSION CLOUD ICE BY SNOW",
        model = :number),
    :NPRA => (kind = :number, long_name = "CHANGE IN N DUE TO DROPLET ACC BY RAIN",
        model = :number),
    :NPRAI => (kind = :number, long_name = "CHANGE N ACCRETION CLOUD ICE", model = :number),
    :NPRACS => (kind = :number, long_name = "CHANGE N RAIN-SNOW COLLECTION",
        model = :number),
    :NPRACG => (kind = :number, long_name = "CHANGE N COLLECTION RAIN BY GRAUPEL",
        model = :number),
    :NPSACWS => (kind = :number, long_name = "CHANGE N DROPLET ACCRETION BY SNOW",
        model = :number),
    :NPSACWI => (kind = :number, long_name = "CHANGE N DROPLET ACCRETION BY CLOUD ICE",
        model = :number),
    :NPSACWG => (kind = :number, long_name = "CHANGE N COLLECTION DROPLETS BY GRAUPEL",
        model = :number),
    :NIACR => (kind = :number, long_name = "CHANGE N, ICE-RAIN COLLECTION",
        model = :number),
    :NIACRS => (kind = :number, long_name = "CHANGE N, ICE RAIN COLLISION, ADDED TO SNOW",
        model = :number),
    :NSAGG => (kind = :number, long_name = "SELF-COLLECTION OF SNOW", model = :number),
    :NRAGG => (kind = :number,
        long_name = "CHANGE IN NR DUE TO SELF-COLLECTION AND BREAKUP OF RAIN",
        model = :number),
    :NSCNG => (kind = :number,
        long_name = "CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW",
        model = :number),
    :NGRACS => (kind = :number,
        long_name = "CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW",
        model = :number),
    :NNUCCC => (kind = :number, long_name = "CHANGE N DUE TO CONTACT FREEZ DROPLETS",
        model = :number),
    :NNUCCD => (kind = :number,
        long_name = "CHANGE N FREEZING AEROSOL (PRIM ICE NUCLEATION)", model = :number),
    :NNUCCI => (kind = :number, long_name = "CHANGE N DUE TO IMMERSION FREEZ DROPLETS",
        model = :number),
    :NNUCCR => (kind = :number, long_name = "CHANGE N DUE TO CONTACT FREEZ RAIN",
        model = :number),
    :NHOMOC => (kind = :number,
        long_name = "CHANGE N DUE TO HOMOGENEOUS FREEZING OF CLOUD WATER TO FORM CLOUD ICE",
        model = :number),
    :NHOMOR => (kind = :number,
        long_name = "CHANGE N DUE TO HOMOGENEOUS FREEZING OF RAIN TO FORM GRAUPEL",
        model = :number),
    :NMULTS => (kind = :number, long_name = "ICE MULT DUE TO RIMING DROPLETS BY SNOW",
        model = :number),
    :NMULTG => (kind = :number, long_name = "ICE MULT DUE TO ACC DROPLETS BY GRAUPEL",
        model = :number),
    :NMULTR => (kind = :number, long_name = "ICE MULT DUE TO RIMING RAIN BY SNOW",
        model = :number),
    :NMULTRG => (kind = :number, long_name = "ICE MULT DUE TO ACC RAIN BY GRAUPEL",
        model = :number),
    :NMELTI => (kind = :number,
        long_name = "CHANGE N DUE TO MELTING OF CLOUD ICE TO FORM RAIN", model = :number),
    :NSMLTS => (kind = :number, long_name = "CHANGE N MELTING SNOW", model = :number),
    :NSMLTR => (kind = :number, long_name = "CHANGE N MELTING SNOW TO RAIN",
        model = :number),
    :NGMLTG => (kind = :number, long_name = "CHANGE N MELTING GRAUPEL", model = :number),
    :NGMLTR => (kind = :number, long_name = "CHANGE N MELTING GRAUPEL TO RAIN",
        model = :number),
    :NSUBC => (kind = :number, long_name = "LOSS OF NC DURING EVAP", model = :number),
    :NSUBR => (kind = :number, long_name = "LOSS OF NR DURING EVAP", model = :number),
    :NSUBI => (kind = :number, long_name = "LOSS OF NI DURING SUB.", model = :number),
    :NSUBS => (kind = :number, long_name = "LOSS OF NS DURING SUB.", model = :number),
    :NSUBG => (kind = :number, long_name = "CHANGE N SUB/DEP OF GRAUPEL", model = :number),
)

"""Mass process rates: `?/kg/day` is a mixing-ratio rate, so the water factor applies."""
const ATLAS_MASS_PROCESS_RATES =
    Tuple(sort!([k for (k, v) in ATLAS_PROCESS_RATES if v.kind === :mass]))

"""Number process rates: `?/kg/day` is a per-mass number rate, carrying no water factor."""
const ATLAS_NUMBER_PROCESS_RATES =
    Tuple(sort!([k for (k, v) in ATLAS_PROCESS_RATES if v.kind === :number]))

_model_terms(::Symbol) = nothing
_model_terms(m::Tuple{Symbol, Int}) = (m,)
_model_terms(m::Tuple{Vararg{Tuple{Symbol, Int}}}) = m

"""
    atlas_model_terms(rate)

The `(term, sign)` pairs [`ATLAS_PROCESS_RATES`](@ref) maps `rate` onto, always as a tuple,
or `nothing` when the entry records a [`NO_MODEL_COUNTERPART`](@ref) reason instead.

A rate carrying one tendency may be written as a bare pair; every reader goes through here,
so only this method sees the two spellings.
"""
function atlas_model_terms(rate::Symbol)
    haskey(ATLAS_PROCESS_RATES, rate) ||
        error("`$rate` is not an Atlas process rate.")
    return _model_terms(ATLAS_PROCESS_RATES[rate].model)
end

"""
    process_rate_tendencies(rates = ATLAS_PROCESS_RATES)

[`ATLAS_PROCESS_RATES`](@ref) grouped the other way round: model tendency => the Atlas rates
summing to it, with signs. A rate mapping onto several tendencies appears under each.
"""
function process_rate_tendencies(rates = ATLAS_PROCESS_RATES)
    out = Dict{Symbol, Vector{Tuple{Symbol, Int}}}()
    for rate in sort!(collect(keys(rates)))
        terms = _model_terms(rates[rate].model)
        isnothing(terms) && continue
        for (term, sign) in terms
            push!(get!(out, term, Tuple{Symbol, Int}[]), (rate, sign))
        end
    end
    return out
end

"""
Transport tendency => the Atlas variables summing to it, with signs. Each is already a
`∂q/∂t` and keeps the archive's sign.

`sed_q_lcl` and `adv_q_lcl` read the derived `QCSED` and `QCADV`: `QTO` is vapour plus cloud
liquid and vapour does not sediment, so `QTOSED` is cloud liquid's sedimentation, and Atlas
reports cloud liquid's vertical advection as a flux rather than a tendency.
"""
const LES_TRANSPORT_TENDENCIES = Dict{Symbol, Vector{Tuple{Symbol, Int}}}(
    :sed_q_lcl => [(:QCSED, +1)],
    :sed_q_icl => [(:QISED, +1)],
    :sed_q_rai => [(:QRSED, +1)],
    :sed_q_sno => [(:QSSED, +1)],
    :adv_q_lcl => [(:QCADV, +1)],
    :adv_q_icl => [(:QIADV, +1)],
    :adv_q_rai => [(:QRADV, +1)],
    :adv_q_sno => [(:QSADV, +1)],
    :dif_q_icl => [(:QIDIFF, +1)],
    :dif_q_rai => [(:QRDIFF, +1)],
    :dif_q_sno => [(:QSDIFF, +1)],
)

"""
Model tendency => why the archive cannot supply it. Together with [`LES_TENDENCIES`](@ref)
this partitions the model budget, so a term nobody mapped cannot pass for a term nobody has.

"""
const LES_TENDENCIES_UNAVAILABLE = Dict{Symbol, String}(
    :dif_q_lcl => "the archive carries no cloud-liquid diffusion",
    :S_accr_melt_lcl_sno => "no droplet-collision melt reported apart from PSMLT",
    :S_accr_melt_rai_sno => "no rain-collision melt reported apart from PSMLT",
)

"""Model tendency => the Atlas variables summing to it, transport and process alike."""
const LES_TENDENCIES = merge(LES_TRANSPORT_TENDENCIES, process_rate_tendencies())



"""Every Atlas variable [`LES_TENDENCIES`](@ref) names."""
les_tendency_variables(tendencies = LES_TENDENCIES) =
    unique(v for terms in values(tendencies) for (v, _) in terms)

"""
    les_tendency_profiles(case; params, variables)

`(; z, time, data)` of the named Atlas variables, as specific-content tendencies
[kg kg^-1 s^-1] with `time` in elapsed seconds.
"""
les_tendency_profiles(
    c::SOCRATESCase;
    params,
    variables = les_tendency_variables(),
) = read_atlas(c, variables; params)

"""
    case_les_tendencies(case, z; params, tendencies)

`(; time, tendencies)` for `case`: `time` in elapsed seconds, and each model tendency name of
`tendencies` as a `(length(z), length(time))` matrix [kg kg^-1 s^-1] on levels `z`.
"""
function case_les_tendencies(
    c::SOCRATESCase,
    z::AbstractVector;
    params,
    tendencies = LES_TENDENCIES,
)
    raw = les_tendency_profiles(c; params, variables = les_tendency_variables(tendencies))
    FT = eltype(raw.z)
    zt = collect(FT, z)

    out = Dict{Symbol, Matrix{FT}}()
    for (term, sources) in tendencies
        total = zeros(FT, length(raw.z), length(raw.time))
        for (name, sign) in sources
            total .+= sign .* raw.data[name]
        end
        resampled = Matrix{FT}(undef, length(zt), length(raw.time))
        for j in axes(total, 2)
            resampled[:, j] .= CA.interp_vertical_prof(zt, raw.z, total[:, j])
        end
        out[term] = resampled
    end
    return (; raw.time, tendencies = out)
end
