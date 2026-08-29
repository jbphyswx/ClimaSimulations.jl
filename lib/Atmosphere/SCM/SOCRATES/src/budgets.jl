"""
    budgets.jl

Which microphysics process rates make up each prognostic tendency, and the diagnostic names that
carry them.
"""

"""
The 1-moment process rates ClimaAtmos registers, available for `NonEquilibriumMicrophysics1M`.

Read out of `CA.Diagnostics.MP1M_SOURCE_TERMS`.
"""
const MP1M_SOURCE_TERMS =
    Tuple(String(field) for (field, _) in CA.Diagnostics.MP1M_SOURCE_TERMS)

"""Where a rate is evaluated. The updraft and environment sets require `PrognosticEDMFX`."""
const MP1M_LOCATIONS = ("mp1m", "mp1mup", "mp1men")

"""The prognostic variables with a budget, in a stable order for figures."""
const MP1M_BUDGET_VARS = ("q_lcl", "q_icl", "q_rai", "q_sno")

"""
Which process rates enter each prognostic tendency, and with what sign, so that summing a variable's
signed terms reproduces its `dq/dt`.

"""
const MP1M_BUDGETS = Dict{String, Vector{Tuple{String, Int}}}(
    "q_lcl" => [
        ("S_phase_change_vap_lcl", +1),
        ("S_acnv_lcl_rai", -1),
        ("S_accr_lcl_rai", -1),
        ("S_accr_lcl_sno_cold", -1),
        ("S_accr_lcl_sno_warm", -1),
        ("S_melt_icl_lcl", +1),
    ],
    "q_icl" => [
        ("S_phase_change_vap_icl", +1),
        ("S_acnv_icl_sno", -1),
        ("S_accr_icl_rai", -1),
        ("S_accr_icl_sno", -1),
        ("S_melt_icl_lcl", -1),
    ],
    "q_rai" => [
        ("S_acnv_lcl_rai", +1),
        ("S_accr_lcl_rai", +1),
        ("S_accr_lcl_sno_warm", +1),
        ("S_accr_melt_lcl_sno", +1),
        ("S_accr_freeze_icl_rai", -1),
        ("S_accr_rai_sno_cold", -1),
        ("S_accr_rai_sno_warm", +1),
        ("S_accr_melt_rai_sno", +1),
        ("S_phase_change_vap_rai", +1),
        ("S_melt_sno_rai", +1),
    ],
    "q_sno" => [
        ("S_acnv_icl_sno", +1),
        ("S_accr_lcl_sno_cold", +1),
        ("S_accr_melt_lcl_sno", -1),
        ("S_accr_icl_rai", +1),
        ("S_accr_freeze_icl_rai", +1),
        ("S_accr_icl_sno", +1),
        ("S_accr_rai_sno_cold", +1),
        ("S_accr_rai_sno_warm", -1),
        ("S_accr_melt_rai_sno", -1),
        ("S_phase_change_vap_sno", +1),
        ("S_melt_sno_rai", -1),
    ],
)

"""Transport terms, each already a `∂q/∂t` and so entering its budget with `+1`."""
const TRANSPORT_PREFIXES = ("sed", "adv", "dif")

const TRANSPORT_BUDGETS = Dict{String, Vector{Tuple{String, Int}}}(
    var => [("$(prefix)_$(var)", +1) for prefix in TRANSPORT_PREFIXES] for
    var in MP1M_BUDGET_VARS
)

"""
[`default_diagnostic_vars`](@ref) plus the state and updraft fields a run's evolution is read
from. `pfull` is the quantity whose sign a thermodynamic blow-up destroys first.
"""
const DEBUG_DIAGNOSTIC_VARS = (
    default_diagnostic_vars...,
    "pfull",
    "rhoa",
    "ta",
    "hus",
    "wa",
    "tke",
    "arup",
    "waup",
    "entr",
    "detr",
)

"""
[`DEBUG_DIAGNOSTIC_VARS`](@ref) plus every process rate at every location and every transport
term — 54 rates
"""
const TENDENCY_DIAGNOSTIC_VARS = (
    DEBUG_DIAGNOSTIC_VARS...,
    ("$(loc)_$(term)" for loc in MP1M_LOCATIONS for term in MP1M_SOURCE_TERMS)...,
    ("$(prefix)_$(var)" for prefix in TRANSPORT_PREFIXES for var in MP1M_BUDGET_VARS)...,
)
