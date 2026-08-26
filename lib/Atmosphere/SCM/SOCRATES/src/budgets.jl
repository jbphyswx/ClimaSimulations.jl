"""
    budgets.jl

Which microphysics process rates make up each prognostic tendency, and the diagnostic names that
carry them. A calibration writes only the scored variables, so reading a finished run for *which
process did what* means asking for these instead.
"""

"""
The 1-moment process rates ClimaAtmos registers, available for `NonEquilibriumMicrophysics1M`.
"""
const MP1M_SOURCE_TERMS = (
    "S_phase_change_vap_lcl",
    "S_phase_change_vap_icl",
    "S_acnv_lcl_rai",
    "S_acnv_icl_sno",
    "S_accr_lcl_rai",
    "S_accr_lcl_sno_cold",
    "S_accr_lcl_sno_warm",
    "S_accr_melt_lcl_sno",
    "S_accr_icl_rai",
    "S_accr_freeze_icl_rai",
    "S_accr_icl_sno",
    "S_accr_rai_sno_cold",
    "S_accr_rai_sno_warm",
    "S_accr_melt_rai_sno",
    "S_phase_change_vap_rai",
    "S_phase_change_vap_sno",
    "S_melt_icl_lcl",
    "S_melt_sno_rai",
)

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
State and updraft fields for reading a run's evolution beyond what is scored. `pfull` is the quantity
whose sign a thermodynamic blow-up destroys first.
"""
const DEBUG_DIAGNOSTIC_VARS = (
    SCORED_VARS...,
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
Every process rate at every location, plus the fields needed to read them: 54 rates, far too much
output for a calibration member, so this is what a postprocessing rerun asks for.
"""
const TENDENCY_DIAGNOSTIC_VARS = (
    DEBUG_DIAGNOSTIC_VARS...,
    ("$(loc)_$(term)" for loc in MP1M_LOCATIONS for term in MP1M_SOURCE_TERMS)...,
    ("$(prefix)_$(var)" for prefix in TRANSPORT_PREFIXES for var in MP1M_BUDGET_VARS)...,
)
