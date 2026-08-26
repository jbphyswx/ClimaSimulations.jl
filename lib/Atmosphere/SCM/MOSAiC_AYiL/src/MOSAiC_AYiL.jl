"""
    MOSAiC_AYiL

Single-column simulation of the MOSAiC "A Year in LES" suite: 190 Arctic days
over the drifting MOSAiC ice camp, forced by ERA5 Lagrangian-trajectory testbed
files and compared against the DALES large-eddy simulations of Schnierstein et
al. (2024).

Part of ClimaSimulations.jl.
"""
module MOSAiC_AYiL

using Artifacts: Artifacts
using ClimaAnalysis: ClimaAnalysis
using ClimaAtmos: ClimaAtmos as CA
using ClimaComms: ClimaComms
using ClimaParams: ClimaParams as CP
using Dates: Dates
using LazyArtifacts: LazyArtifacts
using NCDatasets: NCDatasets as NC
using Statistics: Statistics
using TOML: TOML

# Reached through ClimaAtmos rather than declared as dependencies, so that this
# package never gates a ClimaAtmos upgrade on its own compat bounds.
const CC = CA.CC
const TD = CA.TD
const Insolation = CA.Insolation
const Intp = CA.Intp
const lazy = CA.lazy

include("variable_translations.jl")
include("dales.jl")
include("data.jl")
include("ayil_info.jl")
include("ayil_info_routines.jl") # provenance, not ordinarily needed at runtime
include("cases.jl")
include("params.jl")
include("io.jl")
include("forcing.jl")
include("setup.jl")
include("diagnostics.jl")
include("model.jl")

# the names ClimaAtmos lacks have to be in its registry before a diagnostics config
# can ask for them
__init__() = register_condensate_totals!()

end # module MOSAiC_AYiL
