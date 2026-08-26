"""
    SOCRATES

Single-column simulation and Atlas LES scoring for the SOCRATES campaign.
Part of ClimaSimulations.jl.
"""
module SOCRATES

using ClimaAnalysis: ClimaAnalysis
using ClimaAtmos: ClimaAtmos as CA
using ClimaComms: ClimaComms
using ClimaParams: ClimaParams as CP
using Dates: Dates
using NCDatasets: NCDatasets as NC
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF
using Statistics: Statistics
using TOML: TOML

# Reached through ClimaAtmos rather than declared as dependencies, so that this
# package never gates a ClimaAtmos upgrade on its own compat bounds. `CA.TD` is a
# deliberate alias in ClimaAtmos; `CA.CC` and `CA.Insolation` are imports inside
# `callbacks.jl` and `radiation.jl` that land in its namespace, so they are the
# lines to fix if either ever stops being reachable.
const CC = CA.CC
const TD = CA.TD
const Insolation = CA.Insolation

include("casedata.jl")
include("parameters.jl")
include("atlas_archive.jl")
include("atlas_registry.jl")
include("cases.jl")
include("atlas_reader.jl")
include("forcing.jl")
include("io.jl")
include("radiation.jl")
include("model.jl")
include("budgets.jl")
include("transport_diagnostics.jl")
include("scoring.jl")
include("les_rates.jl")

end # module SOCRATES
