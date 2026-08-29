"""
    BOMEX

Single-column simulation of the BOMEX trade-wind cumulus case. Part of ClimaSimulations.jl.

The vertical profiles come from `AtmosphericProfilesLibrary`; the numbers that are not profiles
come from `cases/bomex.toml`, each carrying its source. 
"""
module BOMEX

using ClimaAnalysis: ClimaAnalysis
using ClimaAtmos: ClimaAtmos as CA
using ClimaComms: ClimaComms
using ClimaParams: ClimaParams as CP
using TOML: TOML

const CC = CA.CC
const TD = CA.TD

include("casedata.jl")
include("cases.jl")
include("setup.jl")
include("io.jl")
include("model.jl")
include("parallel.jl")

end # module BOMEX