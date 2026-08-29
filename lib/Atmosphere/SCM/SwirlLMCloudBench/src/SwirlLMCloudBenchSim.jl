"""
    CloudBench

Single-column simulation of the Swirl-LM CloudBench sites: a GCM-driven LES ensemble of
500 sites × 4 months × 5 experiments, forced as in Shen et al. (2022) and documented by
Chammas et al. (2026).

The case identity, the sounding, the forcing, the initial condition, the insolation and
the surface condition all come from `SwirlLMCloudBench`, which owns the dataset and its
`ClimaAtmos` extension. What lives here is the single-column assembly around them: the
vertical grid, the `AtmosModel`, the `AtmosSimulation`, and reading a case's published
output.

Part of ClimaSimulations.jl.
"""
module SwirlLMCloudBenchSim

using ClimaAnalysis: ClimaAnalysis
using ClimaAtmos: ClimaAtmos as CA
using ClimaComms: ClimaComms
using ClimaParams: ClimaParams as CP
using Dates: Dates
using SwirlLMCloudBench: SwirlLMCloudBench as SW
using SwirlLMCloudBench: Simulation as S

const CC = CA.CC
const TD = CA.TD

include("grid.jl")
include("io.jl")
include("model.jl")
include("parallel.jl")
include("reference.jl")

end # module CloudBench