"""
    SwirlLMCloudBenchCalibration

Comparing a CloudBench single-column run against that case's published output, and the EKI
interface that calibrates against it. Part of ClimaSimulations.jl.
"""
module SwirlLMCloudBenchCalibration

using ClimaAnalysis: ClimaAnalysis
using ClimaAtmos: ClimaAtmos as CA
using ClimaCalibrate: ClimaCalibrate
using EnsembleKalmanProcesses: EnsembleKalmanProcesses as EKP
using JLD2: JLD2
using LinearAlgebra: LinearAlgebra
using Random: Random
using Statistics: Statistics
using SwirlLMCloudBench: SwirlLMCloudBench as SW
using SwirlLMCloudBench: Simulation as S
using SwirlLMCloudBenchSim: SwirlLMCloudBenchSim as CB
using TOML: TOML

include("reference.jl")
include("calibration.jl")
include("postprocess.jl")

end # module SwirlLMCloudBenchCalibration
