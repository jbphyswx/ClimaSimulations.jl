"""
    BOMEXCalibration

EKI calibration of the BOMEX single-column setup. Part of ClimaSimulations.jl.

`scoring.jl` defines what is compared, over what window and on which levels; `calibration.jl`
turns that into the observation EKP takes and runs the ensemble; `postprocess.jl` reruns selected
members with the full diagnostics.
"""
module BOMEXCalibration

using BOMEX: BOMEX
using ClimaAnalysis: ClimaAnalysis
using ClimaAtmos: ClimaAtmos as CA
using ClimaCalibrate: ClimaCalibrate
using EnsembleKalmanProcesses: EnsembleKalmanProcesses as EKP
using JLD2: JLD2
using LinearAlgebra: LinearAlgebra
using NCDatasets: NCDatasets as NC
using Random: Random
using Statistics: Statistics

include("scoring.jl")
include("calibration.jl")
include("postprocess.jl")

end # module BOMEXCalibration