"""
    MOSAiCAYiLCalibration

EKI calibration of the MOSAiC AYiL single-column setup against the DALES reference.

`scoring.jl` defines what is compared; `calibration.jl` turns that into the
observation EKP takes and runs the ensemble; `ayil_info.jl` picks the days.
"""
module MOSAiCAYiLCalibration

using ClimaAnalysis: ClimaAnalysis
using ClimaCalibrate: ClimaCalibrate
using EnsembleKalmanProcesses: EnsembleKalmanProcesses as EKP
using JLD2: JLD2
using LinearAlgebra: LinearAlgebra
using MOSAiC_AYiL: MOSAiC_AYiL
using Statistics: Statistics

include("ayil_info.jl")
include("scoring.jl")
include("calibration.jl")

end # module MOSAiCAYiLCalibration
