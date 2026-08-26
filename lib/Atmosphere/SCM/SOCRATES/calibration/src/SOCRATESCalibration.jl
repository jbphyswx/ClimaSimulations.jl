"""
    SOCRATESCalibration

Ensemble Kalman Inversion (EKI) calibration suite for the SOCRATES campaign.
Part of ClimaSimulations.jl.
"""
module SOCRATESCalibration

using ClimaAnalysis: ClimaAnalysis
using ClimaCalibrate: ClimaCalibrate
using Dates: Dates
using EnsembleKalmanProcesses: EnsembleKalmanProcesses as EKP
using JLD2: JLD2
using LinearAlgebra: LinearAlgebra
using Logging: Logging
using NaNStatistics: NaNStatistics
using Random: Random
using SOCRATES: SOCRATES
using Statistics: Statistics

include("calibration.jl")
include("postprocess.jl")

end # module SOCRATESCalibration
