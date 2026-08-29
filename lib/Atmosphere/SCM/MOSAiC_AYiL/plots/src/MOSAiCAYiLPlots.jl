"""
    MOSAiCAYiLPlots

Figures for MOSAiC AYiL runs and calibrations. Part of ClimaSimulations.jl.

`run_plots.jl` draws one run's own output and needs nothing but the run directory.
`calibration_plots.jl` draws an ensemble and its EKI history.
"""
module MOSAiCAYiLPlots

using CairoMakie: CairoMakie
using ClimaAnalysis: ClimaAnalysis
using ClimaCalibrate: ClimaCalibrate
using EnsembleKalmanProcesses: EnsembleKalmanProcesses as EKP
using MOSAiC_AYiL: MOSAiC_AYiL
using MOSAiCAYiLCalibration: MOSAiCAYiLCalibration as MC
using Statistics: Statistics

"""Colour per ensemble group, in the order the groups are given."""
const GROUP_COLORS = (:steelblue, :firebrick, :seagreen, :darkorange)

"""Write `fig` to `path`, creating its directory."""
function save_figure(fig, path::AbstractString)
    mkpath(dirname(abspath(path)))
    CairoMakie.save(path, fig)
    @info "wrote" path
    return path
end

include("run_plots.jl")
include("calibration_plots.jl")

end # module MOSAiCAYiLPlots
