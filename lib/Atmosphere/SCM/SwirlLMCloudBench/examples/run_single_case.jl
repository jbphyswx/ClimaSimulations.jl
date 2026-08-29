"""
    run_single_case.jl

A standalone forward CloudBench SCM run.

`SwirlLMCloudBenchSim.reference` reads the same case's published output if you want to look
at the two together.
"""

using SwirlLMCloudBench: Simulation as S
using SwirlLMCloudBenchSim: SwirlLMCloudBenchSim as CB

# Select case: site, month, experiment
case = S.CloudBenchInstance(0, 1, :amip)
output_dir = joinpath(@__DIR__, "output", CB.cloudbench_job_id(case))

@info "Running CloudBench SCM simulation" case = CB.cloudbench_job_id(case) output_dir

CB.run_case(
    case;
    FT = Float64,
    output_dir,
    grid = CB.cloudbench_grid(Float64; dz_min = 50.0),
    verbose = true,
)