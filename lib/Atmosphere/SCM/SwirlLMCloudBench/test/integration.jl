# Builds and solves a real column, so this is NOT part of `test/runtests.jl`; CI runs it as a
# separate gated job. It needs the network for the case's sounding and minutes of wall clock.
#
#   julia --project=test test/integration.jl
#
# `t_end` here is a fraction of `cloudbench_t_end()`. A run at the reference length is the
# package's acceptance criterion and belongs in `examples/run_single_case.jl`, never in a test.

using ClimaAtmos: ClimaAtmos as CA
using SwirlLMCloudBench: Simulation as S
using SwirlLMCloudBenchSim: SwirlLMCloudBenchSim as CB
using Test: Test

const CASE = S.CloudBenchInstance(0, 1, :amip)

Test.@testset "a column builds and steps" begin
    mktempdir() do output_dir
        grid = CB.cloudbench_grid(Float64; dz_min = 50.0)

        # every default of `cloudbench_simulation` is evaluated here, which is the only place
        # an undefined default or an unsupported keyword can surface
        sim = CB.cloudbench_simulation(
            Float64, CASE; output_dir, grid, t_end = 1200.0, verbose = false,
        )
        Test.@test sim.job_id == CB.cloudbench_job_id(CASE)

        result = CA.solve_atmos!(sim)
        Test.@test result.ret_code === :success

        # the run wrote the diagnostics a comparison reads
        written = [f for (r, _, fs) in walkdir(output_dir) for f in fs]
        Test.@test !isempty(written)
    end
end
