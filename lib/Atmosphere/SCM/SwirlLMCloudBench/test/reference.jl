# Reads a case's published `data.zarr` over HTTPS, so this is NOT part of `test/runtests.jl`;
# CI runs it as a separate gated job.
#
#   julia --project=test test/reference.jl

using SwirlLMCloudBench: Simulation as S
using SwirlLMCloudBenchSim: SwirlLMCloudBenchSim as CB
using Test: Test

Test.@testset "published output" begin
    inst = S.CloudBenchInstance(0, 1, :amip)
    store = CB.reference_store(inst)

    Test.@testset "the axes line up with a run of this package's own length" begin
        ax = CB.reference_axes(store)
        # `t` is nanoseconds in the file; if the conversion were dropped these would be
        # off by nine orders of magnitude
        Test.@test last(ax.time) == CB.cloudbench_t_end()
        Test.@test ax.time[2] - ax.time[1] == CB.cloudbench_output_interval()
        # the store covers the last simulated day
        Test.@test last(ax.time) - first(ax.time) == 86400
        Test.@test first(ax.z) == 0
        Test.@test last(ax.z) == S.CLOUDBENCH_LES_GRID.lz
    end

    Test.@testset "variables are read by rank" begin
        Test.@test ndims(store.arrays["lwp"]) == 1
        Test.@test ndims(store.arrays["cloud_fraction"]) == 2
        Test.@test ndims(store.arrays["q_c"]) == 4
        r = CB.reference(inst, ("lwp", "cloud_fraction"); store)
        nt = length(r.time)
        Test.@test size(r.data["lwp"]) == (nt,)
        Test.@test size(r.data["cloud_fraction"]) == (length(r.z), nt)
        Test.@test all(isfinite, r.data["lwp"])
        Test.@test all(isfinite, r.data["cloud_fraction"])
        # a profile the column is compared against has to have cloud in it
        Test.@test maximum(r.data["cloud_fraction"]) > 0
        Test.@test_throws KeyError CB.reference(inst, "not_a_variable"; store)
    end
end
