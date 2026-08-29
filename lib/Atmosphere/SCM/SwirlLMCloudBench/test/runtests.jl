using ClimaAtmos: ClimaAtmos as CA
using Dates: Dates
using SwirlLMCloudBench: Simulation as S
using SwirlLMCloudBenchSim: SwirlLMCloudBenchSim as CB
using Test: Test

# Each assertion here has to be able to fail for a reason in this package. Reading a case's
# published output needs the network, so `test/reference.jl` is a separate gated job.

Test.@testset "SwirlLMCloudBenchSim" begin

    Test.@testset "the run length and cadence are the reference's own" begin
        Test.@test CB.cloudbench_t_end() ==
                   S.CLOUDBENCH_LES_GRID.duration_days * 86400.0
        Test.@test CB.cloudbench_output_interval() == S.CLOUDBENCH_OUTPUT_INTERVAL
    end

    Test.@testset "les_faces spans the recorded domain" begin
        faces = CB.les_faces()
        Test.@test first(faces) == 0
        Test.@test last(faces) == S.CLOUDBENCH_LES_GRID.lz
        # one more face than cells, and the cells are the recorded resolution
        Test.@test length(faces) == S.CLOUDBENCH_LES_GRID.nz + 1
        Test.@test all(≈(S.CLOUDBENCH_LES_GRID.dz), diff(faces))
        # a domain that is not a whole number of cells is refused rather than truncated
        Test.@test_throws ErrorException CB.les_faces(;
            les_grid = (; lz = 6000.0, dz = 7.0),
        )
    end

    Test.@testset "coarsening keeps the ends and honours dz_min" begin
        faces = CB.les_faces()
        Test.@test CB.coarsen_faces_to_dz_min(faces, nothing) == faces
        # already coarse enough leaves the grid untouched
        Test.@test CB.coarsen_faces_to_dz_min(faces, 1.0) == faces
        for dz_min in (50.0, 200.0, 1000.0)
            kept = CB.coarsen_faces_to_dz_min(faces, dz_min)
            Test.@test first(kept) == first(faces)
            Test.@test last(kept) == last(faces)
            Test.@test minimum(diff(kept)) >= dz_min
            Test.@test issorted(kept)
        end
        # a dz_min wider than the column cannot leave a cell
        Test.@test_throws ErrorException CB.coarsen_faces_to_dz_min([0.0, 10.0], 100.0)
        Test.@test_throws ErrorException CB.coarsen_faces_to_dz_min([0.0], 1.0)
    end

    Test.@testset "the grid carries the faces it was given" begin
        faces = collect(0.0:500.0:6000.0)
        grid = CB.cloudbench_grid(Float64; faces)
        z = CB.cloudbench_z(grid)
        # centres of the requested cells, so the model column is the grid asked for
        Test.@test length(z) == length(faces) - 1
        Test.@test z ≈ (faces[1:(end - 1)] .+ faces[2:end]) ./ 2
        Test.@test_throws ErrorException CB.cloudbench_grid(
            Float64; faces = [0.0, 2000.0, 1000.0],
        )
    end

    Test.@testset "the sponge covers the reference's damped layer" begin
        sp = CB.cloudbench_sponge(Float64).rayleigh_sponge
        Test.@test sp.zd ==
                   S.CLOUDBENCH_LES_GRID.lz - S.CLOUDBENCH_SPONGE_DEPTH
        # swirl-lm damps at β/dt with β = 1/a_coeff, so the rate is 1/(a_coeff dt)
        Test.@test sp.α_w == 1 / (20 * S.CLOUDBENCH_LES_GRID.dt)
        Test.@test sp.α_w ≈ 0.25
        # the reference damps vertical velocity alone
        Test.@test sp.α_uₕ == 0
        Test.@test sp.α_tracer == 0
        Test.@test_throws ErrorException CB.cloudbench_sponge(Float64; depth = 7000.0)
        Test.@test_throws ErrorException CB.cloudbench_sponge(Float64; depth = 0.0)
    end

    Test.@testset "diagnostics are checked before a model is built" begin
        # a name ClimaAtmos does not register fails here rather than part-way into a run
        Test.@test_throws ErrorException CB.cloudbench_diagnostics(
            ("not_a_diagnostic",); n_levels = 10,
        )
        Test.@test_throws ErrorException CB.cloudbench_diagnostics((); n_levels = 10)
        for name in (CB.default_diagnostic_vars..., CB.STATE_VARS...)
            Test.@test haskey(CA.Diagnostics.ALL_DIAGNOSTICS, name)
        end
        # the cloud radiative effects a comparison needs come from the clear-sky mode,
        # which is what `cloudbench_model` defaults to
        for name in ("rlutcs", "rsutcs")
            Test.@test haskey(CA.Diagnostics.ALL_DIAGNOSTICS, name)
        end
    end

    Test.@testset "the job id names the case" begin
        inst = S.CloudBenchInstance(17, 7, :amip_p4k)
        id = CB.cloudbench_job_id(inst)
        Test.@test occursin("17", id)
        Test.@test occursin("7", id)
        Test.@test occursin("amip_p4k", id)
        # two cases never share a directory
        Test.@test CB.cloudbench_job_id(S.CloudBenchInstance(17, 7, :amip)) != id
        Test.@test CB.cloudbench_job_id(S.CloudBenchInstance(18, 7, :amip_p4k)) != id
    end

    Test.@testset "the start date carries the case's own month" begin
        for month in S.Catalog.CLOUDBENCH_MONTHS
            d = CB.cloudbench_start_date(S.CloudBenchInstance(0, month, :amip))
            Test.@test Dates.month(d) == month
            Test.@test Dates.day(d) == 1
        end
        # a simulation resolves to the same date as the instance it wraps
        inst = S.CloudBenchInstance(0, 4, :amip)
        Test.@test CB.cloudbench_start_date(S.CloudBenchSimulation(inst)) ==
                   CB.cloudbench_start_date(inst)
    end

    include("scheduler.jl")

end
