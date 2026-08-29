using Distributed: Distributed
using MOSAiC_AYiL: MOSAiC_AYiL
using Test

@testset "the choice rule prefers a warm key and never idles" begin
    keys = [:a, :b, :a, :b]

    pending = [1, 2, 3, 4]
    @test MOSAiC_AYiL._take_warm_first!(pending, keys, Set([:b])) == 2
    @test pending == [1, 3, 4]

    pending = [1, 2, 3, 4]
    @test MOSAiC_AYiL._take_warm_first!(pending, keys, Set{Symbol}()) == 1

    pending = [1, 2, 3, 4]
    @test MOSAiC_AYiL._take_warm_first!(pending, keys, Set([:c])) == 1

    pending = [2, 4]
    @test MOSAiC_AYiL._take_warm_first!(pending, keys, Set([:a])) == 2
end

@testset "startup waves" begin
    @test MOSAiC_AYiL._startup_delay.(1:9, 4, 10.0) ==
          [0.0, 10.0, 20.0, 30.0, 0.0, 10.0, 20.0, 30.0, 0.0]
    @test MOSAiC_AYiL._startup_delay(3, 0, 10.0) == 0.0
    @test MOSAiC_AYiL._startup_delay(3, 4, 0) == 0.0
end

@testset "a task list runs once each, in order" begin
    tasks = collect(10:10:80)
    keys = [1, 1, 2, 2, 1, 1, 2, 2]
    serial = MOSAiC_AYiL.run_tasks(x -> 2x, tasks, keys, MOSAiC_AYiL.SerialExecutor())
    @test serial == 2 .* tasks

    @test_throws ErrorException MOSAiC_AYiL.run_tasks(
        identity, tasks, keys[1:2], MOSAiC_AYiL.SerialExecutor(),
    )
end

@testset "a failing task does not stop the others" begin
    results = MOSAiC_AYiL.run_tasks(
        x -> x == 2 ? error("boom") : x,
        [1, 2, 3],
        [1, 1, 1],
        MOSAiC_AYiL.SerialExecutor(),
    )
    @test results == [1, nothing, 3]
end

@testset "a worker pool runs every task exactly once, in order" begin
    added = Distributed.addprocs(2)
    try
        pool = Distributed.WorkerPool(added)
        tasks = collect(1:8)
        keys = [1, 1, 1, 1, 2, 2, 2, 2]
        results = MOSAiC_AYiL.run_tasks(
            x -> (x, Distributed.myid()),
            tasks,
            keys,
            MOSAiC_AYiL.WorkerPoolExecutor(pool),
        )
        @test first.(results) == tasks
        @test all(id -> id in added, last.(results))

        staggered = MOSAiC_AYiL.WorkerPoolExecutor(
            pool; startup_waves = 2, startup_wave_pause = 0.5,
        )
        elapsed = @elapsed staggered_results =
            MOSAiC_AYiL.run_tasks(identity, tasks, keys, staggered)
        @test staggered_results == tasks
        @test elapsed > 0.5
    finally
        Distributed.rmprocs(added; waitfor = 30)
    end
end
