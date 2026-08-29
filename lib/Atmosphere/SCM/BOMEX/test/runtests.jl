using Test: Test

Test.@testset "BOMEX" begin
    include("cases.jl")
    include("setup.jl")
    include("scheduler.jl")
end
