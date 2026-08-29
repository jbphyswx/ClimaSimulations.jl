using Test: Test

Test.@testset "MOSAiC_AYiL" begin
    include("data.jl")
    include("params.jl")
    include("translations.jl")
    include("io.jl")
    include("forcing.jl")
    include("setup.jl")
    include("scheduler.jl")
end
