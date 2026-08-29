using Test: Test

Test.@testset "BOMEXCalibration" begin
    include("scoring.jl")
    include("postprocess.jl")
end
