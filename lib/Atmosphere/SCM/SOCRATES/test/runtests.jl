using Test: Test

Test.@testset "SOCRATES" begin
    include("casedata.jl")
    include("cases.jl")
    include("parameters.jl")
    include("atlas_registry.jl")
    include("forcing.jl")
    include("radiation.jl")
    include("terminal_velocity.jl")
    include("scoring.jl")
end
