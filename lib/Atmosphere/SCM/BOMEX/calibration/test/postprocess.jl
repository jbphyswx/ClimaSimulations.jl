using BOMEXCalibration: BOMEXCalibration as BC
using EnsembleKalmanProcesses: EnsembleKalmanProcesses as EKP
using LinearAlgebra: LinearAlgebra
using Test: Test

struct StubEKP
    y::Vector{Float64}
    Γ::Matrix{Float64}
end
EKP.get_obs(e::StubEKP) = e.y
EKP.get_obs_noise_cov(e::StubEKP) = e.Γ

Test.@testset "member misfit" begin
    y = [1.0, 2.0, 3.0]
    ekp = StubEKP(y, Matrix{Float64}(LinearAlgebra.I, 3, 3))
    g = [1.0 2.0 NaN; 2.0 2.0 3.0; 3.0 2.0 3.0]
    φ = BC.member_misfits(ekp, g)
    Test.@test φ[1] == 0                      # exact match
    Test.@test φ[2] ≈ 1.0                     # d = [1, 0, -1]
    Test.@test isnan(φ[3])                    # a failed member never scores best

    # Γ is inverted, not ignored: with Γ₁₁ = 4 a deviation of 2 costs ½·4/4
    ekp2 = StubEKP(y, [4.0 0 0; 0 1.0 0; 0 0 1.0])
    Test.@test only(BC.member_misfits(ekp2, reshape([3.0, 2.0, 3.0], 3, 1))) ≈ 0.5
end

Test.@testset "G_ensemble location" begin
    # zero-padded to three digits, matching what the observation map writes
    Test.@test BC.g_ensemble_path("/o", 7) ==
               joinpath("/o", "iteration_007", "G_ensemble.jld2")
    Test.@test BC.g_ensemble_path("/o", 123) ==
               joinpath("/o", "iteration_123", "G_ensemble.jld2")
end
