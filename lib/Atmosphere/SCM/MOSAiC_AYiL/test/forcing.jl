using Test: Test
using MOSAiC_AYiL: MOSAiC_AYiL as MA

# The relaxation shape and the inversion search are where a sign or factor error
# would be invisible in a run that merely completes, so they are pinned against
# values worked out by hand.

Test.@testset "nudging ramp" begin
    z_inv, depth = 700.0, 300.0
    z_mid = z_inv + depth
    Test.@test MA.nudge_ramp(0.0, z_inv, z_mid) == 0
    Test.@test MA.nudge_ramp(z_inv, z_inv, z_mid) == 0          # zero *at* the inversion
    Test.@test MA.nudge_ramp(z_inv + 150, z_inv, z_mid) ≈ 0.5
    Test.@test MA.nudge_ramp(z_mid, z_inv, z_mid) ≈ 1
    Test.@test MA.nudge_ramp(z_mid + 5000, z_inv, z_mid) ≈ 1
    # a zero-depth ramp is a step at the inversion
    Test.@test MA.nudge_ramp(z_inv, z_inv, z_inv) == 0
    Test.@test MA.nudge_ramp(z_inv + 1, z_inv, z_inv) == 1
end

Test.@testset "inversion height" begin
    # a θ_l profile with a sharp jump: the centred difference
    # (θ[k+1] - θ[k-1]) / (z[k+1] - z[k-1]) peaks one level *below* the jump,
    # which is a property of the DALES algorithm, not an off-by-one here.
    z = collect(50.0:100.0:1950.0)
    θ = 280.0 .+ 0.001 .* z
    k_jump = 8
    θ[k_jump] += 10.0
    Test.@test MA.inversion_height(θ, z, 100.0, 2000.0) == z[k_jump - 1]

    # the window is honoured: nothing outside it can be selected
    z_above = MA.inversion_height(θ, z, 1200.0, 2000.0)
    Test.@test 1200.0 < z_above < 2000.0

    # an empty window yields zero rather than an arbitrary level
    Test.@test MA.inversion_height(θ, z, 100.0, 120.0) == 0
    Test.@test MA.inversion_height(θ, z, 1900.0, 100.0) == 0

    # the first and last centres are never candidates, so the ramp bottom is
    # always an interior height
    Test.@test MA.inversion_height(θ, z, 0.0, 1.0e5) > first(z)
    Test.@test MA.inversion_height(θ, z, 0.0, 1.0e5) < last(z)

    # a monotonic profile has its steepest centred gradient somewhere inside
    Test.@test MA.inversion_height(280.0 .+ 0.001 .* z, z, 100.0, 2000.0) > 0
end

Test.@testset "forcing from the archive" begin
    c = MA.case("20200503")
    fd = MA.read_scm_in(c.date)
    f = MA.MOSAiCForcing(Float64, c; forcing = fd)

    # τ = 10800 s on every archived day
    Test.@test f.inv_τ ≈ 1 / 10800
    Test.@test f.ramp_depth == 300
    Test.@test (f.z_inv_min, f.z_inv_max) == (100.0, 5000.0)
    Test.@test f.q_tot_threshold == 1.0e-6

    # the profiles are on the ascending scm_in axis
    Test.@test issorted(f.z)
    Test.@test length(f.z) == length(fd.z)
    Test.@test f.ta ≈ fd.ta
    Test.@test f.hus ≈ fd.hus

    # the thermodynamic horizontal advection is identically zero in the archive,
    # while the momentum advection is not
    Test.@test all(iszero, f.dTdt_hadv)
    Test.@test all(iszero, f.dqtdt_hadv)
    Test.@test !all(iszero, f.dudt_hadv)

    # a fixed nudging onset is a different closure and must not be accepted
    Test.@test_throws ErrorException MA.MOSAiCForcing(
        Float64, c;
        forcing = fd,
        nudging = (; timescale = 10800.0, ramp_depth = 300.0, z_min = 500.0),
    )
end

Test.@testset "omega to w conversion" begin
    c = MA.case("20200503")
    fd = MA.read_scm_in(c.date)
    # `wa` is positive upward, so it carries the opposite sign to `omega`
    ω = MA.NC.NCDataset(MA.scm_in_path(c.date)) do ds
        reverse(Array(ds["omega"])[:, 1])
    end
    nz = findall(!iszero, ω)
    Test.@test all(sign.(fd.wa[nz]) .== .-sign.(ω[nz]))
    # and reproduces DALES's own w = -omega / (rho g) with its constants
    (; grav, R_d) = MA.DALES_CONSTANTS
    Tv = fd.ta .* (1 .+ 0.61f0 .* fd.q)
    Test.@test fd.wa ≈ .-ω .* (R_d .* Tv) ./ (fd.p .* grav)
end
