using LinearAlgebra: LinearAlgebra
using Test: Test

# The pipeline built end to end without running a forward model: the observation, the
# EKP object, and the minibatch round trip. A forward run is `test/integration.jl` in
# the package and is not part of this suite.

include(joinpath(@__DIR__, "..", "run_calibration.jl"))

const CASES = [MOSAiC_AYiL.case(d) for d in ("20200503", "20200210")]

Test.@testset "calibrated days" begin
    dates = MC.default_calibration_dates()
    Test.@test length(dates) == 71
    Test.@test issorted(dates)
    Test.@test issubset(dates, MOSAiC_AYiL.available_dates())
    # the five days no canonical rung excludes the anomalous ice on are dropped, and
    # dropping them does not touch the package's own table
    Test.@test !any(in(dates), MC.default_not_calibrated_additions)
    Test.@test length(MOSAiC_AYiL.BEST_SIMULATION_TOP_F) == 76
end

Test.@testset "observation shape" begin
    grid = MC.default_mosaic_ayil_grid(; dz_min = 400)
    interface = MC.MOSAiCInterface(;
        output_dir = mktempdir(), cases = CASES,
        vars = ["ql_all", "qi_all", "lwp", "iwp"], grid,
    )
    obs = MC.observations(interface)
    Test.@test length(obs) == length(CASES)
    for (c, o) in zip(interface.cases, obs)
        n = length(MC.case_levels(interface, c))
        y = EKP.get_obs(o)
        # two profiles over the levels plus two scalar paths
        Test.@test length(y) == 2n + 2
        Test.@test all(isfinite, y)
        Γ = EKP.get_obs_noise_cov(o)
        Test.@test size(Γ) == (2n + 2, 2n + 2)
        Test.@test all(>(0), LinearAlgebra.diag(Γ))
    end
    # a day whose calibration top is below the column is scored on fewer levels than
    # one that runs the whole thing, on the same grid
    Test.@test length(MC.case_levels(interface, CASES[1])) >
               length(MC.case_levels(interface, CASES[2]))
end

Test.@testset "the EKP object and the minibatch round trip" begin
    grid = MC.default_mosaic_ayil_grid(; dz_min = 400)
    interface = MC.MOSAiCInterface(;
        output_dir = mktempdir(), cases = CASES,
        vars = ["ql_all", "qi_all", "lwp", "iwp"], grid,
    )
    series = MC.observation_series(interface; minibatch_size = length(CASES))
    prior = default_ayil_prior()
    Test.@test length(EKP.get_name(prior)) == 4

    rng = Random.MersenneTwister(1234)
    ekp = EKP.EnsembleKalmanProcess(
        EKP.construct_initial_ensemble(rng, prior, 10),
        series,
        EKP.Inversion();
        rng,
        scheduler = EKP.DataMisfitController(; terminate_at = 1000.0),
        localization_method = EKP.SECNice(),
        accelerator = EKP.NesterovAccelerator(),
        failure_handler_method = EKP.SampleSuccGauss(),
        verbose = false,
    )
    Test.@test EKP.get_N_ens(ekp) == 10
    # the whole observation is what EKP sees when nothing is held back
    Test.@test length(EKP.get_obs(ekp)) ==
               sum(length(EKP.get_obs(o)) for o in MC.observations(interface))
    # and the minibatch indices come back as the cases they stand for
    Test.@test [c.date for c in MC.minibatch_cases(interface, ekp)] ==
               [c.date for c in CASES]
end

Test.@testset "a member writes only what is scored" begin
    # the scored profiles, the profiles the paths integrate, and rhoa for that
    Test.@test sort(MC.required_diagnostics(["ql_all", "qi_all", "lwp", "iwp"])) ==
               ["qi_all", "ql_all", "rhoa"]
    Test.@test sort(MC.required_diagnostics(["clw", "clwp"])) == ["clw", "rhoa"]
    # with no paths there is nothing to integrate, so no density either
    Test.@test MC.required_diagnostics(["ta"]) == ["ta"]
    # and it is a strict subset of the package's defaults, which is the point
    Test.@test length(MC.required_diagnostics(collect(MC.default_calibration_vars))) <
               length(MOSAiC_AYiL.DEFAULT_DIAGNOSTIC_VARS)
    Test.@test all(
        in(MOSAiC_AYiL.DEFAULT_DIAGNOSTIC_VARS),
        MC.required_diagnostics(collect(MC.default_calibration_vars)),
    )
end

Test.@testset "a scalar weighs a profile" begin
    n = 40
    q = fill(2.0e-4, n)
    profile = MC.observation_variance("clw", q; dof = n, profile_dof = n)
    scalar = MC.observation_variance("clw", q; dof = 1, profile_dof = n)
    Test.@test scalar ≈ profile / n
    # so both contribute the same total to a sum of squares
    Δ = 1.0e-5
    Test.@test n * Δ^2 / profile ≈ Δ^2 / scalar
end

Test.@testset "a transformed band is the physical band propagated" begin
    t = MC.default_variable_transformations["qi_all"]
    fraction = MC.DEFAULT_UNCERTAINTY_FRACTION
    unweighted = Dict("qi_all" => 1.0)
    # far above the offset, a fractional band on `q` is a fixed width in decades
    σ_high = sqrt(MC.observation_variance(
        "qi_all", fill(1.0e-4, 5); scaling = unweighted, transform = t,
    ))
    Test.@test σ_high ≈ fraction / log(10) rtol = 1.0e-3
    # far below it the stored value is the offset, not the quantity, and the band
    # widens to say so
    σ_low = sqrt(MC.observation_variance(
        "qi_all", fill(1.0e-30, 5); scaling = unweighted, transform = t,
    ))
    Test.@test σ_low > 5 * σ_high
end
