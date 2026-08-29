# Reads a case's published `data.zarr` over HTTPS, so this is NOT part of
# `calibration/test/runtests.jl`; CI runs it as a separate gated job.
#
#   julia --project=calibration/test calibration/test/reference.jl

using SwirlLMCloudBench: Simulation as S
using SwirlLMCloudBenchSim: SwirlLMCloudBenchSim as CB
using SwirlLMCloudBenchCalibration: SwirlLMCloudBenchCalibration as CBC
using Test: Test

Test.@testset "reference on model units" begin
    inst = S.CloudBenchInstance(0, 1, :amip)
    store = CB.reference_store(inst)
    raw = CB.reference(inst, sort(collect(keys(CBC.REFERENCE_PAIRINGS))); store)

    Test.@testset "$name" for name in sort(collect(keys(CBC.REFERENCE_PAIRINGS)))
        got = CBC.reference_comparable(inst; names = (name,), store)
        Test.@test haskey(got.data, name)
        Test.@test size(got.data[name]) == size(raw.data[name])
        Test.@test all(isfinite, got.data[name])
    end

    Test.@testset "the conversions are the ones the pairing declares" begin
        # `lwp` pairs onto ClimaAtmos's own kg m^-2, so nothing is rescaled
        lwp = CBC.reference_comparable(inst; names = ("lwp",), store).data["lwp"]
        Test.@test lwp == raw.data["lwp"]
        # `cl` is a percentage where the store is a fraction
        cf = CBC.reference_comparable(inst; names = ("cloud_fraction",), store)
        Test.@test cf.data["cloud_fraction"] ≈ 100 .* raw.data["cloud_fraction"]
        Test.@test maximum(cf.data["cloud_fraction"]) <= 100
    end

    Test.@testset "an unpaired store variable is refused" begin
        Test.@test_throws ErrorException CBC.reference_comparable(
            inst; names = ("sfc_flux_rad_lw",), store,
        )
    end
end
