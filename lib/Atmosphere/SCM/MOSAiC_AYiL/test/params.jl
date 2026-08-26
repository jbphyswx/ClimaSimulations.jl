using Test: Test
using MOSAiC_AYiL: MOSAiC_AYiL as MA
using ClimaAtmos: ClimaAtmos as CA
using NCDatasets: NCDatasets as NC

# ClimaParams ignores an override whose key no component reads, silently. So the
# check that matters is on the *resolved* values, not on the file's contents: a
# key that moves upstream has to fail here rather than quietly restore a default.

Test.@testset "the resolved thermodynamics is DALES's" begin
    c = MA.case("20200503")
    params = MA.mosaic_params(Float64, c)
    th = CA.Parameters.thermodynamics_params(params)
    TP = CA.TD.Parameters

    Test.@test TP.R_d(th) == 287.04
    Test.@test TP.R_v(th) == 461.5
    Test.@test TP.cp_d(th) == 1004.0
    Test.@test TP.LH_v0(th) == 2.53e6
    Test.@test TP.LH_s0(th) == 2.834e6
    Test.@test TP.T_freeze(th) == 273.16
    Test.@test TP.p_ref_theta(th) == 1.0e5
    Test.@test TP.grav(th) == 9.81

    # derived, so it follows from the two above rather than being set
    Test.@test TP.kappa_d(th) ≈ 287.04 / 1004.0
    Test.@test TP.cv_d(th) ≈ 1004.0 - 287.04

    # and they agree with the constants the readers reconstruct with
    Test.@test TP.R_d(th) == MA.DALES_CONSTANTS.R_d
    Test.@test TP.cp_d(th) == MA.DALES_CONSTANTS.cp_d
    Test.@test TP.LH_v0(th) == MA.DALES_CONSTANTS.L_v
    Test.@test TP.grav(th) == MA.DALES_CONSTANTS.grav
end

Test.@testset "CCN" begin
    c = MA.case("20200503")
    params = MA.mosaic_params(Float64, c)
    tabulated = MA.MOSAiC_AYiL_N_CCNs_default[c.date]

    # the table is the file's own `n_ccn`, which is uniform in height
    fd = MA.read_scm_in(c.date)
    Test.@test extrema(fd.n_ccn) == (tabulated, tabulated)

    # and it is consumed, not silently dropped: ClimaParams tags an override with
    # the component that read it, once something has asked for it
    toml = MA.mosaic_toml_dict(Float64, c)
    MA.mosaic_params(toml)
    entry = toml.data["prescribed_cloud_droplet_number_concentration"]
    Test.@test entry["value"] ≈ tabulated
    Test.@test haskey(entry, "used_in")

    # every archived day has an entry
    Test.@test length(MA.MOSAiC_AYiL_N_CCNs_default) == 190
    Test.@test all(
        haskey(MA.MOSAiC_AYiL_N_CCNs_default, d) for d in MA.available_dates()
    )

    # the table is a cache of the archive, so every entry is checked against it
    disagreeing = filter(MA.available_dates()) do d
        NC.NCDataset(MA.scm_in_path(d), "r") do ds
            lo, hi = extrema(Array(ds["n_ccn"])[:, 1])
            !(lo == hi == MA.MOSAiC_AYiL_N_CCNs_default[d])
        end
    end
    Test.@test isempty(disagreeing)
end

Test.@testset "caller overrides win" begin
    c = MA.case("20200503")
    toml = MA.mosaic_toml_dict(
        Float64, c;
        params = Dict(
            "latent_heat_vaporization_at_reference" =>
                Dict{String, Any}("value" => 2.4e6, "type" => "float"),
        ),
    )
    th = CA.Parameters.thermodynamics_params(MA.mosaic_params(toml))
    Test.@test CA.TD.Parameters.LH_v0(th) == 2.4e6
end
