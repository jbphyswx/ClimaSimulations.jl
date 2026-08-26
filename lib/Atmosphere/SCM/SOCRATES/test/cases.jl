using ClimaAtmos: ClimaAtmos as CA
using SOCRATES: SOCRATES
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF
using Test: Test

const FT = Float64

Test.@testset "cases" begin

    Test.@testset "registry and naming" begin
        cases = SOCRATES.all_cases()
        Test.@test length(cases) == 11
        Test.@test count(c -> c.forcing_type isa SSCF.ObsForcing, cases) == 5
        Test.@test count(c -> c.forcing_type isa SSCF.ERA5Forcing, cases) == 6

        c9 = SOCRATES.case("RF09_Obs")
        Test.@test c9.flight_number == 9
        Test.@test c9.forcing_type isa SSCF.ObsForcing
        Test.@test SOCRATES.case_name(c9) == "RF09_Obs"
        Test.@test SOCRATES.forcing_label(c9) === :Obs
        Test.@test SOCRATES.case(c9) === c9

        c11 = SOCRATES.case("RF11_ERA5")
        Test.@test SOCRATES.case_name(c11) == "RF11_ERA5"

        Test.@test_throws ErrorException SOCRATES.case("RF11_Obs")  # no Obs artifact
        Test.@test_throws ErrorException SOCRATES.case("Invalid_Name")
        Test.@test_throws ErrorException SOCRATES.forcing_type(:Nonsense)
    end

    Test.@testset "clocks" begin
        c9 = SOCRATES.case("RF09_Obs")
        c11 = SOCRATES.case("RF11_ERA5")
        Test.@test SOCRATES.t_end(c9) == 12 * 3600.0
        Test.@test SOCRATES.t_end(c11) == 14 * 3600.0
        Test.@test SOCRATES.n_ccn(c9) == 1.9e8


        # the real campaign date, not a synthetic epoch
        Test.@test SOCRATES.les_start_datetime(c9) ==
              SSCF.get_socrates_initial_time(9)
    end

    Test.@testset "grids" begin
        for c in SOCRATES.all_cases()
            zc = SOCRATES.native_z(c)
            zf = SOCRATES.native_faces(c)
            # `faces_from_centers` errors on non-monotonic input, so ordering is already
            # guaranteed; what it does not guarantee is that the faces it built are the
            # midpoint-consistent ones for these centres
            Test.@test length(zf) == length(zc) + 1
            Test.@test SOCRATES.centers_from_faces(zf) ≈ zc
            Test.@test SOCRATES.z_max(c) == last(zf)
        end

        c = SOCRATES.case("RF09_Obs")
        zf = SOCRATES.native_faces(c)
        coarse = SOCRATES.coarsen_faces_to_dz_min(zf, 50.0)
        Test.@test first(coarse) == first(zf)
        Test.@test last(coarse) == last(zf)
        Test.@test minimum(diff(coarse)) >= 50.0
        Test.@test SOCRATES.coarsen_faces_to_dz_min(zf, nothing) == zf

        # centres that are not cell midpoints must be rejected, not silently used
        Test.@test_throws ErrorException SOCRATES.faces_from_centers([10.0, 5.0])
        Test.@test_throws ErrorException SOCRATES.faces_from_centers(Float64[])

        grid = SOCRATES.socrates_grid(FT, c)
        Test.@test SOCRATES.socrates_z(grid) ≈ SOCRATES.native_z(c)
    end

    Test.@testset "parameters" begin
        c = SOCRATES.case("RF09_Obs")
        dict = SOCRATES.socrates_toml_dict(FT, c)
        Test.@test dict["prescribed_cloud_droplet_number_concentration"] == 1.9e8
        for name in keys(SOCRATES.TERMINAL_VELOCITY_SCALING_PARAMS)
            Test.@test dict[String(name)] == 1.0
        end

        # a caller override wins over both the suite's own set and the case value
        overridden = SOCRATES.socrates_toml_dict(
            FT, c;
            params = Dict(
                "prescribed_cloud_droplet_number_concentration" =>
                    Dict("value" => 5.0e7, "type" => "float"),
            ),
        )
        Test.@test overridden["prescribed_cloud_droplet_number_concentration"] == 5.0e7

        Test.@test_throws ErrorException SOCRATES.socrates_toml_dict(
            FT, c; params = "/nonexistent/path.toml",
        )
        Test.@test_throws ErrorException SOCRATES.socrates_toml_dict(
            FT, c; params = "not_a_toml_file.txt",
        )
    end

    Test.@testset "terminal velocity scaling" begin
        # a scale of 1 must reproduce the diagnostic mode exactly, since that is
        # what the calibration perturbs away from
        mode1 = SOCRATES.ScaledTerminalVelocity{FT}(1.0)
        modeq = SOCRATES.ScaledTerminalVelocity{FT}(0.25)
        mm = CA.NonEquilibriumMicrophysics1M()
        params = SOCRATES.socrates_params(FT, SOCRATES.case("RF09_Obs"))
        cmc = CA.Parameters.microphysics_cloud_params(params)
        cmp = CA.Parameters.microphysics_1m_params(params)
        ρ, q = 1.0, 1.0e-4
        for name in (
            CA.CC.MatrixFields.@name(q_lcl),
            CA.CC.MatrixFields.@name(q_icl),
            CA.CC.MatrixFields.@name(q_rai),
            CA.CC.MatrixFields.@name(q_sno),
        )
            diagnostic = CA.terminal_velocity(
                mm, CA.DiagnosticTerminalVelocity(), name, params, cmc, cmp, ρ, q,
            )
            Test.@test CA.terminal_velocity(mm, mode1, name, params, cmc, cmp, ρ, q) ==
                  diagnostic
            Test.@test CA.terminal_velocity(mm, modeq, name, params, cmc, cmp, ρ, q) ≈
                  0.25 * diagnostic
        end
    end

end