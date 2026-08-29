using BOMEX: BOMEX
using ClimaAtmos: ClimaAtmos as CA
using Test: Test

Test.@testset "cases" begin
    FT = Float64

    Test.@testset "the shipped grid is the one cases/bomex.toml describes" begin
        faces = BOMEX.bomex_faces()
        Test.@test length(faces) == BOMEX.z_elem() + 1
        Test.@test first(faces) == 0
        Test.@test last(faces) == BOMEX.z_max()
        # uniform, which is what z_stretch = false means in the config this came from
        Test.@test all(≈(first(diff(faces))), diff(faces))
    end

    Test.@testset "a column cannot be built from a degenerate domain" begin
        Test.@test_throws ErrorException BOMEX.uniform_faces(3000.0, 0)
        Test.@test_throws ErrorException BOMEX.uniform_faces(0.0, 10)
        Test.@test_throws ErrorException BOMEX.uniform_faces(-1.0, 10)
    end

    Test.@testset "the grid carries the faces it was given" begin
        c = BOMEX.case(FT)
        faces = collect(0.0:500.0:3000.0)
        grid = BOMEX.bomex_grid(FT, c; faces)
        z = BOMEX.bomex_z(grid)
        Test.@test length(z) == length(faces) - 1
        Test.@test z ≈ (faces[1:(end - 1)] .+ faces[2:end]) ./ 2
        Test.@test_throws ErrorException BOMEX.bomex_grid(
            FT, c; faces = [0.0, 2000.0, 1000.0],
        )
    end

    Test.@testset "there is one case, and it names itself" begin
        Test.@test length(BOMEX.all_cases(FT)) == 1
        Test.@test BOMEX.case_name(BOMEX.case(FT)) == "BOMEX"
        Test.@test BOMEX.t_end(BOMEX.case(FT)) == BOMEX.t_end()
    end

    Test.@testset "the profiles are physical over the whole column" begin
        c = BOMEX.case(FT)
        (; θ, q_tot, p, u, tke) = c.profiles
        levels = range(0.0, BOMEX.z_max(); length = 50)
        for z in levels
            Test.@test isfinite(θ(z)) && θ(z) > 0
            Test.@test isfinite(q_tot(z)) && q_tot(z) >= 0
            Test.@test isfinite(p(z)) && p(z) > 0
            Test.@test isfinite(u(z))
            Test.@test isfinite(tke(z)) && tke(z) >= 0
        end
        # hydrostatic: pressure falls monotonically away from the surface value it started at
        Test.@test p(0.0) ≈ BOMEX.reference_pressure()
        pressures = [p(z) for z in levels]
        Test.@test issorted(pressures; rev = true)
    end

    Test.@testset "prognostic_tke selects which TKE the state starts from" begin
        Test.@test BOMEX.case(FT; prognostic_tke = true).prognostic_tke
        Test.@test !BOMEX.case(FT; prognostic_tke = false).prognostic_tke
    end

    Test.@testset "parameters are ClimaAtmos's own file with overrides on top" begin
        # referenced inside ClimaAtmos, not vendored, so the values cannot drift from the
        # configuration they came from
        base_path = BOMEX.reference_parameter_file()
        Test.@test isfile(base_path)
        base = BOMEX._override_dict(base_path)
        Test.@test !isempty(base)

        # a later source wins, and the layers beneath survive
        name = first(sort(collect(keys(base))))
        layered = BOMEX._override_dict([
            base_path,
            Dict{String, Any}(
                "a_name_no_parameter_file_uses" =>
                    Dict{String, Any}("value" => 1.0, "type" => "float"),
            ),
        ])
        Test.@test haskey(layered, "a_name_no_parameter_file_uses")
        Test.@test haskey(layered, name)
        Test.@test length(layered) == length(base) + 1

        overridden = BOMEX._override_dict([
            base_path,
            Dict{String, Any}(name => Dict{String, Any}("value" => -12345.0)),
        ])
        Test.@test overridden[name]["value"] == -12345.0
        Test.@test length(overridden) == length(base)
    end
end
