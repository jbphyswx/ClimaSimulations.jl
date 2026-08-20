"""
    runtests.jl

Unit test suite for SOCRATES.jl.
"""

using Dates: Dates
using LinearAlgebra: LinearAlgebra
using SOCRATES: SOCRATES
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF
using Statistics: Statistics
using Test: @test, @testset, @test_throws

@testset "SOCRATES.jl Suite" begin

    @testset "Cases and Registry" begin
        cases = SOCRATES.all_cases()
        @test length(cases) == 11
        @test count(c -> c.forcing_type isa SSCF.ObsForcing, cases) == 5
        @test count(c -> c.forcing_type isa SSCF.ERA5Forcing, cases) == 6

        # Canonical naming & parsing
        c9 = SOCRATES.case("RF09_Obs")
        @test c9.flight_number == 9
        @test c9.forcing_type isa SSCF.ObsForcing
        @test SOCRATES.case_name(c9) == "RF09_Obs"

        c11 = SOCRATES.case("RF11_ERA5")
        @test c11.flight_number == 11
        @test c11.forcing_type isa SSCF.ERA5Forcing
        @test SOCRATES.case_name(c11) == "RF11_ERA5"

        # Invalid cases
        @test_throws ErrorException SOCRATES.case("RF11_Obs") # Flight 11 has no Obs artifact
        @test_throws ErrorException SOCRATES.case("Invalid_Name")

        # Run duration & droplet number
        @test SOCRATES.t_end(c9) == 12 * 3600.0
        @test SOCRATES.t_end(c11) == 14 * 3600.0
        @test SOCRATES.n_ccn(c9) == 190.0e6
        @test SOCRATES.n_ccn(c11) == 115.0e6

        # Score window
        w9 = SOCRATES.score_window(c9)
        @test w9 == (10 * 3600.0, 12 * 3600.0)
    end

    @testset "Grids and Coordinates" begin
        c = SOCRATES.case("RF09_Obs")
        z_native = SOCRATES.native_z(c)
        faces_native = SOCRATES.native_faces(c)
        @test length(faces_native) == length(z_native) + 1
        @test first(faces_native) == 0.0
        @test issorted(faces_native)

        # Centre to face roundtrip
        calc_faces = SOCRATES.faces_from_centers(z_native)
        @test calc_faces ≈ faces_native
        calc_centers = SOCRATES.centers_from_faces(faces_native)
        @test calc_centers ≈ z_native

        # Coarsening
        coarse = SOCRATES.coarsen_faces_to_dz_min(faces_native, 50.0)
        @test length(coarse) <= length(faces_native)
        @test first(coarse) == first(faces_native)
        @test last(coarse) == last(faces_native)
        @test minimum(diff(coarse)) >= 50.0

        # Grid construction
        grid = SOCRATES.socrates_grid(Float64, c)
        z_grid = SOCRATES.socrates_z(grid)
        @test length(z_grid) == length(z_native)
        @test z_grid ≈ z_native
    end

    @testset "Parameter Composition" begin
        c = SOCRATES.case("RF09_Obs")
        dict = SOCRATES.socrates_toml_dict(Float64, c)
        @test dict["prescribed_cloud_droplet_number_concentration"] == 190.0e6
        @test dict["cloud_liquid_terminal_velocity_scaling_factor"] == 1.0

        params = SOCRATES.socrates_params(Float64, c)
        @test params isa SOCRATES.CA.Parameters.ClimaAtmosParameters
    end

    @testset "Forcing Generation and Memoization" begin
        c = SOCRATES.case("RF09_Obs")
        z = [50.0, 150.0, 300.0, 600.0, 1000.0]
        SOCRATES.empty_forcing_cache!()
        @test SOCRATES.forcing_cache_size() == 0

        arrays1 = SOCRATES.socrates_forcing_arrays(Float64, c; z, dt_sec = 300.0)
        @test SOCRATES.forcing_cache_size() == 1
        @test haskey(arrays1.column, :ta)
        @test haskey(arrays1.column, :hus)
        @test haskey(arrays1.column, :wa)
        @test size(arrays1.column[:ta], 1) == length(z)

        # Cache hit
        arrays2 = SOCRATES.socrates_forcing_arrays(Float64, c; z, dt_sec = 300.0)
        @test arrays1 === arrays2
    end

    @testset "LES Reference and Scoring" begin
        c = SOCRATES.case("RF09_Obs")
        ref = SOCRATES.les_outputvars(c; source = :sscf)
        @test haskey(ref, "clw")
        @test haskey(ref, "lwp")

        clw = ref["clw"]
        @test SOCRATES.ClimaAnalysis.has_altitude(clw)
        @test SOCRATES.ClimaAnalysis.has_time(clw)

        # Scored bounds
        bounds = SOCRATES.z_bounds(c; source = :sscf)
        @test bounds[1] == 0.0
        @test bounds[2] > bounds[1]

        # Score transform
        transform = SOCRATES.ScoreTransform()
        pv = SOCRATES.pool_var(transform, "clw", [1e-4, 2e-4, 1.5e-4])
        @test pv > 0.0
        @test isfinite(pv)
    end

end
