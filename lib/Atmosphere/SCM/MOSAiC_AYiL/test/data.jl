using Test: Test
using MOSAiC_AYiL: MOSAiC_AYiL as MA

Test.@testset "days present" begin
    dates = MA.available_dates()
    Test.@test length(dates) == 190
    Test.@test first(dates) == "20191016"
    Test.@test last(dates) == "20200911"
    Test.@test issorted(dates)
    Test.@test all(d -> length(d) == 8 && all(isdigit, d), dates)
end

Test.@testset "scm_in reader" begin
    c = MA.case("20200503")
    fd = MA.read_scm_in(c.date)

    # native precision: the files are Float32 and nothing widens them
    Test.@test eltype(fd.z) === Float32
    Test.@test eltype(fd.ta) === Float32
    Test.@test eltype(fd.wa) === Float32
    Test.@test fd.surface.ps isa Float32

    # ascending height, and the level count is read rather than assumed
    Test.@test issorted(fd.z)
    Test.@test length(fd.z) == 3040
    for name in (:ta, :hus, :q, :ql, :qi, :ua, :va, :p, :o3, :n_ccn, :wa,
                 :tntha, :tnhusha, :tnua, :tnva, :ug, :vg)
        Test.@test length(getproperty(fd, name)) == length(fd.z)
    end

    # total water is the sum of the three species
    Test.@test fd.hus ≈ fd.q .+ fd.ql .+ fd.qi

    # the surface fluxes are netCDF fill in every file, and must never surface as
    # a number
    Test.@test fd.surface.sensible_heat_flux === missing
    Test.@test fd.surface.latent_heat_flux === missing

    # longitude is wrapped into [-180, 180)
    Test.@test -180 <= fd.surface.trajectory_longitude < 180
end

Test.@testset "the archive's forcing is constant in time" begin
    # each file is one 05:00-11:00 UTC composite written twice; the second record
    # carries no new information, which is what lets the forcing be sampled once
    c = MA.case("20200503")
    Test.@test MA.read_scm_in(c.date; time_index = 1).ta ==
               MA.read_scm_in(c.date; time_index = 2).ta
end

Test.@testset "surface temperature blend" begin
    c = MA.case("20200503")
    fd = MA.read_scm_in(c.date)
    (; sea_ice_fraction, t_skin_ocean, t_skin_seaice, t_skin) = fd.surface
    blended = MA.surface_temperature(fd)
    Test.@test blended ≈
               (1 - sea_ice_fraction) * t_skin_ocean +
               sea_ice_fraction * t_skin_seaice
    # the file's own `t_skin` is that blend, so the two agree
    Test.@test blended ≈ t_skin atol = 1.0e-4
end

Test.@testset "LES output" begin
    c = MA.case("20200503")
    z, ρ = MA.les_density(c.date)
    faces = MA.native_faces(c)

    Test.@test length(z) == 286
    Test.@test length(ρ) == length(z)
    Test.@test all(>(0), ρ)
    Test.@test issorted(z)

    # one more face than centre, ascending, starting at the ground
    Test.@test length(faces) == length(z) + 1
    Test.@test issorted(faces)
    Test.@test first(faces) == 0
    Test.@test last(faces) ≈ 11949.3 atol = 0.1

    # the top face is as far above the last centre as that centre is above the
    # last stored face
    Test.@test last(faces) ≈ 2 * last(z) - faces[end - 1]
end

Test.@testset "namelist" begin
    c = MA.case("20200503")
    Test.@test MA.t_end(c) == 7200
    Test.@test MA.nudging_parameters(c).timescale == 10800
    Test.@test MA.nudging_parameters(c).ramp_depth == 300
    Test.@test MA.nudging_parameters(c).z_min == -1
    Test.@test MA.surface_roughness(c) == (8.0e-4, 8.5e-4)
end

Test.@testset "grid" begin
    c = MA.case("20200503")
    faces = MA.native_faces(c)

    # coarsening keeps the ends and never leaves a cell thinner than dz_min
    for dz_min in (25, 50, 200)
        kept = MA.coarsen_faces_to_dz_min(faces, dz_min)
        Test.@test first(kept) == first(faces)
        Test.@test last(kept) == last(faces)
        Test.@test issorted(kept)
        Test.@test minimum(diff(kept)) >= dz_min
    end
    # no coarsening requested leaves the grid alone
    Test.@test MA.coarsen_faces_to_dz_min(faces, nothing) == faces

    # truncation cuts at the requested top and makes it the new top face
    for top in (2500.0, 5000.0)
        kept = MA.truncate_faces_to_top(faces, top)
        Test.@test last(kept) ≈ top
        Test.@test all(<=(top), kept)
        Test.@test first(kept) == first(faces)
        Test.@test issorted(kept)
    end
    Test.@test MA.truncate_faces_to_top(faces, nothing) == faces

    grid = MA.mosaic_grid(
        Float64, c;
        faces = MA.coarsen_faces_to_dz_min(faces, 200),
    )
    z = MA.mosaic_z(grid)
    Test.@test issorted(z)
    Test.@test first(z) > 0
    Test.@test last(z) < last(faces)
end

Test.@testset "simulation tops" begin
    # every day the table admits is a day the archive has
    Test.@test issubset(keys(MA.BEST_SIMULATION_TOP_F), MA.available_dates())
    # and each is at or below the reference's own column
    Test.@test all(<=(MA.LES_TOP_FACE), values(MA.BEST_SIMULATION_TOP_F))
    Test.@test MA.best_simulation_top(MA.case("20200210")) == 2500.0
    # a day with no entry has no meaningful top, and says so
    Test.@test !haskey(MA.BEST_SIMULATION_TOP_F, "20191025")
    Test.@test_throws ErrorException MA.best_simulation_top(MA.case("20191025"))
end

Test.@testset "a centre becomes the face above it" begin
    faces = MA.LES_FACES
    # idempotent on a face, and never moves a value down
    for f in faces
        Test.@test MA.mosaic_ayil_face_above_center(f) == f
    end
    Test.@test all(
        MA.mosaic_ayil_face_above_center(z) >= z for z in MA.LES_FACES[1:(end - 1)] .+ 1
    )
    # the filter's "nothing flagged" value is the topmost centre, and has to become
    # the top face or a full-column run loses its top cell
    Test.@test MA.mosaic_ayil_face_above_center(11857.2) == MA.LES_TOP_FACE
    Test.@test !any(==(11857.2), values(MA.BEST_SIMULATION_TOP_F))
    Test.@test MA.best_simulation_top(MA.case("20200503")) == MA.LES_TOP_FACE
    # so truncating at it reproduces the reference column exactly
    c = MA.case("20200503")
    Test.@test MA.truncate_faces_to_top(
        MA.native_faces(c), MA.best_simulation_top(c),
    ) == MA.native_faces(c)
    # above the column is an error, not a silent clamp
    Test.@test_throws ErrorException MA.mosaic_ayil_face_above_center(20000.0)
end
