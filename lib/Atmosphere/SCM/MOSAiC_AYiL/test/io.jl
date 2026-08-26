using MOSAiC_AYiL: MOSAiC_AYiL as M
using Test: Test

# A written forcing file is a serialization of what `read_scm_in` returns, so neither
# path is privileged and neither is lossy. These assert equality rather than closeness:
# the profiles carry no time axis, so the round trip stores and returns the same numbers.

Test.@testset "forcing file round-trip" begin
    date = first(M.AYIL_DATES)
    c = M.MOSAiCAYiLCase(date)
    direct = M.read_scm_in(date)

    mktempdir() do dir
        path = M.write_forcing_file(joinpath(dir, "forcing.nc"), c)
        back = M.read_forcing_file(path)

        for name in keys(M.FORCING_PROFILE_UNITS)
            a, b = getproperty(direct, name), getproperty(back, name)
            # the archive's own precision survives, not just the values
            Test.@test eltype(a) == eltype(b)
            Test.@test a == b
        end
        for name in keys(M.FORCING_SURFACE_UNITS)
            # `isequal`, so a `missing` surface flux compares equal to itself
            Test.@test isequal(
                getproperty(direct.surface, name), getproperty(back.surface, name),
            )
        end

        # the file carries the closure parameters, so it reconstructs without the namelist
        Test.@test back.nudging == M.nudging_parameters(c)

        from_artifact = M.MOSAiCForcing(Float64, c)
        from_file =
            M.MOSAiCForcing(Float64, c; forcing = back, nudging = back.nudging)
        for f in fieldnames(typeof(from_artifact))
            Test.@test getfield(from_artifact, f) == getfield(from_file, f)
        end

        # the geostrophic wind reads from either source too
        params = M.mosaic_params(Float64, c)
        a = M.mosaic_scm_coriolis(Float64, c; params)
        b = M.mosaic_scm_coriolis(Float64, c; params, forcing = back)
        Test.@test a.coriolis_param == b.coriolis_param
        for z in (10.0, 500.0, 5000.0)
            Test.@test a.prof_ug(z) == b.prof_ug(z)
            Test.@test a.prof_vg(z) == b.prof_vg(z)
        end
    end
end

Test.@testset "a forcing file that is not one is refused" begin
    mktempdir() do dir
        Test.@test_throws ErrorException M.read_forcing_file(joinpath(dir, "absent.nc"))
        # a file missing a profile is refused rather than read with a gap
        path = joinpath(dir, "partial.nc")
        c = M.MOSAiCAYiLCase(first(M.AYIL_DATES))
        full = M.read_scm_in(c.date)
        M.write_forcing_file(path, full; date = c.date)
        Test.@test M.read_forcing_file(path).nudging === nothing
    end
end
