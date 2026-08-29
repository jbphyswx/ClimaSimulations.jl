using ClimaAtmos: ClimaAtmos as CA
using NCDatasets: NCDatasets as NC
using SOCRATES: SOCRATES
using Test: Test

# A forcing may be built in memory or from a file, and the file is a serialization of the
# in-memory state rather than a second source of truth. So the two have to agree wherever they
# are evaluated, not only at the sample times the file happens to store.

Test.@testset "forcing round-trip" begin
    FT = Float64
    c = SOCRATES.case("RF09_Obs")
    thermo = CA.Parameters.thermodynamics_params(CA.ClimaAtmosParameters(FT))
    # coarsened hard, because this asserts serialization rather than resolution
    z = SOCRATES.socrates_z(SOCRATES.socrates_grid(FT, c; dz_min = 500))

    memory = SOCRATES.sample_forcing(c, z, thermo)
    nodes = SOCRATES.forcing_time_nodes(memory.interpolants)
    Test.@test length(nodes) > 1

    mktempdir() do dir
        path = SOCRATES.write_forcing_file(
            joinpath(dir, "forcing.nc"), c, z; thermo_params = thermo,
        )
        from_file = SOCRATES.sample_forcing(path, z, thermo)

        # the stored sample times, then between them: an interpolant that agreed only at its
        # nodes would still be a different function
        midpoints = [(nodes[i] + nodes[i + 1]) / 2 for i in 1:(length(nodes) - 1)]

        Test.@testset "profiles agree at the stored times" begin
            for name in SOCRATES.SSCF_FORCING_VARS
                mem = memory.interpolants[name]
                fil = from_file.interpolants[name]
                Test.@test length(mem) == length(fil) == length(z)
                for k in eachindex(mem), t in nodes
                    Test.@test mem[k](t) ≈ fil[k](t)
                end
            end
        end

        Test.@testset "profiles agree between the stored times" begin
            for name in SOCRATES.SSCF_FORCING_VARS
                mem = memory.interpolants[name]
                fil = from_file.interpolants[name]
                for k in eachindex(mem), t in midpoints
                    Test.@test mem[k](t) ≈ fil[k](t)
                end
            end
        end

        Test.@testset "the surface series agrees" begin
            for name in SOCRATES.SSCF_SURFACE_VARS
                for t in vcat(nodes, midpoints)
                    Test.@test memory.surface[name](t) ≈ from_file.surface[name](t)
                end
            end
        end

        Test.@testset "a file written for another grid is refused, not regridded" begin
            Test.@test_throws ErrorException SOCRATES.read_forcing_file(
                path, z[1:(end - 1)], thermo,
            )
            Test.@test_throws ErrorException SOCRATES.read_forcing_file(
                path, z .+ 1.0, thermo,
            )
        end

        Test.@testset "a file written with other thermodynamics is refused" begin
            # the profiles depend on the latent heats they were derived with, so the guard is
            # on the recorded value rather than on the file's provenance
            NC.NCDataset(path, "a") do ds
                ds.attrib["latent_heat_vaporization"] =
                    ds.attrib["latent_heat_vaporization"] + 1.0
            end
            Test.@test_throws ErrorException SOCRATES.read_forcing_file(path, z, thermo)
        end
    end
end
