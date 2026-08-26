using Test: Test
using MOSAiC_AYiL: MOSAiC_AYiL as MA
using ClimaAtmos: ClimaAtmos as CA

# The real entry point: the other files test a piece, this tests that the pieces
# assemble into a model that steps. Coarsened and short because the point is the
# path, not the science — a scientific comparison is `run_case` at full resolution
# against the LES output.
#
# Run with `include("test/integration.jl")`, not from `runtests.jl`.

Test.@testset "a case runs" begin
    # the whole testset lives inside the block: the simulation writes into `out`, and the
    # solve at the end needs that directory to still exist
    mktempdir(; prefix = "ayil_test_") do out
        c = MA.case("20200503")
        # truncated at the day's usable top and coarsened, as a comparison run is
        FT = Float64
        faces =
            MA.truncate_faces_to_top(MA.native_faces(c), MA.best_simulation_top(c))
        sim = MA.mosaic_simulation(
            FT,
            c;
            output_dir = out,
            grid = MA.mosaic_grid(
                FT, c; faces = MA.coarsen_faces_to_dz_min(faces, 200),
            ),
            t_end = 120,
            dt = 10,
            verbose = false,
        )

        Y = sim.integrator.u
        p = sim.integrator.p
        ρ = vec(Array(parent(Y.c.ρ)))
        q_tot = vec(Array(parent(Y.c.ρq_tot))) ./ ρ
        Test.@test all(isfinite, ρ) && all(>(0), ρ)
        Test.@test all(isfinite, q_tot) && all(>=(0), q_tot)

        # the forcing cache carries the sampled targets on the model levels
        ef = p.external_forcing
        Test.@test length(ef.z) == length(ρ)
        Test.@test all(isfinite, vec(Array(parent(ef.ᶜT_nudge))))
        Test.@test all(isfinite, vec(Array(parent(ef.ᶜls_subsidence))))

        # the inversion the nudging hangs off is inside the search window, computed the
        # way the tendency computes it, from the cache's own constants
        (; ᶜT, ᶜp, ᶜq_liq) = p.precomputed
        @. ef.ᶜθ_l = (ᶜT - ef.Lv_over_cp * ᶜq_liq) * (ef.p_ref / ᶜp)^ef.Rd_over_cp
        z_inv = MA.inversion_height(ef.ᶜθ_l, ef.z, 100.0, 5000.0)
        Test.@test 100 < z_inv < 5000

        result = CA.solve_atmos!(sim)
        Test.@test result.ret_code === :success
        Test.@test Float64(float(sim.integrator.t)) ≈ 120

        ρ_end = vec(Array(parent(sim.integrator.u.c.ρ)))
        q_end = vec(Array(parent(sim.integrator.u.c.ρq_tot))) ./ ρ_end
        Test.@test all(isfinite, ρ_end) && all(>(0), ρ_end)
        Test.@test all(isfinite, q_end)
    end
end
