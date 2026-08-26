using ClimaAnalysis: ClimaAnalysis
using ClimaAtmos: ClimaAtmos as CA
using SOCRATES: SOCRATES
using Test: Test

# Builds the real simulation, steps it, and checks the output reaches `t_end`. The unit
# tests pass whether or not the model runs, so they cannot stand in for this.
#
# Deliberately not in `runtests.jl`: run with `include("test/integration.jl")`.
Test.@testset "integration" begin
    FT = Float64
    c = SOCRATES.case("RF09_Obs")
    t_end = 120.0
    period_seconds = 60

    # A coarsened grid and a short run keep this affordable; every physics choice is the
    # one a production run uses.
    grid = SOCRATES.socrates_grid(FT, c; dz_min = 200.0)
    z = SOCRATES.socrates_z(grid)

    # the simulation writes here and the run below reads it back, so the whole testset
    # lives inside the block that owns the directory
    mktempdir() do outdir
        sim = SOCRATES.socrates_simulation(
            FT, c;
            grid,
            output_dir = outdir,
            t_end,
            dt = 10,
            diagnostics = SOCRATES.socrates_diagnostics(;
                n_levels = length(z), period_seconds,
            ),
            verbose = false,
        )

        Test.@testset "initial state" begin
            Y = sim.integrator.u
            ρ = vec(Array(parent(Y.c.ρ)))
            q_tot = vec(Array(parent(Y.c.ρq_tot))) ./ ρ
            q_lcl = vec(Array(parent(Y.c.ρq_lcl))) ./ ρ
            q_icl = vec(Array(parent(Y.c.ρq_icl))) ./ ρ
            Test.@test all(isfinite, ρ) && all(>(0), ρ)
            Test.@test all(isfinite, q_tot) && all(>=(0), q_tot)
            # the liquid-only saturation split: condensate where the sounding is
            # supersaturated, none in ice
            Test.@test any(>(0), q_lcl)
            Test.@test all(q_lcl .<= q_tot)
            Test.@test all(iszero, q_icl)
        end

        Test.@testset "forcing tendency at t = 0" begin
            Y = sim.integrator.u
            p = sim.integrator.p
            Yₜ = similar(Y)
            Yₜ .= 0
            CA.external_forcing_tendency!(
                Yₜ, Y, p, sim.integrator.t, p.atmos.external_forcing,
            )
            for name in (:ρe_tot, :ρq_tot)
                Test.@test all(
                    isfinite, vec(Array(parent(getproperty(Yₜ.c, name)))),
                )
            end
            # the initial state is the nudging target, so relaxation contributes nothing
            # to momentum and the wind tendency has to be exactly zero
            Test.@test all(iszero, vec(Array(parent(Yₜ.c.uₕ))))
        end

        Test.@testset "runs to completion and writes the scored variables" begin
            result = CA.solve_atmos!(sim)
            Test.@test result.ret_code === :success

            vars = SOCRATES.run_outputvars(sim.output_dir)
            for name in SOCRATES.SCORED_VARS
                Test.@test haskey(vars, name)
                v = vars[name]
                t = v.dims[ClimaAnalysis.time_name(v)]
                # the record has to reach t_end, not merely exist
                Test.@test maximum(t) >= t_end - period_seconds
                Test.@test all(isfinite, filter(!isnan, v.data))
            end
        end
    end
end
