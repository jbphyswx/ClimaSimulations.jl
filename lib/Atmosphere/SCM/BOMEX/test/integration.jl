using BOMEX: BOMEX
using ClimaAtmos: ClimaAtmos as CA
using Test: Test

# A GATED JOB, deliberately absent from runtests.jl: constructing an AtmosSimulation alone costs
# minutes and gigabytes even on a coarse column, so it is not suite work. Run it directly.
#
#   jlrepl.sh run test/integration.jl <BOMEX>/test
#
# The point is that the machinery works — it constructs, it steps, it writes what it declares.
# How long a physically meaningful BOMEX run needs is a research question, not this one.

const FT = Float64
const COARSE_FACES = BOMEX.uniform_faces(BOMEX.z_max(), 15)

Test.@testset "a case runs" begin
    c = BOMEX.case(FT)
    grid = BOMEX.bomex_grid(FT, c; faces = COARSE_FACES)
    z = BOMEX.bomex_z(grid)

    mktempdir() do dir
        sim = BOMEX.bomex_simulation(
            FT,
            c;
            output_dir = joinpath(dir, "run"),
            grid,
            t_end = 120,
            dt = 20,
            diagnostics = BOMEX.bomex_diagnostics(;
                n_levels = length(z), period_seconds = 60,
            ),
            verbose = false,
        )
        Test.@test sim isa CA.AtmosSimulation

        Y = sim.integrator.u
        ρ = vec(Array(parent(Y.c.ρ)))
        Test.@test all(isfinite, ρ) && all(>(0), ρ)

        result = CA.solve_atmos!(sim)
        Test.@test result.ret_code === :success
        # `integrator.t` is an `ITime` and does not compare against a plain number
        Test.@test Float64(float(sim.integrator.t)) ≈ 120

        ρ_end = vec(Array(parent(sim.integrator.u.c.ρ)))
        Test.@test all(isfinite, ρ_end) && all(>(0), ρ_end)

        written = BOMEX.run_outputvars(sim.output_dir)
        for name in BOMEX.default_diagnostic_vars
            var = written[name]
            Test.@test all(isfinite, Array(var.data))
        end
    end
end
