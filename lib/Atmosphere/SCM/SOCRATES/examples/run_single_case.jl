"""
    run_single_case.jl

Example script demonstrating a standalone forward SOCRATES SCM simulation and LES scoring.
Layer 1 only (zero calibration dependencies required).
"""

using SOCRATES: SOCRATES

# Select case
case = SOCRATES.case("RF09_Obs")
output_dir = joinpath(@__DIR__, "output", SOCRATES.case_name(case))

@info "Running SOCRATES SCM simulation" case = SOCRATES.case_name(case) output_dir

# Build and run simulation
SOCRATES.run_case(
    case;
    FT = Float64,
    output_dir,
    verbose = true,
)

# Compare diagnostics against Atlas LES reference
comparison = SOCRATES.compare_to_les(case, output_dir)
SOCRATES.print_comparison(comparison)
