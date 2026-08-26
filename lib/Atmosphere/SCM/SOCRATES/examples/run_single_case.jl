"""
    run_single_case.jl

A standalone forward SOCRATES SCM run.

`SOCRATES.les_outputvars` and `SOCRATES.run_outputvars` read the reference and the run onto
the same levels and window if you want to look at them.
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