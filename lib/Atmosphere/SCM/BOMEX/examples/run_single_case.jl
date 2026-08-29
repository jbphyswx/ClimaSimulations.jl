"""
    run_single_case.jl

Run BOMEX end to end.

```julia
include("examples/run_single_case.jl")

run_bomex()                                  # the shipped grid and run length
run_bomex(; n_cells = 15, t_end = 600)       # coarse and short
run_bomex(; prognostic_tke = false)          # prescribed TKE instead
```

Reading the output afterwards needs nothing but the directory:

```julia
using BOMEX
vars = BOMEX.run_outputvars(dir)
```
"""

using BOMEX: BOMEX

function run_bomex(;
    FT::Type{<:AbstractFloat} = Float64,
    output_dir::AbstractString = joinpath(@__DIR__, "output", "bomex"),
    n_cells::Integer = BOMEX.z_elem(),
    z_top::Real = BOMEX.z_max(),
    prognostic_tke::Bool = true,
    kwargs...,
)
    c = BOMEX.case(FT; prognostic_tke)
    grid = BOMEX.bomex_grid(FT, c; faces = BOMEX.uniform_faces(z_top, n_cells))
    @info "running BOMEX" levels = n_cells z_top output_dir
    return BOMEX.run_case(c; FT, output_dir, grid, kwargs...)
end
