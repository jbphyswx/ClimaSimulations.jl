"""
    grid.jl

The single-column vertical grid, built from the domain the CloudBench LES ran on.

`SwirlLMCloudBench.CLOUDBENCH_LES_GRID` records that domain — 6000 m deep at a uniform
12.5 m — so a column at the reference resolution is 480 cells. That is affordable for one
case and not for an ensemble, hence `dz_min`.
"""

"""
    les_faces(; les_grid = S.CLOUDBENCH_LES_GRID)

Cell faces [m] of the LES's own uniform vertical grid, from the ground to the domain top.
"""
function les_faces(; les_grid = S.CLOUDBENCH_LES_GRID)
    n = round(Int, les_grid.lz / les_grid.dz)
    n * les_grid.dz ≈ les_grid.lz || error(
        "The recorded CloudBench domain, lz = $(les_grid.lz) m at dz = $(les_grid.dz) m, \
         is not a whole number of cells.",
    )
    return collect(range(0.0, les_grid.lz; length = n + 1))
end

"""
    coarsen_faces_to_dz_min(faces, dz_min)

`faces` with interior faces dropped so every cell is at least `dz_min` thick, keeping the
ground and the domain top. `nothing` leaves the grid alone.
"""
function coarsen_faces_to_dz_min(faces::AbstractVector, dz_min)
    fs = collect(Float64, faces)
    length(fs) >= 2 || error("A grid needs at least two faces, got $(length(fs))")
    isnothing(dz_min) && return fs
    minimum(diff(fs)) >= dz_min && return fs
    kept = [first(fs)]
    for f in fs[2:(end - 1)]
        (f - last(kept)) >= dz_min && push!(kept, f)
    end
    (last(fs) - last(kept)) >= dz_min ? push!(kept, last(fs)) : (kept[end] = last(fs))
    length(kept) >= 2 || error(
        "dz_min = $dz_min m leaves no cells in a column of depth \
         $(last(fs) - first(fs)) m.",
    )
    return kept
end

"""
    cloudbench_grid(FT; faces, dz_min, context)

The column grid for a CloudBench case, from explicit cell faces.

Defaults to the LES's own 480 cells; pass `dz_min` to coarsen, or `faces` to replace the
grid outright.
"""
function cloudbench_grid(
    ::Type{FT};
    faces::AbstractVector = les_faces(),
    dz_min = nothing,
    context = ClimaComms.context(),
) where {FT <: AbstractFloat}
    zf = coarsen_faces_to_dz_min(faces, dz_min)
    issorted(zf) || error(
        "Cell faces must be increasing; got $(length(zf)) faces spanning \
         $(extrema(zf)) m.",
    )
    domain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint(FT(first(zf))),
        CC.Geometry.ZPoint(FT(last(zf)));
        boundary_names = (:bottom, :top),
    )
    z_mesh = CC.Meshes.IntervalMesh(domain, CC.Geometry.ZPoint.(FT.(zf)))
    return CA.ColumnGrid(
        FT; context, z_elem = length(zf) - 1, z_max = FT(last(zf)), z_mesh,
    )
end

"""Centre-level heights [m] of `grid`, ascending."""
function cloudbench_z(grid)
    center_space = CA.get_spaces(grid).center_space
    return vec(Array(parent(CC.Fields.coordinate_field(center_space).z)))
end
