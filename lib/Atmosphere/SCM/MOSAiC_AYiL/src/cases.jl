"""
    cases.jl

The MOSAiC AYiL case type — one Arctic day — its vertical grid, and the DALES
namelist values that define it.
"""

# --- Case ------------------------------------------------------------------- #

"""
    MOSAiCAYiLCase(date)

One AYiL day, identified by its `yyyymmdd` date string.
"""
struct MOSAiCAYiLCase
    date::String
end


Base.show(io::IO, c::MOSAiCAYiLCase) = print(io, "MOSAiCAYiLCase(", c.date, ")")

case_name(c::MOSAiCAYiLCase) = "AYiL_$(c.date)"

"""
    case(date; root = data_root())

The AYiL case for `date`, given as `yyyymmdd` or a `Date`, checked against the
days actually present.
"""
function case(date::AbstractString; root = data_root())
    (length(date) == 8 && all(isdigit, date)) ||
        error("An AYiL date is `yyyymmdd`; got `$date`.")
    day_dir(date; root)   # errors, listing what is available, if absent
    return MOSAiCAYiLCase(date)
end

case(date::Dates.Date; kwargs...) =
    case(Dates.format(date, "yyyymmdd"); kwargs...)
case(c::MOSAiCAYiLCase; kwargs...) = c

"""The calendar date of a case."""
Dates.Date(c::MOSAiCAYiLCase) = Dates.Date(c.date, Dates.dateformat"yyyymmdd")




"""
    best_simulation_top(case)

The height [m] `case` is (subjectively) best simulated to, from [`BEST_SIMULATION_TOP_F`](@ref).

Errors on a day that table has no entry for: those are the days whose reference
ice is not easily reproducible at any height.
"""
best_simulation_top(c::MOSAiCAYiLCase) =
    get(BEST_SIMULATION_TOP_F, c.date) do
        error(
            "AYiL day $(c.date) is not one of the $(length(BEST_SIMULATION_TOP_F)) \
             days a best simulation top is provided for.",
        )
    end

"""
    ayil_cases(; dates = AYIL_DATES, root = data_root())

The AYiL cases
"""
ayil_cases(; dates = AYIL_DATES, root = data_root()) = [case(d; root) for d in dates]

# --- DALES namelist --------------------------------------------------------- #

"""
    namelist(case; root = data_root())

The case's DALES `namoptions` as a flat `Dict` of key to raw string value.

Flat because the keys this package reads are unique across the namelist groups,
and a group-aware parse would be more machinery than that needs.
"""
function namelist(c::MOSAiCAYiLCase; root = data_root())
    path = namoptions_path(c.date; root)
    isfile(path) || error("No namoptions at $path")
    out = Dict{String, String}()
    for line in eachline(path)
        stripped = strip(first(split(line, '!')))          # drop trailing comment
        (isempty(stripped) || startswith(stripped, '&') || stripped == "/") && continue
        parts = split(stripped, '=', limit = 2)
        length(parts) == 2 || continue
        out[strip(parts[1])] = strip(parts[2])
    end
    isempty(out) && error("Parsed no key = value pairs from $path")
    return out
end

function _namelist_number(::Type{T}, nl, key, path_hint) where {T}
    haskey(nl, key) || error("`$key` is not in the AYiL namoptions ($path_hint)")
    raw = replace(nl[key], "d" => "e", "D" => "e")         # Fortran exponents
    value = tryparse(T, raw)
    isnothing(value) && error("`$key = $(nl[key])` is not a $T")
    return value
end

"""
    t_end(case; root = data_root())

The DALES run length [s] of a case, from `runtime` in its namoptions.

The DALES output covers `300` s to `runtime`, so a comparison against it has to
sit inside that window rather than starting at zero.
"""
t_end(c::MOSAiCAYiLCase; root = data_root()) = _namelist_number(
    Float64,
    namelist(c; root),
    "runtime",
    namoptions_path(c.date; root),
)

"""
    site_latitude(case; root) / site_longitude(case; root)

The DALES domain centre in degrees, from `xlat`/`xlon` in the namoptions.
"""
site_latitude(c::MOSAiCAYiLCase; root = data_root()) =
    _namelist_number(Float64, namelist(c; root), "xlat", namoptions_path(c.date; root))
site_longitude(c::MOSAiCAYiLCase; root = data_root()) =
    _namelist_number(Float64, namelist(c; root), "xlon", namoptions_path(c.date; root))

"""
    surface_roughness(case; root)

`(momentum, heat)` roughness lengths [m], from `z0mav`/`z0hav`.
"""
function surface_roughness(c::MOSAiCAYiLCase; root = data_root())
    nl = namelist(c; root)
    hint = namoptions_path(c.date; root)
    return (
        _namelist_number(Float64, nl, "z0mav", hint),
        _namelist_number(Float64, nl, "z0hav", hint),
    )
end

"""
    nudging_parameters(case; root)

The DALES testbed nudging parameters: the relaxation timescale `τ` [s], the ramp
depth `z_mid` [m] above the nudging onset, and `z_min` [m].

`z_min < 0` is DALES's flag for "start at the diagnosed inversion" rather than a
fixed height.
"""
function nudging_parameters(c::MOSAiCAYiLCase; root = data_root())
    nl = namelist(c; root)
    hint = namoptions_path(c.date; root)
    return (;
        timescale = _namelist_number(Float64, nl, "tb_taunudge", hint),
        ramp_depth = _namelist_number(Float64, nl, "tb_zmidnudge", hint),
        z_min = _namelist_number(Float64, nl, "tb_zminnudge", hint),
    )
end

# --- Vertical grid ---------------------------------------------------------- #

"""
    native_faces(case; root = data_root())

The DALES cell faces [m] of a case — the model's default vertical grid.
"""
native_faces(c::MOSAiCAYiLCase; root = data_root()) = les_faces(c.date; root)

"""Domain top [m]: the top face of the DALES grid."""
z_max(c::MOSAiCAYiLCase; root = data_root()) = last(native_faces(c; root))

"""
    truncate_faces_to_top(faces, z_top)

`faces` cut at `z_top`: every face below it, then `z_top` itself as the new top
face. `nothing` leaves the column alone.

A run that will be compared with the reference needs the shorter column
[`best_simulation_top`](@ref) gives, because ice above the trustworthy level
sediments and autoconverts down into the compared region (`docs/design.md` §12).
"""
function truncate_faces_to_top(faces::AbstractVector, z_top)
    isnothing(z_top) && return faces
    top = convert(eltype(faces), z_top)
    top > first(faces) || error(
        "A domain top of $top m is not above the surface face $(first(faces)) m.",
    )
    kept = filter(<(top), faces)
    push!(kept, top)
    length(kept) >= 2 ||
        error("A domain top of $top m leaves no cells above $(first(faces)) m.")
    return kept
end

"""
    coarsen_faces_to_dz_min(faces, dz_min)

`faces` with interior faces dropped so every cell is at least `dz_min` thick,
keeping the bottom and top.
"""
function coarsen_faces_to_dz_min(faces::AbstractVector, dz_min)
    length(faces) >= 2 || error("A grid needs at least two faces, got $(length(faces))")
    isnothing(dz_min) && return faces
    minimum(diff(faces)) >= dz_min && return faces
    kept = [first(faces)]
    for f in faces[2:(end - 1)]
        (f - last(kept)) >= dz_min && push!(kept, f)
    end
    (last(faces) - last(kept)) >= dz_min ? push!(kept, last(faces)) : (kept[end] = last(faces))
    length(kept) >= 2 || error(
        "dz_min = $dz_min m leaves no cells in a column of depth \
         $(last(faces) - first(faces)) m.",
    )
    return kept
end

"""
    mosaic_grid(FT, case; faces, context, root)

The column grid for `case`, from its cell faces — the DALES grid by default.

The face vector is the whole specification. Shorten or thin it by composing the
helpers above before passing it, e.g. for a run that will be compared with the
reference:

```julia
faces = truncate_faces_to_top(native_faces(c), best_simulation_top(c))
grid = mosaic_grid(FT, c; faces = coarsen_faces_to_dz_min(faces, 50))
```
"""
function mosaic_grid(
    ::Type{FT},
    c::MOSAiCAYiLCase;
    root = data_root(),
    faces::AbstractVector = native_faces(c; root),
    context = ClimaComms.context(),
) where {FT <: AbstractFloat}
    zf = FT.(faces)
    issorted(zf) || error(
        "Cell faces must be increasing; got $(length(zf)) faces spanning \
         $(extrema(zf)) m.",
    )
    domain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint(first(zf)),
        CC.Geometry.ZPoint(last(zf));
        boundary_names = (:bottom, :top),
    )
    z_mesh = CC.Meshes.IntervalMesh(domain, CC.Geometry.ZPoint.(zf))
    return CA.ColumnGrid(
        FT;
        context,
        z_elem = length(zf) - 1,
        z_max = FT(last(zf)),
        z_mesh,
    )
end

"""Centre-level heights [m] of `grid`, ascending."""
mosaic_z(grid) =
    vec(Array(parent(CC.Fields.coordinate_field(CA.get_spaces(grid).center_space).z)))