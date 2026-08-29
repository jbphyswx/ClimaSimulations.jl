"""
    reference.jl

Reading a CloudBench case's published `data.zarr`.
"""

"""Variables of the store that carry no horizontal dimensions."""
const REFERENCE_REDUCED_VARS = (
    "asr", "cloud_cover", "cloud_fraction", "cre_lw", "cre_sw", "lwp", "olr",
    "sfc_flux_rad_lw", "sfc_flux_rad_sw", "sfc_heat_flux_latent",
    "sfc_heat_flux_sensible",
)

"""
    reference_store(case; verbose)

The published `data.zarr` for `case`, opened lazily: metadata only until an array is indexed.
"""
reference_store(case; verbose::Union{Nothing, Bool} = nothing) =
    S.open_zarr(case; verbose)

"""
    reference_axes(store) -> (; z, time)

The store's own vertical levels [m] and time axis [s].

`t` is stored in nanoseconds of simulated time and spans the last day of the run.
pass `convert_time = true` to convert
"""
function reference_axes(store, ::Type{FT} = Float32; convert_time::Bool = true) where {FT}


    x = Array{FT}(store.arrays["x"])
    y = Array{FT}(store.arrays["y"])
    z = Array{FT}(store.arrays["z"])

    if convert_time
        time = Array{FT}(store.arrays["t"] ./ 1e9)
    else
        time = Array(store.arrays["t"])
    end
    return (; x, y, z, time)
end


function reference(case, names, ::Type{FT} = Float32; store = reference_store(case)) where {FT}
    names = names isa AbstractString ? (String(names),) : map(String, names)
    (; x, y, z, time) = reference_axes(store)

    return (;
        x,
        y,
        z,
        time,
        data = Dict{String, Array{FT}}(
            name => Array{FT}(store.arrays[name]) for name in names
        ),
    )
end

"""
    reference(case, names; store) -> (; z, time, data)

The arrays `names` from the published `data.zarr`, each reduced to a column.

A four-dimensional `(z, x, y, t)` field is returned as its horizontal mean, shape
`(length(z), length(time))`. A field the store already reduced is returned as stored.
A four-dimensional field is 480 × 124 × 124 × 73 `Float32`, about 2.2 GB over HTTPS.
"""
function reference_profile(case, names, ::Type{FT} = Float32; store = reference_store(case)) where {FT}
    names = names isa AbstractString ? (String(names),) : map(String, names)
    (; z, time) = reference_axes(store)
    return (;
        z,
        time,
        data = Dict{String, Array{FT}}(
            name => _column_field(store, name, FT) for name in names
        ),
    )
end

function _column_field(store, name, ::Type{FT} = Float32) where {FT}
    arr = store.arrays[name]
    n = ndims(arr)
    if n == 4
        nx, ny = size(arr, 2), size(arr, 3)
        return FT.(dropdims(sum(arr; dims = (2, 3)); dims = (2, 3)) ./ (nx * ny))
    elseif n == 1 || n == 2
        return FT.(Array(arr))
    else
        error("`$name` has $n dimensions; expected 1, 2, or 4 (`(z, x, y, t)`).")
    end
end

"""
    reference_derived(case, names, combine; store, FT) -> (z, time, data)

`combine` applied to the `(z, x, y)` slabs of `names` at each time, then averaged horizontally:
`(length(z), length(time))`.

This is the only correct way to form a quantity that is a **product** of store fields — `q_c·f(T)`,
`ρ·q` — because the horizontal mean of a product is not the product of the horizontal means.
Reading the fields with [`reference_profile`](@ref) and multiplying afterwards is biased.

One time slab of each named field is held at once, about 30 MB per field, rather than the 2.2 GB
a whole `(z, x, y, t)` field would take.
"""
function reference_derived(
    case,
    names,
    combine,
    ::Type{FT} = Float32;
    store = reference_store(case),
) where {FT}
    names = names isa AbstractString ? (String(names),) : map(String, names)
    (; z, time) = reference_axes(store)
    arrays = map(n -> store.arrays[n], names)
    for (n, a) in zip(names, arrays)
        ndims(a) == 4 || error(
            "`reference_derived` combines `(z, x, y, t)` fields; `$n` has $(ndims(a)).",
        )
    end
    nx, ny = size(first(arrays), 2), size(first(arrays), 3)
    out = Matrix{FT}(undef, length(z), length(time))
    for j in eachindex(time)
        slabs = map(a -> FT.(a[:, :, :, j]), arrays)
        combined = combine(slabs)
        out[:, j] = dropdims(sum(combined; dims = (2, 3)); dims = (2, 3)) ./ (nx * ny)
    end
    return (; z, time, data = out)
end
