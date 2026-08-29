"""
    reference.jl

Which published CloudBench variable a run's diagnostics are compared against, and how each
side is put on the other's units.
"""

"""
    maximum_overlap_cover(cl; dims = 1)

Column cloud cover from a cloud-fraction profile, under maximum overlap.

`dims` is the vertical axis, the leading one in the published store's `(z, time)` layout.

"""
maximum_overlap_cover(cl::AbstractMatrix; dims::Int = 1) =
    dropdims(maximum(cl; dims); dims)

"""
    column_integral(data, z)

`∫ data dz` over the leading (height) axis of a `(z, time)` array, trapezoidally, giving one
value per time.
"""
function column_integral(data::AbstractMatrix, z::AbstractVector)
    size(data, 1) == length(z) || error(
        "`column_integral` got $(size(data, 1)) levels against $(length(z)) heights.",
    )
    length(z) > 1 || error("`column_integral` needs at least two levels.")
    weights = similar(collect(Float64, z))
    for k in eachindex(z)
        lo = k == 1 ? z[1] : (z[k] + z[k - 1]) / 2
        hi = k == length(z) ? z[end] : (z[k] + z[k + 1]) / 2
        weights[k] = hi - lo
    end
    return vec(sum(data .* weights; dims = 1))
end

"""
Scored quantity => how to form it from each side.

  - `model`          the ClimaAtmos short names it is built from
  - `combine`        applied to those, in order, to form it
  - `reference`      the published store arrays it is built from
  - `from_reference` applied to those, in order, giving it in the model's units
  - `units`          units of the comparable quantity, as ClimaAtmos declares them

Both sides name several inputs and reduce them the same way, so a quantity that is a
combination — `ql` is `q_c` apportioned by temperature — is expressible on the reference side
too. `from_reference` is pointwise, and for a 4-D store field it is applied to `(z, x, y)` slabs
*before* the horizontal average, because the mean of a product is not the product of the means.

`asr` is absorbed shortwave, so it is incident less outgoing. The cloud radiative effects are
clear-sky less all-sky, which needs
`RRTMGPInterface.AllSkyRadiationWithClearSkyDiagnostics`, the default in
`SwirlLMCloudBenchSim.cloudbench_model`.
"""
const REFERENCE_PAIRINGS = Dict{String, NamedTuple}(
    "lwp" => (
        model = ("lwp",), combine = only,
        reference = ("lwp",), from_reference = only, units = "kg m^-2",
    ),
    "olr" => (
        model = ("rlut",), combine = only,
        reference = ("olr",), from_reference = only, units = "W m^-2",
    ),
    "asr" => (
        model = ("rsdt", "rsut"), combine = xs -> xs[1] .- xs[2],
        reference = ("asr",), from_reference = only, units = "W m^-2",
    ),
    "cre_lw" => (
        model = ("rlutcs", "rlut"), combine = xs -> xs[1] .- xs[2],
        reference = ("cre_lw",), from_reference = only, units = "W m^-2",
    ),
    "cre_sw" => (
        model = ("rsutcs", "rsut"), combine = xs -> xs[1] .- xs[2],
        reference = ("cre_sw",), from_reference = only, units = "W m^-2",
    ),
    # ClimaAtmos writes cloud fraction as a percentage where the store is a fraction
    "cloud_fraction" => (
        model = ("cl",), combine = only,
        reference = ("cloud_fraction",), from_reference = xs -> 100 .* only(xs),
        units = "%",
    ),
    "cloud_cover" => (
        model = ("cl",), combine = xs -> maximum_overlap_cover(only(xs)),
        reference = ("cloud_cover",), from_reference = xs -> 100 .* only(xs),
        units = "%",
    ),
    "sfc_heat_flux_latent" => (
        model = ("hfls",), combine = only,
        reference = ("sfc_heat_flux_latent",), from_reference = only, units = "W m^-2",
    ),
    "sfc_heat_flux_sensible" => (
        model = ("hfss",), combine = only,
        reference = ("sfc_heat_flux_sensible",), from_reference = only, units = "W m^-2",
    ),
    # net downward at the surface, so it is down less up
    "sfc_flux_rad_lw" => (
        model = ("rlds", "rlus"), combine = xs -> xs[1] .- xs[2],
        reference = ("sfc_flux_rad_lw",), from_reference = only, units = "W m^-2",
    ),
    "sfc_flux_rad_sw" => (
        model = ("rsds", "rsus"), combine = xs -> xs[1] .- xs[2],
        reference = ("sfc_flux_rad_sw",), from_reference = only, units = "W m^-2",
    ),
    # The store carries total condensate and apportions it by temperature, so there is no
    # separate liquid or ice array to read — see `Simulation.split_q_c`.
    "ql" => (
        model = ("clw",), combine = only,
        reference = ("q_c", "T"),
        from_reference = xs -> first(S.split_q_c(xs[1], xs[2])),
        units = "kg kg^-1",
    ),
    "qi" => (
        model = ("cli",), combine = only,
        reference = ("q_c", "T"),
        from_reference = xs -> last(S.split_q_c(xs[1], xs[2])),
        units = "kg kg^-1",
    ),
    "qr" => (
        model = ("husra",), combine = only,
        reference = ("q_r",), from_reference = only, units = "kg kg^-1",
    ),
    "qs" => (
        model = ("hussn",), combine = only,
        reference = ("q_s",), from_reference = only, units = "kg kg^-1",
    ),
    "iwp" => (
        model = ("iwp",), combine = only,
        reference = ("rho", "q_c", "T"),
        from_reference = xs -> xs[1] .* last(S.split_q_c(xs[2], xs[3])),
        reference_reduce = column_integral,
        units = "kg m^-2",
    ),
    "rwp" => (
        model = ("rwp",), combine = only,
        reference = ("rho", "q_r"),
        from_reference = xs -> xs[1] .* xs[2],
        reference_reduce = column_integral,
        units = "kg m^-2",
    ),
    "swp" => (
        model = ("swp",), combine = only,
        reference = ("rho", "q_s"),
        from_reference = xs -> xs[1] .* xs[2],
        reference_reduce = column_integral,
        units = "kg m^-2",
    ),
)

"""No-op reduction, for a pairing whose reference needs no reduction over height."""
_keep_profile(data, _z) = data

"""Every ClimaAtmos diagnostic a run has to write for [`REFERENCE_PAIRINGS`](@ref)."""
scored_model_vars(pairings = REFERENCE_PAIRINGS) =
    sort!(unique(name for p in values(pairings) for name in p.model))

"""
    assert_pairings_are_diagnostics(; pairings)

Check that every model name a pairing names is a ClimaAtmos diagnostic, and that its declared
units are the ones ClimaAtmos declares.
"""
function assert_pairings_are_diagnostics(; pairings = REFERENCE_PAIRINGS)
    for (reference, pair) in pairings
        for name in pair.model
            entry = get(CA.Diagnostics.ALL_DIAGNOSTICS, name, nothing)
            isnothing(entry) &&
                error("`$reference` is paired with `$name`, which is not a diagnostic.")
            entry.units == pair.units || error(
                "`$reference` declares $(pair.units) but `$name` is $(entry.units).",
            )
        end
    end
    return nothing
end

"""
    model_comparable(run_vars, reference_name)

The quantity of a run that `reference_name` pairs with, built from `run_vars` — a mapping of
ClimaAtmos short name to value.
"""
function model_comparable(run_vars, reference_name::AbstractString)
    pair = get(REFERENCE_PAIRINGS, reference_name, nothing)
    isnothing(pair) && error(
        "`$reference_name` has no model counterpart in `REFERENCE_PAIRINGS`, which pairs \
         $(join(sort(collect(keys(REFERENCE_PAIRINGS))), ", ")).",
    )
    for name in pair.model
        haskey(run_vars, name) ||
            error("The run has no `$name`, which `$reference_name` is built from.")
    end
    return pair.combine([run_vars[name] for name in pair.model])
end

"""
    reference_comparable(case; names, store)

`(; z, time, data)` of the paired reference variables, on the model's units.
"""
function reference_comparable(
    case;
    names = sort(collect(keys(REFERENCE_PAIRINGS))),
    store = CB.reference_store(case),
)
    for name in names
        haskey(REFERENCE_PAIRINGS, name) ||
            error("`$name` has no model counterpart in `REFERENCE_PAIRINGS`.")
    end
    axes_ = CB.reference_axes(store)
    data = Dict{String, Any}()
    for name in names
        pair = REFERENCE_PAIRINGS[name]
        sources = collect(String, pair.reference)
        reduce_over_z = get(pair, :reference_reduce, _keep_profile)
        if any(n -> ndims(store.arrays[n]) == 4, sources)
            # Combine the (z, x, y) slabs first: the horizontal mean of a product is not the
            # product of the horizontal means.
            combined =
                CB.reference_derived(case, sources, pair.from_reference; store).data
            data[name] = reduce_over_z(combined, axes_.z)
        else
            raw = CB.reference_profile(case, sources; store)
            data[name] =
                reduce_over_z(pair.from_reference([raw.data[n] for n in sources]), axes_.z)
        end
    end
    return (; axes_.z, axes_.time, data)
end

