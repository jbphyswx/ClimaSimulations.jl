"""
    atlas_reader.jl

Reading Atlas LES variables through [`atlas_var_specs`](@ref), so every value that comes
out has been converted by its own recipe and nothing is read with unknown units.
"""

"""
    atlas_specs(params)

The full variable table: the raw entries plus the derived ones.

The latent heats and ice density come from `params`, so the conversions and the ice size
distribution follow whatever parameter set is passed. `configs/socrates.toml` defaults them
to the values the Atlas archive was written with, which is what makes the converted values
comparable to it. `kwargs...` reach [`add_derived_specs!`](@ref).
"""
function atlas_specs(params; kwargs...)
    thermo_params = CA.Parameters.thermodynamics_params(params)
    ρ_ice = CA.Parameters.microphysics_1m_params(params).cloud.ice.ρᵢ
    return add_derived_specs!(
        atlas_var_specs(thermo_params), thermo_params, ρ_ice; kwargs...,
    )
end

"""
    atlas_dependencies(specs, names)

`names` together with everything they transitively consume, excluding the seeded
coordinates, ordered so each entry follows its inputs.
"""
function atlas_dependencies(specs, names)
    needed = Set{Symbol}()
    function visit(name)
        name in ATLAS_COORDINATES && return nothing
        name in needed && return nothing
        haskey(specs, name) || error(
            "`$name` has no Atlas spec; it is neither a classified variable nor a \
             seeded coordinate.",
        )
        push!(needed, name)
        spec = specs[name]
        foreach(visit, spec.inputs)
        return nothing
    end
    foreach(visit, names)
    return filter(in(needed), atlas_processing_order(specs))
end

"""
    read_atlas(case, names; params, check_coverage)

`(; z, time, data)` for `names`, each converted by its spec, with `time` in elapsed
seconds and `z` in metres.

Every entry the requested names depend on is computed first. `check_coverage` asserts
that the archive holds no variable the table does not classify; it is the guard that
keeps an unclassified variable from ever being read with unknown units.
"""
function read_atlas(
    c::SOCRATESCase,
    names;
    params,
    specs = atlas_specs(params),
    check_coverage::Bool = true,
)
    les = SSCF.open_atlas_les_output(c.flight_number, c.forcing_type)
    ds = les.data
    check_coverage && assert_specs_cover(specs, keys(ds))

    z = vec(atlas_values(Array(ds["z"])))
    time = _elapsed(days_to_seconds(vec(atlas_values(Array(ds["time"])))))
    seeded = Dict{Symbol, Any}(:z => z, :time => time)

    raw_of(name) = begin
        haskey(ds, String(name)) ||
            error("The Atlas archive has no variable `$name`.")
        atlas_values(Array(ds[String(name)]))
    end

    data = Dict{Symbol, Any}()
    resolved(name) = name in ATLAS_COORDINATES ? seeded[name] : data[name]
    for name in atlas_dependencies(specs, names)
        spec = specs[name]
        data[name] = if spec.kind === :raw
            declared = get(ds[String(name)].attrib, "units", "")
            declared == spec.raw_units || error(
                "Atlas variable `$name` declares units `$declared`, but its spec \
                 converts from `$(spec.raw_units)`.",
            )
            spec.transform(raw_of(name), map(resolved, spec.inputs)...)
        else
            spec.transform(
                map(raw_of, spec.raw_inputs)...,
                map(resolved, spec.inputs)...,
            )
        end
    end
    return (; z, time, data)
end

"""
    atlas_units(name; params)

The units [`read_atlas`](@ref) returns `name` in.
"""
atlas_units(name::Symbol; params) = atlas_specs(params)[name].units
