"""
    io.jl

Serialization of a sampled [`SOCRATESForcing`](@ref): the profiles SSCF produces on
one model grid, written to a file and read back into the same interpolants.

The file is a serialization of the in-memory state, not a different forcing.
`write_forcing_file` stores SSCF's own time nodes rather than resampling, and
`read_forcing_file` rebuilds `SSCF`'s interpolant type with the boundary conditions
SSCF itself uses, so the two paths evaluate identically.
"""

const FORCING_FILE_UNITS = (;
    z = "m",
    time = "s",
    dTdt_hadv = "K s^-1",
    dqtdt_hadv = "kg kg^-1 s^-1",
    T_nudge = "K",
    qt_nudge = "kg kg^-1",
    u_nudge = "m s^-1",
    v_nudge = "m s^-1",
    subsidence = "m s^-1",
    pg = "Pa",
    Tg = "K",
    Tsfc = "K",
    qg = "kg kg^-1",
    qsfc = "kg kg^-1",
)

_build(times, values, bc) = SSCF.Interpolation.build_spline(
    SSCF.Interpolation.FastLinear1DInterpolation,
    times,
    values;
    bc,
    drop_collinear = Val(false),
)

"""
    forcing_time_nodes(interpolants)

The single time axis [s] every level of every variable shares, erroring if they do
not — the file stores one axis, so a ragged set could not round-trip.
"""
function forcing_time_nodes(interpolants)
    reference = collect(Float64, first(first(interpolants)).xp)
    for (name, levels) in pairs(interpolants), (k, itp) in enumerate(levels)
        collect(Float64, itp.xp) == reference || error(
            "`$name` level $k has a different time axis from the first variable, so \
             the forcing cannot be written as one time coordinate.",
        )
    end
    return reference
end

"""
    write_forcing_file(path, case, z; thermo_params)

Sample `case` onto levels `z` and write the profiles to `path`, returning it.

The file belongs to the grid and float type it was written from: `read_forcing_file`
requires `z` to match exactly rather than regridding silently.
"""
function write_forcing_file(
    path::AbstractString,
    c::SOCRATESCase,
    z::AbstractVector;
    thermo_params,
)
    validate(c)
    (; interpolants, surface) = sample_forcing(c, z, thermo_params)
    times = forcing_time_nodes(interpolants)
    nz, nt = length(z), length(times)

    mkpath(dirname(abspath(path)))
    NC.NCDataset(path, "c") do ds
        ds.attrib["title"] = "SOCRATES column forcing, sampled onto one model grid"
        ds.attrib["case"] = case_name(c)
        ds.attrib["flight_number"] = c.flight_number
        ds.attrib["forcing_type"] = String(forcing_label(c))
        # the profiles depend on the thermodynamics they were derived with, so a run
        # whose parameters differ must not silently reuse them
        ds.attrib["latent_heat_vaporization"] =
            Float64(CA.TD.Parameters.LH_v0(thermo_params))
        ds.attrib["latent_heat_sublimation"] =
            Float64(CA.TD.Parameters.LH_s0(thermo_params))

        NC.defDim(ds, "z", nz)
        NC.defDim(ds, "time", nt)
        defvar(name, dims, values) = begin
            v = NC.defVar(ds, String(name), Float64, dims)
            v.attrib["units"] = getproperty(FORCING_FILE_UNITS, name)
            v[:] = values
        end
        defvar(:z, ("z",), collect(Float64, z))
        defvar(:time, ("time",), times)
        for name in SSCF_FORCING_VARS
            levels = interpolants[name]
            profile = Matrix{Float64}(undef, nz, nt)
            for k in 1:nz
                profile[k, :] .= collect(Float64, levels[k].fp)
            end
            defvar(name, ("z", "time"), profile)
        end
        for name in SSCF_SURFACE_VARS
            defvar(name, ("time",), collect(Float64, surface[name].fp))
        end
    end
    return path
end

"""
    read_forcing_file(path, z, thermo)

`(; interpolants, surface)` from a file written by [`write_forcing_file`](@ref),
erroring unless it was written for levels `z` and for thermodynamics matching
`thermo`.
"""
function read_forcing_file(path::AbstractString, z::AbstractVector, thermo)
    isfile(path) || error("SOCRATES forcing file not found: $path")
    return NC.NCDataset(path, "r") do ds
        file_z = Array(ds["z"])
        want_z = collect(Float64, z)
        file_z == want_z || error(
            "$path was written for a $(length(file_z))-level grid spanning \
             $(isempty(file_z) ? "nothing" : extrema(file_z)) m, but this run has \
             $(length(want_z)) levels spanning \
             $(isempty(want_z) ? "nothing" : extrema(want_z)) m. Write a file for \
             this grid and float type rather than regridding one written for another.",
        )
        for (attr, value) in (
            ("latent_heat_vaporization", CA.TD.Parameters.LH_v0(thermo)),
            ("latent_heat_sublimation", CA.TD.Parameters.LH_s0(thermo)),
        )
            ds.attrib[attr] == Float64(value) || error(
                "$path was written with $attr = $(ds.attrib[attr]) J/kg but this run \
                 uses $(Float64(value)) J/kg; the profiles derived from it would not \
                 match the model's thermodynamics.",
            )
        end

        times = Array(ds["time"])
        bc = SSCF.Interpolation.ErrorBoundaryCondition()
        interpolants = NamedTuple{SSCF_FORCING_VARS}(map(SSCF_FORCING_VARS) do name
            profile = Array(ds[String(name)])
            [_build(times, profile[k, :], bc) for k in axes(profile, 1)]
        end)
        sbc = SSCF.Interpolation.ExtrapolateBoundaryCondition()
        surface = NamedTuple{SSCF_SURFACE_VARS}(map(SSCF_SURFACE_VARS) do name
            _build(times, Array(ds[String(name)]), sbc)
        end)
        return (; interpolants, surface)
    end
end