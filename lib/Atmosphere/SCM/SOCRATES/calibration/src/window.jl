"""
    window.jl

Which part of a run is compared against the Atlas LES: the time window and the vertical
region. Both are choices about the comparison, so they live here rather than in the case
package.

The window is chosen by the forcing source, so it dispatches on the case's type parameter
rather than branching at run time.
"""

"""Length [s] of the Obs averaging window: the last two hours of the run, Atlas et al. 2020."""
const OBS_WINDOW_LENGTH = 2 * 3600.0

"""
    reference_window(case, FT = Float64)

The `(t_start, t_end)` window [s] over which model and reference are both averaged.

An Obs run is averaged over its own last [`OBS_WINDOW_LENGTH`](@ref), so the window follows
`t_end` rather than repeating it. An ERA5 run is averaged over `time_bnds` relative to
`reference_time` in the Atlas `SOCRATES_summary.nc` artifact, plus the 12 h offset from the
LES start; those are fixed published metadata, so they are written out here rather than the
artifact being opened on every call. [`era5_score_window`](@ref) is that derivation, and
`test/reference.jl` checks these against it.
"""
@inline reference_window(c::SOCRATES.SOCRATESCase{<:SSCF.ObsForcing}, ::Type{FT} = Float64) where {FT} = (FT(36000), FT(43200))


function reference_window(
    c::SOCRATES.SOCRATESCase{<:SSCF.ERA5Forcing},
    ::Type{FT} = Float64,
) where {FT}
    if c.flight_number == 1
        return FT(40800), FT(44100)
    elseif c.flight_number == 9
        return FT(40200), FT(45000)
    elseif c.flight_number == 10
        return FT(39600), FT(44700)
    elseif c.flight_number == 11
        return FT(39600), FT(46800)
    elseif c.flight_number == 12
        return FT(41100), FT(43800)
    elseif c.flight_number == 13
        return FT(42600), FT(45000)
    else
        error("No ERA5 comparison window for flight $(c.flight_number)")
    end
end

"""
Upper bound [m] of the compared region, by `(flight, forcing)`.

The LES cloud top over the window, plus 1 km, rounded up to 500 m and capped at the domain
top. Flights 1 and 12 Obs sit below that — 3000 rather than 3500, 2000 rather than 2500 —
and are kept as the established compared region for this suite.
"""
const default_calibration_z_top = Dict{SOCRATES.SOCRATESCase, Float64}(
    SOCRATES.SOCRATESCase(1, SSCF.ObsForcing()) => 3000.0,
    SOCRATES.SOCRATESCase(9, SSCF.ObsForcing()) => 4000.0,
    SOCRATES.SOCRATESCase(10, SSCF.ObsForcing()) => 3500.0,
    SOCRATES.SOCRATESCase(12, SSCF.ObsForcing()) => 2000.0,
    SOCRATES.SOCRATESCase(13, SSCF.ObsForcing()) => 2000.0,
    SOCRATES.SOCRATESCase(1, SSCF.ERA5Forcing()) => 3000.0,
    SOCRATES.SOCRATESCase(9, SSCF.ERA5Forcing()) => 4500.0,
    SOCRATES.SOCRATESCase(10, SSCF.ERA5Forcing()) => 2500.0,
    SOCRATES.SOCRATESCase(11, SSCF.ERA5Forcing()) => 4000.0,
    SOCRATES.SOCRATESCase(12, SSCF.ERA5Forcing()) => 2500.0,
    SOCRATES.SOCRATESCase(13, SSCF.ERA5Forcing()) => 2500.0,
)

"""
    score_window(case)

[`reference_window`](@ref), checked against the run it is taken from.
"""
function score_window(c::SOCRATES.SOCRATESCase)
    t0, t1 = reference_window(c, Float64)
    duration = SOCRATES.t_end(c)
    (0 <= t0 < t1 <= duration) || error(
        "Comparison window ($t0, $t1) s for $(SOCRATES.case_name(c)) is not within the \
         run [0, $duration] s.",
    )
    return (t0, t1)
end

"""
    era5_score_window(flight_number)

`(t_start, t_end)` [s] derived from `time_bnds` and `reference_time` in the Atlas
`SOCRATES_summary.nc` artifact.

This opens the artifact, so it is the check on [`reference_window`](@ref) rather than the
way a window is looked up.
"""
function era5_score_window(flight_number::Integer)
    path = SSCF.atlas_socrates_summary_file(Int(flight_number))
    isfile(path) || error("Atlas SOCRATES summary file not found: $path")
    return SOCRATES.NC.NCDataset(path, "r") do ds
        flights = vec(Array(ds["flight_number"]))
        i = findfirst(==(flight_number), flights)
        isnothing(i) && error(
            "Flight $flight_number not present in $path (has $(collect(flights))).",
        )
        bnds = ds["time_bnds"]
        size(bnds) == (2, length(flights)) || error(
            "$path `time_bnds` has size $(size(bnds)); expected \
             (2, $(length(flights))).",
        )
        reference = ds["reference_time"][i]
        offsets =
            (Float64(Dates.value(Dates.Second(b - reference))) for b in bnds[:, i])
        Tuple(o + 12 * 3600.0 for o in offsets)
    end
end

"""
    default_z_bounds(case)

The `(bottom, top)` [m] of the compared region.
"""
@inline default_z_bounds(c::SOCRATES.SOCRATESCase) = (0.0, default_calibration_z_top[c])

"""Model levels inside `bounds`."""
function scored_levels(z_grid, bounds)
    lo, hi = bounds
    levels = filter(z -> lo <= z <= hi, collect(Float64, z_grid))
    isempty(levels) && error("No model levels inside compared region $bounds m.")
    return levels
end