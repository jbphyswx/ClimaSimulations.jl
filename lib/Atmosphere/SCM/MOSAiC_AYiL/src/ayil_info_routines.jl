"""
    ayil_info_routines.jl

The filters behind [`BEST_SIMULATION_TOP_F`](@ref): where a day's reference ice
stops being reproducible, hence how high that day can be simulated, and where its
cloud top sits.

Runnable, so the tables in `ayil_info.jl` can be re-derived rather than trusted.
`docs/design.md` §12 records why the filters exist.

Reads the archive directly — `profiles.001.nc` and `mphysprofiles.001.nc` — with
one open per day shared by every filter. All 190 days take about ten seconds.
"""

"""
SB3 cloud-ice parameters, from DALES `modmicrodata3.f90:138-223`.

Mean particle mass `x` [kg] maps to diameter as `D = a x^b` and to fall speed as
`v = α x^β (ρ_ref/ρ)^γ`, with `x` bounded to `[x_min, x_max]`.

These locate the anomaly: the 60 μm crystals DALES initializes have
`x = 1.7e-12 kg`, essentially the `x_min` floor, and fall at 2.7 mm/s, where a
1 mm crystal falls at 7.6 cm/s.
"""
const SB3_ICE = (;
    a = 0.217,
    b = 0.302115,
    α = 41.9,
    β = 0.36,
    γ = 0.5,
    x_min = 1.0e-12,
    x_max = 1.0e-7,
    ρ_ref = 1.225,
)

"""
    sb3_mean_ice_mass(q_ice, n_ice)

Mean cloud-ice particle mass [kg], `q/n`, clamped to SB3's bounds.

Both arguments are per unit *mass*: `sv007` is a specific number, since
`initiceclouds3` forms it as a mixing ratio over a particle mass
(`modbulkmicro3.f90:2101-2102`). A number carried per unit *volume* would need a
factor of ρ here.
"""
sb3_mean_ice_mass(q_ice, n_ice) =
    clamp(q_ice / n_ice, SB3_ICE.x_min, SB3_ICE.x_max)

"""Mean cloud-ice diameter [m] from the mean particle mass, `D = a x^b`."""
sb3_ice_diameter(x) = SB3_ICE.a * x^SB3_ICE.b

"""Mean cloud-ice fall speed [m/s] from the mean particle mass and air density."""
sb3_ice_fall_speed(x, ρ) =
    SB3_ICE.α * x^SB3_ICE.β * (SB3_ICE.ρ_ref / ρ)^SB3_ICE.γ

# --- One read per day, shared by every filter ------------------------------- #

"""
    ayil_ice_fields(date; root, time_inds, q_min)

The fields the filters need, read once: heights, density, cloud liquid and ice,
ice number, and the mask of levels and times SB3 counts as ice-bearing.

`time_inds` selects the output window. The first record is skipped by default, since some
cases are in principle problematic but start with subsaturated ice which rapidly ameliorates.

We use a window over the beginning by default to speed up analysis, and before significant sedimentation has occured, but
not too small that anomalous ice's transport is ignored in setting the cloud top.
"""
function ayil_ice_fields(
    date;
    root = data_root(),
    time_inds = 2:12,
    q_min::Real = 0.0,
)
    path = les_profiles_path(date; root)
    isfile(path) || error("No DALES output at $path")
    z, ρ, q_liq, q_ice, n_ice = NC.NCDataset(path, "r") do ds
        (
            Float64.(vec(Array(ds["zt"]))),
            Float64.(Array(ds["rhof"])),
            Float64.(Array(ds["sv005"])),
            Float64.(Array(ds["sv008"])),
            Float64.(Array(ds["sv007"])),
        )
    end
    # `ice_rate` is DALES's `sed_q/ρ` — a mixing ratio times a fall speed, despite
    # the raw file labelling it kg/m2 — so dividing out the mixing ratio leaves the
    # mass-weighted fall speed the reference actually realized.
    ice_rate = NC.NCDataset(joinpath(day_dir(date; root), "mphysprofiles.001.nc"), "r") do ds
        var = ds["ice_rate"]
        a = Float64.(Array(var))
        first(NC.dimnames(var)) == "time" ? permutedims(a, (2, 1)) : a
    end
    at(a) = view(a, :, time_inds)
    qi, ni, ql, ρi, rate = at(q_ice), at(n_ice), at(q_liq), at(ρ), at(ice_rate)
    present = (qi .> q_min) .& isfinite.(qi) .& isfinite.(ni)
    return (;
        z, time_inds, ρ = ρi, q_liq = ql, q_ice = qi, n_ice = ni,
        ice_rate = rate, present,
    )
end

"""Highest `z` strictly below the lowest flagged level; `z[end]` if none is flagged."""
function z_max_below_flagged(z, flagged_by_level)
    k = findfirst(flagged_by_level)
    isnothing(k) && return last(z)
    return k == 1 ? 0.0 : z[k - 1]
end

"""
    trim_top_adjacent_cloud(z, z_max, qs; top_tol)

Lower `z_max` below any cloud touching it.

Truncating at `z_max` can leave a cloud sliced against the top; for each field in
`qs`, cloud is more than `top_tol` of that field's maximum inside the region, and
where it reaches the top level the walk continues down until clear air. `0.0` when
nothing is left.
"""
function trim_top_adjacent_cloud(z, z_max, qs; top_tol = eps(Float64))
    for q in qs
        inside = findall(<=(z_max), z)
        isempty(inside) && return 0.0
        window = view(q, inside, :)
        q_max = maximum(window; init = 0.0)
        q_max > 0 || continue
        cloudy = [any(>(top_tol * q_max), view(window, i, :)) for i in axes(window, 1)]
        last(cloudy) || continue
        i = length(cloudy)
        while i > 0 && cloudy[i]      # walk down through the cloud on the boundary
            i -= 1
        end
        z_max = i == 0 ? 0.0 : z[inside[i]]
    end
    return z_max
end

"""Cloud of either phase somewhere at or below `z_max`."""
function has_cloud_below(fields, z_max; cloud_min::Real = 1.0e-8)
    below = fields.z .<= z_max
    any(view(fields.q_liq, below, :) .+ view(fields.q_ice, below, :) .> cloud_min)
end

# --- Filters ---------------------------------------------------------------- #
#
# Each takes the shared fields and returns `(; valid, z_max)`. A driver caps the
# day at the smallest `z_max` any filter allows, so each does its own trim and its
# own cloud check rather than relying on a combined flag.

"""No constraint; the column's own top. The identity of the filter set."""
ayil_identity_filter(fields) = (; valid = true, z_max = last(fields.z))

"""
    ice_size_floor(fields; D_min, top_tol, q_min)

Invalid where the mean ice crystal is smaller than `D_min` while significant ice
is present.

`D_min = 130 μm` sits between the ~60 μm DALES pins at initialization and the
~1 mm of an ice cloud that has evolved.
"""
function ice_size_floor(
    fields;
    D_min::Real = 130.0e-6,
    q_min::Real = 1.0e-8,
    top_tol = eps(Float64),
    cloud_min::Real = 1.0e-8,
)
    (; z, q_liq, q_ice, n_ice, present) = fields
    x_floor = (D_min / SB3_ICE.a)^(1 / SB3_ICE.b)
    # Stage 1: find the valid region from the mean-crystal-size criterion.
    bad = present .& (n_ice .> 0) .& (q_ice .> q_min) .& (q_ice ./ n_ice .< x_floor)
    z_max = z_max_below_flagged(z, vec(any(bad, dims = 2)))
    # Stage 2: cutoff any cloud adjacent to the top of the domain, in liquid or ice
    z_max = trim_top_adjacent_cloud(z, z_max, (q_liq, q_ice); top_tol)
    # Stage 3: make sure there is actually cloud in the surviving region.
    return (; valid = z_max > 0 && has_cloud_below(fields, z_max; cloud_min), z_max)
end

"""
    ice_fall_speed_floor(fields; v_min, top_tol, q_min)

Invalid where the mass-weighted ice fall speed stays below `v_min` while
significant ice is present — ice that neither sediments nor autoconverts.

The speed is the *realized* one, `ice_rate / q_ice`, which carries the large
particles that do most of the falling. [`sb3_ice_fall_speed`](@ref) computes the
speed SB3 predicts from the mean state instead; that is a strictly increasing
function of the same mean particle mass [`ice_size_floor`](@ref) tests, so it adds
no independent constraint and its magnitudes do not share this threshold.
"""
function ice_fall_speed_floor(
    fields;
    v_min::Real = 1.0e-2,
    q_min::Real = 1.0e-8,
    top_tol = eps(Float64),
    cloud_min::Real = 1.0e-8,
)
    (; z, q_liq, q_ice, ice_rate, present) = fields
    v = ice_rate ./ q_ice

    # Stage 1: find the valid region from the ice fall-speed criterion.
    bad = present .& isfinite.(v) .& (v .< v_min) .& (q_ice .> q_min)
    z_max = z_max_below_flagged(z, vec(any(bad, dims = 2)))
    # Stage 2: cutoff any cloud adjacent to the top of the domain, in liquid or ice
    z_max = trim_top_adjacent_cloud(z, z_max, (q_liq, q_ice); top_tol)
    # Stage 3: make sure there is actually cloud in the surviving region.
    return (; valid = z_max > 0 && has_cloud_below(fields, z_max; cloud_min), z_max)
end

"""
    extreme_n_ice(fields; kwargs...)

Invalid where the ice number itself is extreme, independently of the size and
fall-speed criteria. Not implemented; the two filters above already reject the
days this was meant to catch.
"""
extreme_n_ice(fields; kwargs...) =
    error("extreme_n_ice is not implemented; use ice_size_floor and ice_fall_speed_floor")

# --- Drivers ---------------------------------------------------------------- #

"""The filters used in production."""
const CANONICAL_ICE_FILTERS = (ice_size_floor, ice_fall_speed_floor)

"""
    best_z_maxs(dates = AYIL_DATES; filters, z_floor, root, time_inds, kwargs...)

`date => z_max` for every day every filter accepts, capped at the smallest `z_max`
any of them allows.

Days that fail a filter, or that survive only below `z_floor`, are dropped. This
is what [`BEST_SIMULATION_TOP_F`](@ref) derives from, before the per-day hand check
and the collapse onto a few grid heights that the table time_inds.
"""
function best_z_maxs(
    dates = AYIL_DATES;
    filters = CANONICAL_ICE_FILTERS,
    z_floor::Real = 500.0,
    root = data_root(),
    time_inds = 2:12,
    kwargs...,
)
    kept = Dict{String, Float64}()
    for date in dates
        fields = ayil_ice_fields(date; root, time_inds)
        z_max = last(fields.z)
        valid = true
        for filter in filters
            out = filter(fields; kwargs...)
            out.valid || (valid = false; break)
            z_max = min(z_max, out.z_max)
        end
        (valid && z_max > z_floor) && (kept[String(date)] = z_max)
    end
    return kept
end

"""
    get_cloud_tops(dates = keys(BEST_SIMULATION_TOP_F); z_tops, tol, root, time_inds)

`(; cloud_top, failed)`: the highest level carrying at least `tol` of the day's
peak cloud liquid or ice, within the column each day is capped at.

Clouds sitting against the cap are walked down past, so a sliced cloud does not set the top;
a day whose cloud is entirely against the boundary has no top to report and is returned in
`failed` rather than given one derived from the cut.

`z_tops` is the cap per day. It defaults to [`BEST_SIMULATION_TOP_F`](@ref), the
heights days are actually simulated to, but the caps the filters themselves
produced ([`best_z_maxs`](@ref)) are a different and equally sensible choice, so
it is an argument rather than wired in.

`time_inds` spans the whole record here, unlike the ice filters' opening window:
where the cloud reaches is a question about the run, not about its initialization,
and restricting it puts the top hundreds to thousands of metres too low. With the
full record and the filters' own caps this reproduces every one of the 24
hand-checked days, `FAILED` cases included.
"""
function get_cloud_tops(
    dates = keys(BEST_SIMULATION_TOP_F);
    z_tops = BEST_SIMULATION_TOP_F,
    tol::Real = 0.01,
    root = data_root(),
    time_inds = Colon(),
)
    cloud_top = Dict{String, Float64}()
    failed = String[]
    for date in dates
        fields = ayil_ice_fields(date; root, time_inds)
        z_max = get(z_tops, String(date), last(fields.z))
        inside = findall(<=(z_max), fields.z)
        tops = Float64[]
        ok = true
        for (q, phase) in ((fields.q_liq, :liquid), (fields.q_ice, :ice))
            window = view(q, inside, :)
            q_max = maximum(window; init = 0.0)
            iszero(q_max) && continue      # this phase has no cloud; not a failure
            result = _cloud_top_below_boundary(fields.z, inside, window, q_max, tol)
            if result.success
                push!(tops, result.z_top)
            else
                @error "cloud only against the domain top" date phase z_max q_max
                ok = false
            end
        end
        if ok
            cloud_top[String(date)] = isempty(tops) ? z_max : maximum(tops)
        else
            push!(failed, String(date))
        end
    end
    return (; cloud_top, failed)
end

# The highest cloudy level that is not part of a cloud touching the boundary.
function _cloud_top_below_boundary(z, inside, window, q_max, tol)
    cloudy = [any(>=(tol * q_max), view(window, i, :)) for i in axes(window, 1)]
    i_top = findlast(cloudy)
    isnothing(i_top) && return (; z_top = 0.0, success = true)
    i_top < length(inside) && return (; z_top = z[inside[i_top]], success = true)

    i = length(inside)
    while i > 0 && cloudy[i]
        i -= 1
    end
    i == 0 && return (; z_top = 0.0, success = false)
    i_cloud = findlast(view(cloudy, 1:i))
    isnothing(i_cloud) && return (; z_top = 0.0, success = false)
    return (; z_top = z[inside[i_cloud]], success = true)
end