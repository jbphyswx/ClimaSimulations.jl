#=

    Machinery to convert between the raw DALES variables and the variables
    and ClimaAtmos comparable fields

    The archive carries ~700 distinct fields, so the naming is handled by rule
    rather than by list: the SB3 scalar families, the microphysics species x
    process grid, and the tendency-budget process x sample grid each follow a
    scheme, and every name in them is covered.

    Units are corrected in one place, because the files get three things wrong:
    the number scalars are per unit mass but labelled with the mass family's
    units, `precep_*`/`*_rate` claim kg/m2 but are a mixing ratio times a fall
    speed, and the potential-temperature tendencies claim K/kg/s but are K/s.
    `docs/design.md` section 3 records all three.

=#

# --- SB3 scalars ------------------------------------------------------------ #

#=
    One target vocabulary: ClimaAtmos diagnostic short names.

    Raw DALES names (`sv008`, `wsv008t`, `dq_i_sed`) are the source side and stay
    inside the readers. Everything a caller asks for is a ClimaAtmos short name,
    so the two sides of a comparison are addressed identically.

    ClimaAtmos ships names for only part of what the reference carries — it has no
    counterpart for the per-process microphysics tendencies, the conditional-sample
    budgets, the resolved/SFS flux splits or the variances. The answer is to
    register those with `add_diagnostic_variable!` rather than to invent a parallel
    naming scheme, and that is forced rather than optional: comparing a quantity
    means the *model* has to write it, which means it has to be a registered
    diagnostic in the first place. Any name this package adds is defined on types
    it owns, so ClimaAtmos itself is untouched.

    The parsers below exist for the archive's own naming schemes — the SB3 scalar
    families, the microphysics species x process grid, the budget process x sample
    grid — and what they return is a readable description for an output attribute
    or a plot axis, not a third set of identifiers.
=#

"""DALES SB3 scalar index → what that scalar holds, as a description."""
const SB3_TO_PHYSICAL = Dict{String, String}(
    "sv001" => "n_rain",
    "sv002" => "q_rain",
    "sv003" => "n_cloud_liquid",
    "sv004" => "n_inp",
    "sv005" => "q_cloud_liquid",
    "sv006" => "n_ccn",
    "sv007" => "n_cloud_ice",
    "sv008" => "q_cloud_ice",
    "sv009" => "n_snow",
    "sv010" => "q_snow",
    "sv011" => "n_graupel",
    "sv012" => "q_graupel",
)

"""
What each SB3 scalar is in prose, replacing the `scalar NNN` the archive puts in
every `sv`-family long name.
"""
const SB3_DESCRIPTION = Dict{String, String}(
    "sv001" => "rain number",
    "sv002" => "rain mass",
    "sv003" => "cloud droplet number",
    "sv004" => "ice nucleating particle number",
    "sv005" => "cloud liquid mass",
    "sv006" => "CCN number",
    "sv007" => "cloud ice number",
    "sv008" => "cloud ice mass",
    "sv009" => "snow number",
    "sv010" => "snow mass",
    "sv011" => "graupel number",
    "sv012" => "graupel mass",
)

_sb3_physical(index) =
    get(SB3_TO_PHYSICAL, "sv" * index) do
        error("sv$index is not one of DALES's twelve SB3 scalars")
    end

# --- Microphysics process rates --------------------------------------------- #

"""SB3 species code in an `mphysprofiles` name → the species."""
const SB3_SPECIES = Dict{String, String}(
    "c" => "cloud_liquid",
    "r" => "rain",
    "i" => "ice",
    "s" => "snow",
    "g" => "graupel",
    "ccn" => "ccn",
)

"""
SB3 process code → the process.

Keyed on the code in the variable name, never on the archive's long name:
`dq_i_col_rig` is announced as a graupel tendency but is filled from the *ice*
loss in that collision.
"""
const MPHYS_PROCESS = Dict{String, String}(
    "nuc" => "nucleation",
    "inuc" => "nucleation",
    "dep" => "deposition",
    "sub" => "sublimation",
    "hom" => "freezing_homogeneous",
    "het" => "freezing_heterogeneous",
    "cond" => "condensation",
    "ev" => "evaporation",
    "sc" => "self_collection",
    "br" => "breakup",
    "au" => "autoconversion",
    "ac" => "accretion",
    "agg" => "aggregation",
    "rime" => "riming",
    "rime_ic" => "riming_by_cloud_liquid",
    "rime_sc" => "riming_by_cloud_liquid",
    "rime_gc" => "riming_by_cloud_liquid",
    "rime_gr" => "riming_by_rain",
    "col_rig" => "collision_rain_ice_to_graupel",
    "col_rsg" => "collision_rain_snow_to_graupel",
    "col_si" => "collection_by_snow",
    "col_gsg" => "collection_by_graupel",
    "eme_ic" => "enhanced_melting_by_cloud_liquid",
    "eme_sc" => "enhanced_melting_by_cloud_liquid",
    "eme_gc" => "enhanced_melting_by_cloud_liquid",
    "eme_ri" => "enhanced_melting_by_rain",
    "eme_rs" => "enhanced_melting_by_rain",
    "eme_gr" => "enhanced_melting_by_rain",
    "me" => "melting",
    "cv_ig" => "conversion_to_graupel",
    "cv_sg" => "conversion_to_graupel",
    "mul" => "ice_multiplication",
    "sadj" => "saturation_adjustment",
    "sed" => "sedimentation",
    "tot" => "microphysics",
    "totf" => "all_processes",
)

"""Potential-temperature process code in a `dth*` name → the process."""
const THETA_PROCESS = Dict{String, String}(
    "mphys" => "microphysics",
    "freeze" => "freezing",
    "melt" => "melting",
    "cond" => "condensation",
    "ev" => "evaporation",
    "dep" => "deposition",
    "sub" => "sublimation",
)

"""
    BULKMICROSTAT3

The 23 variables `modbulkmicrostat3` contributes to `profiles.001.nc`, as
`raw => (description, units)`.

They need their own table because two things about them are wrong as stored.

**They are one record late, and the first and last samples are gone.**
`writebulkmicrostat3` writes at modgenstat's record counter
(`modbulkmicrostat3.f90:2301`), which profile writes take as `intent(in)` while
only genstat's time write increments it (`modstat_nc.f90:425,437`) — but the module
runs from `microsources` (`modbulkmicro3.f90:1944`, `program.f90:235`), *before*
`genstat` (`program.f90:272`) in the same iteration. So the k-th sample is written
to record k−1: the first goes to record 0 and is lost, the value at `time[t]` is
the sample for `time[t+1]`, and the last record is never written. Verified against
the data: `qrmn[t]` matches genstat's own `sv002[t+1]` to 5e-4 relative against
0.6 at every other offset, and the final record is fill on 190/190 days.
[`dales_field`](@ref) puts them back on the times they are averages for.

**Five state the wrong units**, sharing the label `W/m^2` with the two that really
are fluxes: `qrmn` is `sum(q_hr)` (`:1126`); `nrrain` is `rhof*sum(n_hr)` (`:1122`),
so per volume already; `raincount` and `preccount` are column fractions
(`:1116,1118`); `dvrmn` divides a diameter sum by a count (`:1129`). Only
`rainrate` and `precmn` carry the `rhof*rlv` that makes a flux (`:2273,2277`). The
tendencies' own `#/m3/s` is right — `Npav = rhof*avfield` (`:1416`).

`rainrate`, `precmn` and `preccount` come from `precep_l`, so they are liquid
precipitation and not total. `dvrmn` sums `modmicrodata`'s `Dvr`, which SB3 never
fills — its `imicro` guard is commented out (`:1128`) — so it carries no
information here.
"""
const BULKMICROSTAT3 = Dict{String, Tuple{String, String}}(
    "cfrac" => ("cloud_fraction", "1"),
    "rainrate" => ("rain_rate_precipitating_columns", "W/m^2"),
    "preccount" => ("precipitating_fraction", "1"),
    "nrrain" => ("n_rain", "m^-3"),
    "raincount" => ("rain_fraction", "1"),
    "precmn" => ("rain_rate", "W/m^2"),
    "dvrmn" => ("rain_mean_volume_diameter", "m"),
    "qrmn" => ("q_rain", "kg/kg"),
    "npauto" => ("n_rain_tendency_autoconversion", "m^-3 s^-1"),
    "npaccr" => ("n_rain_tendency_accretion", "m^-3 s^-1"),
    "npsed" => ("n_rain_tendency_sedimentation", "m^-3 s^-1"),
    "npevap" => ("n_rain_tendency_evaporation", "m^-3 s^-1"),
    "nptot" => ("n_rain_tendency_microphysics", "m^-3 s^-1"),
    "qrpauto" => ("q_rain_tendency_autoconversion", "kg/kg/s"),
    "qrpaccr" => ("q_rain_tendency_accretion", "kg/kg/s"),
    "qrpsed" => ("q_rain_tendency_sedimentation", "kg/kg/s"),
    "qrpevap" => ("q_rain_tendency_evaporation", "kg/kg/s"),
    "qrptot" => ("q_rain_tendency_microphysics", "kg/kg/s"),
    "qtpauto" => ("q_tot_tendency_autoconversion", "kg/kg/s"),
    "qtpaccr" => ("q_tot_tendency_accretion", "kg/kg/s"),
    "qtpsed" => ("q_tot_tendency_sedimentation", "kg/kg/s"),
    "qtpevap" => ("q_tot_tendency_evaporation", "kg/kg/s"),
    "qtptot" => ("q_tot_tendency_microphysics", "kg/kg/s"),
)

"""`mphysprofiles` names that are not tendencies."""
const MPHYS_DIAGNOSTIC = Dict{String, String}(
    "cfrac_l" => "cloud_fraction_liquid",
    "cfrac_i" => "cloud_fraction_ice",
    "cfrac_tot" => "cloud_fraction_total",
    "ice_rate" => "ice_fall_rate",
    "snow_rate" => "snow_fall_rate",
    "rain_rate" => "rain_fall_rate",
    "graupel_rate" => "graupel_fall_rate",
)

# --- Tendency budgets ------------------------------------------------------- #

"""
`samptend` field code → the field. `qr`/`nr` index `iqr`/`inr`, which SB3 aliases
onto rain.
"""
const SAMPTEND_FIELD = Dict{String, String}(
    "u" => "u", "v" => "v", "w" => "w", "thl" => "thl", "qt" => "qt",
    "qr" => "q_rain", "nr" => "n_rain",
)

"""`samptend` process code → the process."""
const SAMPTEND_PROCESS = Dict{String, String}(
    "adv" => "advective",
    "dif" => "diffusive",
    "for" => "forces",
    "cor" => "coriolis",
    "ls" => "large_scale",
    "top" => "top_boundary",
    "pois" => "pressure_gradient",
    "addon" => "addons",
    "rad" => "radiative",
    "micro" => "microphysics",
    "tot" => "total",
    "leib" => "total_leibniz",
)

"""Conditional-sample suffixes, longest first so `cldcr` is not read as `cld`."""
const SAMPTEND_SAMPLES = ("cldcr", "cldup", "buup", "cld", "upd", "all")

const SAMPTEND_SAMPLE_NAME = Dict{String, String}(
    "all" => "all",
    "upd" => "updraft",
    "buup" => "buoyant_updraft",
    "cld" => "cloud",
    "cldcr" => "cloud_core",
    "cldup" => "cloud_updraft",
)

"""
`samptend` names carrying no sample suffix.

`modsamptend.f90:184` registers `utendcor` without a sample suffix — every sample
writes to the one variable — so which sample survives is not recoverable from the
file, and the name says so rather than claiming `all`.
"""
const SAMPTEND_IRREGULAR =
    Dict{String, String}("utendcor" => "u_tendency_coriolis_unknown_sample")

# --- The general translator ------------------------------------------------- #

"""
    mphys_name(raw) -> String or nothing

The name for an `mphysprofiles` variable, or `nothing` when it is not one. A name
that *is* of the family but carries an unknown species or process errors rather
than passing through.
"""
function mphys_name(raw::AbstractString)
    haskey(MPHYS_DIAGNOSTIC, raw) && return MPHYS_DIAGNOSTIC[raw]
    m = match(r"^d(q|n)_([a-z]+)_(.+)$", raw)
    if m !== nothing
        kind, species_code, process_code = m.captures
        species = get(SB3_SPECIES, species_code) do
            error("$raw: species code `$species_code` is not an SB3 species")
        end
        process = get(MPHYS_PROCESS, process_code) do
            error("$raw: process code `$process_code` is not an SB3 process")
        end
        return string(kind, "_", species, "_tendency_", process)
    end
    m = match(r"^dth(l?)_(.+)$", raw)
    if m !== nothing
        liquid, process_code = m.captures
        process = get(THETA_PROCESS, process_code) do
            error("$raw: process code `$process_code` is not a theta process")
        end
        return string(isempty(liquid) ? "theta" : "thetal", "_tendency_", process)
    end
    return nothing
end

"""
    samptend_name(raw) -> String or nothing

The name for a `samptend` variable, `<field>_tendency_<process>_<sample>`, or
`nothing` when it is not one.
"""
function samptend_name(raw::AbstractString)
    haskey(SAMPTEND_IRREGULAR, raw) && return SAMPTEND_IRREGULAR[raw]
    m = match(r"^(thl|qt|qr|nr|u|v|w)tend([a-z]+)$", raw)
    m === nothing && return nothing
    field_code, rest = m.captures
    for sample_code in SAMPTEND_SAMPLES
        endswith(rest, sample_code) || continue
        process_code = rest[1:(end - length(sample_code))]
        isempty(process_code) && continue
        process = get(SAMPTEND_PROCESS, process_code) do
            error("$raw: process code `$process_code` is not a samptend process")
        end
        return string(
            SAMPTEND_FIELD[field_code], "_tendency_", process, "_",
            SAMPTEND_SAMPLE_NAME[sample_code],
        )
    end
    error("$raw: no known sample suffix, and not one of the irregular names")
end

"""
    dales_description(raw) -> String

A readable name for any variable the archive carries, or the name unchanged when
it belongs to no renamed family.

`svNNN` alone is unreadable downstream, and every form the archive uses is
covered:

| archive | becomes | from |
|---|---|---|
| `svNNN` | `q_ice` | mixing ratio or number |
| `svNNN2r` | `q_ice_variance_resolved` | resolved variance |
| `svpNNN` | `q_ice_tendency` | scalar tendency |
| `svptNNN` | `q_ice_tendency_turbulence` | turbulence tendency |
| `wsvNNNr/s/t` | `q_ice_flux_resolved/sfs/total` | fluxes |

An `svNNN` index that is not one of the twelve errors: a raw `sv004` downstream is
indistinguishable from a variable deliberately left alone.
"""
function dales_description(raw::AbstractString)
    haskey(BULKMICROSTAT3, raw) && return first(BULKMICROSTAT3[raw])
    m = match(r"^sv(\d{3})$", raw)
    m === nothing || return _sb3_physical(m.captures[1])
    m = match(r"^sv(\d{3})2r$", raw)
    m === nothing || return _sb3_physical(m.captures[1]) * "_variance_resolved"
    m = match(r"^svp(\d{3})$", raw)
    m === nothing || return _sb3_physical(m.captures[1]) * "_tendency"
    m = match(r"^svpt(\d{3})$", raw)
    m === nothing || return _sb3_physical(m.captures[1]) * "_tendency_turbulence"
    m = match(r"^wsv(\d{3})([rst])$", raw)
    if m !== nothing
        kind = Dict("r" => "resolved", "s" => "sfs", "t" => "total")
        return _sb3_physical(m.captures[1]) * "_flux_" * kind[m.captures[2]]
    end
    name = mphys_name(raw)
    name === nothing || return name
    name = samptend_name(raw)
    name === nothing || return name
    return String(raw)
end

# --- Units ------------------------------------------------------------------ #

"""
Per-volume units of a number quantity, and the power of density that converts it,
by every spelling the archive reaches it with.

DALES carries the `sv` scalars as specific quantities — number per unit *mass* —
and labels all twelve with the mass units of the family, so a number arrives as
`(kg/kg)`. Multiplying a variance by `ρ^2` treats the density as having no
fluctuation of its own; on a slab mean that is a good approximation, not an
identity.
"""
const NUMBER_UNITS = Dict{String, Tuple{String, Int}}(
    "kg/kg" => ("m^-3", 1),
    "(kg/kg)" => ("m^-3", 1),
    "#/kg" => ("m^-3", 1),
    "/kg" => ("m^-3", 1),
    "(kg/kg/s)" => ("m^-3 s^-1", 1),
    "#/kg/s" => ("m^-3 s^-1", 1),
    "/kg/s" => ("m^-3 s^-1", 1),
    "(kg/kg)^2" => ("m^-6", 2),
    "kg/kg m/s" => ("m^-3 m/s", 1),
)

"""
Units the archive states wrongly, and what they are.

`precep_*`/`*_rate` is `sed_q/ρ`, a mixing ratio times a fall speed, not the
`kg/m2` claimed. The `K/kg/s` on the potential-temperature tendencies is likewise
wrong: they are formed as `(L_v/(c_p Π)) dq/dt`, and `L_v/c_p` is kelvin per unit
mixing ratio, so the product is `K/s`.
"""
const MISLABELLED_UNITS =
    Dict{String, String}("kg/m2" => "kg/kg m/s", "K/kg/s" => "K/s")

"""One spelling for each unit the archive writes more than one way."""
const UNIT_SPELLINGS = Dict{String, String}(
    "(kg/kg)" => "kg/kg",
    "(kg/kg/s)" => "kg/kg/s",
    "kg/m^2 /s" => "kg/m^2/s",
    "Km/s" => "K m/s",
    "#/m3/s" => "m^-3 s^-1",
    "-" => "1",
)

"""
    dales_variable_attributes(raw, name, units, long_name)

`(units, long_name, ρ_power)` for a translated variable: the units it is really
in, a long name with `scalar NNN` replaced by what that scalar is, and the power
of air density that converts its values.

`ρ_power` is 0 for a relabelling, 1 for a number, 2 for a number variance. A
number carrying units in no known spelling errors rather than being mislabelled.
The [`BULKMICROSTAT3`](@ref) group states its own units and needs no conversion —
its numbers were multiplied by density where they were formed.
"""
function dales_variable_attributes(
    raw::AbstractString,
    name::AbstractString,
    units::AbstractString,
    long_name::AbstractString,
)
    haskey(BULKMICROSTAT3, raw) && return (last(BULKMICROSTAT3[raw]), long_name, 0)
    out = String(long_name)
    m = match(r"[Ss]calar (\d{3})", out)
    if m !== nothing
        description = get(SB3_DESCRIPTION, "sv" * m.captures[1]) do
            error("$raw: sv$(m.captures[1]) has no description")
        end
        out = replace(out, m.match => description)
        startswith(m.match, "S") && (out = uppercasefirst(out))
    end
    if startswith(name, "n_")
        converted, ρ_power = get(NUMBER_UNITS, units) do
            error("$name is a number but is labelled `$units`")
        end
        return (converted, out, ρ_power)
    end
    return (spelled_units(units), out, 0)
end

"""
    spelled_units(units)

The archive's own units, mislabellings corrected and one spelling per unit.

A relabelling only: nothing here implies a change to the values, which is what lets
it stand in for the converted units when a caller asked for untranslated data.
"""
function spelled_units(units::AbstractString)
    corrected = get(MISLABELLED_UNITS, units, String(units))
    return get(UNIT_SPELLINGS, corrected, corrected)
end

# --- The ClimaAtmos bridge -------------------------------------------------- #

"""
    dales_presf(p_face, ρ, z_center, z_face)

Cell-centre pressure [Pa], which the archive does not carry, as one hydrostatic
step up from the face below.

`profiles.001.nc`'s `presh` is DALES's **half**-level pressure: `modgenstat.f90:371`
names the slot that `:1523` fills with `presh`, `fromztop` sets `presh(1) = ps` and
steps it across whole cells, and its "Pressure at cell center" label and `zt`
dimension are both wrong. Using it as a centre pressure is a half-cell error worth
1.5 % in pressure and 0.91 K in a temperature formed from it **[verified]**.

Centres are face midpoints, so `p − ρ g (z_c − z_f)` is second-order accurate in the
cell thickness. Measured against two routes sharing none of its arithmetic — `p^κ`
interpolation of the faces, which uses no density, and `fromztop` replicated from
`ps` upward, which uses no face pressure above the surface — it agrees to 9e-5
relative at worst, 170 times smaller than the error it removes, and it lands
strictly between its bracketing faces at every level and time.

Not from the archived `thv`, which is the slab mean of the 3D `θ_v` field
(`modgenstat.f90:676`) and so differs from `θ_v` of the slab means by the θ–`q_l`
covariance — largest in cloud, where this matters most.
"""
dales_presf(p_face, ρ, z_center, z_face) =
    @. p_face - ρ * DALES_CONSTANTS.grav * (z_center - z_face)

"""
    dales_temperature(θ_l, q_l, p_face, ρ, z_center, z_face)

Slab-mean temperature [K], `θ_l Π + (L_v/c_p) q_l` on [`dales_presf`](@ref).

`q_l` and not the SB3 cloud-liquid scalar `sv005`: the thermodynamic liquid DALES's
saturation adjustment produced is what its θ_l carries.
"""
dales_temperature(θ_l, q_l, p_face, ρ, z_center, z_face) =
    @. dales_exner(dales_presf(p_face, ρ, z_center, z_face)) * θ_l +
       (DALES_CONSTANTS.L_v / DALES_CONSTANTS.cp_d) * q_l

"""
    CLIMAATMOS_FROM_DALES

ClimaAtmos diagnostic short name → how it is built from the archive, and the units
that produces.

Each entry is a function of a reader — `raw` → that variable's values in the units
the archive stores them in — so every conversion is written where it happens: the
`× ρ` that takes a `sv` number out of DALES's per-mass convention into ClimaAtmos's
per-volume one, the `× 100` from a fraction to a percentage, and the sum over
phases that makes `hus` total water rather than `qt`.

Only names ClimaAtmos has a diagnostic for appear; everything the archive carries
besides is reachable through [`dales_description`](@ref) and [`dales_field`](@ref)
under its own name. `hussn` sums snow and graupel because ClimaAtmos has one
frozen-precipitation class where DALES splits it in two. `ta` inverts DALES's θ_l
with DALES's own constants, using `ql` — the thermodynamic liquid its saturation
adjustment produced, which is what its θ_l carries — and not the SB3 cloud-liquid
scalar `sv005`.
"""
const CLIMAATMOS_FROM_DALES =
    Dict{String, @NamedTuple{f::Function, units::String}}(
        "clw" => (f = read -> read("sv005"), units = "kg kg^-1"),
        "cli" => (f = read -> read("sv008"), units = "kg kg^-1"),
        "husra" => (f = read -> read("sv002"), units = "kg kg^-1"),
        "hussn" => (f = read -> read("sv010") .+ read("sv012"), units = "kg kg^-1"),
        "hus" => (
            f = read -> read("qt") .+ read("sv008") .+ read("sv002") .+
                        read("sv010") .+ read("sv012"),
            units = "kg kg^-1",
        ),
        "ql_all" => (f = read -> read("sv005") .+ read("sv002"), units = "kg kg^-1"),
        "qi_all" => (
            f = read -> read("sv008") .+ read("sv010") .+ read("sv012"),
            units = "kg kg^-1",
        ),
        "ua" => (f = read -> read("u"), units = "m s^-1"),
        "va" => (f = read -> read("v"), units = "m s^-1"),
        "rhoa" => (f = read -> read("rhof"), units = "kg m^-3"),
        "pfull" => (
            f = read -> dales_presf(
                read("presh"), read("rhof"), read("zt"), read("zm"),
            ),
            units = "Pa",
        ),
        "cl" => (f = read -> 100 .* read("cfrac"), units = "%"),
        # radiative flux profiles, on the 286 lower faces. DALES signs them along
        # z and CF names carry the direction, so the magnitude is what is wanted —
        # as DALES itself takes wherever it uses them (`modradstat.f90:256`).
        "rld" => (f = read -> abs.(read("lwd")), units = "W m^-2"),
        "rlu" => (f = read -> abs.(read("lwu")), units = "W m^-2"),
        "rsd" => (f = read -> abs.(read("swd")), units = "W m^-2"),
        "rsu" => (f = read -> abs.(read("swu")), units = "W m^-2"),
        "rldcs" => (f = read -> abs.(read("lwdca")), units = "W m^-2"),
        "rlucs" => (f = read -> abs.(read("lwuca")), units = "W m^-2"),
        "rsdcs" => (f = read -> abs.(read("swdca")), units = "W m^-2"),
        "rsucs" => (f = read -> abs.(read("swuca")), units = "W m^-2"),
        "cdnc" => (f = read -> read("sv003") .* read("rhof"), units = "m^-3"),
        "ncra" => (f = read -> read("sv001") .* read("rhof"), units = "m^-3"),
        "ta" => (
            f = read -> dales_temperature(
                read("thl"), read("ql"), read("presh"), read("rhof"),
                read("zt"), read("zm"),
            ),
            units = "K",
        ),
    )

"""
Radiative fluxes, which ClimaAtmos reports as surface or top scalars and DALES
writes as face profiles.

DALES's faces are the 286 *lower* faces, so `:top` is the highest it has —
11765 m, about 184 m below the column's own top face — not the true TOA.

DALES signs these along z, so its downward fluxes are negative, while every CF name
here carries the direction in the name and the value is a magnitude. The magnitude
is taken, which is what DALES itself does wherever it uses them together
(`modradstat.f90:256`).
"""
const CLIMAATMOS_FLUX_FROM_DALES = Dict(
    "rlds" => ("lwd", :surface),
    "rlus" => ("lwu", :surface),
    "rsds" => ("swd", :surface),
    "rsus" => ("swu", :surface),
    "rlut" => ("lwu", :top),
    "rsut" => ("swu", :top),
    "rsdt" => ("swd", :top),
)

_archive_path(::Val{:profiles}, date, root) = les_profiles_path(date; root)
_archive_path(::Val{:mphys}, date, root) =
    joinpath(day_dir(date; root), "mphysprofiles.001.nc")
_archive_path(::Val{:samptend}, date, root) =
    joinpath(day_dir(date; root), "samptend.001.nc")
_archive_path(::Val{:tmser}, date, root) =
    joinpath(day_dir(date; root), "tmser.001.nc")

"""Which archive file a raw variable lives in."""
function dales_archive_file(raw::AbstractString)
    samptend_name(raw) === nothing || return :samptend
    mphys_name(raw) === nothing || return :mphys
    return :profiles
end

# 2D fields are stored vertical-axis-first, but take the orientation from the
# declared dimensions rather than from which axis is longer: that is a property of
# the run, not of the layout.
function _read_oriented(ds, name)
    var = ds[name]
    a = Array(var)
    ndims(a) == 2 || return a
    return first(NC.dimnames(var)) == "time" ? permutedims(a, (2, 1)) : a
end

"""
    dales_field(raw, date; root, translate_units)

One raw archive variable, as `(; name, z, time, data, units, long_name)`.

`name` is its [`dales_description`](@ref). With `translate_units` the values are
converted to the units [`dales_variable_attributes`](@ref) reports — a number
becomes per volume — so a caller never handles the archive's mislabelling itself.

A [`BULKMICROSTAT3`](@ref) variable comes back on the times it is an average for,
which is one record on from where it is stored, so it has one sample fewer than
the rest of the file and none of them missing.
"""
function dales_field(
    raw::AbstractString,
    date;
    root = data_root(),
    translate_units::Bool = true,
)
    file = dales_archive_file(raw)
    path = _archive_path(Val(file), date, root)
    isfile(path) || error("No archive file at $path")
    return NC.NCDataset(path, "r") do ds
        haskey(ds, raw) || error("`$raw` is not in $path")
        var = ds[raw]
        data = _read_oriented(ds, raw)
        description = dales_description(raw)
        raw_units = get(var.attrib, "units", "")
        long_name =
            get(var.attrib, "longname", get(var.attrib, "long_name", raw))
        units, long_name, ρ_power =
            dales_variable_attributes(raw, description, raw_units, long_name)
        if ρ_power != 0
            if translate_units
                ρ = _read_oriented(ds, "rhof")
                data = data .* ρ .^ ρ_power
            else
                # the values were left as the archive holds them, so the units must
                # say so too rather than reporting the conversion that did not happen
                units = spelled_units(raw_units)
            end
        end
        vertical = first(NC.dimnames(var)) == "time" ?
                   NC.dimnames(var)[2] : first(NC.dimnames(var))
        z = haskey(ds, vertical) ? vec(Array(ds[vertical])) : Float64[]
        time = Float64.(vec(Array(NC.variable(ds, "time"))))
        if haskey(BULKMICROSTAT3, raw)
            data = identity.(data[:, 1:(end - 1)])   # the dropped record was the only fill
            time = time[2:end]
        end
        return (; raw, description, z, time, data, units, long_name)
    end
end

"""
    climaatmos_field(short_name, date; root)

The DALES equivalent of a ClimaAtmos diagnostic, in ClimaAtmos's units.

`(; z, time, data, units)`, with `data` a `(z, time)` array for a profile and a
time series for a surface or top flux.
"""
function climaatmos_field(short_name::AbstractString, date; root = data_root())
    if haskey(CLIMAATMOS_FLUX_FROM_DALES, short_name)
        raw, level = CLIMAATMOS_FLUX_FROM_DALES[short_name]
        f = dales_field(raw, date; root, translate_units = false)
        k = level === :surface ? 1 : size(f.data, 1)
        return (; z = f.z[k], f.time, data = abs.(f.data[k, :]), units = "W m^-2")
    end
    spec = get(CLIMAATMOS_FROM_DALES, short_name) do
        error(
            "No DALES translation for ClimaAtmos `$short_name`; known names are " *
            join(sort(climaatmos_translated_names()), ", "),
        )
    end
    sources = NamedTuple[]
    function read(raw)
        f = dales_field(raw, date; root, translate_units = false)
        push!(sources, f)
        return f.data
    end
    data = spec.f(read)
    axes = first(sources)
    return (; axes.z, axes.time, data, spec.units)
end

"""Every ClimaAtmos short name this can produce from the archive."""
climaatmos_translated_names() = vcat(
    collect(keys(CLIMAATMOS_FROM_DALES)),
    collect(keys(CLIMAATMOS_FLUX_FROM_DALES)),
)

# --- Derived quantities ----------------------------------------------------- #

"""
    dales_radiative_heating(date; root, band = :total)

Radiative heating rate [K/s] of the reference column, as `(; z, time, data)`, for
`:longwave`, `:shortwave` or `:total`.

The archive's `thllwtend`/`thlswtend`/`thltend` are labelled K/s but are **potential**
temperature tendencies: `modradstat.f90:256` forms each as a face-flux divergence
over `rhof*exnf*cp*dzf`, so the `exnf` has to be put back to get a temperature
tendency. Π comes from [`dales_presf`](@ref), the archive's centre pressure being
absent — using its `presh` instead would add a further 1.5 %.

This is the quantity a radiation scheme is judged on, `thltend` being the sum of the
two bands and `thlradls` a prescribed large-scale term that is separate from both.
"""
function dales_radiative_heating(date; root = data_root(), band::Symbol = :total)
    raw = get(
        Dict(:longwave => "thllwtend", :shortwave => "thlswtend", :total => "thltend"),
        band,
    ) do
        error("No radiative heating for `$band`; try :longwave, :shortwave or :total")
    end
    f = dales_field(raw, date; root, translate_units = false)
    p = climaatmos_field("pfull", date; root)
    return (; f.z, f.time, data = f.data .* dales_exner.(p.data), units = "K s^-1")
end

"""
    dales_fall_speed(species, date; root, q_min)

Mass-weighted fall speed [m/s] of `:ice`, `:snow`, `:rain` or `:graupel`, as the
reference realized it.

`<species>_rate` is DALES's `sed_q/ρ` — a mixing ratio times a fall speed, whatever
the archive's `kg/m2` label says — so dividing out the mixing ratio leaves the
speed. `NaN` where there is too little of the species for the ratio to mean
anything.
"""
function dales_fall_speed(
    species::Symbol,
    date;
    root = data_root(),
    q_min::Real = 1.0e-8,
)
    raw_rate = Dict(
        :ice => "ice_rate", :snow => "snow_rate",
        :rain => "rain_rate", :graupel => "graupel_rate",
    )
    raw_mass = Dict(
        :ice => "sv008", :snow => "sv010", :rain => "sv002", :graupel => "sv012",
    )
    haskey(raw_rate, species) ||
        error("No fall speed for `$species`; try :ice, :snow, :rain or :graupel")
    rate = dales_field(raw_rate[species], date; root, translate_units = false)
    q = dales_field(raw_mass[species], date; root, translate_units = false).data
    return (;
        rate.z, rate.time,
        data = ifelse.(q .> q_min, rate.data ./ q, NaN),
        units = "m s^-1",
    )
end

"""
    surface_heat_fluxes(date; root)

`(; time, hfss, hfls)` [W m^-2], upward positive, from the reference's own surface
kinematic fluxes.

`wthls` and `wqts` are the subfilter θ_l and total-water fluxes at the lowest
face, so the energy fluxes are `ρ c_p wθ_l` and `ρ L_v wq_t`. DALES computes these
itself at `isurf = 2`; `scm_in`'s `sfc_*_flx` are netCDF fill and are not them
(`docs/design.md` section 10).
"""
function surface_heat_fluxes(date; root = data_root())
    ρ = dales_field("rhof", date; root, translate_units = false).data[1, :]
    wθ = dales_field("wthls", date; root, translate_units = false)
    wq = dales_field("wqts", date; root, translate_units = false).data[1, :]
    return (;
        wθ.time,
        hfss = ρ .* DALES_CONSTANTS.cp_d .* wθ.data[1, :],
        hfls = ρ .* DALES_CONSTANTS.L_v .* wq,
    )
end

"""
    column_water_path(q, ρ, faces)

`∫ ρ q dz` [kg m^-2] for each time.

A water path is a *derived* quantity on both sides of a comparison, never a stored
one: the reference's own bars were integrated by DALES over its 286 levels with
its `ρ` and `Δz`, a model's `lwp` by ClimaAtmos over the model's levels, so
differencing them mixes a physical difference with an integration one. Build both
with this, from profiles on one grid (`docs/design.md` section 12).
"""
function column_water_path(q::AbstractMatrix, ρ::AbstractMatrix, faces::AbstractVector)
    nz = size(q, 1)
    length(faces) == nz + 1 || error(
        "A column of $nz cells needs $(nz + 1) faces, got $(length(faces)).",
    )
    dz = diff(faces)
    return [sum(q[k, t] * ρ[k, t] * dz[k] for k in 1:nz) for t in axes(q, 2)]
end
