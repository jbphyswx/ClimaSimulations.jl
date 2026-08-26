"""
    data.jl

Access to the MOSAiC AYiL day directories, and the readers for the two files a
case needs from each: the ERA5 testbed forcing `scm_in.*.nc` and the DALES output
`profiles.001.nc`.
"""

const ARTIFACT_NAME = "ayil_config_input_results"

"""Environment variable that overrides the artifact with a local directory."""
const DATA_ROOT_ENV = "MOSAIC_AYIL_DATA_ROOT"

"""
    data_root()

Directory holding one subdirectory per AYiL day.

Resolves `\$$DATA_ROOT_ENV` if set, otherwise the lazily-downloaded
`$ARTIFACT_NAME` artifact (Zenodo 10.5281/zenodo.10491362, ~911 MB).
"""
function data_root()
    if haskey(ENV, DATA_ROOT_ENV)
        root = ENV[DATA_ROOT_ENV]
        isdir(root) || error("$DATA_ROOT_ENV points at $root, which is not a directory")
        return root
    end
    toml = joinpath(dirname(@__DIR__), "Artifacts.toml")
    hash = Artifacts.artifact_hash(ARTIFACT_NAME, toml)
    isnothing(hash) && error(
        "The `$ARTIFACT_NAME` artifact is not bound in $toml. Either run \
         `gen/build_data_artifact.jl` to create it from Zenodo, or set \
         `ENV[\"$DATA_ROOT_ENV\"]` to a directory of AYiL day subdirectories.",
    )
    return Artifacts.artifact_path(hash)
end

"""
    available_dates(; root = data_root())

The AYiL days present under `root`, as `yyyymmdd` strings, ascending.

Read from the directory rather than tabulated, so the list cannot drift from the
data.
"""
function available_dates(; root = data_root())
    dates = filter(readdir(root)) do name
        length(name) == 8 && all(isdigit, name) && isdir(joinpath(root, name))
    end
    isempty(dates) && error(
        "No AYiL day directories (yyyymmdd) under $root.",
    )
    return sort!(dates)
end

"""
    day_dir(date; root = data_root())

The directory of one AYiL day, given `date` as `yyyymmdd` or a `Date`.
"""
function day_dir(date::AbstractString; root = data_root())
    dir = joinpath(root, date)
    isdir(dir) || error(
        "No AYiL day $date under $root; it has $(length(available_dates(; root))) \
         days from $(first(available_dates(; root))) to \
         $(last(available_dates(; root))).",
    )
    return dir
end

day_dir(date::Dates.Date; kwargs...) =
    day_dir(Dates.format(date, "yyyymmdd"); kwargs...)

scm_in_path(date; kwargs...) =
    joinpath(day_dir(date; kwargs...), "scm_in.a_year_in_les.$date.nc")

les_profiles_path(date; kwargs...) =
    joinpath(day_dir(date; kwargs...), "profiles.001.nc")

namoptions_path(date; kwargs...) =
    joinpath(day_dir(date; kwargs...), "namoptions")

# --- ERA5 testbed forcing --------------------------------------------------- #

# `sfc_sens_flx` and `sfc_lat_flx` are absent from every testbed file, and the
# files declare no `_FillValue`, so NCDatasets hands back the netCDF default fill
# rather than `missing`. That is what turns an absent surface flux into 1e37 W/m².
const NC_FILL_FLOAT = 9.9692099683868690f36
const NC_FILL_DOUBLE = 9.9692099683868690e36

_is_fill(::Missing) = true
_is_fill(v::Float32) = v === NC_FILL_FLOAT
_is_fill(v::Float64) = v === NC_FILL_DOUBLE
_is_fill(::Any) = false

"""
    read_scm_in(date; root, time_index)

The ERA5 testbed forcing for one AYiL day, on an ascending height axis.

`scm_in` stores its levels top-down from about 85 km — 3037 to 3042 of them,
depending on the day — so every profile is reversed here; ClimaAtmos wants
ascending `z` throughout.

`time_index` selects which time record to read. Each file is a single 05:00–11:00
UTC composite average written twice, and the two records are bitwise identical on
every shipped day, so the choice is immaterial for the archive; the keyword exists
so a file that is *not* constant can still be read deliberately rather than by
accident.

# Returns

A `NamedTuple` with the height axis `z` [m], the state profiles `ta`, `hus`, `q`,
`ql`, `qi`, `ua`, `va`, `p`, `o3` and `n_ccn`, the tendencies `tntha`, `tnhusha`,
`tnua`, `tnva`, the subsidence `wa` [m/s, positive up], the geostrophic wind
`ug`/`vg`, and the scalars in `surface`.

`hus` is total water, `q + ql + qi`, matching the quantity ClimaAtmos relaxes
`ρq_tot` toward; DALES instead relaxed `q + ql` and left its ice unnudged, a
deviation `docs/design.md` §4 quantifies. `tnhusha` likewise sums all three
advective tendencies where DALES applies the ice part to a separate scalar — moot
in the archive, where all four are identically zero.

The `_local` variants are used: they are at the trajectory point rather than
averaged over the ERA5 box.
"""
function read_scm_in(date; root = data_root(), time_index::Int = 1)
    path = scm_in_path(date; root)
    isfile(path) || error("No testbed forcing at $path")
    return NC.NCDataset(path, "r") do ds
        # Native precision throughout: the files are Float32, and the caller
        # converts once to the float type its simulation runs in.
        prof(name) = reverse(Array(ds[name])[:, time_index])
        scalar(name) = Array(ds[name])[time_index]
        function optional(name)
            v = scalar(name)
            return _is_fill(v) ? missing : v
        end

        z = reverse(Array(ds["height_f"])[:, time_index])

        ta = prof("t_local")
        q = prof("q_local")
        ql = prof("ql_local")
        qi = prof("qi_local")
        hus = q .+ ql .+ qi
        p = prof("pressure_f")

        # omega is a pressure velocity [Pa/s], positive downward; the geometric
        # velocity is w = -omega / (rho g), positive up, which is ClimaAtmos's
        # `wa` convention. Reconstructed exactly as `modtestbed.f90:739-741` does
        # it — DALES's constants, its 0.61 literal, and vapour alone in `Tv`.
        FT = eltype(ta)
        grav = FT(DALES_CONSTANTS.grav)
        R_d = FT(DALES_CONSTANTS.R_d)
        Tv = ta .* (1 .+ FT(0.61) .* q)
        wa = .-prof("omega") .* (R_d .* Tv) ./ (p .* grav)

        surface = (;
            ps = scalar("ps"),
            # The camp's position on this day. `namoptions` carries a single
            # nominal `xlat`/`xlon` for all 190 days, which does not track the
            # drift; these do.
            trajectory_latitude = scalar("lat"),
            trajectory_longitude = mod(scalar("lon") + 180, 360) - 180,
            albedo = scalar("albedo"),
            albedo_snow = scalar("albedo_snow"),
            snow = scalar("snow"),
            z0_momentum = scalar("mom_rough"),
            z0_heat = scalar("heat_rough"),
            sea_ice_fraction = scalar("sea_ice_frct"),
            t_skin = scalar("t_skin"),
            t_skin_ocean = scalar("t_skin_ocean"),
            t_skin_seaice = scalar("t_skin_seaice"),
            open_sst = scalar("open_sst"),
            land_sea_mask = scalar("lsm"),
            # Absent from the testbed files: DALES reads them but does not use
            # them as the boundary flux at `isurf = 2`. `missing`, never a number,
            # so a consumer cannot mistake the fill value for a flux.
            sensible_heat_flux = optional("sfc_sens_flx"),
            latent_heat_flux = optional("sfc_lat_flx"),
        )

        (;
            z,
            ta,
            hus,
            q,
            ql,
            qi,
            ua = prof("u_local"),
            va = prof("v_local"),
            p,
            o3 = prof("o3"),
            n_ccn = prof("n_ccn"),
            wa,
            tntha = prof("tadv"),
            tnhusha = prof("qadv") .+ prof("ladv") .+ prof("iadv"),
            tnua = prof("uadv"),
            tnva = prof("vadv"),
            ug = prof("ug"),
            vg = prof("vg"),
            surface,
        )
    end
end

"""
    surface_temperature(forcing)

The skin temperature [K] of a blended sea-ice and open-ocean surface,
`(1 - f) T_ocean + f T_seaice` with `f` the sea-ice fraction, as DALES does at
`isurf = 2` with `l_surficefrac`.
"""
surface_temperature(forcing) =
    (1 - forcing.surface.sea_ice_fraction) * forcing.surface.t_skin_ocean +
    forcing.surface.sea_ice_fraction * forcing.surface.t_skin_seaice

# --- DALES output ----------------------------------------------------------- #

"""
    les_density(date; root = data_root(), time_index = 1)

`(z, ρ)`: the DALES cell-centre heights [m] and the slab-mean air density
[kg/m³] the reference run had, from `profiles.001.nc` `rhof`.

This is the density of the simulation being reproduced: DALES integrates its own
pressure from `ps` at startup and forms `rhof = presf/(R_d θ_v Π)`
(`modstartup.f90:750`). `rhof` evolves with the state, and the output starts at
t = 300 s, so `time_index = 1` is the earliest sample rather than t = 0 exactly.

Not `rhobf`/`rhobh`: those are the anelastic base state, which differs from the
air density by 4–9%.
"""
function les_density(date; root = data_root(), time_index::Int = 1)
    path = les_profiles_path(date; root)
    isfile(path) || error("No DALES output at $path")
    return NC.NCDataset(path, "r") do ds
        return (vec(Array(ds["zt"])), Array(ds["rhof"])[:, time_index])
    end
end




"""
The 287 cell faces [m] of the DALES grid: uniform 10 m to 1220 m, geometric stretch,
then uniform ~185 m above 7139 m.

Written out because it is the same grid on all 190 days **[verified: read from every
day's `profiles.001.nc` and compared]**, and it is what every domain top and every
model grid is resolved against — there is nothing per-day to look up. These are the
file's own `Float32` values, so the literal reproduces it exactly rather than nearly.

`zm` holds the *lower* faces, so the top one is `zh(287) = 2 zt(286) - zm(286)`. Not
`zt[end] + 185/2`: the face recursion is unstable and leaves the thicknesses
alternating 184.147 / 185.853, so the last cell is 184.149 m thick and the nominal
half-spacing misses by 0.43 m (`docs/design.md` §2).
"""
const LES_FACES = Float32[
    0.0f0, 10.0f0, 20.0f0, 30.0f0, 40.0f0, 50.0f0, 60.0f0, 70.0f0,
    80.0f0, 90.0f0, 100.0f0, 110.0f0, 120.0f0, 130.0f0, 140.0f0, 150.0f0,
    160.0f0, 170.0f0, 180.0f0, 190.0f0, 200.0f0, 210.0f0, 220.0f0, 230.0f0,
    240.0f0, 250.0f0, 260.0f0, 270.0f0, 280.0f0, 290.0f0, 300.0f0, 310.0f0,
    320.0f0, 330.0f0, 340.0f0, 350.0f0, 360.0f0, 370.0f0, 380.0f0, 390.0f0,
    400.0f0, 410.0f0, 420.0f0, 430.0f0, 440.0f0, 450.0f0, 460.0f0, 470.0f0,
    480.0f0, 490.0f0, 500.0f0, 510.0f0, 520.0f0, 530.0f0, 540.0f0, 550.0f0,
    560.0f0, 570.0f0, 580.0f0, 590.0f0, 600.0f0, 610.0f0, 620.0f0, 630.0f0,
    640.0f0, 650.0f0, 660.0f0, 670.0f0, 680.0f0, 690.0f0, 700.0f0, 710.0f0,
    720.0f0, 730.0f0, 740.0f0, 750.0f0, 760.0f0, 770.0f0, 780.0f0, 790.0f0,
    800.0f0, 810.0f0, 820.0f0, 830.0f0, 840.0f0, 850.0f0, 860.0f0, 870.0f0,
    880.0f0, 890.0f0, 900.0f0, 910.0f0, 920.0f0, 930.0f0, 940.0f0, 950.0f0,
    960.0f0, 970.0f0, 980.0f0, 990.0f0, 1000.0f0, 1010.0f0, 1020.0f0, 1030.0f0,
    1040.0f0, 1050.0f0, 1060.0f0, 1070.0f0, 1080.0f0, 1090.0f0, 1100.0f0, 1110.0f0,
    1120.0f0, 1130.0f0, 1140.0f0, 1150.0f0, 1160.0f0, 1170.0f0, 1180.0f0, 1190.0f0,
    1200.0f0, 1210.0f0, 1220.0f0, 1230.1648f0, 1240.3326f0, 1250.6676f0, 1261.0084f0, 1271.5194f0,
    1282.039f0, 1292.732f0, 1303.4366f0, 1314.3176f0, 1325.2136f0, 1336.2894f0, 1347.3834f0, 1358.6608f0,
    1369.96f0, 1381.4462f0, 1392.9583f0, 1404.6608f0, 1416.3932f0, 1428.3206f0, 1440.2819f0, 1452.4423f0,
    1464.6407f0, 1477.0436f0, 1489.489f0, 1502.1432f0, 1514.8455f0, 1527.7614f0, 1540.7306f0, 1553.9188f0,
    1567.1663f0, 1580.6382f0, 1594.1752f0, 1607.943f0, 1621.782f0, 1635.8586f0, 1650.0128f0, 1664.4116f0,
    1678.8951f0, 1693.6309f0, 1708.4587f0, 1723.5466f0, 1738.7347f0, 1754.1917f0, 1769.7572f0, 1785.6006f0,
    1801.5618f0, 1817.811f0, 1834.1874f0, 1850.8624f0, 1867.6757f0, 1884.798f0, 1902.0706f0, 1919.6642f0,
    1937.4204f0, 1955.5106f0, 1973.7769f0, 1992.391f0, 2011.1958f0, 2030.3634f0, 2049.7373f0, 2069.4907f0,
    2089.4675f0, 2109.8413f0, 2130.4575f0, 2151.4895f0, 2172.784f0, 2194.5159f0, 2216.532f0, 2239.0088f0,
    2261.7937f0, 2285.0647f0, 2308.6702f0, 2332.7896f0, 2357.2727f0, 2382.2998f0, 2407.7227f0, 2433.7236f0,
    2460.1553f0, 2487.2026f0, 2514.7195f0, 2542.8936f0, 2571.5808f0, 2600.9707f0, 2630.9229f0, 2661.6287f0,
    2692.9504f0, 2725.0835f0, 2757.8933f0, 2791.5786f0, 2826.0088f0, 2861.3867f0, 2897.5864f0, 2934.816f0,
    2972.9534f0, 3012.2136f0, 3052.4802f0, 3093.9746f0, 3136.5872f0, 3180.5476f0, 3225.754f0, 3272.444f0,
    3320.5269f0, 3370.2495f0, 3421.5315f0, 3474.6328f0, 3529.485f0, 3586.3616f0, 3645.2085f0, 3706.315f0,
    3769.6426f0, 3835.498f0, 3903.8604f0, 3975.0554f0, 4049.08f0, 4126.28f0, 4206.671f0, 4290.615f0,
    4378.146f0, 4469.64f0, 4565.141f0, 4665.0317f0, 4769.3564f0, 4878.4907f0, 4992.4624f0, 5111.62f0,
    5235.9517f0, 5365.7534f0, 5500.947f0, 5641.746f0, 5787.9805f0, 5939.7627f0, 6096.813f0, 6259.134f0,
    6426.3413f0, 6598.34f0, 6774.6636f0, 6955.1543f0, 7139.301f0, 7325.1543f0, 7509.301f0, 7695.1543f0,
    7879.301f0, 8065.1543f0, 8249.301f0, 8435.154f0, 8619.301f0, 8805.154f0, 8989.301f0, 9175.154f0,
    9359.301f0, 9545.154f0, 9729.301f0, 9915.154f0, 10099.301f0, 10285.154f0, 10469.301f0, 10655.154f0,
    10839.301f0, 11025.154f0, 11209.301f0, 11395.154f0, 11579.301f0, 11765.154f0, 11949.301f0,
]

"""Top face [m] of the DALES grid, the last of [`LES_FACES`](@ref)."""
const LES_TOP_FACE = last(LES_FACES)

"""
    MOSAiC_AYiL_face_above_center(z)

The [`LES_FACES`](@ref) face at or above `z` [m].

    The first face above a center in a grid, defaulting tothe default MOSAiC AYiL grid
"""
function MOSAiC_AYiL_face_above_center(z, faces::AbstractVector = LES_FACES)
    k = findfirst(>=(z), faces)
    isnothing(k) &&
        error("$z m is above the faces top face $(last(faces)) m.")
    return faces[k]
end

"""
    les_faces(date; root = data_root())

The cell faces [m] of the DALES grid for one AYiL day.

`profiles.001.nc` stores `zt` cell centres and `zm` cell *lower* faces, so the
face set is `zm` plus the top face implied by the last centre.
"""
function les_faces(date; LES_FACES = LES_FACES, root = data_root())
    if !isnothing(LES_FACES)
        return LES_FACES
    end
    path = les_profiles_path(date; root)
    isfile(path) || error("No DALES output at $path")
    return NC.NCDataset(path, "r") do ds
        zt = vec(Array(ds["zt"]))
        zm = vec(Array(ds["zm"]))
        # the top face sits as far above the last centre as that centre is above
        # the last stored face
        return vcat(zm, 2 * last(zt) - last(zm))
    end
end
