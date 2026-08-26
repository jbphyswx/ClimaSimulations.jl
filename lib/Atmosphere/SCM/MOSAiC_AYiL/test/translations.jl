using Statistics: mean
using Test: Test
using MOSAiC_AYiL: MOSAiC_AYiL as MA

# The archive is static published data, so these are exact statements about it rather
# than tolerances chosen to pass. 20200503 carries liquid, rain and a little ice.

const DATE = "20200503"

Test.@testset "the namer covers the archive" begin
    named = Dict{String, Vector{String}}()
    errored = String[]
    for path in (
        MA.les_profiles_path(DATE),
        joinpath(MA.day_dir(DATE), "mphysprofiles.001.nc"),
        joinpath(MA.day_dir(DATE), "samptend.001.nc"),
    )
        MA.NC.NCDataset(path) do ds
            for k in keys(ds)
                k in ("time", "zt", "zm") && continue
                name = try
                    MA.dales_description(k)
                catch
                    push!(errored, k)
                    continue
                end
                push!(get!(named, name, String[]), k)
            end
        end
    end

    # a name the archive contains must never throw, whatever the rules make of it
    Test.@test isempty(errored)
    Test.@test !any(isempty, keys(named))

    # Two names colliding means one description stands for two different quantities.
    # The archive genuinely reports these eight twice — `bulkmicrostat3` recomputes
    # the rain scalars and their tendencies that SB3 already carries — so the
    # collisions are the data's, not the rules'. Any other pair is a rule that has
    # become ambiguous.
    collisions = Dict(n => sort(v) for (n, v) in named if length(v) > 1)
    Test.@test sort(collect(keys(collisions))) == [
        "n_rain",
        "n_rain_tendency_evaporation",
        "n_rain_tendency_microphysics",
        "n_rain_tendency_sedimentation",
        "q_rain",
        "q_rain_tendency_evaporation",
        "q_rain_tendency_microphysics",
        "q_rain_tendency_sedimentation",
    ]
    Test.@test collisions["q_rain"] == ["qrmn", "sv002"]
    # and the two really are the same quantity, so the collision is benign. `qrmn` is
    # a BULKMICROSTAT3 name, whose leading fill record is dropped on read, so it is
    # the later times of `sv002` that it lines up with
    a = MA.dales_field("qrmn", DATE; translate_units = false)
    b = MA.dales_field("sv002", DATE; translate_units = false)
    Test.@test a.time == b.time[2:end]
    Test.@test maximum(abs.(a.data .- b.data[:, 2:end])) / maximum(abs.(b.data)) < 0.01
end

Test.@testset "the naming rules" begin
    Test.@test MA.dales_description("sv005") == "q_cloud_liquid"
    Test.@test MA.dales_description("sv006") == "n_ccn"
    Test.@test MA.dales_description("svp008") == "q_cloud_ice_tendency"
    Test.@test MA.dales_description("wsv002r") == "q_rain_flux_resolved"
    Test.@test MA.dales_description("dq_i_dep") == "q_ice_tendency_deposition"
    Test.@test MA.dales_description("utendcor") ==
               "u_tendency_coriolis_unknown_sample"
    Test.@test MA.dales_description("thltendmicroall") ==
               "thl_tendency_microphysics_all"
    # a name of a known family carrying an unknown species or process is an error,
    # not a passthrough
    Test.@test_throws ErrorException MA.dales_description("dq_zz_dep")
    Test.@test_throws ErrorException MA.dales_description("thltendmicro")
end

Test.@testset "moisture is specific, and no w->q conversion belongs anywhere" begin
    # The archive's labels disagree — `profiles` says "specific humidity" for `qt`
    # and "specific mixing ratio" for the `sv` scalars, `scm_in` says "mixing ratio"
    # — so the data decides (`docs/design.md` section 4). Both checks must prefer the
    # as-stored reading; a `w -> q` conversion added anywhere would flip them.

    # 1. DALES consumed `scm_in` unconverted, so the IC and the nudge targets do too
    fd = MA.read_scm_in(DATE)
    zt, _ = MA.les_density(DATE)
    total = Float64.(fd.q .+ fd.ql .+ fd.qi)
    on_les(v) = MA.Intp.extrapolate(
        MA.Intp.interpolate((Float64.(fd.z),), v, MA.Intp.Gridded(MA.Intp.Linear())),
        MA.Intp.Flat(),
    )(
        zt,
    )
    les = MA.climaatmos_field("hus", DATE).data[:, 1]
    k = findall(<(2000), zt)
    err(v) = sum(abs.(on_les(v)[k] .- les[k])) / length(k)
    Test.@test err(total) < err(total ./ (1 .+ total))

    # 2. the `sv` scalars pair with the moist slab density, not a dry one
    faces = MA.native_faces(MA.case(DATE))
    ρ = MA.dales_field("rhof", DATE; translate_units = false).data
    qt = MA.dales_field("qt", DATE; translate_units = false).data
    q = MA.dales_field("sv005", DATE; translate_units = false).data
    bar = MA.NC.NCDataset(joinpath(MA.day_dir(DATE), "tmser.001.nc")) do ds
        Float64.(vec(Array(ds["clwp_bar"])))
    end
    windowed = [mean(bar[(5j - 4):(5j)]) for j in axes(q, 2)]
    moist = maximum(abs.(MA.column_water_path(q, ρ, faces) .- windowed))
    dry = maximum(abs.(MA.column_water_path(q, ρ .* (1 .- qt), faces) .- windowed))
    Test.@test moist < dry / 5
end

Test.@testset "numbers leave DALES's per-mass convention" begin
    per_mass = MA.dales_field("sv007", DATE; translate_units = false)
    per_volume = MA.dales_field("sv007", DATE)
    ρ = MA.dales_field("rhof", DATE; translate_units = false).data
    Test.@test per_mass.units == "kg/kg"
    Test.@test per_volume.units == "m^-3"
    Test.@test per_volume.data ≈ per_mass.data .* ρ
    # a mass is left alone
    Test.@test MA.dales_field("sv005", DATE).units == "kg/kg"
end

Test.@testset "the ClimaAtmos bridge" begin
    names = MA.climaatmos_translated_names()
    Test.@test length(names) == length(unique(names))
    for n in names
        v = MA.climaatmos_field(n, DATE)
        Test.@test all(isfinite, v.data)
        Test.@test !isempty(v.units)
    end

    # `hus` is total water over every phase, so it is at least `qt`, which is only
    # vapour plus cloud liquid
    hus = MA.climaatmos_field("hus", DATE).data
    qt = MA.dales_field("qt", DATE; translate_units = false).data
    Test.@test all(hus .>= qt .- 1.0e-12)

    # cloud fraction becomes a percentage
    cl = MA.climaatmos_field("cl", DATE)
    Test.@test cl.units == "%"
    Test.@test maximum(cl.data) <= 100
    Test.@test maximum(cl.data) > 1   # a fraction would never exceed 1 here

    # the numbers arrive per volume
    Test.@test MA.climaatmos_field("cdnc", DATE).units == "m^-3"
end

Test.@testset "`ta` inverts θ_l, checked against the gas law" begin
    # nothing in the θ_l inversion touches ρ, so forming the virtual temperature from
    # `presf` and `rhof` is an independent statement about the same T
    ta = MA.climaatmos_field("ta", DATE).data
    p = MA.climaatmos_field("pfull", DATE).data
    ρ = MA.dales_field("rhof", DATE; translate_units = false).data
    qt = MA.dales_field("qt", DATE; translate_units = false).data
    ql = MA.dales_field("ql", DATE; translate_units = false).data
    T_v = p ./ (ρ .* MA.DALES_CONSTANTS.R_d)
    # the two sides share no arithmetic, so what is left is the Float32 storage and
    # the centre-pressure reconstruction, not a difference in the inversion
    Test.@test maximum(abs.(ta .* (1 .+ 0.6078 .* (qt .- ql) .- ql) .- T_v)) < 0.25
end

Test.@testset "water paths match DALES's own bars" begin
    # the bars are 60 s snapshots and the profiles 300 s averages, so each profile
    # sample takes the five bars it covers
    faces = MA.native_faces(MA.case(DATE))
    ρ = MA.dales_field("rhof", DATE; translate_units = false).data
    bars = MA.NC.NCDataset(joinpath(MA.day_dir(DATE), "tmser.001.nc")) do ds
        Dict(
            k => Float64.(vec(Array(ds[k]))) for
            k in ("clwp_bar", "icwp_bar", "rlwp_bar", "siwp_bar", "giwp_bar")
        )
    end
    for (raws, bar_names) in (
        (("sv005",), ("clwp_bar",)),
        (("sv008",), ("icwp_bar",)),
        (("sv002",), ("rlwp_bar",)),
        (("sv010", "sv012"), ("siwp_bar", "giwp_bar")),
    )
        q = sum(MA.dales_field(r, DATE; translate_units = false).data for r in raws)
        mine = MA.column_water_path(q, ρ, faces)
        bar = sum(bars[b] for b in bar_names)
        windowed = [mean(bar[(5j - 4):(5j)]) for j in eachindex(mine)]
        scale = max(maximum(windowed), 1.0e-12)
        Test.@test maximum(abs.(mine .- windowed)) / scale < 0.05
    end
end

Test.@testset "a fall speed is the realized flux, not a mean-state estimate" begin
    v = MA.dales_fall_speed(:ice, DATE)
    finite = filter(isfinite, v.data)
    Test.@test !isempty(finite)
    Test.@test all(>=(0), finite)
    Test.@test v.units == "m s^-1"
    # where there is too little of the species the ratio means nothing and says so
    Test.@test any(isnan, v.data)
end

Test.@testset "surface fluxes are the reference's own" begin
    f = MA.surface_heat_fluxes(DATE)
    Test.@test length(f.hfss) == length(f.time) == length(f.hfls)
    Test.@test all(isfinite, f.hfss) && all(isfinite, f.hfls)
end

Test.@testset "a water path needs one face per cell plus one" begin
    q = fill(1.0e-4, 3, 2)
    ρ = fill(1.2, 3, 2)
    Test.@test MA.column_water_path(q, ρ, [0.0, 10.0, 20.0, 30.0]) ≈
               fill(1.2e-4 * 30, 2)
    Test.@test_throws ErrorException MA.column_water_path(q, ρ, [0.0, 10.0, 20.0])
end
