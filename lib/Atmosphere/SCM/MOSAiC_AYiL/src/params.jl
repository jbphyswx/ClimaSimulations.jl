"""
    params.jl

The parameter set for an AYiL case: the reference runs' thermodynamic constants,
the day's CCN number, and any caller overrides.
"""

""" Fast no IO cache for the MOSAiC N_CCNs, native data resolution (todo, figure out how to make this immutable w/o stressing the compiler with a 190 element NamedTuple (plus the keys are numbers). Base.ImmuableDict does O(N) lookups so nogo, maybe Base.PersistentDict but that can fail string lookups at runtime due to `===` comparisons) """
const MOSAiC_AYiL_N_CCNs_default = Dict{String, Float32}( "20191016" => 2.1808302f7, "20191022" => 2.5889902f7, "20191024" => 1.1388299f7, "20191025" => 1.0518072f7, "20191026" => 5.5844325f6, "20191027" => 7.0263755f6, "20191028" => 3.8039692f7, "20191029" => 1.0481482f7, "20191030" => 2.9511376f7, "20191031" => 1.1005444f7, "20191101" => 6.338704f6, "20191102" => 1.2887198f7, "20191103" => 2.915645f7, "20191105" => 9.312496f7, "20191106" => 3.8504636f7, "20191107" => 3.865106f7, "20191108" => 4.585421f7, "20191109" => 1.9233083f9, "20191110" => 3.901741f7, "20191111" => 9.407267f7, "20191112" => 8.254113f7, "20191113" => 8.3525624f7, "20191114" => 7.121436f7, "20191115" => 3.0735566f7, "20191125" => 6.2868936f7, "20191126" => 3.4734508f7, "20191129" => 6.796808f7, "20191130" => 6.796808f7, "20191201" => 6.796808f7, "20191203" => 1.196908f8, "20191204" => 1.3052354f8, "20191205" => 8.1746024f7, "20191206" => 5.610962f7, "20191209" => 2.0445376f8, "20191210" => 1.5398891f8, "20191211" => 1.1855675f8, "20191212" => 1.3028856f8, "20191213" => 4.6548452f7, "20191214" => 6.4998756f7, "20191215" => 1.01733384f8, "20191217" => 6.718056f7, "20191218" => 5.1007964f7, "20191219" => 2.9257304f7, "20191221" => 4.625443f7, "20191222" => 7.022983f7, "20191223" => 3.73152f7, "20191225" => 2.864112f7, "20191227" => 3.506839f7, "20191229" => 6.429136f7, "20191230" => 4.976068f7, "20191231" => 3.9297136f7, "20200101" => 2.0695853f8, "20200102" => 4.0131174f8, "20200103" => 1.4058053f8, "20200105" => 1.1207502f8, "20200106" => 1.3522086f8, "20200109" => 1.5685293f8, "20200110" => 1.5685293f8, "20200111" => 6.668138f7, "20200112" => 5.7951868f7, "20200114" => 8.5749576f7, "20200115" => 1.5993214f8, "20200116" => 2.3534558f8, "20200117" => 2.3534558f8, "20200118" => 7.8987384f7, "20200119" => 4.6855604f7, "20200120" => 1.4341648f8, "20200121" => 1.7496646f8, "20200122" => 1.7672224f8, "20200123" => 1.7672224f8, "20200124" => 1.7672224f8, "20200126" => 5.5122696f7, "20200127" => 1.064669f8, "20200128" => 1.0224261f8, "20200130" => 8.614944f7, "20200131" => 9.12797f7, "20200202" => 1.1428027f8, "20200203" => 1.4039206f8, "20200204" => 1.7365226f8, "20200205" => 1.1284066f8, "20200206" => 1.9205264f8, "20200207" => 2.1335366f8, "20200208" => 2.375227f8, "20200209" => 8.310831f7, "20200210" => 1.131489f8, "20200211" => 5.0745236f7, "20200212" => 5.272561f7, "20200214" => 2.0723115f8, "20200215" => 9.848197f7, "20200216" => 8.530302f7, "20200217" => 8.884924f7, "20200218" => 1.0717726f8, "20200219" => 3.8811923f8, "20200220" => 5.0022336f8, "20200221" => 2.4481613f8, "20200222" => 2.1923862f8, "20200223" => 9.9804264f7, "20200224" => 2.9005715f8, "20200225" => 1.7121126f8, "20200226" => 1.9910726f8, "20200227" => 7.787144f7, "20200228" => 8.403689f7, "20200301" => 8.4471096f7, "20200302" => 9.700246f7, "20200303" => 1.06670376f8, "20200304" => 1.0282744f8, "20200305" => 1.0169768f8, "20200306" => 9.946501f7, "20200307" => 9.630498f7, "20200308" => 9.7667896f7, "20200310" => 5.0640256f8, "20200311" => 1.5998003f8, "20200316" => 1.74925f8, "20200317" => 1.4482454f8, "20200318" => 1.3992355f8, "20200319" => 3.927191f8, "20200321" => 2.2376366f8, "20200323" => 1.1727019f8, "20200324" => 9.642301f7, "20200325" => 8.88796f7, "20200326" => 1.593104f8, "20200327" => 1.6972246f8, "20200329" => 1.4891424f8, "20200330" => 1.3902528f8, "20200402" => 1.8176342f8, "20200403" => 1.7722173f8, "20200404" => 1.7252898f8, "20200405" => 2.4548243f8, "20200406" => 2.4548243f8, "20200407" => 1.375808f8, "20200408" => 1.6784349f8, "20200409" => 1.5517746f8, "20200410" => 1.2898544f8, "20200411" => 1.00136424f8, "20200412" => 1.4903608f8, "20200413" => 1.1635323f8, "20200414" => 1.514105f8, "20200415" => 3.244332f8, "20200416" => 4.9236563f8, "20200417" => 1.5908976f8, "20200418" => 1.5158534f8, "20200419" => 9.179611f7, "20200420" => 1.0756734f8, "20200422" => 1.8942224f8, "20200423" => 1.6711186f8, "20200424" => 1.5597688f8, "20200425" => 1.3677205f8, "20200426" => 1.7712192f8, "20200427" => 1.6728906f8, "20200429" => 1.1647699f8, "20200430" => 1.1647699f8, "20200501" => 2.1472355f8, "20200502" => 1.9990091f8, "20200503" => 1.6630994f8, "20200504" => 1.32955416f8, "20200505" => 1.3233245f8, "20200506" => 1.5006672f8, "20200701" => 1.0942572f7, "20200702" => 6.3689896f7, "20200706" => 1.14227144f8, "20200707" => 2.0938348f7, "20200708" => 2.0938348f7, "20200709" => 1.6263541f8, "20200710" => 1.3222135f8, "20200713" => 2.2485108f7, "20200714" => 5.885514f7, "20200715" => 6.0800136f7, "20200717" => 5.4780244f7, "20200720" => 4.5565536f7, "20200721" => 5.4446172f7, "20200722" => 1.5000542f8, "20200723" => 5.5026172f7, "20200724" => 4.415541f7, "20200725" => 6.495131f7, "20200726" => 5.1083782f8, "20200826" => 1.4043322f8, "20200827" => 1.3228816f8, "20200828" => 1.4409632f8, "20200829" => 9.064879f7, "20200830" => 3.0284972f7, "20200831" => 1.2496875f7, "20200901" => 5.419664f7, "20200902" => 6.923689f6, "20200903" => 8.743732f6, "20200905" => 3.3025138f7, "20200906" => 4.3755224f7, "20200907" => 1.5069349f7, "20200909" => 2.325295f7, "20200910" => 1.7169656f7, "20200911" => 9.736504f6, )


"""Path the shipped copy of [`DALES_THERMODYNAMICS`](@ref) is written to."""
dales_thermodynamics_path() =
    joinpath(dirname(@__DIR__), "configs", "dales_thermodynamics.toml")

"""
Thermodynamic constants of the reference DALES runs, from `modglobal.f90:73-101` and
`modmicrodata3.f90:163-164`, offered as the default layer of [`mosaic_toml_dict`](@ref).

Every DALES constant with a counterpart is here, including those whose values coincide
with today's ClimaParams defaults: an upstream default may change, and a case definition
that moved with it silently would not be a case definition.

The key is the one each consumer actually reads, which is not obvious and not stable
across versions. `Thermodynamics/ext/CreateParametersExt.jl` maps `gas_constant_dry_air`
to `R_d` and `isobaric_specific_heat_dry_air` to `cp_d`, while `molar_mass_dry_air` and
`adiabatic_exponent_dry_air` reach RRTMGP's gas optics instead, so both pairs are set and
kept consistent. Setting only the molar masses would leave the thermodynamics untouched.

Two differences no parameter can remove, and so part of the comparison's error budget:
DALES holds the latent heat of vaporization constant where Thermodynamics uses
`L_v(T) = LH_v0 + (c_pv - c_pl)(T - T_0)`, so matching at `T_0` leaves CliMA 1.2% high at
260 K, 2.1% at 250 K and 4.0% at 230 K; and DALES's saturation vapour pressure is
Tetens/Murray rather than the Clausius-Clapeyron integration Thermodynamics performs.
"""
const DALES_THERMODYNAMICS = (;
    # read by Thermodynamics
    gas_constant_dry_air =
        (value = 287.04, description = "DALES `rd`; ClimaParams default 287.0 (+0.014%)"),
    gas_constant_vapor = (
        value = 461.5,
        description = "DALES `rv`; equals the ClimaParams default, set so it cannot drift",
    ),
    isobaric_specific_heat_dry_air = (
        value = 1004.0,
        description = "DALES `cp`; ClimaParams default 1004.5 (-0.050%)",
    ),
    latent_heat_vaporization_at_reference = (
        value = 2.53e6,
        description = "DALES `rlv`; ClimaParams default 2.5008e6 (+1.168%)",
    ),
    latent_heat_sublimation_at_reference = (
        value = 2.834e6,
        description = "DALES SB3 `rlvi`, the value its microphysics used; `riv = 2.84e6` \
                       in modglobal is declared but never referenced. ClimaParams default \
                       2.8344e6 (-0.014%)",
    ),
    temperature_water_freeze = (
        value = 273.16, description = "DALES `tmelt`; ClimaParams default 273.15",
    ),
    potential_temperature_reference_pressure = (
        value = 1.0e5,
        description = "DALES `pref0`, the reference pressure of its Exner function; \
                       equals the ClimaParams default, set so it cannot drift",
    ),
    gravitational_acceleration = (
        value = 9.81,
        description = "DALES `grav`; equals the ClimaParams default, set so it cannot \
                       drift",
    ),
    # read by RRTMGP's gas optics, kept consistent with R_d and R_v
    molar_mass_dry_air = (
        value = 0.028966206103678928,
        description = "gas_constant / 287.04, so the molar mass matches \
                       `gas_constant_dry_air`",
    ),
    molar_mass_water = (
        value = 0.018016164247020586,
        description = "gas_constant / 461.5, so the molar mass matches \
                       `gas_constant_vapor`",
    ),
    adiabatic_exponent_dry_air = (
        value = 0.28589641434262952,
        description = "287.04 / 1004, matching `gas_constant_dry_air` / \
                       `isobaric_specific_heat_dry_air`",
    ),
    # other DALES constants with a counterpart
    von_karman_constant = (
        value = 0.4,
        description = "DALES `fkar`, used by its Monin-Obukhov surface layer; equals the \
                       ClimaParams default, set so it cannot drift",
    ),
    stefan_boltzmann_constant = (
        value = 5.67e-8,
        description = "DALES `boltz`; equals the ClimaParams default, set so it cannot \
                       drift",
    ),
    density_liquid_water =
        (value = 998.0, description = "DALES `rhow`; ClimaParams default 1000 (-0.2%)"),
)

"""
    parameter_overrides(set = DALES_THERMODYNAMICS)

`set` in the `ClimaParams` override form, carrying each entry's description.
"""
parameter_overrides(set = DALES_THERMODYNAMICS) = Dict{String, Any}(
    String(name) => Dict{String, Any}(
        "value" => entry.value, "type" => "float", "description" => entry.description,
    ) for (name, entry) in pairs(set)
)

"""
    write_parameter_file(path = dales_thermodynamics_path(); set = DALES_THERMODYNAMICS)

Write `set` as a `ClimaParams` override TOML, returning `path`. The shipped
`configs/dales_thermodynamics.toml` is this file; nothing reads it to build a run.
"""
function write_parameter_file(
    path::AbstractString = dales_thermodynamics_path();
    set = DALES_THERMODYNAMICS,
)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(
            io,
            "# Generated by `MOSAiC_AYiL.write_parameter_file`; edit `src/params.jl`.",
        )
        println(io)
        TOML.print(io, parameter_overrides(set))
    end
    return path
end

_override_dict(source::AbstractDict) =
    Dict{String, Any}(string(k) => v for (k, v) in source)
function _override_dict(source::AbstractString)
    isfile(source) || error("Parameter TOML not found: $source")
    endswith(source, ".toml") ||
        error("Parameter source is not a .toml file: $source")
    return TOML.parsefile(source)
end

"""
    mosaic_toml_dict(FT, case; ccn, params)

The `ClimaParams.ParamDict{FT}` for `case`: [`DALES_THERMODYNAMICS`](@ref), then the day's
CCN number, then `params` — a TOML path or a dictionary, as a calibration supplies. No
file is read unless `params` names one, and `params` can override any value beneath it.

`ccn` defaults to the day's entry in `MOSAiC_AYiL_N_CCNs_default` [m⁻³], the
value DALES read from `scm_in n_ccn` through `inittbccn3`. The file's profile is
uniform in height, so one number is the whole of it.
"""
function mosaic_toml_dict(
    ::Type{FT},
    c::MOSAiCAYiLCase;
    ccn::Real = MOSAiC_AYiL_N_CCNs_default[c.date],
    params = nothing,
) where {FT <: AbstractFloat}
    overrides = parameter_overrides()
    overrides["prescribed_cloud_droplet_number_concentration"] =
        Dict{String, Any}("value" => FT(ccn), "type" => "float")
    if !isnothing(params)
        for (name, entry) in _override_dict(params)
            overrides[name] = entry
        end
    end
    return CP.create_toml_dict(FT; override_file = overrides)
end

"""`ClimaAtmosParameters` for `case`."""
mosaic_params(
    ::Type{FT},
    c::MOSAiCAYiLCase;
    microphysics_model = CA.NonEquilibriumMicrophysics1M(),
    kwargs...,
) where {FT <: AbstractFloat} =
    mosaic_params(mosaic_toml_dict(FT, c; kwargs...); microphysics_model)

mosaic_params(
    toml_dict::CP.ParamDict;
    microphysics_model = CA.NonEquilibriumMicrophysics1M(),
) = CA.ClimaAtmosParameters(toml_dict; microphysics_model)
