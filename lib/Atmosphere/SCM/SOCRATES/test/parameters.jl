using ClimaAtmos: ClimaAtmos as CA
using ClimaParams: ClimaParams as CP
using SOCRATES: SOCRATES
using Test: Test
using TOML: TOML

# ClimaParams ignores an unknown override key silently, so every value the suite sets is
# asserted where the model reads it, not where it is written.

const CASE = SOCRATES.case("RF09_Obs")

resolved(; params = nothing) =
    SOCRATES.socrates_params(SOCRATES.socrates_toml_dict(Float64, CASE; params))

Test.@testset "Atlas LES constants" begin
    tp = CA.Parameters.thermodynamics_params(resolved())
    # recovered from the Atlas output; not the ClimaParams defaults 2.5008e6 / 2.8344e6
    Test.@test CA.TD.Parameters.LH_v0(tp) == 2_510_400
    Test.@test CA.TD.Parameters.LH_s0(tp) == 2_844_000
    # fusion is derived as the difference, and is unchanged by the two above
    Test.@test CA.TD.Parameters.LH_f0(tp) == 333_600

    # `density_ice_water` does not reach cloud ice; `cloud_ice_apparent_density` does
    mp = CA.Parameters.microphysics_1m_params(resolved())
    Test.@test mp.cloud.ice.ρᵢ == 500
end

"""Config key => the `TurbulenceConvectionParameters` field the model reads it as."""
const EDMFX_CLOSURES = (
    EDMF_max_surface_area = :max_surface_area,
    EDMF_min_area = :min_area,
    EDMF_max_area = :max_area,
    turb_entr_param_vec = :turb_entr_param_vec,
    entr_inv_tau = :entr_inv_tau,
    entr_coeff = :entr_coeff,
    min_area_limiter_scale = :min_area_limiter_scale,
    min_area_limiter_power = :min_area_limiter_power,
    detr_coeff = :detr_coeff,
    detr_buoy_coeff = :detr_buoy_coeff,
    detr_vertdiv_coeff = :detr_vertdiv_coeff,
    detr_massflux_vertdiv_coeff = :detr_massflux_vertdiv_coeff,
    max_area_limiter_scale = :max_area_limiter_scale,
    max_area_limiter_power = :max_area_limiter_power,
    pressure_normalmode_drag_coeff = :pressure_normalmode_drag_coeff,
)

Test.@testset "EDMFX closures reach the model" begin
    tcp = CA.Parameters.turbconv_params(resolved())
    for (key, field) in pairs(EDMFX_CLOSURES)
        Test.@test haskey(SOCRATES.EDMFX_CLOSURE_PARAMS, key)
        Test.@test getproperty(tcp, field) ==
                   getproperty(SOCRATES.EDMFX_CLOSURE_PARAMS, key)
    end
    # the values that would otherwise silently come from the ClimaParams defaults
    Test.@test tcp.entr_coeff == 0.1
    Test.@test tcp.detr_buoy_coeff == 1
    Test.@test tcp.detr_massflux_vertdiv_coeff == 0.3
    Test.@test tcp.turb_entr_param_vec == [0.001, 1.0e5]
end

Test.@testset "the run needs no file, and the shipped one matches it" begin
    # `SOCRATES_PARAMETERS` is the source of truth and a run reads no TOML, so this is
    # the only thing keeping the generated `configs/socrates.toml` from drifting
    Test.@test TOML.parsefile(SOCRATES.config_path()) == SOCRATES.parameter_overrides()
    Test.@test SOCRATES.parameter_overrides()["latent_heat_vaporization_at_reference"]["value"] ==
               2_510_400
    # write it somewhere else and it round-trips
    mktempdir() do dir
        path = SOCRATES.write_parameter_file(joinpath(dir, "p.toml"))
        Test.@test TOML.parsefile(path) == SOCRATES.parameter_overrides()
    end
end

Test.@testset "terminal velocity scaling factors" begin
    td = SOCRATES.socrates_toml_dict(Float64, CASE)
    scaling = CP.get_parameter_values(
        td,
        SOCRATES.TERMINAL_VELOCITY_SCALING_PARAMS,
        "SOCRATES",
    )
    for species in (:liquid, :ice, :rain, :snow)
        Test.@test getproperty(scaling, species) == 1
    end
end

Test.@testset "droplet number is per case" begin
    md = Dict("prescribed_cloud_droplet_number_concentration" => "n")
    for c in SOCRATES.all_cases()
        td = SOCRATES.socrates_toml_dict(Float64, c)
        Test.@test CP.get_parameter_values(td, md, "SOCRATES").n == SOCRATES.n_ccn(c)
    end
    # pinned against the Atlas per-flight droplet numbers
    Test.@test SOCRATES.n_ccn(SOCRATES.case("RF09_Obs")) == 1.9e8
    Test.@test SOCRATES.n_ccn(SOCRATES.case("RF01_Obs")) == 7.5e7
end

Test.@testset "override layering" begin
    entry(v) = Dict{String, Any}("value" => v, "type" => "float")
    key = "entr_coeff"
    md = Dict(key => "v")
    value(params) =
        CP.get_parameter_values(
            SOCRATES.socrates_toml_dict(Float64, CASE; params),
            md,
            "SOCRATES",
        ).v

    Test.@test value(nothing) == 0.1                       # the suite's own config
    Test.@test value(Dict(key => entry(0.5))) == 0.5       # a dictionary overrides it
    # a vector is applied in order, so the last source wins
    Test.@test value([Dict(key => entry(0.5)), Dict(key => entry(0.7))]) == 0.7

    # a path and a dictionary are interchangeable, and mix
    mktempdir() do dir
        path = joinpath(dir, "p.toml")
        open(io -> TOML.print(io, Dict(key => entry(0.9))), path, "w")
        Test.@test value(path) == 0.9
        Test.@test value([path, Dict(key => entry(0.3))]) == 0.3
        Test.@test value([Dict(key => entry(0.3)), path]) == 0.9
    end

    mktempdir() do dir
        Test.@test_throws ErrorException value(joinpath(dir, "absent.toml"))
        Test.@test_throws ErrorException value(@__FILE__)      # not a .toml
    end
end
