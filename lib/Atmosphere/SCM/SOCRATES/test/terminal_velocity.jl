using Test: Test
using ClimaAtmos: ClimaAtmos as CA
using SOCRATES: SOCRATES

const FT = Float64
const MODEL = CA.NonEquilibriumMicrophysics1M()
const SPECIES = (
    CA.MatrixFields.FieldName{(:q_lcl,)}(),
    CA.MatrixFields.FieldName{(:q_icl,)}(),
    CA.MatrixFields.FieldName{(:q_rai,)}(),
    CA.MatrixFields.FieldName{(:q_sno,)}(),
)

Test.@testset "the grid-scale velocity does not reapply the scale" begin
    # `gs_terminal_velocity` forms `ρwχ/ρχ` from subdomain fluxes the scale has already acted on; scaling again
    # squares it, which no test of the scaling itself would catch
    params = CA.ClimaAtmosParameters(FT)
    ρχ = FT(1.0e-4)
    ρwχ = FT(5.0e-6)          # a sedimentation flux, so the ratio is 0.05 m/s
    args = (params, ρwχ, ρχ)
    for name in SPECIES, scale in (FT(0.25), FT(2.5))
        Test.@test CA.gs_terminal_velocity(
            MODEL, SOCRATES.ScaledTerminalVelocity{FT}(scale), name, args...,
        ) === CA.gs_terminal_velocity(
            MODEL, CA.DiagnosticTerminalVelocity(), name, args...,
        )
    end
end

Test.@testset "all four factors exist and default to 1" begin
    # These four are SOCRATES's own parameters — ClimaParams has no such entry — so they arrive only through
    # `configs/socrates.toml`. An unscaled column must be the stock model, so a missing or non-unit default would
    # perturb every run.
    modes = SOCRATES.terminal_velocity_modes(
        FT, SOCRATES.socrates_toml_dict(FT, SOCRATES.case("RF09_Obs")),
    )
    for species in values(SOCRATES.TERMINAL_VELOCITY_SCALING_PARAMS)
        Test.@test getproperty(modes, Symbol(:terminal_velocity_, species)).scale == 1
    end
end

Test.@testset "a scaling factor in the TOML reaches the mode" begin
    # ClimaParams ignores an override it does not recognise, so the only way to know these four arrived is to read
    # them back off the constructed modes
    wanted = Dict(
        :cloud_liquid_terminal_velocity_scaling_factor => 0.4,
        :cloud_ice_terminal_velocity_scaling_factor => 0.5,
        :rain_terminal_velocity_scaling_factor => 0.6,
        :snow_terminal_velocity_scaling_factor => 0.7,
    )
    overrides = Dict{String, Any}(
        String(k) => Dict{String, Any}("value" => v, "type" => "float") for
        (k, v) in wanted
    )
    modes = SOCRATES.terminal_velocity_modes(
        FT,
        SOCRATES.socrates_toml_dict(FT, SOCRATES.case("RF09_Obs"); params = overrides),
    )
    for (param, species) in pairs(SOCRATES.TERMINAL_VELOCITY_SCALING_PARAMS)
        mode = getproperty(modes, Symbol(:terminal_velocity_, species))
        Test.@test mode isa SOCRATES.ScaledTerminalVelocity{FT}
        Test.@test mode.scale == FT(wanted[param])
    end
end
