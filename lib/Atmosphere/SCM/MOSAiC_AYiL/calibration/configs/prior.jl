
#=

    The parameters an AYiL calibration varies, and their priors

    Configuration: the names are ClimaParams keys, so a sampled value reaches the
    model through the parameter TOML `ClimaCalibrate` writes per member and nothing
    here has to know how the model reads it.

    Bounds are the constraint, not the distribution — `constrained_gaussian` puts the
    named mean and width in physical space and transforms to the unconstrained space
    EKP samples in.

=#

"""
    default_ayil_prior()

The prior over the mixed-phase microphysics parameters an AYiL calibration varies.

These four are what set whether a supercooled liquid layer survives: the two
non-equilibrium relaxation timescales decide how fast vapour is moved onto liquid
versus ice, and the two snow-autoconversion parameters decide how fast cloud ice is
converted to snow and falls out. The reference's SB3 scheme partitions condensate
differently from a one-moment scheme, so these are the terms that have to absorb the
difference (`docs/design.md` §13).

Defaults in ClimaParams are 10 s, 10 s, 1e-6 and 100 s.
"""
default_ayil_prior() = EKP.combine_distributions([
    EKP.constrained_gaussian(
        "condensation_evaporation_timescale", 10.0, 8.0, 0.1, Inf,
    ),
    EKP.constrained_gaussian(
        "sublimation_deposition_timescale", 10000.0, 8.0, 0.1, Inf,
    ),
    EKP.constrained_gaussian(
        "cloud_ice_specific_humidity_autoconversion_threshold",
        1.0e-6, 1.0e-6, 0.0, Inf,
    ),
    EKP.constrained_gaussian(
        "snow_autoconversion_timescale", 100.0, 80.0, 1.0, Inf,
    ),
])