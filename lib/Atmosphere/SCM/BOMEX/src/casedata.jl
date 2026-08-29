"""Path to the shipped per-case table."""
case_data_path() = joinpath(dirname(@__DIR__), "cases", "bomex.toml")

const CASE_DATA = TOML.parsefile(case_data_path())

"""Run length [s] of the case."""
t_end() = CASE_DATA["run_duration"]["t_end"]

"""Surface pressure [Pa] the hydrostatic profile is integrated from."""
reference_pressure() = CASE_DATA["reference_pressure"]["p_0"]

"""Surface temperature [K]."""
surface_temperature() = CASE_DATA["surface"]["temperature"]

"""Aerodynamic roughness length [m]."""
surface_roughness() = CASE_DATA["surface"]["z0"]

"""Prescribed surface `θ` flux [K m s^-1]."""
surface_theta_flux() = CASE_DATA["surface"]["theta_flux"]

"""Prescribed surface specific-humidity flux [m s^-1]."""
surface_q_flux() = CASE_DATA["surface"]["q_flux"]

"""Prescribed friction velocity [m s^-1]."""
surface_ustar() = CASE_DATA["surface"]["ustar"]

"""Surface pressure [Pa] imposed on the boundary state."""
surface_override_pressure() = CASE_DATA["surface_overrides"]["p"]

"""Surface vapour specific humidity [kg kg^-1] imposed on the boundary state."""
surface_override_q_vap() = CASE_DATA["surface_overrides"]["q_vap"]

"""Coriolis parameter [s^-1]."""
coriolis_parameter() = CASE_DATA["coriolis"]["f"]

"""Domain top [m]."""
z_max() = CASE_DATA["grid"]["z_max"]

"""Number of cells the shipped grid uses."""
z_elem() = CASE_DATA["grid"]["z_elem"]
