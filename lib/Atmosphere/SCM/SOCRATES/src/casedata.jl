"""
    casedata.jl

Readers for the declarative case data in `cases/` and the parameter overrides in
`configs/`. Values are returned as the file holds them; callers convert to their
own float type.
"""

case_data_path() = joinpath(dirname(@__DIR__), "cases", "socrates.toml")
config_path() = joinpath(dirname(@__DIR__), "configs", "socrates.toml")

const CASE_DATA = TOML.parsefile(case_data_path())

function flight_entry(flight_number::Integer)
    flights = CASE_DATA["flight"]
    i = findfirst(e -> e["number"] == flight_number, flights)
    isnothing(i) && error(
        "No SOCRATES flight data for flight $flight_number; \
         $(case_data_path()) has flights $([e["number"] for e in flights]).",
    )
    return flights[i]
end

"""Prescribed cloud droplet number concentration [m^-3] for a flight."""
droplet_number(flight_number::Integer) = flight_entry(flight_number)["droplet_number"]

"""Run length [s] for a forcing source, `:Obs` or `:ERA5`."""
run_duration(label::Symbol) = CASE_DATA["run_duration"][String(label)]

"""Relaxation timescale [s] of temperature and total water."""
scalar_nudge_timescale() = CASE_DATA["nudging"]["scalar_timescale"]

"""Relaxation timescale [s] of horizontal wind, by forcing source."""
wind_nudge_timescale(label::Symbol) = CASE_DATA["nudging"]["wind_timescale_$label"]

"""Aerodynamic roughness length [m] of the ocean surface."""
surface_roughness() = CASE_DATA["surface_roughness"]["z0"]
