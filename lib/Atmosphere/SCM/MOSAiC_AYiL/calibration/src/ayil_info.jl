
#=

    Calibration setup defaults

=#

# a copy: `delete!` on the package's own table would take these days out of it
default_calibration_tops = copy(MOSAiC_AYiL.BEST_SIMULATION_TOP_F)




#=
    Calibration defaults make a few adjustments to the case list, 
    including dropping a few who would be valid but who don't neatly fit into the
    chosen canonical grid heights, or are by manual judgement deemed unworthy of calibration

    It is possible that after future tuning (tbd, we still have some cases that need further trimming for example),
    we can whittle down to 3-4 canonical grid heights that capture the full range of plausible cases.


        Days whose filter top is below the lowest canonical rung, so no rung excludes
        their anomalous ice: 665, 2273, 875, 2100 and 1983 m. Dropped rather than run,
        and [`best_simulation_top`] errors on them.
            ("20200429", "20200702", "20200717", "20200721", "20200906")

=#
const default_not_calibrated_additions =
    ("20200429", "20200702", "20200717", "20200721", "20200906")

for date in default_not_calibrated_additions
    delete!(default_calibration_tops, date)
end