"""
    atlas_registry.jl

One entry per variable of the Atlas LES archive, each carrying its complete recipe.
Conversions come from [`atlas_archive.jl`](@ref).

Each spec carries

  - `kind`        `:raw` (read from the file) or `:derived` (computed from other entries)
  - `inputs`      entries the transform consumes, in order, already converted
  - `raw_inputs`  raw variables the transform consumes unconverted (derived entries only)
  - `transform`   `:raw` → `f(x, inputs...)`; `:derived` → `f(raw_inputs..., inputs...)`
  - `units`       units of the converted variable; never empty
  - `raw_units`   the `units` attribute expected on the file, asserted before converting
  - `quantity`    what the variable is, independent of how the file spells it
  - `provenance`  which authority fixed the classification:
      * `:units_attr`         the `units` attribute, corroborated by magnitude
      * `:magnitude_verified` the attribute is wrong; magnitude against a sibling settles it
      * `:long_name`          the attribute is ambiguous; `long_name` settles it
      * `:measured`           recovered numerically from the Atlas output itself
      * `:derived`            computed here rather than read

[`assert_specs_cover`](@ref) refuses a file holding a variable with no spec, so nothing is
ever copied through with unknown units.
"""

# Atlas reports `?/kg/day` for three unrelated kinds of quantity, distinguished only by
# `long_name`: mass mixing-ratio rates, number rates, and inverse lengths.
const ATLAS_MASS_PROCESS_RATES = (
    :EPRD, :EPRDG, :EPRDS, :EVPMG, :EVPMS, :MNUCCC, :MNUCCD, :MNUCCI, :MNUCCR, :PCC,
    :PCCN, :PGMLT, :PGRACS, :PGSACW, :PIACR, :PIACRS, :PITOSN, :PRA, :PRACG, :PRACI,
    :PRACIS, :PRACS, :PRAI, :PRC, :PRCI, :PRD, :PRDG, :PRDS, :PRE, :PSACR, :PSACWG,
    :PSACWI, :PSACWS, :PSMLT, :QHOMOC, :QHOMOR, :QMELTI, :QMULTG, :QMULTR, :QMULTRG,
    :QMULTS,
)

const ATLAS_NUMBER_PROCESS_RATES = (
    :NGMLTG, :NGMLTR, :NGRACS, :NHOMOC, :NHOMOR, :NIACR, :NIACRS, :NMELTI, :NMULTG,
    :NMULTR, :NMULTRG, :NMULTS, :NNUCCC, :NNUCCD, :NNUCCI, :NNUCCR, :NPRA, :NPRACG,
    :NPRACS, :NPRAI, :NPRC, :NPRC1, :NPRCI, :NPSACWG, :NPSACWI, :NPSACWS, :NRAGG,
    :NSAGG, :NSCNG, :NSMLTR, :NSMLTS, :NSUBC, :NSUBG, :NSUBI, :NSUBR, :NSUBS,
)

const ATLAS_SLOPE_PARAMETERS = (:LAMC, :LAMG, :LAMI, :LAMR, :LAMS)

"""Coordinates seeded rather than specified, and so exempt from coverage."""
const ATLAS_COORDINATES = (:time, :z, :p)

"""
    atlas_var_specs(thermo_params)

Every variable of the raw Atlas file, keyed by name.

The W/m² water fluxes are divided by the latent heats `thermo_params` carries, so those
must be the ones the archive was written with.
"""
function atlas_var_specs(thermo_params)
    L_cond = TD.Parameters.LH_v0(thermo_params)
    L_sub = TD.Parameters.LH_s0(thermo_params)
    cp_d = TD.Parameters.cp_d(thermo_params)
    specs = Dict{Symbol, NamedTuple}()

    function register!(
        names::Tuple{Vararg{Symbol}};
        kind::Symbol = :raw,
        inputs::Tuple = (),
        raw_inputs::Tuple = (),
        transform,
        units::AbstractString,
        raw_units::AbstractString = "",
        quantity::Symbol,
        provenance::Symbol,
    )
        isempty(units) &&
            error("empty units for $names; every converted variable needs units")
        for name in names
            haskey(specs, name) && error("duplicate spec for $name")
            specs[name] = (;
                kind, inputs, raw_inputs, transform, units, raw_units, quantity,
                provenance,
            )
        end
        return nothing
    end

    # ---- the water reference state ---------------------------------------------------
    # `QT` is vapour + cloud liquid + cloud ice (QT == QV + QN to Float32); rain and snow
    # sit in `QP`, so the self-consistent total-water specific humidity needs both.
    register!((:q_tot,);
        kind = :derived, raw_inputs = (:QT, :QP),
        transform = (qt_g, qp_g) -> wt_to_qt(g_to_kg(qt_g .+ qp_g)),
        units = "kg/kg", quantity = :water_total_specific, provenance = :measured)

    # ---- pressure --------------------------------------------------------------------
    register!((:PRES, :ISCCPPTOP, :MODISPTOP); transform = hPa_to_Pa,
        units = "Pa", raw_units = "mb", quantity = :pressure, provenance = :units_attr)
    # `Ps` is 984.205 hPa under an empty units attribute.
    register!((:Ps,); transform = hPa_to_Pa,
        units = "Pa", raw_units = "", quantity = :pressure,
        provenance = :magnitude_verified)

    # ---- water mass mixing ratios [g/kg] ---------------------------------------------
    register!((
        :QC, :QCI, :QCL, :QCOND, :QG, :QGCLD, :QI, :QICLD, :QN, :QNCLD, :QP, :QPCLD,
        :QPI, :QPL, :QR, :QRCLD, :QS, :QSAT, :QSCLD, :QT, :QTCLD, :QTO, :QTOCLD, :QV,
        :QVOBS,
        ); inputs = (:q_tot,), transform = ComposeFirst(wc_to_qc, g_to_kg),
        units = "kg/kg", raw_units = "g/kg", quantity = :water_mass,
        provenance = :units_attr)
    # `QCCLD`/`QVCLD` are labelled kg/kg but span the same range as `QC`/`QV` in g/kg.
    register!((:QCCLD, :QVCLD); inputs = (:q_tot,),
        transform = ComposeFirst(wc_to_qc, g_to_kg),
        units = "kg/kg", raw_units = "kg/kg", quantity = :water_mass,
        provenance = :magnitude_verified)

    # ---- water mass moments ----------------------------------------------------------
    register!((:QC2, :QI2, :QS2, :QT2); inputs = (:q_tot,),
        transform = wc2_g2kg2_to_kg2kg2,
        units = "(kg/kg)^2", raw_units = "g2/kg2", quantity = :water_mass_moment2,
        provenance = :units_attr)
    # QTO2 == QT2 + QT^2, a second moment rather than the variance its long_name claims;
    # either way it carries two water factors.
    register!((:QTO2,); inputs = (:q_tot,), transform = wc2_g2kg2_to_kg2kg2,
        units = "(kg/kg)^2", raw_units = "(g/kg)^2", quantity = :water_mass_moment2,
        provenance = :measured)
    # Covariance of HL [K] and QT: TQ/sqrt(TL2*QT2) reaches ±1 only once QT2 is scaled out
    # of g2/kg2, so TQ is already K*kg/kg and the attribute `K2` is wrong; one water factor.
    register!((:TQ,); inputs = (:q_tot,), transform = wc_to_qc,
        units = "K kg/kg", raw_units = "K2", quantity = :water_temperature_covariance,
        provenance = :magnitude_verified)

    # ---- water mass rates ------------------------------------------------------------
    register!((
        :QGADV, :QGDIFF, :QGLSADV, :QGMPHY, :QGSED, :QHTEND, :QIADV, :QIDIFF, :QILSADV,
        :QIMPHY, :QISED, :QNUDGE, :QRADV, :QRDIFF, :QRLSADV, :QRMPHY, :QRSED, :QSADV,
        :QSDIFF, :QSLSADV, :QSMPHY, :QSSED, :QTEND, :QTOADV, :QTODIFF, :QTOLSADV,
        :QTOMPHY, :QTOSED, :QVTEND,
        ); inputs = (:q_tot,),
        transform = ComposeFirst(dwcdt_to_dqcdt, perday_to_persec ∘ g_to_kg),
        units = "kg/kg/s", raw_units = "g/kg/day", quantity = :water_mass_rate,
        provenance = :units_attr)
    # M2005 carries these as kg/kg/s (`micro_pumas_v1.F90`); Atlas prints the numerator as
    # `?`. Signs are Atlas's own and uniform across the family: growth and collection
    # positive, evaporation, sublimation and melting negative.
    register!(ATLAS_MASS_PROCESS_RATES; inputs = (:q_tot,),
        transform = ComposeFirst(dwcdt_to_dqcdt, perday_to_persec),
        units = "kg/kg/s", raw_units = "?/kg/day", quantity = :water_mass_rate,
        provenance = :long_name)
    # Yanai apparent moisture sink and total-water storage, both reported as K/day:
    # Q2 = -(L/cp) dq/dt.
    register!((:Q2, :QTSTOR); inputs = (:q_tot,),
        transform = ComposeFirst(
            dwcdt_to_dqcdt, Base.Fix1(*, cp_d / L_cond) ∘ perday_to_persec,
        ),
        units = "kg/kg/s", raw_units = "K/day", quantity = :water_mass_rate,
        provenance = :long_name)

    # ---- water mass fluxes and mass-weighted covariances ------------------------------
    register!((:QCFLUX, :QRFLXR, :QRFLXS, :QRSDFLX, :QTFLUX, :QTOFLXR, :QTOFLXS,
        :QTOSDFLX); transform = x -> wm2_to_kgm2s_liq(x, L_cond),
        units = "kg/m2/s", raw_units = "W/m2", quantity = :water_mass_flux,
        provenance = :measured)
    register!((:QGFLXR, :QGFLXS, :QGSDFLX, :QIFLUX, :QIFLXR, :QIFLXS, :QISDFLX, :QSFLXR,
        :QSFLXS, :QSSDFLX); transform = x -> wm2_to_kgm2s_ice(x, L_sub),
        units = "kg/m2/s", raw_units = "W/m2", quantity = :water_mass_flux,
        provenance = :measured)
    register!((:MFQTCLD, :MFQTCLDA); inputs = (:q_tot,),
        transform = ComposeFirst(wc_to_qc, g_to_kg),
        units = "kg/m2/s", raw_units = "g/m2/s", quantity = :water_mass_flux,
        provenance = :units_attr)
    register!((:QCWCLD, :QIWCLD, :QTWCLD); inputs = (:q_tot,),
        transform = ComposeFirst(wc_to_qc, g_to_kg),
        units = "kg/kg m/s", raw_units = "g/kg m/s",
        quantity = :water_velocity_covariance, provenance = :units_attr)
    # Budget of total-water variance: QT2/|Q2GRAD| is a physical timescale (6.5e3-9.1e3 s),
    # so these carry the variance units the attribute omits.
    register!((:Q2ADVTR, :Q2DIFTR, :Q2DISSIP, :Q2GRAD, :Q2PREC); inputs = (:q_tot,),
        transform = wc2_g2kg2_to_kg2kg2,
        units = "(kg/kg)^2/s", raw_units = "1/s", quantity = :water_mass_moment2_rate,
        provenance = :magnitude_verified)

    # ---- water paths and precipitation ------------------------------------------------
    register!((:CWP, :GWP, :IWP, :LWP, :MODISIWP, :MODISLWP, :RWP, :SWP);
        transform = g_to_kg,
        units = "kg/m2", raw_units = "g/m2", quantity = :water_path,
        provenance = :units_attr)
    register!((:LWP2,); transform = g_to_kg ∘ g_to_kg,
        units = "(kg/m2)^2", raw_units = "(g/m2)^2", quantity = :water_path_moment2,
        provenance = :units_attr)
    # Depth of liquid water; 1 mm over 1 m² is 1 kg for ρ_w = 1000 kg/m³.
    register!((:PW, :PWOBS); transform = identity,
        units = "kg/m2", raw_units = "mm", quantity = :water_path,
        provenance = :units_attr)
    register!((:PREC, :PRECIP); transform = perday_to_persec,
        units = "kg/m2/s", raw_units = "mm/day", quantity = :precipitation_rate,
        provenance = :units_attr)
    register!((:PRECMAX, :PRECMN); transform = perday_to_persec,
        units = "kg/m2/s", raw_units = "mm/d", quantity = :precipitation_rate,
        provenance = :units_attr)
    register!((:PREC2,); transform = perday_to_persec ∘ perday_to_persec,
        units = "(kg/m2/s)^2", raw_units = "(mm/d)^2",
        quantity = :precipitation_rate_moment2, provenance = :units_attr)

    # ---- number concentrations, rates and fluxes --------------------------------------
    register!((:NCMN, :NG, :NGCLD, :NI, :NICLD, :NR, :NRCLD, :NRMN, :NS, :NSCLD);
        transform = percm3_to_perm3,
        units = "1/m3", raw_units = "#/cm3", quantity = :number_concentration,
        provenance = :units_attr)
    register!((
        :NGADV, :NGDIFF, :NGLSADV, :NGMPHY, :NGSED, :NIADV, :NIDIFF, :NILSADV, :NIMPHY,
        :NISED, :NRADV, :NRDIFF, :NRLSADV, :NRMPHY, :NRSED, :NSADV, :NSDIFF, :NSLSADV,
        :NSMPHY, :NSSED,
        ); transform = x -> percm3_to_perm3(perday_to_persec(x)),
        units = "1/m3/s", raw_units = "#/cm3/day", quantity = :number_concentration_rate,
        provenance = :units_attr)
    register!(ATLAS_NUMBER_PROCESS_RATES; transform = perday_to_persec,
        units = "1/kg/s", raw_units = "?/kg/day", quantity = :number_mass_rate,
        provenance = :long_name)
    register!((:NGFLXR, :NGFLXS, :NGSDFLX, :NIFLXR, :NIFLXS, :NISDFLX, :NRFLXR, :NRFLXS,
        :NRSDFLX, :NSFLXR, :NSFLXS, :NSSDFLX); transform = identity,
        units = "1/m2/s", raw_units = "#/m2/s", quantity = :number_flux,
        provenance = :units_attr)

    # ---- inverse lengths --------------------------------------------------------------
    # `long_name` gives "<X> DIAMETER = 1/(2*LAM<X>)", so these are slope parameters.
    register!(ATLAS_SLOPE_PARAMETERS; transform = identity,
        units = "1/m", raw_units = "?/kg/day", quantity = :inverse_length,
        provenance = :long_name)

    # ---- temperatures and their rates -------------------------------------------------
    register!((:ISCCPTB, :ISCCPTBCLR, :SST, :SSTOBS, :TABS, :TABSOBS, :TACLD, :THETA,
        :THETAE, :THETAL, :THETAV, :TVCLD, :TVCLDA); transform = identity,
        units = "K", raw_units = "K", quantity = :temperature, provenance = :units_attr)
    register!((:RADQR, :RADQRC, :RADQRLW, :RADQRS, :RADQRSW, :THTEND, :TNUDGE, :TTEND,
        :TVTEND); transform = perday_to_persec,
        units = "K/s", raw_units = "K/day", quantity = :temperature_rate,
        provenance = :units_attr)
    register!((:RADQRCLW, :RADQRCSW); transform = perday_to_persec,
        units = "K/s", raw_units = "K/d", quantity = :temperature_rate,
        provenance = :units_attr)
    register!((:TVWCLD,); transform = identity,
        units = "K m/s", raw_units = "Km/s",
        quantity = :temperature_velocity_covariance, provenance = :units_attr)
    register!((:MFTVCLD, :MFTVCLDA); transform = identity,
        units = "K kg/m2/s", raw_units = "K kg/m2/s", quantity = :temperature_mass_flux,
        provenance = :units_attr)

    # ---- static energies --------------------------------------------------------------
    # SAM writes these divided by c_pd, so their `K` attribute is a temperature only by
    # that normalisation and the `long_name` gives the quantity away.
    register!((:DSE, :DSECLD, :HFCLD, :HFCLDA, :MSE, :MSECLD, :SSE);
        transform = x -> K_to_J_kg(x, cp_d),
        units = "J/kg", raw_units = "K", quantity = :static_energy,
        provenance = :long_name)
    register!((:MFHCLD, :MFHCLDA); transform = x -> K_to_J_kg(x, cp_d),
        units = "W/m2", raw_units = "K kg/m2/s", quantity = :static_energy_mass_flux,
        provenance = :long_name)

    # ---- the TL / HL family: θ_li, in K -----------------------------------------------
    # Hvar / HQTcov / QTvar are formed from θ_liq_ice and TL2 / TQ / QT2 are compared
    # against them, so this family stays in the c_pd-normalised units it is stored in.
    register!((:TL, :TLCLD); transform = identity,
        units = "K", raw_units = "K", quantity = :temperature, provenance = :units_attr)
    register!((:TL2,); transform = identity,
        units = "K2", raw_units = "K2", quantity = :temperature_moment2,
        provenance = :units_attr)
    register!((:T2ADVTR, :T2DIFTR, :T2DISSIP, :T2GRAD, :T2PREC); transform = identity,
        units = "K2/s", raw_units = "K2/s", quantity = :temperature_moment2_rate,
        provenance = :units_attr)
    register!((:HLADV, :HLDIFF, :HLLAT, :HLRAD, :HLSTOR, :Q1C);
        transform = perday_to_persec,
        units = "K/s", raw_units = "K/day", quantity = :temperature_rate,
        provenance = :units_attr)
    register!((:TLWCLD,); transform = identity,
        units = "K m/s", raw_units = "Km/s",
        quantity = :temperature_velocity_covariance, provenance = :units_attr)
    register!((:MFTLCLD, :MFTLCLDA); transform = identity,
        units = "K kg/m2/s", raw_units = "K kg/m2/s", quantity = :temperature_mass_flux,
        provenance = :units_attr)
    register!((:TWADV, :TWBUOY, :TWDIFF, :TWGRAD, :TWPREC, :TWPRES); transform = identity,
        units = "K m/s2", raw_units = "m s-2 K", quantity = :temperature_flux_budget,
        provenance = :units_attr)

    # ---- energy fluxes ----------------------------------------------------------------
    register!((:LHF, :LWDS, :LWNS, :LWNSC, :LWNT, :LWNTOA, :LWNTOAC, :RADLWDN, :RADLWUP,
        :RADSWDN, :RADSWUP, :SHF, :SOLIN, :SWDS, :SWNS, :SWNSC, :SWNT, :SWNTOA,
        :SWNTOAC, :TLFLUX, :TLFLUXS, :TVFLUX); transform = identity,
        units = "W/m2", raw_units = "W/m2", quantity = :energy_flux,
        provenance = :units_attr)
    register!((:LHFOBS,); transform = mask_undeclared_sentinel,
        units = "W/m2", raw_units = "W/m2", quantity = :energy_flux,
        provenance = :magnitude_verified)
    # The units attribute is the literal string "SHFOBS".
    register!((:SHFOBS,); transform = mask_undeclared_sentinel,
        units = "W/m2", raw_units = "SHFOBS", quantity = :energy_flux,
        provenance = :long_name)

    # ---- lengths ----------------------------------------------------------------------
    register!((:MISRZTOP, :ZCB, :ZCBMIN, :ZCT, :ZCTMAX, :ZINV); transform = km_to_m,
        units = "m", raw_units = "km", quantity = :height, provenance = :units_attr)
    register!((:ZCT2, :ZINV2); transform = km2_to_m2,
        units = "m2", raw_units = "km2", quantity = :height_moment2,
        provenance = :units_attr)
    # Labelled `km`, but a variance of a height cannot be a length; its siblings say km2.
    register!((:ZCB2,); transform = km2_to_m2,
        units = "m2", raw_units = "km", quantity = :height_moment2,
        provenance = :long_name)
    # "mkm" is micrometres: MODISREL spans 8.9-10.3, a cloud-droplet effective radius.
    register!((:MODISREI, :MODISREL); transform = micron_to_m,
        units = "m", raw_units = "mkm", quantity = :length,
        provenance = :magnitude_verified)

    # ---- momentum, turbulence and remaining SI diagnostics ----------------------------
    register!((:U, :UCLD, :UCLDA, :UMAX, :UOBS, :V, :VCLD, :VCLDA, :VOBS, :WCLD, :WCLDA,
        :WMAX, :WOBS); transform = identity,
        units = "m/s", raw_units = "m/s", quantity = :velocity, provenance = :units_attr)
    register!((:TKE, :TKES, :U2, :UW, :UWCLD, :UWSB, :UWSBCLD, :V2, :VW, :VWCLD, :VWSB,
        :VWSBCLD, :W2); transform = identity,
        units = "m2/s2", raw_units = "m2/s2", quantity = :velocity_moment2,
        provenance = :units_attr)
    register!((:W3, :WSTAR3); transform = identity,
        units = "m3/s3", raw_units = "m3/s3", quantity = :velocity_moment3,
        provenance = :units_attr)
    register!((
        :ADVTR, :ADVTRS, :BUOYA, :BUOYAS, :DIFTR, :DISSIP, :DISSIPS, :PRESSTR, :SHEAR,
        :SHEARS, :W2ADV, :W2BUOY, :W2DIFF, :W2PRES, :W2REDIS, :WUADV, :WUANIZ, :WUBUOY,
        :WUDIFF, :WUPRES, :WUSHEAR, :WVADV, :WVANIZ, :WVBUOY, :WVDIFF, :WVPRES, :WVSHEAR,
        ); transform = identity,
        units = "m2/s3", raw_units = "m2/s3", quantity = :velocity_moment2_rate,
        provenance = :units_attr)
    register!((
        :UADV, :UDIFF, :ULSADVV, :UNUDGE, :URESID, :USTOR, :UTENDCOR, :VADV, :VDIFF,
        :VLSADVV, :VNUDGE, :VRESID, :VSTOR, :VTENDCOR,
        ); transform = perday_to_persec,
        units = "m/s2", raw_units = "m/s/day", quantity = :velocity_rate,
        provenance = :units_attr)
    register!((:UPGFCLD, :VPGFCLD, :WPGFCLD); transform = identity,
        units = "m/s2", raw_units = "m/s2", quantity = :velocity_rate,
        provenance = :units_attr)
    # All six are terms in the same ⟨w'q_t'⟩ budget, but only three carry the `K` — and
    # they sit within a factor of 10 of the other three on RF01/RF09/RF12, not 1000, so
    # the water is a mixing ratio in both halves and the `K` is a stray label from TW*.
    register!((:QWADV, :QWDIFF, :QWGRAD); inputs = (:q_tot,), transform = wc_to_qc,
        units = "m/s2 kg/kg", raw_units = "m s-2 K", quantity = :water_flux_budget,
        provenance = :magnitude_verified)
    register!((:QWBUOY, :QWPREC, :QWPRES); inputs = (:q_tot,), transform = wc_to_qc,
        units = "m/s2 kg/kg", raw_units = "m s-2", quantity = :water_flux_budget,
        provenance = :magnitude_verified)
    register!((:RUWCLD, :RVWCLD, :RWWCLD); transform = identity,
        units = "kg/m/s2", raw_units = "kg/m/s2", quantity = :momentum_flux,
        provenance = :units_attr)
    register!((:MC, :MCDNS, :MCDNU, :MCR, :MCRDNS, :MCRDNU, :MCRUP, :MCUP, :MFCLD);
        transform = identity,
        units = "kg/m2/s", raw_units = "kg/m2/s", quantity = :mass_flux,
        provenance = :units_attr)
    register!((:TK, :TKH); transform = identity,
        units = "m2/s", raw_units = "m2/s", quantity = :diffusivity,
        provenance = :units_attr)
    register!((:RHO,); transform = identity,
        units = "kg/m3", raw_units = "kg/m3", quantity = :density,
        provenance = :units_attr)
    register!((:CAPE, :CAPEOBS, :CIN, :CINOBS); transform = identity,
        units = "J/kg", raw_units = "J/kg", quantity = :specific_energy,
        provenance = :units_attr)
    register!((:RELH,); transform = percent_to_fraction,
        units = "1", raw_units = "per cent", quantity = :relative_humidity,
        provenance = :units_attr)
    register!((:QCOEFFR, :QGOEFFR, :QIOEFFR, :QROEFFR, :QSOEFFR); inputs = (:q_tot,),
        transform = ComposeFirst(wc_to_qc, permicron_to_perm ∘ g_to_kg),
        units = "kg/kg/m", raw_units = "g/kg/micro", quantity = :water_mass_per_length,
        provenance = :units_attr)

    # ---- dimensionless ----------------------------------------------------------------
    register!((
        :AREAPREC, :AREAPRTHR, :AUP, :CLD, :CLD245, :CLDCRM30, :CLDCRM40, :CLDCUMDN,
        :CLDCUMUP, :CLDHI, :CLDLOW, :CLDMID, :CLDSHD, :CLRAD00, :CLRAD05, :CLRAD10,
        :CLRAD15, :CLRAD20, :CLRADHI, :CLRADLO, :CLRADM05, :CLRADM10, :CLRADM15,
        :CLRADM20, :CLRADM25, :CLRADM30, :CLRADM35, :CORECL, :COREDNCL, :HYDRO,
        :ISCCPALB, :ISCCPHGH, :ISCCPLOW, :ISCCPMID, :ISCCPTAU, :ISCCPTOT, :MISRTOT,
        :MODISHGH, :MODISLOW, :MODISMID, :MODISTAU, :MODISTAUI, :MODISTAUL, :MODISTOT,
        :MODISTOTI, :MODISTOTL, :WSKEW,
        ); transform = identity,
        units = "1", raw_units = "", quantity = :dimensionless, provenance = :long_name)
    register!((:QGFRAC, :QIFRAC, :QRFRAC, :QSFRAC, :TAUQC, :TAUQG, :TAUQI, :TAUQR,
        :TAUQS); transform = identity,
        units = "1", raw_units = "1", quantity = :dimensionless,
        provenance = :units_attr)

    return specs
end

"""Atlas's `QSMALL`: the specific humidity below which a species counts as absent."""
const ATLAS_Q_SMALL = 1.0e-14

"""Terminal-velocity caps [m/s] for velocities diagnosed from a sedimentation flux."""
const ATLAS_SED_W_CAP = (; cloud = -0.2, ice = -1.0, snow = -2.0, rain = -10.0)

"""
    sedimentation_velocity(F, q, ρ, cap)

Mass-weighted fall speed [m/s] implied by a sedimentation flux `F` [kg m^-2 s^-1] of a
species with specific content `q` and air density `ρ`, positive downward.

Where the species is absent `F/(qρ)` is `0/0` and is reported as zero fall speed. Where
`q` is tiny but non-zero the quotient is unbounded — it reaches 1e6 m/s for ice on the raw
files — so it is clamped into `[cap, 0]`, which binds on under 0.6% of points.
"""
function sedimentation_velocity(F, q, ρ, cap)
    w = F ./ (q .* ρ)
    w = ifelse.(isfinite.(w), w, zero(eltype(w)))
    return .-clamp.(w, convert(eltype(w), cap), zero(eltype(w)))
end

"""
Threshold for cloud-ice autoconversion, as a radius [m].

Atlas splits deposition at the *diameter* `DCS = 125 µm`; `λ` here is radius-conjugate
where Morrison's is diameter-conjugate (`λ_r = 2λ_D`), so `λ_r r_th = λ_D DCS`.
"""
const ATLAS_ICE_THRESHOLD_RADIUS = 62.5e-6

"""Smallest particle radius resolved [m]; regularizes the slope as `q` goes to zero."""
const ATLAS_PARTICLE_MIN_RADIUS = 0.2e-6

"""
    ice_distribution_slope(q, N, ρ, ρ_ice)

Slope `λ` [m^-1] of the gamma ice size distribution implied by specific content `q`,
number per volume `N` and air density `ρ`, with shape parameter `μ`.

For `n(r) = n₀ r^μ e^{-λr}` and `m(r) = (4/3)πρ_ice r³`,
`N = n₀ Γ(μ+1)/λ^{μ+1}` and `ρq = (4/3)πρ_ice n₀ Γ(μ+4)/λ^{μ+4}`, so

    λ³ = (4/3)πρ_ice (μ+1)(μ+2)(μ+3) N/(ρq)

using `Γ(μ+4)/Γ(μ+1) = (μ+1)(μ+2)(μ+3)`, which is exact for any real `μ` and needs no
gamma function. `μ = 0` recovers the exponential case, `λ³ = 8πρ_ice N/(ρq)`.

`N` counts particles per unit volume. `NaN` where there is no ice.
"""
function ice_distribution_slope(
    q::FT,
    N::FT,
    ρ::FT,
    ρ_ice::FT,
    μ::FT = FT(0);
    r_min::FT = FT(ATLAS_PARTICLE_MIN_RADIUS),
) where {FT}
    (isfinite(q) && isfinite(N) && isfinite(ρ) && q > 0 && N > 0) || return FT(NaN)
    m_min = FT(4 // 3) * FT(π) * ρ_ice * r_min^3
    q_eff = q + N * m_min / ρ
    shape = (μ + 1) * (μ + 2) * (μ + 3)
    return cbrt(FT(4 // 3) * FT(π) * ρ_ice * shape * N / (ρ * q_eff))
end

"""
    ice_mean_radius(q, N, ρ, ρ_ice, μ = 0)

Number-weighted mean ice radius [m]: `⟨r⟩ = ∫r n dr / ∫n dr = (μ+1)/λ`.
"""
ice_mean_radius(q::FT, N::FT, ρ::FT, ρ_ice::FT, μ::FT = zero(FT); r_min::FT = FT(ATLAS_PARTICLE_MIN_RADIUS)) where {FT} =
    (μ + one(FT)) / ice_distribution_slope(q, N, ρ, ρ_ice, μ; r_min = r_min)

"""
    ice_process_threshold_weight(q, N, ρ, ρ_ice)

Deposition-weighted fraction of the ice distribution below
[`ATLAS_ICE_THRESHOLD_RADIUS`](@ref): `1 - e^{-λ r_th}(1 + λ r_th)`, Morrison's `DUM`.

The radius weighting is the one deposition takes, since `dm/dt = S/τ` with `τ ∝ r`. This
closed form holds only for the exponential distribution, which is the Atlas LES ice
distribution (`N0I = NI·λᵢ`), so `μ` is fixed at zero here rather than exposed.
"""
function ice_process_threshold_weight(
    q::FT,
    N::FT,
    ρ::FT,
    ρ_ice::FT;
    r_threshold::FT = FT(ATLAS_ICE_THRESHOLD_RADIUS),
    kwargs...,
) where {FT}
    λ = ice_distribution_slope(q, N, ρ, ρ_ice; kwargs...)
    x = λ * r_threshold
    return one(FT) - exp(-x) * (one(FT) + x)
end

"""Water-vapour diffusivity in air [m^2 s^-1] (Pruppacher and Klett 1997)."""
water_vapor_diffusivity_in_air(T, p) =
    oftype(T, 2.11e-5) * (T / oftype(T, 273.15))^oftype(T, 1.94) * (oftype(T, 101325) / p)

"""
    ice_theory_timescale(q, N, ρ, T, p, ρ_ice, μ = 0)

Supersaturation relaxation timescale for ice [s], `1/(4π D N r ρ)`, with the mean radius
implied by `(q, N, ρ)` and the vapour diffusivity at `(T, p)`.
"""
function ice_theory_timescale(
    q::FT,
    N::FT,
    ρ::FT,
    T::FT,
    p::FT,
    ρ_ice::FT,
    μ::FT = zero(FT);
    r_min::FT = FT(ATLAS_PARTICLE_MIN_RADIUS),
) where {FT}
    r = ice_mean_radius(q, N, ρ, ρ_ice, μ; r_min = r_min)
    D = water_vapor_diffusivity_in_air(T, p)
    return inv(4 * FT(π) * D * N * r * ρ)
end

"""Cooper (1986) ice-nucleus concentration [m^-3] at temperature `T`."""
get_N_i_Cooper_curve(T) =
    oftype(T, 0.005) * exp(oftype(T, 0.304) * (oftype(T, 273.15) - T)) * oftype(T, 1000)

"""
    add_derived_specs!(specs, thermo_params, ρ_ice)

Add the derived entries to `specs`, in place. Their inputs are entries of `specs`, so
they arrive already converted.

`ρ_ice` sets the ice size distribution, so it should be the density the archive was
written with; `configs/socrates.toml` defaults `cloud_ice_apparent_density` to the Atlas
value. `q_small` and `sed_w_cap` default to Atlas's own thresholds and may be replaced.
"""
function add_derived_specs!(
    specs::AbstractDict{Symbol, <:NamedTuple},
    thermo_params,
    ρ_ice;
    q_small = ATLAS_Q_SMALL,
    sed_w_cap = ATLAS_SED_W_CAP,
)
    function register_derived!(
        name::Symbol;
        inputs::Tuple,
        transform,
        units::AbstractString,
        quantity::Symbol,
    )
        haskey(specs, name) && error("duplicate spec for $name")
        isempty(units) && error("empty units for $name")
        specs[name] = (;
            kind = :derived, inputs, raw_inputs = (), transform, units, raw_units = "",
            quantity, provenance = :derived,
        )
        return nothing
    end

    # ---- standard deviations of the water moments ------------------------------------
    for (std_name, var_name) in
        ((:QT2_STD, :QT2), (:QC2_STD, :QC2), (:QI2_STD, :QI2), (:QS2_STD, :QS2))
        register_derived!(std_name; inputs = (var_name,), transform = x -> sqrt.(x),
            units = "kg/kg", quantity = :water_mass)
    end

    # ---- aliases ----------------------------------------------------------------------
    # `QTO` is vapour + cloud liquid and vapour does not sediment, so the sedimentation of
    # QTO is that of cloud liquid alone.
    register_derived!(:QCSED; inputs = (:QTOSED,), transform = identity,
        units = "kg/kg/s", quantity = :water_mass_rate)
    register_derived!(:QCSDFLX; inputs = (:QTOSDFLX,), transform = identity,
        units = "kg/m2/s", quantity = :water_mass_flux)
    register_derived!(:QCFLXR; inputs = (:QCFLUX,), transform = identity,
        units = "kg/m2/s", quantity = :water_mass_flux)

    # ---- resolved vertical advection of cloud liquid -----------------------------------
    # A flux divergence of a true mass flux gives a mixing-ratio rate, so it takes one
    # water factor.
    register_derived!(:QCADV; inputs = (:QCFLUX, :z, :RHO, :q_tot),
        transform = (F, z, ρ, q) ->
            wc_to_qc(atlas_forward_flux_divergence_tendency(F, z, ρ), q),
        units = "kg/kg/s", quantity = :water_mass_rate)

    # ---- net deposition rates ----------------------------------------------------------
    register_derived!(:PRDT; inputs = (:PRD, :EPRD), transform = +,
        units = "kg/kg/s", quantity = :water_mass_rate)
    register_derived!(:PRDST; inputs = (:PRDS, :EPRDS), transform = +,
        units = "kg/kg/s", quantity = :water_mass_rate)

    # ---- relative humidity -------------------------------------------------------------
    for (name, phase) in ((:RHL, TD.Liquid()), (:RHI, TD.Ice()))
        register_derived!(name; inputs = (:QV, :TABS, :RHO),
            transform = (qv, T, ρ) ->
                qv ./ TD.q_vap_saturation.(Ref(thermo_params), T, ρ, Ref(phase)),
            units = "1", quantity = :relative_humidity)
        # bracketing values one total-water standard deviation either side of the vapour
        for (suffix, sign) in ((:_low, -1), (:_high, +1))
            register_derived!(Symbol(name, suffix); inputs = (:QV, :TABS, :RHO, :QT2_STD),
                transform = (qv, T, ρ, σ) ->
                    (qv .+ sign .* σ) ./
                    TD.q_vap_saturation.(Ref(thermo_params), T, ρ, Ref(phase)),
                units = "1", quantity = :relative_humidity)
        end
    end

    # ---- sedimentation velocities ------------------------------------------------------
    for (name, flux, species, cap) in (
        (:WI, :QISDFLX, :QI, sed_w_cap.ice),
        (:WS, :QSSDFLX, :QS, sed_w_cap.snow),
        (:WR, :QRSDFLX, :QR, sed_w_cap.rain),
        (:WC, :QCSDFLX, :QC, sed_w_cap.cloud),
    )
        register_derived!(name; inputs = (flux, species, :RHO),
            transform = (F, q, ρ) -> sedimentation_velocity(F, q, ρ, cap),
            units = "m/s", quantity = :velocity)
    end

    # ---- ice number --------------------------------------------------------------------
    register_derived!(:NI_ALL; inputs = (:NI, :NS), transform = +,
        units = "1/m3", quantity = :number_concentration)
    register_derived!(:NINP; inputs = (:TABS,),
        transform = T -> get_N_i_Cooper_curve.(T),
        units = "1/m3", quantity = :number_concentration)

    # ---- ice particle radius and relaxation timescales ---------------------------------
    register_derived!(:RI; inputs = (:QI, :NI, :RHO),
        transform = (q, N, ρ) -> ice_mean_radius.(q, N, ρ, convert(eltype(q), ρ_ice)),
        units = "m", quantity = :length)

    # `_full` is the whole ice distribution; the bare name is the part below the 125 µm
    # threshold, which is the part Atlas reports a deposition rate for. A smaller
    # population takes up less vapour, hence `τ_below = τ_full / DUM >= τ_full`.
    register_derived!(:τ_ice_theory_full; inputs = (:QI, :NI, :RHO, :TABS, :PRES),
        transform = (q, N, ρ, T, p) ->
            ice_theory_timescale.(q, N, ρ, T, p, convert(eltype(q), ρ_ice)),
        units = "s", quantity = :timescale)
    register_derived!(:τ_ice_theory; inputs = (:τ_ice_theory_full, :QI, :NI, :RHO),
        transform = (τ, q, N, ρ) -> τ ./ ice_process_threshold_weight.(q, N, ρ, convert(eltype(q), ρ_ice)),
        units = "s", quantity = :timescale)

    # `dq/dt = (q_vap - q_sat)/τ`, so the numerator is the supersaturation excess rather
    # than a humidity ratio. `PRDT`, not `PRD`: exactly one of deposition and sublimation
    # is non-zero at a point and both are the same vapour exchange. Atlas books the tail
    # above the threshold to snow, except with no snow to book it to, when `PRDT` is
    # already the whole rate. Negative where the slab mean is subsaturated while its
    # saturated fraction still deposits — a mean-state diagnostic of a rate that is not
    # linear in the mean state.
    register_derived!(:τ_ice_full; inputs = (:QV, :TABS, :RHO, :PRDT, :QS, :QI, :NI),
        transform = function (qv, T, ρ, prdt, q_sno, q, N)
            excess =
                qv .-
                TD.q_vap_saturation.(Ref(thermo_params), T, ρ, Ref(TD.Ice()))
            below = ice_process_threshold_weight.(q, N, ρ, convert(eltype(q), ρ_ice))
            rate = ifelse.(q_sno .>= oftype.(prdt, q_small), prdt ./ below, prdt)
            return excess ./ ifelse.(rate .== 0, oftype.(rate, NaN), rate)
        end,
        units = "s", quantity = :timescale)
    register_derived!(:τ_ice; inputs = (:τ_ice_full, :QI, :NI, :RHO),
        transform = (τ, q, N, ρ) -> τ ./ ice_process_threshold_weight.(q, N, ρ, convert(eltype(q), ρ_ice)),
        units = "s", quantity = :timescale)

    return specs
end

"""
    atlas_processing_order(specs; seeds = ATLAS_COORDINATES)

Entries ordered so every `inputs` dependency precedes its consumer, by Kahn's algorithm.
Ties break by name, so the order and any failure are reproducible run to run.
"""
function atlas_processing_order(
    specs::AbstractDict{Symbol, <:NamedTuple};
    seeds = ATLAS_COORDINATES,
)
    for (name, spec) in specs, dep in spec.inputs
        (haskey(specs, dep) || dep in seeds) ||
            error("input $dep of $name is neither a spec nor a seeded coordinate")
    end

    indegree =
        Dict{Symbol, Int}(name => count(!in(seeds), spec.inputs) for (name, spec) in specs)
    dependents = Dict{Symbol, Vector{Symbol}}(name => Symbol[] for name in keys(specs))
    for (name, spec) in specs, dep in spec.inputs
        dep in seeds || push!(dependents[dep], name)
    end

    order = Vector{Symbol}(undef, length(specs))
    ready = sort!([name for (name, deg) in indegree if deg == 0])
    filled = 0
    head = 1
    while head <= length(ready)
        name = ready[head]
        head += 1
        filled += 1
        order[filled] = name
        for dependent in sort!(dependents[name])
            indegree[dependent] -= 1
            indegree[dependent] == 0 && push!(ready, dependent)
        end
    end
    filled == length(specs) || error(
        "cyclic or unsatisfiable `inputs` among: " *
        join(sort!([n for (n, d) in indegree if d > 0]), ", "),
    )
    return order
end

"""
    assert_specs_cover(specs, raw_names)

Error unless every variable of the raw Atlas file has a spec, so an unclassified variable
can never be read with unknown units, and unless every `:raw` spec exists in the file.
"""
function assert_specs_cover(specs::AbstractDict{Symbol, <:NamedTuple}, raw_names)
    raw = Set(Symbol(n) for n in raw_names if Symbol(n) ∉ ATLAS_COORDINATES)
    unclassified = sort!(collect(setdiff(raw, keys(specs))))
    isempty(unclassified) || error(
        "$(length(unclassified)) Atlas variable(s) have no spec; classify them before \
         reading: " * join(unclassified, ", "),
    )
    absent = sort!([name for (name, spec) in specs if spec.kind === :raw && name ∉ raw])
    isempty(absent) ||
        error(":raw spec(s) absent from the raw Atlas file: " * join(absent, ", "))
    return nothing
end
