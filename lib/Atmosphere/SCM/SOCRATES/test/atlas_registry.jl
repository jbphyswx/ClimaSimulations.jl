using ClimaAtmos: ClimaAtmos as CA
using SOCRATES: SOCRATES
using SOCRATESSingleColumnForcings: SOCRATESSingleColumnForcings as SSCF
using Test: Test

# Each assertion here has to be able to fail for a reason in this package. Facts about the
# published archive — that a Southern Ocean column is cold, that ice supersaturates — are
# properties of the data, not of the code, and are recorded in docstrings instead.

const RF09 = SOCRATES.case("RF09_Obs")
const PARAMS = SOCRATES.socrates_params(Float64, RF09)
const SPECS = SOCRATES.atlas_specs(PARAMS)

Test.@testset "no archive variable is read unclassified" begin
    for c in SOCRATES.all_cases()
        les = SSCF.open_atlas_les_output(c.flight_number, c.forcing_type)
        Test.@test SOCRATES.assert_specs_cover(SPECS, keys(les.data)) === nothing
    end
    Test.@test_throws ErrorException SOCRATES.assert_specs_cover(SPECS, ("NOT_A_VAR",))
    Test.@test_throws ErrorException SOCRATES.atlas_dependencies(SPECS, (:NOT_A_VAR,))

    # a consumer must never be evaluated before an input it reads
    order = SOCRATES.atlas_processing_order(SPECS)
    position = Dict(name => i for (i, name) in enumerate(order))
    Test.@test length(order) == length(SPECS)
    for (name, spec) in SPECS, dep in spec.inputs
        dep in SOCRATES.ATLAS_COORDINATES && continue
        Test.@test position[dep] < position[name]
    end
end

Test.@testset "a spec that disagrees with the file is refused" begin
    Test.@test haskey(
        SOCRATES.read_atlas(RF09, (:TABS,); params = PARAMS, specs = SPECS).data, :TABS,
    )
    tampered = copy(SPECS)
    tampered[:TABS] = merge(SPECS[:TABS], (; raw_units = "degC"))
    Test.@test_throws ErrorException SOCRATES.read_atlas(
        RF09, (:TABS,); params = PARAMS, specs = tampered,
    )
end

Test.@testset "the sedimentation identity holds through the conversions" begin
    # Q*SED = -(1/rho) d(Q*SDFLX/L)/dz is a relation between two independently written
    # archive variables, so it tests this package's flux-divergence operator, its unit
    # scalings and the latent heats together. Nothing here is a fact about the weather.
    raw = SOCRATES.read_atlas(
        RF09, (:QCSED, :QCSDFLX, :RHO, :q_tot); params = PARAMS, specs = SPECS,
    )
    d = raw.data
    div = SOCRATES.atlas_forward_flux_divergence_tendency(d[:QCSDFLX], raw.z, d[:RHO])
    rhs = SOCRATES.wc_to_qc.(div, d[:q_tot])
    lhs = d[:QCSED]
    idx = [(i, j) for i in 1:(size(lhs, 1) - 1), j in axes(lhs, 2)
           if isfinite(lhs[i, j]) && isfinite(rhs[i, j]) && abs(lhs[i, j]) > 1e-12]
    rel = sort([abs(rhs[i, j] - lhs[i, j]) / abs(lhs[i, j]) for (i, j) in idx])
    # the archive is Float32, so agreement is bounded by its stored precision
    Test.@test rel[cld(length(rel), 2)] < 1.0e-5
end

Test.@testset "sedimentation_velocity" begin
    # exercised with inputs chosen to hit each branch rather than with whatever the
    # archive happens to contain
    cap = SOCRATES.ATLAS_SED_W_CAP.ice           # -1.0 m/s
    F = [-5.0e-6, -1.0e-4, 0.0, -1.0e-4]
    q = [1.0e-5, 1.0e-30, 1.0e-5, 0.0]           # normal, tiny, no flux, absent
    ρ = [1.0, 1.0, 1.0, 1.0]
    w = SOCRATES.sedimentation_velocity(F, q, ρ, cap)
    Test.@test w[1] ≈ 0.5                        # -F/(q rho), positive downward, uncapped
    Test.@test w[2] == -cap                      # unbounded quotient meets the cap
    Test.@test w[3] == 0                         # no flux, no fall speed
    Test.@test w[4] == 0                         # 0/0 reported as zero, not NaN
end

Test.@testset "ice_distribution_slope inverts the mass relation" begin
    # lambda^3 = (4/3) pi rho_ice (mu+1)(mu+2)(mu+3) N / (rho q): recovering q from lambda
    # is an independent check of the algebra rather than a restatement of it
    ρ_ice, N, ρ, q = 500.0, 1.0e5, 1.0, 1.0e-5
    λ = SOCRATES.ice_distribution_slope(q, N, ρ, ρ_ice)
    Test.@test 8 * π * ρ_ice * N / (ρ * λ^3) ≈ q rtol = 1.0e-6
    # the mean radius of an exponential distribution is 1/lambda
    Test.@test SOCRATES.ice_mean_radius(q, N, ρ, ρ_ice) ≈ inv(λ)
    # no ice gives no radius rather than a number
    Test.@test isnan(SOCRATES.ice_distribution_slope(0.0, N, ρ, ρ_ice))
    Test.@test isnan(SOCRATES.ice_distribution_slope(q, 0.0, ρ, ρ_ice))
    # the threshold weight is a fraction, and a coarser distribution keeps less below it
    fine = SOCRATES.ice_process_threshold_weight(1.0e-7, N, ρ, ρ_ice)
    coarse = SOCRATES.ice_process_threshold_weight(1.0e-4, N, ρ, ρ_ice)
    Test.@test 0 < coarse < fine < 1
end

Test.@testset "every tendency source is a rate" begin
    for (tendency, sources) in SOCRATES.LES_TENDENCIES
        name = String(tendency)
        Test.@test name in SOCRATES.MP1M_SOURCE_TERMS ||
                   any(startswith(name, p) for p in SOCRATES.TRANSPORT_PREFIXES)
        for (source, _) in sources
            Test.@test haskey(SPECS, source)
            # a tendency mapped to a variable the registry converts to something else
            Test.@test SPECS[source].units == "kg/kg/s"
        end
        name in SOCRATES.MP1M_SOURCE_TERMS || continue
        Test.@test any(
            any(term == name for (term, _) in terms) for
            terms in values(SOCRATES.MP1M_BUDGETS)
        )
    end
end

Test.@testset "mapped and unavailable partition the model budget" begin
    budget = Set(
        Symbol.((
            SOCRATES.MP1M_SOURCE_TERMS...,
            ("$(p)_$(v)" for p in SOCRATES.TRANSPORT_PREFIXES for
             v in SOCRATES.MP1M_BUDGET_VARS)...,
        )),
    )
    mapped = Set(keys(SOCRATES.LES_TENDENCIES))
    absent = Set(keys(SOCRATES.LES_TENDENCIES_UNAVAILABLE))
    Test.@test isempty(intersect(mapped, absent))
    Test.@test isempty(setdiff(budget, union(mapped, absent)))
    Test.@test isempty(setdiff(union(mapped, absent), budget))
end

Test.@testset "each process rate records the archive's own description" begin
    les = SSCF.open_atlas_les_output(RF09.flight_number, RF09.forcing_type)
    for (rate, entry) in SOCRATES.ATLAS_PROCESS_RATES
        Test.@test haskey(les.data, String(rate))
        long_name = strip(get(les.data[String(rate)].attrib, "long_name", ""))
        Test.@test endswith(long_name, entry.long_name)
        terms = SOCRATES.atlas_model_terms(rate)
        if isnothing(terms)
            Test.@test entry.model in SOCRATES.NO_MODEL_COUNTERPART
        else
            for (term, sign) in terms
                Test.@test String(term) in SOCRATES.MP1M_SOURCE_TERMS
                Test.@test abs(sign) == 1
            end
        end
    end
    Test.@test_throws ErrorException SOCRATES.atlas_model_terms(:NOT_A_RATE)
end

Test.@testset "the rate lists are the table's own mass and number entries" begin
    Test.@test length(SOCRATES.ATLAS_PROCESS_RATES) ==
               length(SOCRATES.ATLAS_MASS_PROCESS_RATES) +
               length(SOCRATES.ATLAS_NUMBER_PROCESS_RATES)
    for rate in SOCRATES.ATLAS_MASS_PROCESS_RATES
        Test.@test SOCRATES.ATLAS_PROCESS_RATES[rate].kind === :mass
    end
    for rate in SOCRATES.ATLAS_NUMBER_PROCESS_RATES
        Test.@test SOCRATES.ATLAS_PROCESS_RATES[rate].kind === :number
    end
end
