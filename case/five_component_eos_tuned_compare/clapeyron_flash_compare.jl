#!/usr/bin/env julia

using Clapeyron
using Printf

const COMPONENTS = ["water", "carbon dioxide", "methane", "ethane", "n-butane"]
const COMPONENT_LABELS = ["H2O", "CO2", "CH4", "C2H6", "nC4H10"]
const P = 6.0e6
const T = 305.0
const ZFEED = [0.30, 0.10, 0.15, 0.15, 0.30]
const RGAS = 8.31446261815324

const TC = [647.096, 304.128200003, 190.564002651, 305.322, 425.125]
const PC = [22.064000000e6, 7.377298373e6, 4.599200474e6,
            4.872200000e6, 3.796000017e6]
const MW_G_PER_MOL = 1000.0 .* [0.018015268, 0.0440098, 0.0160428,
                                0.03006904, 0.0581222]
const OMEGA = [0.3442920843, 0.22394, 0.01142, 0.099, 0.200810094644]

const PR_K = [
     0.0                 -0.0391803910093841  -0.2                 0.5                 0.5;
    -0.0391803910093841   0.0                  0.110007640625       0.174325125         0.1277;
    -0.2                  0.110007640625       0.0                 0.002281521088      0.01088232117888;
     0.5                  0.174325125          0.002281521088      0.0                 0.0;
     0.5                  0.1277               0.01088232117888    0.0                 0.0
]

const CPA_K = [
    0.0                 0.158412686728       0.0204118168704    0.04415             0.08750;
    0.158412686728      0.0                  0.0096233516       0.0677109300736     0.11220;
    0.0204118168704     0.0096233516         0.0               -0.004599425       -0.02189179616;
    0.04415             0.0677109300736     -0.004599425        0.0                 0.0;
    0.08750             0.11220             -0.02189179616      0.0                 0.0
]

const CPA_B = [14.52e-6, 27.2e-6, 29.1e-6, 42.9e-6, 72.081e-6]
const CPA_GAMMA = [1017.3, 1551.22, 959.02, 1544.54, 2193.08]
const CPA_A0 = CPA_B .* RGAS .* CPA_GAMMA
const CPA_C1 = [0.6736, 0.7602, 0.4472, 0.5846, 0.7077]

value(row, name) = row[name]

function read_our_flash(path)
    lines = readlines(path)
    length(lines) >= 2 || error("Missing data row in $path")
    header = split(lines[1], ',')
    fields = split(lines[2], ',')
    length(header) == length(fields) || error("Malformed CSV row in $path")
    row = Dict(header[i] => parse(Float64, fields[i]) for i in eachindex(header))
    beta = [value(row, "beta_oil"), value(row, "beta_gas"),
            value(row, "beta_water")]
    kg = [value(row, "K_gas_over_oil_$(i)_$(COMPONENT_LABELS[i + 1])")
          for i in 0:4]
    kw = [value(row, "K_water_over_oil_$(i)_$(COMPONENT_LABELS[i + 1])")
          for i in 0:4]
    oil = ZFEED ./ (beta[1] .+ beta[2] .* kg .+ beta[3] .* kw)
    gas = kg .* oil
    water = kw .* oil
    compositions = [oil ./ sum(oil), gas ./ sum(gas), water ./ sum(water)]
    zfactor = [value(row, "Z_oil"), value(row, "Z_gas"), value(row, "Z_water")]
    volumes = zfactor .* RGAS .* T ./ P
    saturation = beta .* volumes ./ sum(beta .* volumes)
    return (; beta, saturation, zfactor, compositions)
end

function make_pr()
    return PR78(
        COMPONENTS;
        userlocations=(; Tc=TC, Pc=PC, Mw=MW_G_PER_MOL, k=PR_K),
        alpha_userlocations=(; acentricfactor=OMEGA),
    )
end

function make_cpa()
    epsilon_assoc = Dict(
        (("water", "e"), ("water", "H")) => 2003.2,
        (("carbon dioxide", "e"), ("water", "H")) => 14200.0 / RGAS,
    )
    bondvol = Dict(
        (("water", "e"), ("water", "H")) => 0.0692,
        (("carbon dioxide", "e"), ("water", "H")) => 0.0162,
    )
    return CPA(
        COMPONENTS;
        radial_dist=:KG,
        assoc_options=Clapeyron.AssocOptions(implicit_ad=true),
        userlocations=(;
            Mw=MW_G_PER_MOL,
            Tc=TC,
            Pc=PC,
            a=CPA_A0,
            b=CPA_B,
            c1=CPA_C1,
            k=CPA_K,
            n_H=[2, 0, 0, 0, 0],
            n_e=[2, 1, 0, 0, 0],
            epsilon_assoc=epsilon_assoc,
            bondvol=bondvol,
        ),
    )
end

function canonical_clapeyron_result(model, result)
    composition_matrix, mole_matrix, _ = result
    nphase = size(composition_matrix, 1)
    nphase == 3 || error("Expected three phases, Clapeyron returned $nphase")

    compositions_raw = [collect(composition_matrix[i, :]) for i in 1:nphase]
    beta_raw = [sum(mole_matrix[i, :]) for i in 1:nphase]
    water = argmax([x[1] for x in compositions_raw])
    hydrocarbon_phases = filter(i -> i != water, 1:nphase)
    vapour = hydrocarbon_phases[
        argmax([compositions_raw[i][2] + compositions_raw[i][3]
                for i in hydrocarbon_phases])]
    oil = only(filter(i -> i != vapour, hydrocarbon_phases))
    order = [oil, vapour, water]

    beta = beta_raw[order]
    beta ./= sum(beta)
    volumes = [
        volume(model, P, T, compositions_raw[oil]; phase=:liquid),
        volume(model, P, T, compositions_raw[vapour]; phase=:vapour),
        volume(model, P, T, compositions_raw[water]; phase=:liquid),
    ]
    saturation = beta .* volumes ./ sum(beta .* volumes)
    zfactor = P .* volumes ./ (RGAS * T)
    compositions = [compositions_raw[i] ./ sum(compositions_raw[i])
                    for i in order]
    all(isfinite, beta) || error("Non-finite phase fractions")
    all(isfinite, saturation) || error("Non-finite phase saturations")
    all(isfinite, zfactor) || error("Non-finite compressibility factors")
    all(all(isfinite, x) for x in compositions) ||
        error("Non-finite phase compositions")
    return (; beta, saturation, zfactor, compositions)
end

function valid_flash_tuple(result)
    composition_matrix, mole_matrix, gibbs = result
    return size(composition_matrix, 1) == 3 &&
           all(isfinite, composition_matrix) && all(isfinite, mole_matrix) &&
           isfinite(gibbs) && all(vec(sum(mole_matrix; dims=2)) .> 1.0e-10)
end

function run_flash(model, seed; fallbacks=false)
    methods = Tuple{String,Any}[]
    push!(methods,
          ("seeded multiphase Newton",
           MultiPhaseTPFlash(n0=seed, K_tol=1.0e-10, ss_iters=12,
                             phase_iters=30, second_order=true)))
    if fallbacks
        push!(methods,
              ("full-TPD multiphase L-BFGS",
               MultiPhaseTPFlash(full_tpd=true, max_phases=3, K_tol=1.0e-10,
                                 ss_iters=20, phase_iters=40,
                                 second_order=false)))
        push!(methods,
              ("three-phase global SASS",
               DETPFlash(numphases=3, max_steps=12000, population_size=80,
                         time_limit=180.0, seed=20260825,
                         stagnation_evals=2000, stagnation_tol=1.0e-10,
                         logspace=true)))
    end

    total_elapsed = 0.0
    for (label, method) in methods
        result = nothing
        elapsed = @elapsed result = tp_flash(model, P, T, ZFEED, method)
        total_elapsed += elapsed
        if valid_flash_tuple(result)
            try
                state = canonical_clapeyron_result(model, result)
                return state, total_elapsed, label
            catch exception
                @printf("  %s produced unusable derived properties after %.6f s: %s\n",
                        label, elapsed, sprint(showerror, exception))
                flush(stdout)
                continue
            end
        end
        @printf("  %s returned an invalid three-phase result after %.6f s; trying fallback.\n",
                label, elapsed)
        flush(stdout)
    end
    error("All configured Clapeyron flash methods failed to return a finite three-phase state")
end

function run_full_tpd(model)
    method = MultiPhaseTPFlash(
        full_tpd=true, max_phases=3, K_tol=1.0e-10,
        ss_iters=20, phase_iters=40, second_order=false)
    result = nothing
    elapsed = @elapsed result = tp_flash(model, P, T, ZFEED, method)
    if !valid_flash_tuple(result)
        return nothing, elapsed, "invalid_three_phase_state"
    end
    try
        return canonical_clapeyron_result(model, result), elapsed, "converged"
    catch
        return nothing, elapsed, "invalid_derived_properties"
    end
end

function write_summary(path, records)
    open(path, "w") do io
        println(io, "eos,software,method,nphase,beta_oil,beta_gas,beta_water,",
                "saturation_oil,saturation_gas,saturation_water,",
                "Z_oil,Z_gas,Z_water,elapsed_seconds")
        for record in records
            state = record.state
            @printf(io, "%s,%s,%s,3,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,",
                    record.eos, record.software, record.method, state.beta...,
                    state.saturation...)
            @printf(io, "%.16e,%.16e,%.16e,%.9f\n", state.zfactor...,
                    record.elapsed)
        end
    end
end

function write_full_tpd(path, records)
    open(path, "w") do io
        println(io, "eos,status,elapsed_seconds,max_abs_beta,",
                "max_abs_saturation,max_abs_Z,max_abs_phase_composition")
        for record in records
            if isnothing(record.state)
                @printf(io, "%s,%s,%.9f,,,,\n",
                        record.eos, record.status, record.elapsed)
                continue
            end
            state, ours = record.state, record.ours
            max_beta = maximum(abs.(state.beta .- ours.beta))
            max_sat = maximum(abs.(state.saturation .- ours.saturation))
            max_z = maximum(abs.(state.zfactor .- ours.zfactor))
            max_comp = maximum(maximum(abs.(state.compositions[i] .-
                                            ours.compositions[i])) for i in 1:3)
            @printf(io, "%s,%s,%.9f,%.16e,%.16e,%.16e,%.16e\n",
                    record.eos, record.status, record.elapsed,
                    max_beta, max_sat, max_z, max_comp)
        end
    end
end

function write_compositions(path, records)
    phases = ["oil", "gas", "water"]
    open(path, "w") do io
        println(io, "eos,software,phase,", join(COMPONENT_LABELS, ","))
        for record in records, phase in 1:3
            x = record.state.compositions[phase]
            @printf(io, "%s,%s,%s,%.16e,%.16e,%.16e,%.16e,%.16e\n",
                    record.eos, record.software, phases[phase], x...)
        end
    end
end

function write_differences(path, pairs)
    open(path, "w") do io
        println(io, "eos,max_abs_beta,max_abs_saturation,max_abs_Z,",
                "max_abs_phase_composition")
        for pair in pairs
            ours, clap = pair.ours, pair.clap
            max_beta = maximum(abs.(clap.beta .- ours.beta))
            max_sat = maximum(abs.(clap.saturation .- ours.saturation))
            max_z = maximum(abs.(clap.zfactor .- ours.zfactor))
            max_comp = maximum(maximum(abs.(clap.compositions[i] .-
                                            ours.compositions[i])) for i in 1:3)
            @printf(io, "%s,%.16e,%.16e,%.16e,%.16e\n",
                    pair.eos, max_beta, max_sat, max_z, max_comp)
        end
    end
end

if get(ENV, "CLAPEYRON_SKIP_MAIN", "0") != "1"
let
    case_dir = @__DIR__
    input_root = length(ARGS) >= 1 ? abspath(ARGS[1]) :
        joinpath(case_dir, "results", "clapeyron_compare_20260825")
    output_dir = length(ARGS) >= 2 ? abspath(ARGS[2]) :
        joinpath(case_dir, "clapeyron_comparison")

    our_pr = read_our_flash(joinpath(input_root, "our_pr", "phase_state_step_0.csv"))
    our_cpa = read_our_flash(joinpath(input_root, "our_cpa", "phase_state_step_0.csv"))

    println("Running Clapeyron PR78 multiphase flash...")
    flush(stdout)
    clap_pr, pr_time, pr_method = run_flash(make_pr(), our_pr.compositions)
    @printf("PR78 completed in %.6f s\n", pr_time)
    flush(stdout)
    println("Running Clapeyron sCPA multiphase flash...")
    flush(stdout)
    clap_cpa, cpa_time, cpa_method = run_flash(
        make_cpa(), our_cpa.compositions; fallbacks=true)

    records = [
        (; eos="PR78", software="Our", method="production_flash",
           state=our_pr, elapsed=NaN),
        (; eos="PR78", software="Clapeyron", method=pr_method,
           state=clap_pr, elapsed=pr_time),
        (; eos="sCPA", software="Our", method="production_flash",
           state=our_cpa, elapsed=NaN),
        (; eos="sCPA", software="Clapeyron", method=cpa_method,
           state=clap_cpa, elapsed=cpa_time),
    ]
    pairs = [
        (; eos="PR78", ours=our_pr, clap=clap_pr),
        (; eos="sCPA", ours=our_cpa, clap=clap_cpa),
    ]

    mkpath(output_dir)
    write_summary(joinpath(output_dir, "initial_flash_summary.csv"), records)
    write_compositions(joinpath(output_dir, "phase_compositions.csv"), records)
    write_differences(joinpath(output_dir, "comparison_metrics.csv"), pairs)

    if get(ENV, "CLAPEYRON_FULL_TPD", "0") == "1"
        println("Running optional full-TPD searches without phase seeds...")
        full_records = []
        for (eos, model, ours) in (
                ("PR78", make_pr(), our_pr), ("sCPA", make_cpa(), our_cpa))
            state, elapsed, status = run_full_tpd(model)
            push!(full_records, (; eos, state, elapsed, status, ours))
            @printf("%s full TPD: %s (%.6f s)\n", eos, status, elapsed)
        end
        write_full_tpd(joinpath(output_dir, "full_tpd_check.csv"), full_records)
    end

    println("Clapeyron version: ", pkgversion(Clapeyron))
    println("PR78 accepted method: ", pr_method)
    println("sCPA accepted method: ", cpa_method)
    @printf("PR78 flash wall time: %.6f s\n", pr_time)
    @printf("sCPA flash wall time: %.6f s\n", cpa_time)
    for pair in pairs
        max_beta = maximum(abs.(pair.clap.beta .- pair.ours.beta))
        max_sat = maximum(abs.(pair.clap.saturation .- pair.ours.saturation))
        max_comp = maximum(maximum(abs.(pair.clap.compositions[i] .-
                                        pair.ours.compositions[i])) for i in 1:3)
        @printf("%s: max |Δbeta|=%.6e, max |ΔS|=%.6e, max |Δx|=%.6e\n",
                pair.eos, max_beta, max_sat, max_comp)
    end
    println("Results: ", output_dir)
end
end
