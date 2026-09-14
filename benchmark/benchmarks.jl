#!/usr/bin/env julia
# benchmark/benchmarks.jl — BenchmarkTools suite for ERGMCount.jl's hot loops.
#
# Locks in (panel 2026-09, item 25) the compressed MPLE design, the
# allocation-free derivative closure and the de-abstracted Gibbs sweep:
#   fit/zach_default      fit_ergm_count(zach, sum + nonzero) on the default
#                         (error-controlled) support path — the fixture model
#   fit/zach_maxval60     the same at a fixed, wide support
#   derivatives/zach      one (ll, grad, hess) evaluation of the closure
#   sweep/n34, sweep/n100 one Gibbs sweep (every dyad redrawn) on directed
#                         count networks with the SAME mean degree, and the
#                         per-dyad scaling between them asserted: the profile
#                         work per dyad is O(degree), so the sweep must scale
#                         with the number of dyads, not faster.
#
# Defines the standard `SUITE::BenchmarkGroup`. Run standalone with
#     julia --project=benchmark benchmark/benchmarks.jl
# which tunes + runs the suite, prints one tab-separated `BENCHJL` line per
# benchmark (consumed by the site repo's tools/run_benchmarks.jl), and exits
# non-zero if the scaling assertion fails.

using BenchmarkTools
using ERGMCount
using Networks
using Random

include(joinpath(@__DIR__, "fixtures.jl"))

const N_SMALL = 34
const N_LARGE = 100
const SCALING_LIMIT = 2.0    # tolerated (t_large / t_small) / (dyads_large / dyads_small)

const ZACH = zach_network()
const ZACH_MODEL = CountERGMModel((SumTerm(), NonzeroTerm()), ZACH, PoissonReference())
const ZACH_DESIGN = ERGMCount._count_design(ZACH_MODEL, 0:30)
const ZACH_DERIV = ERGMCount._count_derivatives(ZACH_DESIGN, [1, 2],
                                                fill(true, 31, length(ZACH_DESIGN.n_tot)))
const ZACH_THETA = [-0.5, 1.0]

const STATES = Dict(n => sweep_state(n) for n in (N_SMALL, N_LARGE))

# ---------------------------------------------------------------------------
# Suite
# ---------------------------------------------------------------------------

const SUITE = BenchmarkGroup()

let g = addgroup!(SUITE, "fit")
    g["zach_default"] = @benchmarkable fit_ergm_count($ZACH, [SumTerm(), NonzeroTerm()];
                                                      warn=false)
    g["zach_maxval60"] = @benchmarkable fit_ergm_count($ZACH, [SumTerm(), NonzeroTerm()];
                                                       max_val=60, warn=false)
end

let g = addgroup!(SUITE, "derivatives")
    g["zach"] = @benchmarkable $ZACH_DERIV($ZACH_THETA)
end

let g = addgroup!(SUITE, "sweep")
    for n in (N_SMALL, N_LARGE)
        g["n$(n)"] = @benchmarkable one_sweep!($(STATES[n]))
    end
end

# ---------------------------------------------------------------------------
# Standalone entry point
# ---------------------------------------------------------------------------

function print_benchjl(results::BenchmarkGroup)
    for (path, trial) in BenchmarkTools.leaves(results)
        est = median(trial)
        println("BENCHJL\t", join(path, "/"), "\t",
                BenchmarkTools.time(est), "\t",
                BenchmarkTools.allocs(est), "\t",
                BenchmarkTools.memory(est))
    end
end

"Assert that a sweep scales with the number of dyads (O(degree) per dyad, not O(n))."
function assert_scaling(results::BenchmarkGroup)
    t_small = BenchmarkTools.time(median(results["sweep"]["n$(N_SMALL)"]))
    t_large = BenchmarkTools.time(median(results["sweep"]["n$(N_LARGE)"]))
    dyad_ratio = n_dyads(N_LARGE) / n_dyads(N_SMALL)
    ratio = (t_large / t_small) / dyad_ratio
    println("SCALING\tsweep\tn", N_LARGE, "/n", N_SMALL, " per dyad\t", round(ratio, digits=2))
    if ratio > SCALING_LIMIT
        println(stderr, "SCALING FAILURE: a Gibbs sweep at n=$(N_LARGE) costs ",
                round(ratio, digits=2), "x more PER DYAD than at n=$(N_SMALL) ",
                "(same mean degree; limit $(SCALING_LIMIT)x). The per-dyad ",
                "conditional is no longer O(degree).")
        return false
    end
    return true
end

function main()
    tune!(SUITE)
    results = run(SUITE; verbose=false, seconds=1)
    print_benchjl(results)
    assert_scaling(results) || exit(1)
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
