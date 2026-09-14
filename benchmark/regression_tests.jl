#!/usr/bin/env julia
# benchmark/regression_tests.jl — allocation- and complexity-regression
# assertions for ERGMCount.jl's hot loops. Standalone; run with
#     julia --project=benchmark benchmark/regression_tests.jl
# (CI runs it on one matrix cell; `Pkg.test()` guards the same pins from the
# testsets "Gibbs dyad update: 0 B per dyad, pinned", "MPLE design build
# allocates O(unique slabs), not O(dyads)" and "MPLE derivative evaluations
# allocate O(p²)".)
#
# What is asserted (panel 2026-09, item 25):
#   1. 0 B per Gibbs dyad update — the conditional folded from the support
#      profiles, the inverse-CDF draw and the typed weight snapshot, measured
#      as the WORST over every dyad of a warmed chain, mutation included:
#      exactly 0 for every update that does not insert an edge, and an
#      amortised per-insertion bound (Base's adjacency-vector growth) for
#      the insertions, over five whole sweeps at n=34 and at n=100;
#   2. the MPLE design build allocates O(unique slabs), not O(dyads): a
#      dyad-independent model builds ONE row at 12 and at 40 nodes for the
#      same bytes;
#   3. ≤ 512 B per derivative evaluation, independent of rows and support;
#   4. SCALING: a Gibbs sweep at n=100 costs at most 2× more per dyad than at
#      n=34 with the same mean degree (the profiles are O(degree)).

using ERGMCount
using Networks
using Random
using Test

include(joinpath(@__DIR__, "fixtures.jl"))

"""
Allocation of `_gibbs_update_dyad!` over `sweeps` sweeps of a warmed chain,
every update measured, as `(worst over updates that did not insert an edge,
bytes over insertions, number of insertions)`. Of Networks.jl's mutations only
an edge insertion can allocate (Base's `_growat!` reallocating a sorted
adjacency vector once its front slack is used up; amortised, independent of
the model), so the first must be exactly 0 and the second bounded per
insertion.
"""
function update_allocs(st; sweeps::Int=1)
    n = Int(nv(st.cur))
    ERGMCount._gibbs_update_dyad!(st.rng, st.cur, st.weights, SWEEP_TERMS, st.θ, 1, 2,
                                  st.support, st.log_h, st.η, st.buf)
    worst_noninsert = 0; insert_bytes = 0; n_insert = 0
    for _ in 1:sweeps, i in 1:n, j in 1:n
        i == j && continue
        old = ERGMCount.dyad_value(st.cur, st.weights, i, j)
        bytes = @allocated new = ERGMCount._gibbs_update_dyad!(
            st.rng, st.cur, st.weights, SWEEP_TERMS, st.θ, i, j, st.support, st.log_h,
            st.η, st.buf)
        if old == 0 && new > 0
            insert_bytes += bytes; n_insert += 1
        else
            worst_noninsert = max(worst_noninsert, bytes)
        end
    end
    return worst_noninsert, insert_bytes, n_insert
end

function design_rows_allocs(n, terms)
    net = sparse_count_network(Random.Xoshiro(21), n; mean_degree=4)
    weights = get_edge_attribute(net, :weight, Int)
    support = 0:10
    slab = zeros(length(terms) * length(support)); buf = zeros(length(support))
    rows = Dict{Vector{Float64}, Int}(); counts = Vector{Vector{Float64}}()
    ERGMCount._count_design_rows!(rows, counts, slab, buf, net, weights, terms,
                                  PoissonReference(), support)
    R = length(rows)
    empty!(rows); empty!(counts)
    a = @allocated ERGMCount._count_design_rows!(rows, counts, slab, buf, net, weights,
                                                 terms, PoissonReference(), support)
    return a, R
end

function derivative_allocs(net, max_val, terms)
    model = CountERGMModel(terms, net, PoissonReference())
    D = ERGMCount._count_design(model, 0:max_val)
    d = ERGMCount._count_derivatives(D, collect(1:length(terms)),
                                     fill(true, length(D.support), length(D.n_tot)))
    β = fill(0.1, length(terms))
    d(β)
    return @allocated d(β)
end

"Minimum over `reps` of the wall time of one sweep (the chain keeps moving)."
function min_sweep_time(st; reps::Int=5)
    one_sweep!(st)
    return minimum(@elapsed(one_sweep!(st)) for _ in 1:reps)
end

@testset "ERGMCount allocation and scaling regressions" begin
    @testset "0 B per Gibbs dyad update (warmed chain, mutation included)" begin
        for n in (34, 100)
            worst_noninsert, insert_bytes, n_insert = update_allocs(sweep_state(n); sweeps=5)
            @test n_insert > 0
            @test worst_noninsert == 0
            @test insert_bytes <= 128 * n_insert + 1024
        end
    end

    @testset "MPLE design build is O(unique slabs)" begin
        small, R_small = design_rows_allocs(12, (SumTerm(), NonzeroTerm()))
        big, R_big = design_rows_allocs(40, (SumTerm(), NonzeroTerm()))
        @test R_small == 1 && R_big == 1
        @test abs(big - small) <= 64
        dep_small, Rd_small = design_rows_allocs(12, (SumTerm(), NodeOSumTerm()))
        dep_big, Rd_big = design_rows_allocs(40, (SumTerm(), NodeOSumTerm()))
        @test Rd_big < n_dyads(40) && Rd_small < n_dyads(12)
        @test dep_big / Rd_big <= 2 * (dep_small / Rd_small) + 512
    end

    @testset "≤ 512 B per derivative evaluation" begin
        zach = zach_network()
        @test derivative_allocs(zach, 30, (SumTerm(), NonzeroTerm())) <= 512
        @test derivative_allocs(zach, 60, (SumTerm(), NonzeroTerm())) <= 512
        @test derivative_allocs(sparse_count_network(Random.Xoshiro(2), 40), 30,
                                (SumTerm(), NodeOSumTerm())) <= 512
    end

    @testset "SCALING: per-dyad sweep cost within 2x between n=34 and n=100" begin
        t34 = min_sweep_time(sweep_state(34))
        t100 = min_sweep_time(sweep_state(100))
        ratio = (t100 / t34) / (n_dyads(100) / n_dyads(34))
        println("SCALING\tsweep\tn100/n34 per dyad\t", round(ratio, digits=2))
        @test ratio <= 2.0
    end
end
