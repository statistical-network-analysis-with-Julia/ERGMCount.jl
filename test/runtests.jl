using ERGMCount
using ERGM
using Networks
using Graphs
using Distributions
using Random
using Statistics
using LinearAlgebra: norm
using Test

# Build a small count-valued network from a (i, j, w) list
function count_net(n, ties; directed=true)
    net = network(n; directed=directed)
    for (i, j, w) in ties
        add_edge!(net, i, j)
        set_edge_attribute!(net, :weight, i, j, w)
    end
    return net
end

# Set dyad (i,j) to count value y (0 removes the edge; a NEGATIVE value —
# admissible under DiscUnif2Reference(a < 0, b) — is an edge with that weight)
function set_dyad!(net, i, j, y)
    if y != 0
        has_edge(net, i, j) || add_edge!(net, i, j)
        set_edge_attribute!(net, :weight, i, j, y)
    elseif has_edge(net, i, j)
        rem_edge!(net, i, j)
    end
    return net
end

# Brute-force change statistic by actually editing a copy of the network
function brute_change_count(term, net, i, j, old, new)
    test = deepcopy(net)
    set_dyad!(test, i, j, old)
    s0 = compute(term, test)
    set_dyad!(test, i, j, new)
    s1 = compute(term, test)
    return s1 - s0
end

function random_count_net(n; directed=true, p=0.3, maxw=4, seed=1)
    rng = Random.Xoshiro(seed)
    net = network(n; directed=directed)
    for i in 1:n
        for j in (directed ? (1:n) : (i+1:n))
            i == j && continue
            if rand(rng) < p
                add_edge!(net, i, j)
                set_edge_attribute!(net, :weight, i, j, rand(rng, 1:maxw))
            end
        end
    end
    return net
end

const ALL_TERMS = [SumTerm(), NonzeroTerm(), GreaterthannTerm(2),
                   CountAtleastnTerm(2), CountMutualTerm(),
                   TransitiveTiesTerm(), CyclicalTiesTerm(),
                   NodeOSumTerm(), NodeISumTerm(), NodeSumTerm(),
                   # R-parity terms (WP4): the (min, max, min) weights, the
                   # other valued `mutual` forms and the threshold terms
                   TransitiveWeightsTerm(), CyclicalWeightsTerm(),
                   CountMutualTerm(:nabsdiff), CountMutualTerm(:geometric),
                   CountMutualTerm(:product), CountMutualTerm(:threshold; threshold=2),
                   SmallerthanTerm(2), EqualToTerm(3), InIntervalTerm(1, 3),
                   InIntervalTerm(1, 3; open=(false, false))]

@testset "ERGMCount.jl" begin
    @testset "Reference measures" begin
        pois = PoissonReference(2.0)
        # h(y) = λ^y / y!
        @test log_reference(pois, 0) ≈ 0.0
        @test log_reference(pois, 3) ≈ 3 * log(2.0) - log(6)

        # Geometric reference is the counting measure h(y) = 1
        geo = GeometricReference()
        @test log_reference(geo, 0) == 0.0
        @test log_reference(geo, 7) == 0.0

        # Binomial reference is h(y) = C(n, y) — no probability parameter
        bin = BinomialReference(5)
        @test log_reference(bin, 2) ≈ log(binomial(5, 2))
        @test log_reference(bin, 6) == -Inf

        du = DiscUnifReference(4)
        @test log_reference(du, 3) ≈ -log(5)

        du2 = DiscUnif2Reference(1, 3)
        @test log_reference(du2, 2) ≈ -log(3)

        @test_throws ArgumentError BinomialReference(0)
        @test_throws ArgumentError DiscUnifReference(-1)
        @test_throws ArgumentError DiscUnif2Reference(3, 1)
        # `PoissonReference(λ ≤ 0)` used to construct and fail later as a bare
        # `DomainError` from `log` (λ < 0) or newton_fit's "objective is not
        # finite" (λ = 0); the constructor names the reference and the rule
        for bad in (-1.0, 0.0, 0, NaN, Inf)
            err = try; PoissonReference(bad); nothing catch e; e end
            @test err isa ArgumentError
            @test occursin("PoissonReference: lambda must be positive", sprint(showerror, err))
        end
        @test PoissonReference(2).lambda === 2.0          # an Integer rate is accepted
        @test PoissonReference().lambda === 1.0
    end

    @testset "Term compute on known networks" begin
        # Directed triangle with counts: 1→2 (3), 2→3 (2), 1→3 (1), plus
        # reciprocal 2→1 (2)
        net = count_net(3, [(1, 2, 3), (2, 3, 2), (1, 3, 1), (2, 1, 2)])

        @test compute(SumTerm(), net) == 8.0
        @test compute(NonzeroTerm(), net) == 4.0
        @test compute(GreaterthannTerm(1), net) == 3.0
        @test compute(CountAtleastnTerm(2), net) == 3.0
        @test compute(CountMutualTerm(), net) == 2.0  # min(3,2) for dyad {1,2}

        # Transitive triple 1→2→3 with shortcut 1→3:
        # ordered distinct triples contribute min terms; brute check below
        # verifies exact value — here check a hand value:
        # triple (1,2,3): min(y_12, y_23, y_13) = min(3,2,1) = 1
        # triple (2,1,3): min(y_21, y_13, y_23) = min(2,1,2) = 1
        # all other ordered triples have a zero dyad → 0
        @test compute(TransitiveTiesTerm(), net) == 2.0

        # No directed 3-cycle present
        @test compute(CyclicalTiesTerm(), net) == 0.0

        # Add 3→1 (2) to close the cycle 1→2→3→1: min(3,2,2)=2, once
        set_dyad!(net, 3, 1, 2)
        @test compute(CyclicalTiesTerm(), net) == 2.0

        # Node strengths: out: 1: 3+1=4, 2: 2+2=4, 3: 2 → 16+16+4=36
        @test compute(NodeOSumTerm(), net) == 36.0
        # in: 1: 2+2=4, 2: 3, 3: 2+1=3 → 16+9+9=34
        @test compute(NodeISumTerm(), net) == 34.0

        # Edges without a weight attribute count as 1
        bare = network(3)
        add_edge!(bare, 1, 2)
        @test compute(SumTerm(), bare) == 1.0
    end

    @testset "change_stat_count matches brute force" begin
        for directed in (true, false)
            net = random_count_net(7; directed=directed, seed=directed ? 3 : 4)
            weights = get_edge_attribute(net, :weight)
            n = nv(net)

            for term in ALL_TERMS
                for i in 1:n, j in 1:n
                    i == j && continue
                    !directed && j < i && continue
                    old = ERGMCount.dyad_value(net, weights, i, j)
                    for new in (0, 1, 3)
                        expected = brute_change_count(term, net, i, j, old, new)
                        actual = change_stat_count(term, net, weights, i, j, old, new)
                        @test actual ≈ expected atol = 1e-9
                    end
                end
            end
        end
    end

    @testset "MPLE recovers Poisson rate (Sum-only model)" begin
        # With Poisson(λ) reference and only a SumTerm, dyads are iid
        # Poisson(λ·e^θ); MPLE should give θ̂ ≈ log(ȳ/λ)
        rng = Random.Xoshiro(2026)
        n = 14
        μ = 2.0
        net = network(n)
        for i in 1:n, j in 1:n
            i == j && continue
            y = rand(rng, Poisson(μ))
            y > 0 || continue
            add_edge!(net, i, j)
            set_edge_attribute!(net, :weight, i, j, y)
        end

        result = ergm_count(net, [SumTerm()]; reference=PoissonReference(1.0),
                            max_val=30)
        ȳ = compute(SumTerm(), net) / (n * (n - 1))

        @test result.converged
        @test result.coefficients[1] ≈ log(ȳ) atol = 1e-3
        @test isfinite(result.loglik)
        @test result.std_errors[1] > 0
    end

    @testset "Multi-term estimation runs" begin
        net = random_count_net(8; seed=9)
        # This exact model formerly threw MethodError
        result = ergm_count(net, [SumTerm(), NonzeroTerm(), CountMutualTerm()])
        @test result isa CountERGMResult
        @test length(result.coefficients) == 3
        @test isfinite(result.loglik)
        @test all(isfinite, result.coefficients)
    end

    # ------------------------------------------------------------------
    # Allocation regressions on the count MPLE (panel 2026-09, item 25).
    #
    # The design is COMPRESSED: dyads with identical (terms × support) change-
    # statistic slabs share one row, so a dyad-independent model has ONE row
    # whatever the network size, and the per-dyad sweep writes into a reused
    # buffer (only a NEW slab is copied). The derivative closure fills its
    # conditional moments in place on workspaces allocated once, so an
    # evaluation allocates only the gradient and Hessian it hands to
    # `newton_fit`, independent of dyads and support.
    # ------------------------------------------------------------------
    @testset "MPLE design build allocates O(unique slabs), not O(dyads)" begin
        function rows_allocs(n, terms)
            net = random_count_net(n; seed=21)
            weights = get_edge_attribute(net, :weight, Int)
            support = 0:10
            slab = zeros(length(terms) * length(support))
            buf = zeros(length(support))
            rows = Dict{Vector{Float64}, Int}()
            counts = Vector{Vector{Float64}}()
            ERGMCount._count_design_rows!(rows, counts, slab, buf, net, weights, terms,
                                          PoissonReference(), support)
            R = length(rows)
            empty!(rows); empty!(counts)
            a = @allocated ERGMCount._count_design_rows!(rows, counts, slab, buf, net,
                                                         weights, terms,
                                                         PoissonReference(), support)
            return a, R
        end
        indep = (SumTerm(), NonzeroTerm())
        small, R_small = rows_allocs(12, indep)    # 132 dyads
        big, R_big = rows_allocs(40, indep)        # 1560 dyads
        @test R_small == 1 && R_big == 1           # one slab: dyad-independent
        @test small <= 1024
        @test abs(big - small) <= 64                # 12x the dyads, same bytes (± Dict bookkeeping)
        # A dyad-dependent term gives one slab per distinct neighbourhood, and
        # the allocation is proportional to the slabs, not to the dyads
        dep_small, Rd_small = rows_allocs(12, (SumTerm(), NodeOSumTerm()))
        dep_big, Rd_big = rows_allocs(40, (SumTerm(), NodeOSumTerm()))
        # one slab per distinct (out-strength minus own value); far fewer than dyads
        @test Rd_big < 1560 && Rd_small < 132
        @test dep_big / Rd_big <= 2 * (dep_small / Rd_small) + 512

        # The per-dyad slab fill itself — every term's change statistic over
        # the whole support, written into the reused buffer — is allocation-free
        # for every term, dyad-dependent ones included
        net = random_count_net(10; seed=21)
        weights = get_edge_attribute(net, :weight, Int)
        # measured inside a function specialised on the term tuple (a dynamic
        # call from a loop over differently-typed tuples would box `0:10`)
        function fill_allocs(terms, net, weights)
            slab = zeros(length(terms) * 11)
            buf = zeros(11)
            ERGMCount._fill_slab!(slab, buf, terms, net, weights, 2, 3, 0:10)
            return (@allocated ERGMCount._fill_slab!(slab, buf, terms, net, weights, 2, 3, 0:10)), slab
        end
        for terms in ((SumTerm(), NonzeroTerm()), (ALL_TERMS...,))
            bytes, slab = fill_allocs(terms, net, weights)
            @test bytes == 0
            # ... and it agrees with the change statistics it is made of
            for (sidx, yv) in enumerate(0:10), (k, t) in enumerate(terms)
                @test slab[(sidx - 1) * length(terms) + k] ==
                      change_stat_count(t, net, weights, 2, 3, 0, yv)
            end
        end
    end

    @testset "MPLE derivative evaluations allocate O(p²), not O(rows · |support| · p²)" begin
        function evaluation_allocs(n, max_val, terms)
            net = random_count_net(n; seed=21)
            model = CountERGMModel(terms, net, PoissonReference())
            D = ERGMCount._count_design(model, 0:max_val)
            mask = fill(true, length(D.support), length(D.n_tot))
            cols = collect(1:length(terms))
            d = ERGMCount._count_derivatives(D, cols, mask)
            β = fill(0.1, length(terms))
            d(β)                    # warm up: @allocated on a first call
            return @allocated d(β)  # would measure compilation
        end
        small = evaluation_allocs(6, 5, (SumTerm(), NodeOSumTerm()))    # 30 rows, 6 values
        big = evaluation_allocs(20, 30, (SumTerm(), NodeOSumTerm()))    # ≤380 rows, 31 values
        @test small <= 512
        @test big <= 512
        @test big <= small + 64
    end

    @testset "Support profiles equal the per-value change statistics" begin
        # `change_stats_support!` is the O(degree + |support|) profile both hot
        # paths consume; `change_stat_count` (brute-force-tested above) stays
        # the definition. Every term, both directednesses, every support value,
        # every starting value, supports that do and do not start at 0.
        for directed in (true, false), seed in (1, 2), support in (0:6, 0:1, 1:5, 3:9)
            net = random_count_net(9; directed=directed, seed=seed, maxw=7)
            weights = get_edge_attribute(net, :weight, Int)
            dest = zeros(length(support))
            for t in ALL_TERMS, i in 1:9, j in 1:9
                (i == j || (!directed && i > j)) && continue
                for old in 0:7
                    ERGMCount.change_stats_support!(dest, t, net, weights, i, j, old, support)
                    @test dest == [change_stat_count(t, net, weights, i, j, old, y)
                                   for y in support]
                end
            end
        end
        # ... and they are allocation-free for every term (the specialised
        # ones walk neighbour lists, never a temporary)
        function profile_allocs(t, net, weights, dest)
            ERGMCount.change_stats_support!(dest, t, net, weights, 2, 3, 1, 0:10)
            return @allocated ERGMCount.change_stats_support!(dest, t, net, weights, 2, 3, 1, 0:10)
        end
        for directed in (true, false)
            net = random_count_net(30; directed=directed, seed=4)
            weights = get_edge_attribute(net, :weight, Int)
            dest = zeros(11)
            for t in ALL_TERMS
                @test profile_allocs(t, net, weights, dest) == 0
            end
        end
        # The profile is what the static fold accumulates into the conditional
        net = random_count_net(6; seed=5)
        weights = get_edge_attribute(net, :weight, Int)
        terms = Tuple(ALL_TERMS)
        θ = 0.1 .* (1:length(ALL_TERMS))
        η = zeros(4); buf = zeros(4)
        ERGMCount._accumulate_conditional!(η, buf, terms, θ, net, weights, 1, 2, 0, 0:3)
        @test η ≈ [sum(θ[k] * change_stat_count(t, net, weights, 1, 2, 0, y)
                       for (k, t) in enumerate(terms)) for y in 0:3]
    end

    @testset "Gibbs dyad update: 0 B per dyad, pinned" begin
        # The per-dyad kernel `_gibbs_update_dyad!` — the conditional folded
        # statically over the term tuple from the support profiles, the
        # inverse-CDF draw on the un-normalised weights, the typed weight
        # snapshot — allocates nothing. Only an actual change of a dyad's
        # value touches the network, and of Networks.jl's own mutations only
        # an edge INSERTION can allocate: Base's `_growat!` reallocates a
        # sorted adjacency vector on a first-half `insert!` once its front
        # slack is used up (a few hundred bytes, amortised over many
        # insertions, independent of the model). Changing an existing edge's
        # value and removing an edge never allocate. So: every update that is
        # not an insertion is exactly 0 B, and insertions cost a bounded
        # amortised amount, at 12 and at 40 nodes alike.
        # Function barrier: `Tuple(ALL_TERMS)` has a concrete type at run time
        # but not at inference time in the caller, exactly as `model.terms`
        # does once it reaches the sampler through dispatch.
        support = 0:10
        log_h = [log_reference(PoissonReference(), y) for y in support]
        # The state is warmed first — two sweeps at a dense specification
        # (sum coefficient 2, mean count ≈ 7) fill both weight dictionaries
        # with every key they will need (containers never shrink), then three
        # sweeps at the target coefficients. The warmed chain AND its snapshot
        # are what the pins run on: a `copy` would hand back exact-capacity
        # containers and a snapshot without the zeroed keys
        function grown(net, terms, θ)
            cur = copy(net); weights = get_edge_attribute(cur, :weight, Int)
            η = zeros(11); buf = zeros(11); rng = Xoshiro(2)
            dense = [k == 1 ? 2.0 : 0.0 for k in eachindex(θ)]
            for _ in 1:2
                ERGMCount._gibbs_sweep!(rng, cur, weights, terms, dense, support, log_h, η, buf)
            end
            for _ in 1:3
                ERGMCount._gibbs_sweep!(rng, cur, weights, terms, θ, support, log_h, η, buf)
            end
            return cur, weights
        end
        # `sweeps` sweeps of per-dyad updates on a warmed chain, every update
        # measured: (worst over the updates that did not insert an edge, bytes
        # over insertions, number of insertions, conditional, draw)
        function chain_allocs(terms::Tuple, (cur, weights), θ; sweeps::Int=1)
            η = zeros(length(support)); buf = zeros(length(support))
            rng = Xoshiro(1)
            ERGMCount._gibbs_update_dyad!(rng, cur, weights, terms, θ, 1, 2, support,
                                          log_h, η, buf)
            worst_noninsert = 0; insert_bytes = 0; n_insert = 0
            for _ in 1:sweeps, i in 1:nv(cur), j in 1:nv(cur)
                (i == j || (!is_directed(cur) && i > j)) && continue
                old = ERGMCount.dyad_value(cur, weights, i, j)
                bytes = @allocated new = ERGMCount._gibbs_update_dyad!(
                    rng, cur, weights, terms, θ, i, j, support, log_h, η, buf)
                if old == 0 && new > 0
                    insert_bytes += bytes; n_insert += 1
                else
                    worst_noninsert = max(worst_noninsert, bytes)
                end
            end
            ERGMCount._dyad_conditional!(η, buf, terms, θ, log_h, support, cur, weights, 1, 2, 0)
            a = @allocated ERGMCount._dyad_conditional!(η, buf, terms, θ, log_h, support,
                                                        cur, weights, 1, 2, 0)
            total = sum(η)
            ERGMCount._draw_index(rng, η, total)
            b = @allocated ERGMCount._draw_index(rng, η, total)
            return worst_noninsert, insert_bytes, n_insert, a, b
        end
        four_d = (SumTerm(), NonzeroTerm(), CountMutualTerm(), NodeOSumTerm())
        four_u = (SumTerm(), NonzeroTerm(), TransitiveTiesTerm(), NodeSumTerm())
        θ4 = [0.2, -0.3, 0.4, -0.02]
        # Every term of the package, both directednesses: the kernel is 0 B
        for (terms, θ, directed) in ((four_d, θ4, true), (four_u, θ4, false),
                                     (Tuple(ALL_TERMS), fill(0.05, length(ALL_TERMS)), true))
            state = grown(random_count_net(10; directed=directed, seed=5), terms, θ)
            worst_noninsert, insert_bytes, n_insert, a, b = chain_allocs(terms, state, θ)
            @test a == 0                                # the conditional
            @test b == 0                                # the draw
            @test worst_noninsert == 0                  # unchanged / changed / removed
            @test insert_bytes <= 128 * n_insert + 1024 # insertions: amortised growth
        end

        # A whole sweep does not grow with the number of dyads: on a warmed
        # sparse chain (mean out-degree ≈ 4, the regime of `benchmark/`), ten
        # sweeps at 12 nodes (132 dyads) and at 40 nodes (1560 dyads) — with
        # hundreds of insertions each — allocate nothing outside insertions,
        # and the insertions stay under the same per-insertion bound
        # (measured 0–35 B per insertion)
        function sparse_seed(n)
            rng = Xoshiro(21); net = network(n; directed=true)
            for i in 1:n, j in 1:n
                i != j && rand(rng) < 4 / n || continue
                add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, rand(rng, 1:3))
            end
            return net
        end
        for n in (12, 40)
            θ = [log(4 / n), 0.0, 0.3, -0.01]           # mean out-degree ≈ 4 at any n
            state = grown(sparse_seed(n), four_d, θ)
            worst_noninsert, insert_bytes, n_insert, a, b = chain_allocs(four_d, state, θ;
                                                                          sweeps=10)
            @test n_insert > 0
            @test worst_noninsert == 0
            @test insert_bytes <= 128 * n_insert + 1024
        end
        # A sweep in which no dyad changes value is exactly 0 B at any size,
        # without any warm-up
        function frozen_sweep_bytes(n)
            cur = network(n; directed=true)
            weights = get_edge_attribute(cur, :weight, Int)
            η = zeros(11); buf = zeros(11); rng = Xoshiro(3)
            θ = [-40.0, 0.0, 0.0, 0.0]              # every draw is 0 on an empty network
            ERGMCount._gibbs_sweep!(rng, cur, weights, four_d, θ, support, log_h, η, buf)
            return @allocated ERGMCount._gibbs_sweep!(rng, cur, weights, four_d, θ, support,
                                                      log_h, η, buf)
        end
        @test frozen_sweep_bytes(12) == 0
        @test frozen_sweep_bytes(40) == 0

        # The inverse-CDF draw is a correct sampler on un-normalised weights:
        # frequencies match the probabilities it is handed
        # (distinct names from the closures' locals: a name assigned in both
        # this scope and a closure would be captured and boxed there, and the
        # boxed argument would turn the measured call dynamic)
        draw_rng = Xoshiro(11)
        draw_w = [0.4, 1.0, 0.6]                     # sums to 2: p = 0.2, 0.5, 0.3
        draws = [ERGMCount._draw_index(draw_rng, draw_w, 2.0) for _ in 1:20_000]
        @test count(==(2), draws) / 20_000 ≈ 0.5 atol = 0.02
        @test count(==(3), draws) / 20_000 ≈ 0.3 atol = 0.02
        # A degenerate conditional (one weight) always returns its index, and
        # the last index absorbs rounding
        @test ERGMCount._draw_index(draw_rng, [0.0, 0.0, 1.0], 1.0) == 3
        @test ERGMCount._draw_index(draw_rng, [1.0, 0.0, 0.0], 1.0) == 1
        @test ERGMCount._draw_index(draw_rng, [0.3, 0.3], 0.6 + 1e-12) in (1, 2)
    end

    @testset "Simulation targets the model" begin
        Random.seed!(7)
        n = 8
        seed_net = network(n)

        # Sum-only model with Poisson(1) reference and θ = log 2:
        # dyads iid Poisson(2) (truncated); check the mean
        sims = simulate_count_ergm(seed_net, [SumTerm()], [log(2.0)];
                                   reference=PoissonReference(1.0),
                                   n_sim=30, burnin=20, interval=2, max_val=15)
        @test length(sims) == 30
        mean_sum = mean(compute(SumTerm(), s) for s in sims)
        @test mean_sum ≈ 2.0 * n * (n - 1) rtol = 0.15

        # Structural term influences draws: positive mutual coefficient
        # yields more reciprocity than the zero-coefficient model
        sims_mut = simulate_count_ergm(seed_net, [SumTerm(), CountMutualTerm()],
                                       [log(0.5), 1.5];
                                       reference=PoissonReference(1.0),
                                       n_sim=20, burnin=20, interval=2, max_val=10)
        sims_ind = simulate_count_ergm(seed_net, [SumTerm(), CountMutualTerm()],
                                       [log(0.5), 0.0];
                                       reference=PoissonReference(1.0),
                                       n_sim=20, burnin=20, interval=2, max_val=10)
        mut_pos = mean(compute(CountMutualTerm(), s) for s in sims_mut)
        mut_zero = mean(compute(CountMutualTerm(), s) for s in sims_ind)
        @test mut_pos > mut_zero
    end

    @testset "rng keyword gives reproducible draws" begin
        seed_net = network(6)
        sim = rng -> simulate_count_ergm(seed_net, [SumTerm(), CountMutualTerm()],
                                         [log(1.5), 0.5];
                                         reference=PoissonReference(1.0),
                                         n_sim=3, burnin=10, interval=2,
                                         max_val=10, rng=rng)
        sims1 = sim(Random.Xoshiro(99))
        sims2 = sim(Random.Xoshiro(99))
        @test [compute(SumTerm(), s) for s in sims1] ==
              [compute(SumTerm(), s) for s in sims2]
        @test all(get_edge_attribute(a, :weight) == get_edge_attribute(b, :weight)
                  for (a, b) in zip(sims1, sims2))

        # The fitted-result method accepts rng too
        net = random_count_net(6; seed=17)
        result = ergm_count(net, [SumTerm()])
        r1 = simulate_count_ergm(result; n_sim=2, burnin=5, interval=2,
                                 rng=Random.Xoshiro(4))
        r2 = simulate_count_ergm(result; n_sim=2, burnin=5, interval=2,
                                 rng=Random.Xoshiro(4))
        @test [compute(SumTerm(), s) for s in r1] ==
              [compute(SumTerm(), s) for s in r2]

        # Terms as a Tuple or a single term simulate identically
        t1 = simulate_count_ergm(seed_net, (SumTerm(), CountMutualTerm()), [log(1.5), 0.5];
                                 reference=PoissonReference(1.0), n_sim=3, burnin=10,
                                 interval=2, max_val=10, rng=Random.Xoshiro(99))
        @test [compute(SumTerm(), s) for s in t1] == [compute(SumTerm(), s) for s in sims1]
        s1 = simulate_count_ergm(seed_net, SumTerm(), [log(1.5)]; n_sim=2, burnin=5,
                                 interval=1, max_val=10, rng=Random.Xoshiro(1))
        s2 = simulate_count_ergm(seed_net, [SumTerm()], [log(1.5)]; n_sim=2, burnin=5,
                                 interval=1, max_val=10, rng=Random.Xoshiro(1))
        @test [compute(SumTerm(), s) for s in s1] == [compute(SumTerm(), s) for s in s2]
        # Retained draws are independent copies of the chain, not views of it
        @test sims1[1] !== sims1[2]
        add_edge!(sims1[1], 1, 2); set_edge_attribute!(sims1[1], :weight, 1, 2, 99)
        @test compute(SumTerm(), sims1[2]) == compute(SumTerm(), sims2[2])
        # Mismatched coefficients and non-finite ones are refused
        @test_throws ArgumentError simulate_count_ergm(seed_net, [SumTerm()], [0.1, 0.2])
        @test_throws ArgumentError simulate_count_ergm(seed_net, [SumTerm()], [-Inf])

        # sample_reference draws flow through the rng keyword
        @test sample_reference(PoissonReference(2.0); rng=Random.Xoshiro(3)) ==
              sample_reference(PoissonReference(2.0); rng=Random.Xoshiro(3))
        @test sample_reference(DiscUnif2Reference(1, 5); rng=Random.Xoshiro(8)) ==
              sample_reference(DiscUnif2Reference(1, 5); rng=Random.Xoshiro(8))
    end

    @testset "Estimation-simulation round trip" begin
        Random.seed!(21)
        n = 8
        seed_net = network(n)
        θ_true = log(1.5)
        sims = simulate_count_ergm(seed_net, [SumTerm()], [θ_true];
                                   reference=PoissonReference(1.0),
                                   n_sim=1, burnin=50, interval=1, max_val=15)
        result = ergm_count(sims[1], [SumTerm()];
                            reference=PoissonReference(1.0), max_val=15)
        @test result.converged
        @test result.coefficients[1] ≈ θ_true atol = 0.35
    end

    @testset "StatsAPI surface (item 15)" begin
        net = random_count_net(6; seed=11)
        result = ergm_count(net, [SumTerm(), NonzeroTerm()])
        @test coef(result) == result.coefficients
        @test stderror(result) == result.std_errors
        @test vcov(result) == result.vcov
        @test size(vcov(result)) == (2, 2)
        @test all(stderror(result) .≈
                  sqrt.(abs.([vcov(result)[k, k] for k in 1:2])))
        @test loglikelihood(result) == result.loglik
        @test nobs(result) == 6 * 5  # directed dyads
        @test dof(result) == 2

        # The ONE checker every model package calls: the full ten-verb surface,
        # on a Poisson fit, a bounded-reference fit and a bootstrap fit
        @test Networks.check_statsapi(result; strict=true) !== nothing
        bin = ergm_count(net, [SumTerm(), NonzeroTerm()]; reference=BinomialReference(4))
        @test Networks.check_statsapi(bin; strict=true) !== nothing
        boot = ergm_count(net, [SumTerm(), NonzeroTerm()]; se=:bootstrap, n_boot=10,
                          rng=Xoshiro(5))
        @test Networks.check_statsapi(boot; strict=true) !== nothing
        @test all(Networks.check_statsapi(result))

        # Pseudo-likelihood AIC/BIC from the same loglik/dof/nobs
        @test aic(result) ≈ -2 * loglikelihood(result) + 2 * dof(result)
        @test bic(result) ≈ -2 * loglikelihood(result) + dof(result) * log(nobs(result))
        # Normal-theory intervals from the reported SEs
        ci = confint(result)
        @test size(ci) == (2, 2)
        @test all(ci[:, 1] .< coef(result) .< ci[:, 2])
        @test ci ≈ hcat(coef(result) .- 1.959963984540054 .* stderror(result),
                        coef(result) .+ 1.959963984540054 .* stderror(result))
        ci90 = confint(result; level=0.9)
        @test all(ci90[:, 2] .- ci90[:, 1] .< ci[:, 2] .- ci[:, 1])
        @test_throws ArgumentError confint(result; level=1.5)
        # The coefficient table, labelled with the shared `name(term, net)`
        tbl = coeftable(result)
        @test tbl isa Networks.CoefficientTable
        @test tbl.names == [name(t, net) for t in result.model.terms] == ["sum", "nonzero"]
        @test tbl["sum"].estimate == coef(result)[1]
        @test tbl["nonzero"].std_error == stderror(result)[2]
        @test tbl[1].z_value == result.z_values[1]
        @test tbl[2].p_value == result.p_values[2]
        @test tbl.p_values == Networks.z_pvalues(coef(result), stderror(result)).p
        # Co-loading: the verbs are the same bindings as ERGM.jl's / StatsAPI's
        @test ERGMCount.coef === ERGM.coef
        @test ERGMCount.aic === ERGM.aic
        @test ERGMCount.coeftable === ERGM.coeftable
        @test ERGMCount.gof === ERGM.gof
    end

    @testset "CountERGMModel{T,D}: concrete, tuple-backed, no directed field" begin
        for directed in (true, false)
            net = random_count_net(6; directed=directed, seed=5)
            terms = directed ? [SumTerm(), NonzeroTerm(), CountMutualTerm()] :
                               [SumTerm(), NonzeroTerm(), NodeSumTerm()]
            model = CountERGMModel(terms, net, PoissonReference())
            @test model isa CountERGMModel{Int, directed}
            @test is_directed(model) == is_directed(net) == directed
            @test isconcretetype(fieldtype(typeof(model), :network))
            @test isconcretetype(fieldtype(typeof(model), :terms))
            @test isconcretetype(fieldtype(typeof(model), :reference))
            @test model.terms isa Tuple
            @test model.terms == Tuple(terms)
            @test !hasfield(typeof(model), :directed)
            # Tuple, Vector and single-term constructors agree; the reference
            # defaults to Poisson
            @test CountERGMModel(Tuple(terms), net, PoissonReference()).terms == model.terms
            @test CountERGMModel(SumTerm(), net).terms == (SumTerm(),)
            @test CountERGMModel(terms, net).reference == PoissonReference()
            # The result is parameterised on the model type (mirrors ERGMResult)
            fit = ERGMCount.count_mple(model)
            @test fit isa CountERGMResult{typeof(model)}
            @test fit.model === model
        end
        # The static folds over the tuple are @generated per term count, so a
        # model wider than Base's 32-element `map` unrolling limit still
        # allocates nothing per dyad (36 terms here)
        net = random_count_net(6; seed=5)
        weights = get_edge_attribute(net, :weight, Int)
        wide = Tuple(repeat(ALL_TERMS, 4)[1:36])
        @test length(wide) == 36
        function wide_allocs(terms, net, weights)
            θ = fill(0.01, length(terms))
            η = zeros(4); buf = zeros(4); slab = zeros(4 * length(terms))
            ERGMCount._accumulate_conditional!(η, buf, terms, θ, net, weights, 1, 2, 0, 0:3)
            ERGMCount._fill_slab!(slab, buf, terms, net, weights, 1, 2, 0:3)
            return (@allocated ERGMCount._accumulate_conditional!(η, buf, terms, θ, net,
                                                                  weights, 1, 2, 0, 0:3)),
                   (@allocated ERGMCount._fill_slab!(slab, buf, terms, net, weights, 1, 2, 0:3))
        end
        @test wide_allocs(wide, net, weights) == (0, 0)
        @test !isdefined(ERGMCount, :_weighted_change)   # replaced by the profile fold
        # Validation
        @test_throws ArgumentError CountERGMModel((), net)
        @test_throws ArgumentError CountERGMModel((SumTerm(), 1), net)
        @test_throws ArgumentError CountERGMModel(AbstractERGMTerm[], net)
        @test !isdefined(ERGMCount, :_dyads)     # the dyad-pair vector is gone
    end

    @testset "compute is invariant under Base.copy(::Network)" begin
        # Base.copy preserves edge attributes, so every count statistic
        # must agree between a network and its copy
        for directed in (true, false)
            net = random_count_net(6; directed=directed, seed=5)
            for term in ALL_TERMS
                @test compute(term, copy(net)) == compute(term, net)
            end
        end
    end

    @testset "fit alias and API" begin
        @test fit_count_ergm === ergm_count
        @test fit_ergm_count === ergm_count
        net = random_count_net(5; seed=13)

        # A single term and a Tuple of terms are accepted and give the same fit
        fv = fit_ergm_count(net, [SumTerm(), NonzeroTerm()])
        ft = fit_ergm_count(net, (SumTerm(), NonzeroTerm()))
        @test coef(ft) == coef(fv)
        @test coef(fit_ergm_count(net, SumTerm())) == coef(fit_ergm_count(net, [SumTerm()]))
        @test fv.model.terms isa Tuple

        # method=:mcmle (or anything but :mple) is refused with the reason
        for m in (:mcmle, :mcmc)
            err = try; ergm_count(net, [SumTerm()]; method=m); nothing
                  catch e; e end
            @test err isa ArgumentError
            msg = sprint(showerror, err)
            @test occursin("maximum pseudo-likelihood only", msg)
            @test occursin("exact MLE for a dyad-independent model", msg)
            # ... naming EVERY dyad-independent term, not the first four
            @test occursin("sum, nonzero, greaterthan, atleast, smallerthan, equalto, ininterval", msg)
            @test occursin("se=:bootstrap", msg)
            @test occursin("not implemented", msg)
        end

        # Swapped arguments name the right order instead of a MethodError
        err = try; fit_ergm_count([SumTerm()], net); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("the network comes first", sprint(showerror, err))
        err = try; fit_ergm_count(SumTerm(), net); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("fit_ergm_count(net, [SumTerm()])", sprint(showerror, err))

        # An empty term list is refused
        @test_throws ArgumentError fit_ergm_count(net, AbstractERGMTerm[])

        # Observed count outside the support: the message branches on whether
        # the bound is the MODEL's (bounded reference: say so, name the
        # reference and the dyad, do not advise `max_val`, which it ignores)
        # or OURS (unbounded reference: advise `max_val`).
        bnet = count_net(4, [(1, 2, 5), (2, 3, 1), (3, 4, 2)])
        err = try; fit_ergm_count(bnet, [SumTerm()]; reference=BinomialReference(3)); nothing
              catch e; e end
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("Observed count 5 at dyad (1,2)", msg)
        @test occursin("part of the model", msg)
        @test occursin("BinomialReference", msg)
        @test occursin("BinomialReference(5)", msg)      # the fix it suggests
        @test !occursin("max_val", msg)
        err = try; fit_ergm_count(bnet, [SumTerm()]; reference=DiscUnifReference(1)); nothing
              catch e; e end
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("part of the model", msg)
        @test occursin("DiscUnifReference", msg)
        @test occursin("dyad (1,2)", msg)
        @test !occursin("max_val", msg)
        err = try; fit_ergm_count(bnet, [SumTerm()]; reference=DiscUnif2Reference(0, 1)); nothing
              catch e; e end
        @test occursin("DiscUnif2Reference(0, 5)", sprint(showerror, err))
        # ... and the truncating branch advises the bound we chose
        err = try; fit_ergm_count(bnet, [SumTerm()]; max_val=2); nothing catch e; e end
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("truncation you chose", msg)
        @test occursin("`max_val` ≥ 5", msg)

        # A directed-only term on an undirected network is refused with a hint
        unet = count_net(4, [(1, 2, 2), (2, 3, 1)]; directed=false)
        for (t, hint) in ((CountMutualTerm(), "drop the term"),
                          (CyclicalTiesTerm(), "drop the term"),
                          (NodeOSumTerm(), "NodeSumTerm()"),
                          (NodeISumTerm(), "NodeSumTerm()"))
            err = try; fit_ergm_count(unet, [SumTerm(), t]); nothing catch e; e end
            @test err isa ArgumentError
            msg = sprint(showerror, err)
            @test occursin("requires a directed network", msg)
            @test occursin(hint, msg)
        end
        @test fit_ergm_count(unet, [SumTerm(), NodeSumTerm()]) isa CountERGMResult
        @test fit_ergm_count(unet, [SumTerm(), TransitiveTiesTerm()]) isa CountERGMResult
    end

    @testset "Shared contracts are imported, not re-implemented" begin
        # ONE z → p helper, ONE Newton kernel, ONE se validator, ONE dependence
        # predicate (panel 2026-09, items 13/14/28)
        @test ERGMCount.z_pvalues === Networks.z_pvalues
        @test ERGMCount.newton_fit === Networks.newton_fit
        @test ERGMCount.check_se === Networks.check_se
        @test ERGMCount.has_dyad_dependent === ERGM.has_dyad_dependent
        @test ERGMCount.CoefficientTable === Networks.CoefficientTable
        @test ERGMCount.coeftable === Networks.coeftable
        @test !isdefined(ERGMCount, :_z_pvalues)
        @test !isdefined(ERGMCount, :_has_dyad_dependent)
        @test hasmethod(ERGM.has_dyad_dependent, Tuple{CountERGMModel})
        @test hasmethod(ERGM.requires_directed, Tuple{NodeOSumTerm})
        @test ERGMCount.count_mple isa Function
        # No StatsBase: the Gibbs draw is an inline inverse-CDF loop
        @test !isdefined(ERGMCount, :StatsBase)
        @test !isdefined(ERGMCount, :Weights)

        # Every cross-package `Pkg._name` reach-in targets a PUBLIC binding
        # (panel item 13): the one such reference in the source is the
        # sampler-defaults rule, and ERGM.jl declares it `public` — dropping
        # that declaration would go red here, not silently keep working
        src = read(joinpath(dirname(@__DIR__), "src", "ERGMCount.jl"), String)
        code = join(filter(l -> !startswith(strip(l), "#"), split(src, '\n')), '\n')
        reach = unique(String(m.match) for m in eachmatch(r"\b(?:ERGM|Networks)\._\w+", code))
        @test reach == ["ERGM._mcmc_defaults"]
        @test Base.ispublic(ERGM, :_mcmc_defaults)
        @test all(Base.ispublic(ERGM, Symbol(r[6:end])) for r in reach if startswith(r, "ERGM."))
        # ... and this package's own documented names are public bindings:
        # `count_mple` is exported as ERGM.jl exports `mple` (a warning tells
        # the user to call `count_mple(model; max_val = …)`, so it must resolve
        # as written), `dyad_value` too, `change_stats_support!` is `public`
        @test Base.isexported(ERGMCount, :count_mple)
        @test Base.isexported(ERGMCount, :dyad_value)
        @test Base.ispublic(ERGMCount, :count_mple)
        @test Base.ispublic(ERGMCount, :dyad_value)
        @test Base.ispublic(ERGMCount, :change_stats_support!)
        @test !Base.isexported(ERGMCount, :change_stats_support!)
        @test count_mple === ERGMCount.count_mple
        # Co-loading leaves the name free of conflicts: neither ERGM nor
        # Networks exports a `count_mple` or `dyad_value`
        @test !Base.isexported(ERGM, :count_mple) && !Base.isexported(Networks, :dyad_value)

        # The shared `se=` validator's message shape
        net = random_count_net(5; seed=13)
        err = try; fit_ergm_count(net, [SumTerm()]; se=:sandwich); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("count_mple: se must be one of", sprint(showerror, err))
        @test occursin(":sandwich", sprint(showerror, err))
    end

    @testset "show renders the shared coefficient table" begin
        net = random_count_net(6; seed=11)
        result = fit_ergm_count(net, [SumTerm(), NonzeroTerm()])
        out = sprint(show, result)
        @test occursin("Count ERGM Results", out)
        @test occursin("Estimate", out)
        @test occursin("Pr(>|z|)", out)
        @test occursin("Signif. codes", out)
        @test occursin("sum", out)
        @test occursin("nonzero", out)
        @test occursin("AIC:", out)
        @test occursin("pseudo-likelihood", out)
        @test occursin("Coefficients:", out)
        # The printed table IS `coeftable(result)`
        @test occursin(sprint(show, coeftable(result)), out)
    end

    @testset "p-values are floored, never exactly 0.0" begin
        # Two coefficients with |z| ≈ 14 and 18 (the zach fixture) underflow the
        # naive 2(1 − Φ(|z|)); the shared `z_pvalues` floors them at floatmin
        # and the shared printer renders "<1e-16", never "0.0".
        g = load_golden(joinpath(@__DIR__, "fixtures", "zach_poisson.toml"))
        n = Int(g.values["n_actors"])
        net = network(n; directed=false)
        for (a, b, w) in zip(Int.(g.values["edge_src"]), Int.(g.values["edge_dst"]),
                             Int.(g.values["edge_weight"]))
            add_edge!(net, a, b); set_edge_attribute!(net, :weight, a, b, w)
        end
        fit = fit_ergm_count(net, [SumTerm(), NonzeroTerm()])
        @test all(abs.(fit.z_values) .> 8.3)
        @test all(p -> 0 < p < 1e-16, fit.p_values)   # far below the print floor
        @test fit.p_values == Networks.z_pvalues(fit.z_values)
        out = sprint(show, fit)
        @test occursin("<1e-16", out)
        @test !occursin("0.0000 ***", out)
    end

    @testset "gof extends the shared Networks.gof generic" begin
        # One generic across the ecosystem: the method is added to
        # Networks.gof, not a package-local function
        @test ERGMCount.gof === Networks.gof

        net = random_count_net(6; seed=19)
        result = fit_ergm_count(net, [SumTerm(), NonzeroTerm()])
        g = ERGMCount.gof(result; n_sim=8, burnin=10, interval=2,
                          rng=Random.Xoshiro(31))
        @test g isa Networks.GOFResult
        @test Networks.n_simulations(g) == 8
        @test length(g.statistics) == 2
        @test g.statistics[1].name == "model statistics"
        @test g.statistics[1].labels == ["sum", "nonzero"]
        @test g.statistics[1].observed ==
              [compute(SumTerm(), net), compute(NonzeroTerm(), net)]
        @test g.statistics[2].name == "dyad count values"
        # Dyad-value counts sum to the number of dyads in every row
        n_dyads = 6 * 5
        @test sum(g.statistics[2].observed) == n_dyads
        @test all(sum(g.statistics[2].simulated, dims=2) .== n_dyads)
        @test all(p -> 0 < p <= 1, g.statistics[1].p_values)
        # Formatted display renders the shared GOF table
        out = sprint(show, g)
        @test occursin("Goodness-of-fit assessment: Count ERGM", out)
        @test occursin("MC p-value", out)
    end

    @testset "Truncation is explicit" begin
        # The support `0:max_val` is part of the ESTIMAND for an unbounded
        # reference, not an implementation detail. It must be visible, and a
        # fit that leans on the bound must say so.

        @testset "is_truncating trait" begin
            # Unbounded references: the enumeration truncates them
            @test is_truncating(PoissonReference())
            @test is_truncating(GeometricReference())
            # Genuinely bounded references: the support IS the model
            @test !is_truncating(BinomialReference(5))
            @test !is_truncating(DiscUnifReference(4))
            @test !is_truncating(DiscUnif2Reference(1, 3))
        end

        @testset "result records support and boundary mass" begin
            net = network(6; directed=true)
            for (i, j, w) in [(1, 2, 2), (2, 3, 1), (3, 1, 3), (4, 5, 1)]
                add_edge!(net, i, j)
                set_edge_attribute!(net, :weight, i, j, w)
            end

            res = fit_ergm_count(net, [SumTerm()]; reference=PoissonReference())
            @test res.truncated
            @test res.max_val >= 10
            @test 0.0 <= res.boundary_mass <= 1.0

            out = sprint(show, res)
            @test occursin("Support:", out)
            @test occursin("TRUNCATED", out)
            @test occursin("Boundary mass", out)

            # A bounded reference is not a truncation, and says so
            resb = fit_ergm_count(net, [SumTerm()]; reference=BinomialReference(5))
            @test !resb.truncated
            @test resb.boundary_mass == 0.0
            outb = sprint(show, resb)
            @test occursin("bounded reference", outb)
            @test !occursin("TRUNCATED", outb)
        end

        @testset "an inadequate bound warns" begin
            # Squeeze the support until the fitted conditionals must pile up on
            # the boundary; the fit is then a different (truncated) family and
            # the user has to be told.
            net = network(5; directed=true)
            for (i, j) in [(1, 2), (2, 3), (3, 4), (4, 5), (5, 1)]
                add_edge!(net, i, j)
                set_edge_attribute!(net, :weight, i, j, 2)
            end

            res = @test_logs (:warn, r"TRUNCATED|truncated|max_val"i) match_mode = :any begin
                fit_ergm_count(net, [SumTerm()]; reference=PoissonReference(),
                               max_val=2)
            end
            @test res.truncated
            @test res.boundary_mass > BOUNDARY_MASS_TOL
        end

        @testset "an adequate bound is quiet" begin
            # Sparse counts under a generous bound: no boundary warning
            net = network(6; directed=true)
            add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 1)
            res = fit_ergm_count(net, [SumTerm()]; reference=PoissonReference(),
                                 max_val=30)
            @test res.boundary_mass < BOUNDARY_MASS_TOL
        end
    end

    @testset "Missing-data contract: masked dyads are refused everywhere" begin
        # Count MPLE enumerates every dyad as observed, and the Gibbs sampler
        # conditions every dyad on every other's face value: a masked dyad
        # would enter both at a value that was never observed. No `missing=`
        # keyword exists, so the contract's declared vocabulary is `:error` only.
        @test Networks.missing_policies(fit_ergm_count) == (:error,)
        @test Networks.missing_policies(ergm_count) == (:error,)
        @test Networks.missing_policies(simulate_count_ergm) == (:error,)
        @test Networks.missing_policies(ERGMCount.count_mple) == (:error,)
        @test Networks.supports_missing(fit_ergm_count) == false
        @test Networks.supports_missing(simulate_count_ergm) == false

        net = network(5; directed=true)
        add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 2)
        add_edge!(net, 2, 3); set_edge_attribute!(net, :weight, 2, 3, 1)

        set_missing_dyad!(net, 3, 4)          # absent face
        err = try; fit_ergm_count(net, [SumTerm()]); nothing catch e; e end
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("fit_ergm_count", msg)
        @test occursin("clear_missing_dyads!", msg)
        @test !occursin(":face", msg)         # a keyword that does not exist

        # The unfitted simulator is guarded too (it used to Gibbs-resample from
        # the face values of a masked network)
        err = try
            simulate_count_ergm(net, [SumTerm()], [0.1]; n_sim=1, burnin=1, interval=1)
            nothing
        catch e; e end
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("simulate_count_ergm", msg)
        @test occursin("clear_missing_dyads!", msg)
        @test !occursin(":face", msg)

        clear_missing_dyads!(net)
        set_missing_dyad!(net, 1, 2)          # present face
        @test_throws ArgumentError fit_ergm_count(net, [SumTerm()])

        # `count_mple(model)` — the documented estimator entry point — is
        # guarded too: it used to enumerate the masked dyad's face value as an
        # observed zero (coef [0.179, -2.025], nobs 30 on a 6-node network the
        # round-2 grader ran) while `fit_ergm_count` on the same network refused
        set_missing_dyad!(net, 3, 4)
        model = CountERGMModel([SumTerm(), NonzeroTerm()], net)   # a model can be built ...
        err = try; count_mple(model); nothing catch e; e end        # ... but not fit
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("count_mple", msg)
        @test occursin("clear_missing_dyads!", msg)
        @test !occursin(":face", msg)
        @test_throws ArgumentError count_mple(model; max_val=20)
        @test_throws ArgumentError count_mple(model; se=:bootstrap, n_boot=2)
        # ... and masking after construction is caught the same way
        clear_missing_dyads!(net)
        model2 = CountERGMModel([SumTerm()], net)
        set_missing_dyad!(model2.network, 4, 5)
        @test_throws ArgumentError count_mple(model2)
        clear_missing_dyads!(net)
        @test count_mple(model2) isa CountERGMResult

        fit = fit_ergm_count(net, [SumTerm()])
        @test fit isa CountERGMResult
        @test missing_method(fit) == :rejected

        # `gof` and the fitted-result simulator inherit the guard through the
        # result's network: mask a dyad after the fit and both refuse
        set_missing_dyad!(fit.model.network, 4, 5)
        @test_throws ArgumentError simulate_count_ergm(fit; n_sim=1, burnin=1, interval=1)
        @test_throws ArgumentError ERGMCount.gof(fit; n_sim=2, burnin=1, interval=1)
        clear_missing_dyads!(fit.model.network)
        @test length(simulate_count_ergm(fit; n_sim=1, burnin=1, interval=1)) == 1
    end

    @testset "Result metadata protocol" begin
        net = network(6; directed=true)
        for (i, j, w) in [(1, 2, 2), (2, 3, 1), (3, 1, 3), (4, 5, 1), (5, 6, 2)]
            add_edge!(net, i, j)
            set_edge_attribute!(net, :weight, i, j, w)
        end

        # Poisson: unbounded reference, so the enumerated support is a
        # TRUNCATION — the approximation is named, with the boundary mass
        pois = fit_ergm_count(net, [SumTerm()]; reference=PoissonReference())
        md = fit_metadata(pois)
        @test md.estimand == :count_ergm
        @test md.objective == :pseudolikelihood
        @test md.se_method == :hessian
        @test md.missing_method == :rejected
        @test !md.is_exact                       # truncated
        @test any(occursin("truncated at 0:$(pois.max_val)", a)
                  for a in md.approximations)
        @test any(occursin("boundary mass", a) for a in md.approximations)
        # The prose the `show` method prints and the protocol agree
        @test occursin("TRUNCATED", sprint(show, pois))

        # Bounded reference + a dyad-independent term: the support IS the model
        # and the dyad conditionals ARE the model's, so the pseudo-likelihood is
        # the likelihood. Same estimator, exact fit.
        bin = fit_ergm_count(net, [SumTerm()]; reference=BinomialReference(5))
        @test is_exact(bin)
        @test isempty(approximations(bin))
        @test objective(bin) == :pseudolikelihood     # same estimator as above

        # One dyad-dependent term (squared out-strengths) and it is not exact
        # any more, with the anticonservative-SE caveat attached
        dep = fit_ergm_count(net, [SumTerm(), NodeOSumTerm()];
                             reference=BinomialReference(5))
        @test !is_exact(dep)
        @test any(occursin("anticonservative", a) for a in approximations(dep))

        # The dependence classification the above reads
        @test !is_dyad_dependent(SumTerm())
        @test !is_dyad_dependent(NonzeroTerm())
        @test is_dyad_dependent(NodeOSumTerm())
        @test is_dyad_dependent(CountMutualTerm())
    end
    @testset "Robust standard errors: se=:bootstrap" begin
        # Issue #9 / ERGMCount#2: the inverse-pseudo-Hessian SEs of a
        # dyad-dependent count model are anticonservative and there was no
        # alternative. `se=:bootstrap` adds a parametric bootstrap (Gibbs-simulate
        # at θ̂ with `simulate_count_ergm`, refit, empirical covariance) on the ONE
        # shared `Networks.bootstrap_cov` loop — same API as `ERGM.mple`.
        # The network carries reciprocated pairs: on one without any,
        # `mutual.min` sits at its minimum attainable value, has no finite
        # MPLE, is fixed at -Inf (see "Boundary statistics") and cannot be
        # bootstrapped — 0.2 fixed both of those, which silently produced
        # θ̂ ≈ −21 with a 1.7e4 standard error here before.
        net = count_net(8, [(1,2,2), (2,1,1), (2,3,1), (3,2,2), (3,1,3), (1,4,1),
                            (4,5,2), (5,4,1), (5,6,1), (6,1,2), (2,7,3), (7,8,1),
                            (8,2,2), (8,7,1)])
        terms = [SumTerm(), CountMutualTerm()]      # dyad-DEPENDENT (mutual)
        ref = BinomialReference(3)
        @test compute(CountMutualTerm(), net) > 0   # not a boundary statistic

        hess = fit_ergm_count(net, terms; reference=ref)
        boot = fit_ergm_count(net, terms; reference=ref, se=:bootstrap,
                              n_boot=100, rng=MersenneTwister(7))

        # The bootstrap replaces the COVARIANCE, not the point estimate
        @test coef(boot) == coef(hess)
        @test loglikelihood(boot) == loglikelihood(hess)
        @test stderror(boot) != stderror(hess)
        @test vcov(boot) != vcov(hess)
        @test all(isfinite, stderror(boot))

        # Reproducible under a fixed rng
        boot2 = fit_ergm_count(net, terms; reference=ref, se=:bootstrap,
                               n_boot=100, rng=MersenneTwister(7))
        @test stderror(boot2) == stderror(boot)
        @test vcov(boot2) == vcov(boot)
        @test stderror(fit_ergm_count(net, terms; reference=ref, se=:bootstrap,
                                      n_boot=100, rng=MersenneTwister(8))) !=
              stderror(boot)

        # The dependent model's robust SEs EXCEED the Hessian ones — by ~40% for
        # `sum` and ~60% for `mutual.min` here — and that gap IS the
        # anticonservatism of the pseudo-likelihood Hessian, the whole point of
        # the option. The direction of the correction is a property of the fit,
        # not a law; it is pinned for this fixture.
        @test all(stderror(boot) .> stderror(hess))
        @test stderror(boot)[1] / stderror(hess)[1] > 1.1
        @test boot.boot_replicates isa Matrix{Float64}
        @test size(boot.boot_replicates) == (100, 2)
        @test all(isfinite, boot.boot_replicates)
        @test hess.boot_replicates === nothing

        # `se_method` reports what was ACTUALLY used, in both directions
        @test se_method(hess) === :hessian
        @test se_method(boot) === :bootstrap
        @test fit_metadata(hess).se_method === :hessian
        @test fit_metadata(boot).se_method === :bootstrap

        # ... and so does the printed output: the anticonservatism caveat is a
        # claim about the inverse Hessian, so it must NOT be made of a bootstrap
        out_h = sprint(show, hess)
        out_b = sprint(show, boot)
        @test occursin("inverse pseudo-Hessian", out_h)
        @test occursin("anticonservative", out_h)
        @test occursin("parametric bootstrap", out_b)
        @test !occursin("anticonservative", out_b)

        # The approximations list agrees with the printed prose (one predicate)
        @test any(occursin("anticonservative", a) for a in approximations(hess))
        @test !any(occursin("anticonservative", a) for a in approximations(boot))
        @test any(occursin("parametric bootstrap", a) for a in approximations(boot))
        # The POINT ESTIMATE is a pseudo-likelihood estimate either way, and both
        # fits still say so
        @test any(occursin("biased in finite samples", a)
                  for a in approximations(boot))

        # A dyad-INDEPENDENT model: the pseudo-likelihood is the likelihood, so
        # the Hessian SEs are correct there and the bootstrap is merely optional
        indep = fit_ergm_count(net, [SumTerm()]; reference=ref)
        indep_b = fit_ergm_count(net, [SumTerm()]; reference=ref, se=:bootstrap,
                                 n_boot=20, rng=MersenneTwister(3))
        @test coef(indep_b) == coef(indep)
        @test is_exact(indep) && is_exact(indep_b)
        @test !occursin("anticonservative", sprint(show, indep))
        @test se_method(indep_b) === :bootstrap

        # Unknown se symbols are rejected, not silently ignored (through the
        # shared validator, whose message shape is pinned elsewhere)
        @test_throws ArgumentError fit_ergm_count(net, terms; reference=ref,
                                                  se=:sandwich)
        @test_throws ArgumentError fit_ergm_count(net, terms; reference=ref,
                                                  se=:bootstrap, n_boot=1)
        # The bootstrap's Gibbs controls default through the shared rule too
        @test_throws ArgumentError fit_ergm_count(net, terms; reference=ref,
                                                  se=:bootstrap, n_boot=4,
                                                  boot_interval=0)
    end

    @testset "Bootstrap is thread-count independent (fresh process)" begin
        # The replicates are drawn serially from the caller's `rng`; only the
        # refits run under `Threads.@threads` (in `Networks.bootstrap_cov`), and
        # each refit is deterministic — so the standard errors are bit-identical
        # whatever the thread count. Prove it in a fresh process with a
        # DIFFERENT thread count, as ERGM.jl does for `mcmle`.
        net = count_net(8, [(1,2,2), (2,1,1), (2,3,1), (3,2,2), (3,1,3), (1,4,1),
                            (4,5,2), (5,4,1), (5,6,1), (6,1,2), (2,7,3), (7,8,1),
                            (8,2,2), (8,7,1)])
        boot = fit_ergm_count(net, [SumTerm(), CountMutualTerm()];
                              reference=BinomialReference(3), se=:bootstrap,
                              n_boot=20, rng=Xoshiro(7))
        other_threads = Threads.nthreads() == 1 ? 4 : 1
        script = """
            using ERGMCount, Networks, Graphs, Random
            net = network(8; directed=true)
            for (i, j, w) in [(1,2,2), (2,1,1), (2,3,1), (3,2,2), (3,1,3), (1,4,1),
                              (4,5,2), (5,4,1), (5,6,1), (6,1,2), (2,7,3), (7,8,1),
                              (8,2,2), (8,7,1)]
                add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
            end
            f = fit_ergm_count(net, [SumTerm(), CountMutualTerm()];
                               reference=BinomialReference(3), se=:bootstrap,
                               n_boot=20, rng=Xoshiro(7))
            println(Threads.nthreads()); println(repr(coef(f))); println(repr(stderror(f)))
            """
        cmd = `$(Base.julia_cmd()) --startup-file=no --threads=$other_threads --project=$(dirname(@__DIR__)) -e $script`
        lines = split(strip(read(pipeline(cmd; stderr=devnull), String)), '\n')
        @test length(lines) == 3
        @test parse(Int, lines[1]) == other_threads
        @test lines[2] == repr(coef(boot))
        @test lines[3] == repr(stderror(boot))
    end

    @testset "Bootstrap replicates without a finite MPLE are excluded, loudly" begin
        # A sparse `greaterthan.3` statistic: many simulated replicates contain no
        # count above 3, so the statistic sits at its boundary there and the
        # refit has no finite MPLE. Those replicates must not enter the
        # covariance as NaN rows (the NaN-SE hazard); they are excluded and
        # counted, in a warning, in `approximations` and in `show`.
        net = count_net(7, [(1,2,1), (2,3,1), (3,4,1), (4,5,1), (5,6,1), (6,7,1),
                            (7,1,1), (1,3,4), (2,4,1), (3,5,1)])
        terms = [SumTerm(), GreaterthannTerm(3)]
        # EXACTLY one warning: the refits run quietly (a boundary in a
        # simulated replicate is not a fact about the data) and the exclusion
        # is reported once, in aggregate — as `ERGM.mple` does
        boot = @test_logs (:warn, r"refits had no finite, converged count MPLE") begin
            fit_ergm_count(net, terms; reference=BinomialReference(4),
                           se=:bootstrap, n_boot=40, rng=Xoshiro(1))
        end
        reps = boot.boot_replicates
        @test size(reps) == (40, 2)
        ok = [all(isfinite, view(reps, b, :)) for b in 1:40]
        @test 2 <= count(ok) < 40
        @test any(isnan, reps)                     # the excluded rows are NaN
        @test all(isfinite, stderror(boot))
        @test all(isfinite, vcov(boot))
        @test se_method(boot) === :bootstrap
        @test fit_metadata(boot).se_method === :bootstrap
        # The covariance is over the surviving replicates only
        @test vcov(boot) ≈ cov(reps[ok, :])
        @test any(occursin("were excluded", a) for a in approximations(boot))
        @test occursin("were excluded", sprint(show, boot))
    end

    @testset "Boundary statistics: R's drop semantics, not a silent number" begin
        # `mutual.min` is 0 on a network with no reciprocated pair — its
        # smallest attainable value on every dyad's conditional support — so the
        # pseudo-likelihood increases monotonically as its coefficient → −∞ and
        # no finite MPLE exists. 0.2 "converged" to θ̂ ≈ −21 with a 1.7e4
        # standard error and three significance stars. Now, as in R ergm: the
        # coefficient is fixed at -Inf with SE 0 and p 0, the rest is fit on
        # the restricted supports (the exact limit), `dof` counts only the
        # finite coefficients, and the user is warned.
        net = count_net(8, [(1,2,2), (2,3,1), (3,1,3), (1,4,1), (4,5,2),
                            (5,6,1), (6,1,2), (2,7,3), (7,8,1), (8,2,2)])
        @test compute(CountMutualTerm(), net) == 0
        fit = @test_logs (:warn, r"mutual.min are at their smallest attainable value") match_mode = :any begin
            fit_ergm_count(net, [SumTerm(), CountMutualTerm()]; reference=BinomialReference(3))
        end
        @test fit.converged
        @test coef(fit)[2] == -Inf
        @test stderror(fit)[2] == 0.0
        @test fit.z_values[2] == -Inf
        @test fit.p_values[2] == 0.0
        @test isfinite(coef(fit)[1]) && stderror(fit)[1] > 0
        @test dof(fit) == 1
        @test aic(fit) ≈ -2 * loglikelihood(fit) + 2
        @test all(vcov(fit)[2, :] .== 0) && all(vcov(fit)[:, 2] .== 0)
        @test confint(fit)[2, :] == [-Inf, -Inf]
        # A limit is not a maximizer: the fit is not exact, whatever the terms
        @test !is_exact(fit)
        # `warn=false` silences R's sentence (the bootstrap refits use it); the
        # result is identical and still records the boundary
        quiet = @test_logs min_level = Base.CoreLogging.Warn begin
            fit_ergm_count(net, [SumTerm(), CountMutualTerm()];
                           reference=BinomialReference(3), warn=false)
        end
        @test coef(quiet) == coef(fit) && stderror(quiet) == stderror(fit)
        @test any(occursin("fixed by a boundary statistic", a) for a in approximations(quiet))
        # The limit θ_mutual → −∞ forbids reciprocation: the 10 dyads whose
        # reverse carries a count are restricted to 0 (min(y, y_ji) would grow
        # with y there) and contribute nothing; the other 46 dyads keep the
        # full support 0:3 and, under Binomial(3) with a `sum` coefficient θ,
        # are iid Binomial(3, logistic(θ)) — so θ̂ is the logit of their mean
        # y/3 = 18/138, and the loglik is that binomial log-likelihood.
        p̂ = 18 / 138
        @test coef(fit)[1] ≈ log(p̂ / (1 - p̂)) atol = 1e-8
        @test loglikelihood(fit) ≈ 8 * log(3) + 18 * log(p̂) + 120 * log(1 - p̂) atol = 1e-8
        # Reported everywhere a user looks
        @test any(occursin("fixed by a boundary statistic", a) for a in approximations(fit))
        out = sprint(show, fit)
        @test occursin("-Inf", out)
        @test occursin("Note: coefficient fixed by a boundary statistic (mutual.min at -Inf)", out)
        @test coeftable(fit)["mutual.min"].estimate == -Inf
        # No finite model to simulate from: the bootstrap, the simulator and
        # gof all refuse rather than produce NaN
        err = try
            fit_ergm_count(net, [SumTerm(), CountMutualTerm()]; reference=BinomialReference(3),
                           se=:bootstrap, n_boot=10)
            nothing
        catch e; e end
        @test err isa ArgumentError
        @test occursin("se=:bootstrap is not available", sprint(showerror, err))
        @test occursin("mutual.min", sprint(showerror, err))
        @test_throws ArgumentError simulate_count_ergm(fit; n_sim=1, burnin=1, interval=1)
        @test_throws ArgumentError ERGMCount.gof(fit; n_sim=2, burnin=1, interval=1)

        # The other side: on a complete network every dyad is nonzero, so
        # `nonzero` is at its LARGEST attainable value → +Inf, and `sum` is the
        # zero-truncated Poisson MLE: λ/(1 − e^{-λ}) = ȳ = 9/6
        cn = count_net(3, [(1,2,1), (2,1,2), (1,3,1), (3,1,1), (2,3,3), (3,2,1)])
        fc = @test_logs (:warn, r"nonzero are at their largest attainable value") match_mode = :any begin
            fit_ergm_count(cn, [SumTerm(), NonzeroTerm()])
        end
        @test coef(fc)[2] == Inf
        @test fc.p_values[2] == 0.0
        λ = exp(coef(fc)[1])
        @test λ / (1 - exp(-λ)) ≈ 9 / 6 atol = 1e-6
        @test dof(fc) == 1
        @test !is_exact(fc)

        # `greaterthan.1` when every observed count exceeds 1: at its largest
        # attainable value (one per dyad) → +Inf, and the restricted supports
        # are 2:max_val, so `sum` is the MLE of a Poisson truncated below 2
        gt = count_net(4, [(1,2,3), (2,1,2), (1,3,5), (3,1,2), (2,3,4), (3,2,2),
                           (1,4,2), (4,1,3), (2,4,2), (4,2,6), (3,4,2), (4,3,3)])
        fg = @test_logs (:warn, r"greaterthan.1 are at their largest attainable value") match_mode = :any begin
            fit_ergm_count(gt, [SumTerm(), GreaterthannTerm(1)])
        end
        @test coef(fg)[2] == Inf
        @test stderror(fg)[2] == 0.0
        @test isfinite(coef(fg)[1])
        λg = exp(coef(fg)[1])
        # E[y | y ≥ 2] = ȳ = 36/12 for Poisson(λ): (λ − λe^{-λ}) / (1 − e^{-λ} − λe^{-λ})
        @test (λg - λg * exp(-λg)) / (1 - exp(-λg) - λg * exp(-λg)) ≈ 36 / 12 atol = 1e-6

        # An empty network: every statistic at its minimum, nothing free
        fe = @test_logs (:warn, r"sum, nonzero are at their smallest") match_mode = :any begin
            fit_ergm_count(network(5), [SumTerm(), NonzeroTerm()])
        end
        @test coef(fe) == [-Inf, -Inf]
        @test dof(fe) == 0
        @test loglikelihood(fe) == 0.0
        @test fe.converged
        @test fe.iterations == 0 && fe.gradient_norm == 0.0
        @test !is_exact(fe)
    end

    @testset "Non-convergence is loud, recorded, and never exact" begin
        # One Newton iteration on a single-rung support (max_val = 10 is the
        # ladder's base here, so there is no warm start to lean on): the fit
        # cannot converge, and must say so everywhere a user looks — a warning
        # naming maxiter, tol and the pseudo-score, `converged = false`,
        # `approximations`, `show`, and `is_exact = false`.
        net = random_count_net(8; seed=9)
        terms = [SumTerm(), NonzeroTerm(), CountMutualTerm()]
        unc = @test_logs (:warn, r"did not converge in 1 iteration") match_mode = :any begin
            fit_ergm_count(net, terms; maxiter=1, max_val=10)
        end
        @test !unc.converged
        @test unc.iterations == 1
        @test unc.gradient_norm > 1e-4
        @test !is_exact(unc)
        @test !fit_metadata(unc).is_exact
        @test any(occursin("did not converge in 1 iteration", a) for a in approximations(unc))
        @test any(occursin("raise `maxiter`", a) for a in approximations(unc))
        @test occursin("Converged: false", sprint(show, unc))
        @test occursin("did not converge", sprint(show, unc))
        # The warning names the budget, the tolerance and the score
        rec = Test.collect_test_logs() do
            fit_ergm_count(net, terms; maxiter=1, max_val=10)
        end
        warns = [r.message for r in rec[1] if r.level == Base.CoreLogging.Warn]
        @test any(occursin("maxiter = 1", w) && occursin("tol = 1.0e-8", w) &&
                  occursin("pseudo-score norm", w) && occursin("boundary or non-identified", w)
                  for w in warns)
        # `warn=false` silences the warning; the record is unchanged
        quiet = @test_logs min_level = Base.CoreLogging.Warn begin
            fit_ergm_count(net, terms; maxiter=1, max_val=10, warn=false)
        end
        @test !quiet.converged && coef(quiet) == coef(unc)
        # The budget is the caller's: running out of it is not rescued by the
        # coordinate-ascent start (that exists for a Newton that cannot MOVE)
        conv = fit_ergm_count(net, terms; max_val=10)
        @test conv.converged && conv.iterations > 1
        @test conv.gradient_norm < 1e-4
        @test loglikelihood(conv) > loglikelihood(unc)
        @test is_exact(fit_ergm_count(net, [SumTerm()]; reference=BinomialReference(4)))
    end

    @testset "Error-controlled count support (N6, item 17)" begin
        # The default `max_val` used to be `max(10, 2·max count)` — on zach that
        # is 14, which moves θ̂ by 1.4e-6/3.9e-6, ABOVE the fixture's own 1e-6
        # tolerance (it passed only because the testset pinned max_val=30). The
        # support is now doubled until the estimates stop moving; the default
        # path is pinned against the fixture below ("Golden fixture" testset).
        net = network(6; directed=true)
        for (i, j, w) in [(1, 2, 2), (2, 3, 1), (3, 1, 3), (4, 5, 1), (5, 6, 2)]
            add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
        end

        # Adaptive default: doubled from 10 to 20, converged
        fit = fit_ergm_count(net, [SumTerm()]; reference=PoissonReference())
        @test fit.support_control === :converged
        @test fit.support_stable
        @test fit.max_val == 20
        @test fit.support_tol == 1e-3
        # The stopping rule is a conjunction: the estimates settled AND the
        # previous bound omitted negligible mass AND the new bound is not
        # shaping the fit
        @test 0 <= fit.support_delta <= 1e-3
        @test 0 <= fit.omitted_tail <= 1e-3
        @test fit.boundary_mass <= BOUNDARY_MASS_TOL
        @test isfinite(fit.support_delta) && isfinite(fit.omitted_tail)
        @test any(occursin("error-controlled doubling", a) for a in approximations(fit))
        @test occursin("chosen by doubling", sprint(show, fit))
        # ... and the fit at the larger bound IS the reported one
        fixed20 = fit_ergm_count(net, [SumTerm()]; reference=PoissonReference(), max_val=20)
        @test coef(fit) == coef(fixed20)
        @test fixed20.support_control === :fixed
        @test fixed20.support_stable
        @test isnan(fixed20.support_delta) && isnan(fixed20.omitted_tail)
        @test occursin("max_val fixed by the caller", sprint(show, fixed20))
        @test !any(occursin("doubling", a) for a in approximations(fixed20))

        # A tighter tolerance doubles further; the default stopping rule is a
        # bound on the truncation error in SE units, so tightening it can only
        # move the estimate by less than the previous delta
        tight = fit_ergm_count(net, [SumTerm()]; reference=PoissonReference(),
                               support_tol=1e-8)
        @test tight.max_val >= fit.max_val
        @test abs(coef(tight)[1] - coef(fit)[1]) <= fit.support_delta * stderror(fit)[1] + 1e-12

        # Bounded references have nothing to control: nothing was truncated,
        # so the truncation error is exactly zero, and the fit never loops —
        # the ladder has a single rung and `max_val` is the reference's own
        bin = fit_ergm_count(net, [SumTerm()]; reference=BinomialReference(4))
        @test bin.support_control === :bounded
        @test bin.support_stable
        @test !bin.truncated
        @test bin.support_delta == 0.0 && bin.omitted_tail == 0.0
        @test bin.max_val == 4
        @test ERGMCount._support_ladder(ERGMCount._support(BinomialReference(4), 0), 10) == [0:4]
        for ref in (BinomialReference(4), DiscUnifReference(6), DiscUnif2Reference(0, 5))
            b = fit_ergm_count(net, [SumTerm()]; reference=ref,
                               support_tol=1e-300, max_doublings=1)   # would loop if it could
            @test b.support_control === :bounded && b.support_stable
            @test b.support_delta == 0.0
            @test b.max_val == last(ERGMCount._support(ref, 0))
        end
        @test coef(fit_ergm_count(net, [SumTerm()]; reference=BinomialReference(4),
                                  max_val=99)) == coef(bin)   # ignored, as documented

        # The doubling cap: a tolerance no doubling can meet within the budget
        # is reported as :unconverged, with a warning, in `approximations` and
        # in `show` — never as a quiet fit
        unc = @test_logs (:warn, r"error-controlled support did not converge") match_mode = :any begin
            fit_ergm_count(net, [SumTerm()]; reference=PoissonReference(),
                           support_tol=1e-300, max_doublings=1)
        end
        @test unc.support_control === :unconverged
        @test !unc.support_stable
        @test unc.max_val == 20
        @test any(occursin("did NOT converge", a) for a in approximations(unc))
        @test occursin("did NOT converge", sprint(show, unc))
        @test !is_exact(unc)
        @test ERGMCount._MAX_SUPPORT_DOUBLINGS == 8
        # ... and the default-path twin of "an inadequate bound warns" does
        # NOT warn: the doubling moves past the bound that bites
        tight_net = network(5; directed=true)
        for (i, j) in [(1, 2), (2, 3), (3, 4), (4, 5), (5, 1)]
            add_edge!(tight_net, i, j); set_edge_attribute!(tight_net, :weight, i, j, 2)
        end
        ok = @test_logs min_level = Base.CoreLogging.Warn begin
            fit_ergm_count(tight_net, [SumTerm()]; reference=PoissonReference())
        end
        @test ok.support_control === :converged && ok.boundary_mass <= BOUNDARY_MASS_TOL
        @test ok.max_val >= 20

        # Continuation in the support: a cold Newton start on a wide support
        # overshoots (GeometricReference + sum + nonzero on 0:40 from θ = 0
        # gave up after two iterations); every fit now climbs the ladder
        # base, 2·base, ..., top with warm starts, so a wide fixed `max_val`,
        # a wide bounded reference and the adaptive path all converge and agree
        gnet = network(12; directed=true)
        grng = Xoshiro(1)
        for i in 1:12, j in 1:12
            if i != j && rand(grng) < 0.2
                add_edge!(gnet, i, j); set_edge_attribute!(gnet, :weight, i, j, rand(grng, 1:5))
            end
        end
        gterms = [SumTerm(), NonzeroTerm()]
        g40 = fit_ergm_count(gnet, gterms; reference=GeometricReference(), max_val=40)
        @test g40.converged
        g200 = fit_ergm_count(gnet, gterms; reference=GeometricReference(), max_val=200)
        @test g200.converged
        @test coef(g200) ≈ coef(g40) atol = 1e-4       # P(y > 40) ≈ 0.68^40 per dyad
        gad = fit_ergm_count(gnet, gterms; reference=GeometricReference())
        @test gad.converged && gad.support_control === :converged
        @test coef(gad) ≈ coef(g200) atol = 1e-3
        # BinomialReference(100): at θ = 0 the conditional ∝ C(100, y) is a point
        # mass at the top of the first rung, the Hessian is singular and Newton
        # has no direction; the coordinate-ascent rescue start recovers it, and
        # the fit is the binomial-logit MLE (each dyad Binomial(100, logistic(θ_sum))
        # once nonzero is accounted for) — check the score is zero there
        b100 = fit_ergm_count(gnet, gterms; reference=BinomialReference(100))
        @test b100.converged && b100.support_control === :bounded
        @test all(isfinite, stderror(b100))
        Db = ERGMCount._count_design(b100.model, 0:100)
        maskb = fill(true, 101, length(Db.n_tot))
        @test norm(ERGMCount._count_derivatives(Db, [1, 2], maskb)(coef(b100))[2]) < 1e-6
        # The rescue itself: from zeros (where Newton has no direction) the
        # coordinate start climbs the objective, and Newton from there reaches
        # the optimum
        db = ERGMCount._count_derivatives(Db, [1, 2], maskb)
        θr = ERGMCount._coordinate_ascent_start(db, zeros(2))
        @test db(θr)[1] > db(zeros(2))[1]
        @test !ERGMCount._quiet_newton(db, zeros(2); maxiter=100, tol=1e-8).converged
        nr = ERGMCount._quiet_newton(db, θr; maxiter=100, tol=1e-8)
        @test nr.converged
        @test nr.θ ≈ coef(b100) atol = 1e-6
        @test ERGMCount._support_ladder(0:40, 10) == [0:10, 0:20, 0:40]
        @test ERGMCount._support_ladder(0:100, 14) == [0:14, 0:28, 0:56, 0:100]
        @test ERGMCount._support_ladder(0:3, 10) == [0:3]
        @test ERGMCount._support_ladder(1:5, 10) == [1:5]

        # Validation
        @test_throws ArgumentError fit_ergm_count(net, [SumTerm()]; support_tol=0.0)
        @test_throws ArgumentError fit_ergm_count(net, [SumTerm()]; max_doublings=0)
        @test_throws ArgumentError fit_ergm_count(net, [SumTerm()]; max_val=0)

        # The Gibbs sampler's defaults come from THE dyad-scaled rule shared with
        # ERGM.jl (`ERGM._mcmc_defaults`), converted from toggles to sweeps
        for nd in (30, 100, 870, 5000)
            d = ERGM._mcmc_defaults(nd)
            @test ERGMCount._gibbs_defaults(nd) ==
                  (burnin=cld(d.burnin, nd), interval=cld(d.interval, nd))
        end
        @test ERGMCount._gibbs_defaults(30) == (burnin=20, interval=4)
        @test ERGMCount._gibbs_defaults(870) == (burnin=20, interval=1)
        sims = simulate_count_ergm(fit; n_sim=2, rng=Xoshiro(2))   # defaults resolve
        @test length(sims) == 2

        # ... and they are the sweeps that actually run: on a 34-node directed
        # network (1122 dyads; 561 undirected) the default is 20 burn-in
        # sweeps and 1 between draws, so a default-path draw equals the one
        # made with those numbers spelled out, and differs from a 19-sweep one
        @test ERGMCount._gibbs_defaults(561) == (burnin=20, interval=1)
        @test ERGMCount._gibbs_defaults(1122) == (burnin=20, interval=1)
        net34 = random_count_net(34; seed=34, p=0.1)
        fit34 = fit_ergm_count(net34, [SumTerm(), NonzeroTerm()]; max_val=10)
        weights_of(s) = get_edge_attribute(s, :weight)
        default_draw = simulate_count_ergm(fit34; n_sim=1, rng=Xoshiro(5))[1]
        spelled_out = simulate_count_ergm(fit34; n_sim=1, burnin=20, interval=1,
                                          rng=Xoshiro(5))[1]
        shorter = simulate_count_ergm(fit34; n_sim=1, burnin=19, interval=1,
                                      rng=Xoshiro(5))[1]
        @test weights_of(default_draw) == weights_of(spelled_out)
        @test weights_of(default_draw) != weights_of(shorter)
        # gof and the bootstrap inherit the same resolution (nothing → rule)
        g_default = gof(fit34; n_sim=2, rng=Xoshiro(6))
        g_spelled = gof(fit34; n_sim=2, burnin=20, interval=1, rng=Xoshiro(6))
        @test g_default.statistics[1].simulated == g_spelled.statistics[1].simulated
        b_default = fit_ergm_count(net34, [SumTerm(), NonzeroTerm()]; max_val=10,
                                   se=:bootstrap, n_boot=3, rng=Xoshiro(7))
        b_spelled = fit_ergm_count(net34, [SumTerm(), NonzeroTerm()]; max_val=10,
                                   se=:bootstrap, n_boot=3, boot_burnin=20,
                                   boot_interval=1, rng=Xoshiro(7))
        @test b_default.boot_replicates == b_spelled.boot_replicates
        # The explicit-specification form keeps its literal `max_val=20` but
        # resolves burnin/interval through the same rule
        e_default = simulate_count_ergm(net34, [SumTerm()], [-0.5]; n_sim=1, rng=Xoshiro(8))[1]
        e_spelled = simulate_count_ergm(net34, [SumTerm()], [-0.5]; n_sim=1, burnin=20,
                                        interval=1, max_val=20, rng=Xoshiro(8))[1]
        @test weights_of(e_default) == weights_of(e_spelled)
    end

    @testset "Two-mode (bipartite) networks are refused everywhere" begin
        # The estimator enumerates and the sampler resamples EVERY off-diagonal
        # dyad: on a two-mode network the within-mode dyads are structurally
        # impossible, and round 1 fitted them as observed zeros (nobs = 15 on
        # `network(6; bipartite=3)`, wrong pseudo-likelihood, BIC and
        # conditionals) while the README claimed there was no entry point.
        # Now every entry point refuses, as ERGM.jl's `_refuse_two_mode` does.
        b = network(6; bipartite=3)
        add_edge!(b, 1, 4); set_edge_attribute!(b, :weight, 1, 4, 2)
        add_edge!(b, 2, 5); set_edge_attribute!(b, :weight, 2, 5, 1)
        @test is_two_mode(b)
        bp = BipartiteNetwork(2, 3)
        for (f, ctx) in ((() -> fit_ergm_count(b, [SumTerm(), NonzeroTerm()]), "CountERGMModel"),
                         (() -> ergm_count(b, SumTerm()), "CountERGMModel"),
                         (() -> CountERGMModel([SumTerm()], b), "CountERGMModel"),
                         (() -> fit_ergm_count(bp, [SumTerm()]), "fit_ergm_count"),
                         (() -> CountERGMModel([SumTerm()], bp), "CountERGMModel"),
                         (() -> simulate_count_ergm(b, [SumTerm()], [0.1]; n_sim=1, burnin=1, interval=1),
                          "simulate_count_ergm"),
                         (() -> simulate_count_ergm(bp, [SumTerm()], [0.1]; n_sim=1), "simulate_count_ergm"))
            err = try; f(); nothing catch e; e end
            @test err isa ArgumentError
            msg = sprint(showerror, err)
            @test occursin("$ctx: ERGMCount.jl fits one-mode networks only", msg)
            @test occursin("two-mode (bipartite)", msg)
            @test occursin("within-mode dyads as observed zeros", msg)
            @test occursin("Not implemented", msg)
        end
        # `gof` and the fitted-result simulator inherit the refusal through the
        # result's network (flag it two-mode after the fit, as the missing-data
        # test masks a dyad after the fit)
        net = random_count_net(6; seed=3)
        fit = fit_ergm_count(net, [SumTerm(), NonzeroTerm()])
        @test !is_two_mode(fit.model.network)
        fit.model.network.bipartite = 3
        @test is_two_mode(fit.model.network)
        @test_throws ArgumentError gof(fit; n_sim=2, burnin=1, interval=1)
        @test_throws ArgumentError simulate_count_ergm(fit; n_sim=1, burnin=1, interval=1)
        fit.model.network.bipartite = nothing
        @test length(simulate_count_ergm(fit; n_sim=1, burnin=1, interval=1)) == 1
    end

    @testset "Negative counts under DiscUnif2Reference(a < 0, b)" begin
        # R's DiscUnif(a, b) admits a negative `a`. The estimator was already
        # exact on negative data, but the Gibbs sampler stored a draw of -1 or
        # -2 as an ABSENT edge (`y > 0` decided "edge present"), so every
        # negative value collapsed to 0: gof compared the exact fit with a
        # sampler of a different law, and `se=:bootstrap` was silently wrong.
        # An edge now exists for every non-zero count, negative ones included.
        ref = DiscUnif2Reference(-2, 2)
        seed = network(3; directed=true)
        # θ = 0 on `sum`: every dyad's conditional is exactly uniform on -2:2,
        # so each retained sweep is an iid draw and the per-dyad frequencies
        # are binomial (sd 0.007 at 3000 draws)
        sims = simulate_count_ergm(seed, (SumTerm(),), [0.0]; reference=ref,
                                   n_sim=3000, burnin=5, interval=1, rng=Xoshiro(1))
        for i in 1:3, j in 1:3
            i == j && continue
            vals = [ERGMCount.dyad_value(s, get_edge_attribute(s, :weight, Int), i, j)
                    for s in sims]
            for v in -2:2
                @test count(==(v), vals) / 3000 ≈ 0.2 atol = 0.03
            end
        end
        # A negative count is an edge carrying that weight, read back as such
        s = sims[findfirst(s -> any(v < 0 for v in values(get_edge_attribute(s, :weight))), sims)]
        w = get_edge_attribute(s, :weight, Int)
        @test any(v < 0 for v in values(w))
        @test all(has_edge(s, k...) for (k, v) in w if v != 0)
        @test compute(NonzeroTerm(), s) == count(v -> v != 0, [ERGMCount.dyad_value(s, w, i, j)
                                                            for i in 1:3 for j in 1:3 if i != j])
        # The exact fit on negative data, and the sampler it now agrees with:
        # `sum + nonzero` is dyad-independent, so the MLE sets E[g] = g_obs
        # and 200 iid Gibbs draws must be centred on the observed statistics
        neg = simulate_count_ergm(network(8; directed=true), (SumTerm(), NonzeroTerm()),
                                  [0.3, 0.2]; reference=ref, n_sim=1, burnin=20,
                                  interval=1, rng=Xoshiro(2))[1]
        @test any(v < 0 for v in values(get_edge_attribute(neg, :weight)))
        fit = fit_ergm_count(neg, [SumTerm(), NonzeroTerm()]; reference=ref)
        @test is_exact(fit) && fit.converged && !fit.separated
        @test fit.support_control === :bounded
        @test occursin("Support:   -2:2  (bounded reference)", sprint(show, fit))
        g = gof(fit; n_sim=200, rng=Xoshiro(3))
        st = g.statistics[1]
        for k in 1:2
            μ = mean(st.simulated[:, k]); σ = std(st.simulated[:, k])
            @test abs(μ - st.observed[k]) < 4 * σ / sqrt(200)
            @test st.p_values[k] > 0.05
        end
        @test g.statistics[2].labels == ["-2", "-1", "0", "1", "2"]
        @test sum(g.statistics[2].observed) == 56
        @test all(sum(g.statistics[2].simulated, dims=2) .== 56)
        # ... and the parametric bootstrap of an exact fit reproduces the
        # Hessian standard errors (they are the information-based SEs of the
        # same likelihood), with replicates centred on θ̂
        boot = fit_ergm_count(neg, [SumTerm(), NonzeroTerm()]; reference=ref,
                              se=:bootstrap, n_boot=150, rng=Xoshiro(4))
        ok = [all(isfinite, view(boot.boot_replicates, b, :)) for b in 1:150]
        @test count(ok) >= 120
        for k in 1:2
            @test stderror(boot)[k] ≈ stderror(fit)[k] rtol = 0.3
            @test abs(mean(boot.boot_replicates[ok, k]) - coef(fit)[k]) < 0.5 * stderror(fit)[k]
        end
        # The positive-`a` reference: 0 is impossible, so an absent edge is
        # refused with the wider-reference hint, and the support prints as R's
        cn = network(4; directed=true)
        for i in 1:4, j in 1:4
            i != j && (add_edge!(cn, i, j); set_edge_attribute!(cn, :weight, i, j, 1 + (i + j) % 3))
        end
        pos = fit_ergm_count(cn, [SumTerm()]; reference=DiscUnif2Reference(1, 5))
        @test occursin("Support:   1:5  (bounded reference)", sprint(show, pos))
        @test pos.max_val == 5
        rem_edge!(cn, 1, 2)
        err = try; fit_ergm_count(cn, [SumTerm()]; reference=DiscUnif2Reference(1, 5)); nothing
              catch e; e end
        @test err isa ArgumentError
        @test occursin("DiscUnif2Reference(0, 5)", sprint(showerror, err))
    end

    @testset "Separation: the MPLE does not exist along a combination of statistics" begin
        # On a network whose every count is 0 or 1, `sum − nonzero` (= Σ (y−1)
        # I(y>0)) is at its minimum on every dyad: the pseudo-likelihood is
        # flat along θ_sum → −∞, θ_nonzero → +∞ and no finite MPLE exists.
        # `_count_boundary_columns` sees single columns only; round 1 reported
        # coef [−21.4, 21.3], SE 2.4e4, `converged = true`, `iterations = 1`,
        # no warning — the very number pattern the boundary fix removed. The
        # count analogue of `ERGM._separated` now catches it: an unobserved
        # support value pushed > 18.42 nats down by the tilt alone, AND a
        # Newton step still O(1) (or an uninvertible Hessian).
        s = count_net(6, [(1,2,1), (2,3,1), (3,4,1), (4,5,1), (5,6,1), (1,3,1), (2,5,1)];
                      directed=false)
        @test compute(SumTerm(), s) == compute(NonzeroTerm(), s) == 7
        fs = @test_logs (:warn, r"the MPLE does not exist \(perfect separation\)") match_mode = :any begin
            fit_ergm_count(s, [SumTerm(), NonzeroTerm()])
        end
        @test !fs.converged
        @test fs.separated
        @test !is_exact(fs)
        @test coef(fs)[1] < -15 && coef(fs)[2] > 15          # the asymptote Newton stopped on
        @test fs.hessian_cond > ERGMCount._HESSIAN_COND_TOL   # and it is numerically singular there
        @test fs.collinear == ["sum", "nonzero"]
        @test any(occursin("MPLE does not exist", a) for a in approximations(fs))
        out = sprint(show, fs)
        @test occursin("Converged: false", out)
        @test occursin("MPLE does not exist (perfect separation)", out)
        @test occursin("R ergm warns", out)
        # The separation caveat subsumes the conditioning one (one message
        # about meaningless standard errors, not two)
        @test !occursin("numerically singular", out)
        @test !any(occursin("numerically singular", a) for a in approximations(fs))
        rec = Test.collect_test_logs() do
            fit_ergm_count(s, [SumTerm(), NonzeroTerm()])
        end
        warns = [r.message for r in rec[1] if r.level == Base.CoreLogging.Warn]
        @test any(occursin("The MPLE does not exist!", w) && occursin("sum − nonzero", w) for w in warns)
        @test !any(occursin("numerically singular", w) for w in warns)
        # `warn=false` silences it; the verdict is unchanged
        quiet = @test_logs min_level = Base.CoreLogging.Warn begin
            fit_ergm_count(s, [SumTerm(), NonzeroTerm()]; warn=false)
        end
        @test quiet.separated && !quiet.converged
        # A mixed direction: every count 0 or 2 with `sum + atleast(2)` —
        # neither column is at a boundary (sum is observed at 0 and 2, atleast.2
        # at 0 and 1), but (−1, 2)'(y, I(y ≥ 2)) is maximal exactly at y ∈ {0, 2},
        # so the objective is flat along θ_sum → −∞, θ_atleast → +∞ at 1:2
        m2 = count_net(6, [(1,2,2), (2,3,2), (3,4,2), (4,5,2), (5,6,2), (1,3,2), (2,5,2)];
                       directed=false)
        fm = @test_logs (:warn, r"MPLE does not exist") match_mode = :any begin
            fit_ergm_count(m2, [SumTerm(), CountAtleastnTerm(2)])
        end
        @test fm.separated && !fm.converged
        @test coef(fm)[2] > 2 * abs(coef(fm)[1]) * 0.9      # the (−1, 2) direction
        # ... and a bootstrap replicate on such a design is excluded, never a
        # row of the covariance: a separated refit is `converged = false`
        @test_throws ArgumentError fit_ergm_count(s, [SumTerm(), NonzeroTerm()];
                                                  reference=BinomialReference(1),
                                                  se=:bootstrap, n_boot=4, warn=false)
        # Well-posed fits are NOT flagged: the fixture model on zach (step
        # 3e-12, condition number 54), a wide truncated support whose tilt
        # alone spans > 18 nats (the first signature fires, the second does
        # not), a dyad-dependent model, and the bounded references
        g = load_golden(joinpath(@__DIR__, "fixtures", "zach_poisson.toml"))
        zach = network(Int(g.values["n_actors"]); directed=false)
        for (a, b, w) in zip(Int.(g.values["edge_src"]), Int.(g.values["edge_dst"]),
                             Int.(g.values["edge_weight"]))
            add_edge!(zach, a, b); set_edge_attribute!(zach, :weight, a, b, w)
        end
        for kw in ((;), (; max_val=60), (; reference=GeometricReference()),
                   (; reference=BinomialReference(8)))
            f = fit_ergm_count(zach, [SumTerm(), NonzeroTerm()]; kw...)
            @test f.converged && !f.separated
            @test f.hessian_cond < 1e4
            @test isempty(f.collinear)
        end
        fd = fit_ergm_count(zach, [SumTerm(), NonzeroTerm(), TransitiveWeightsTerm()]; max_val=14)
        @test fd.converged && !fd.separated && fd.hessian_cond < 1e4
        # The detector itself, on the design: separated at the asymptote,
        # not at zach's maximum
        D = ERGMCount._count_design(fs.model, 0:fs.max_val)
        mask = fill(true, length(D.support), length(D.n_tot))
        d = ERGMCount._count_derivatives(D, [1, 2], mask)
        @test ERGMCount._count_separated(d, D, [1, 2], mask, coef(fs), stderror(fs))
        Dz = ERGMCount._count_design(fd.model, 0:14)
        mz = fill(true, 15, length(Dz.n_tot))
        dz = ERGMCount._count_derivatives(Dz, [1, 2, 3], mz)
        @test !ERGMCount._count_separated(dz, Dz, [1, 2, 3], mz, coef(fd), stderror(fd))
    end

    @testset "Collinear statistics: a numerically singular pseudo-Hessian is loud" begin
        # `greaterthan(2)` and `atleast(3)` are the SAME statistic on integer
        # counts: the design has two identical columns, the pseudo-Hessian is
        # exactly singular, and the fit used to say nothing beyond NaN standard
        # errors. Now the condition number is recorded, the two statistics are
        # named, and `show`/`approximations` carry the caveat.
        g = load_golden(joinpath(@__DIR__, "fixtures", "zach_poisson.toml"))
        zach = network(Int(g.values["n_actors"]); directed=false)
        for (a, b, w) in zip(Int.(g.values["edge_src"]), Int.(g.values["edge_dst"]),
                             Int.(g.values["edge_weight"]))
            add_edge!(zach, a, b); set_edge_attribute!(zach, :weight, a, b, w)
        end
        terms = [SumTerm(), NonzeroTerm(), GreaterthannTerm(2), CountAtleastnTerm(3)]
        @test compute(GreaterthannTerm(2), zach) == compute(CountAtleastnTerm(3), zach)
        # SVD roundoff can report a huge finite condition number for an
        # exactly singular design; the numerical-rank policy is the contract.
        fc = @test_logs (:warn, r"numerically singular \(condition number") match_mode = :any begin
            fit_ergm_count(zach, terms; max_val=14)
        end
        @test fc.hessian_cond > ERGMCount._HESSIAN_COND_TOL
        @test fc.collinear == ["greaterthan.2", "atleast.3"]
        @test all(isnan, stderror(fc))
        @test !fc.separated
        rec = Test.collect_test_logs() do
            fit_ergm_count(zach, terms; max_val=14)
        end
        warns = [r.message for r in rec[1] if r.level == Base.CoreLogging.Warn]
        @test any(occursin("loading on greaterthan.2, atleast.3", w) &&
                  occursin("reported as NaN", w) for w in warns)
        out = sprint(show, fc)
        @test occursin("Converged: false", out)
        @test occursin("numerically singular", out)
        @test occursin("loading on greaterthan.2, atleast.3", out)
        @test any(occursin("numerically singular", a) for a in approximations(fc))
        @test !is_exact(fc)
        # The professor's case: a 0/1 network fit as counts (`:weight` = 1 on
        # every edge) — `sum ≡ nonzero` on the data, which round 1 printed as
        # `Converged: true` with SE 2.5e4 and nothing else. It is separated
        # (the two columns differ on the unobserved values), numerically
        # singular, and named
        bn = network(8; directed=true)
        rng = Xoshiro(1)
        for i in 1:8, j in 1:8
            i != j && rand(rng) < 0.3 || continue
            add_edge!(bn, i, j); set_edge_attribute!(bn, :weight, i, j, 1)
        end
        fb = @test_logs (:warn, r"MPLE does not exist") match_mode = :any begin
            fit_ergm_count(bn, [SumTerm(), NonzeroTerm(), CountMutualTerm()])
        end
        @test !fb.converged && fb.separated
        @test fb.hessian_cond > ERGMCount._HESSIAN_COND_TOL
        @test fb.collinear == ["sum", "nonzero"]
        @test !occursin("Converged: true", sprint(show, fb))
        @test ERGMCount._HESSIAN_COND_TOL == 1e8
        @test ERGMCount._hessian_cond(zeros(0, 0)) == 1.0
        @test ERGMCount._hessian_cond([NaN 0.0; 0.0 -1.0]) == Inf
        @test ERGMCount._hessian_cond([-1.0 0.0; 0.0 -1.0]) == 1.0
    end

    @testset "Counts must be integers under :weight (the response= migration)" begin
        # `dyad_value` counts an edge without `:weight` as 1, so a network whose
        # counts live under another name (R's `response="w"`) — or were never
        # attached — fitted silently as a 0/1 network, on which `sum` and
        # `nonzero` coincide and the printed fit is non-identified garbage.
        # Non-integer weights surfaced as a bare InexactError / MethodError
        # from the typed attribute read. The model constructor validates.
        bn = network(8; directed=true)
        rng = Xoshiro(2)
        for i in 1:8, j in 1:8
            i != j && rand(rng) < 0.3 && add_edge!(bn, i, j)
        end
        @test ne(bn) > 0
        err = try; fit_ergm_count(bn, [SumTerm(), NonzeroTerm(), CountMutualTerm()]); nothing
              catch e; e end
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("$(ne(bn)) edges but no `:weight` edge attribute", msg)
        @test occursin("response=\"w\"", msg)
        @test occursin("weight=:w", msg)
        @test occursin("set_edge_attribute!(net, :weight, i, j, w)", msg)
        @test_throws ArgumentError CountERGMModel([SumTerm()], bn)
        # An empty network has no counts to validate (its statistics are at
        # their boundary, which is a different, reported, story)
        @test coef(fit_ergm_count(network(4), [SumTerm()]; warn=false)) == [-Inf]
        # `weight=:w` maps R's `response="w"`: the fit runs on a copy carrying
        # the counts under `:weight`; the caller's network is untouched
        set_edge_attribute!(bn, :w, Dict((src(e), dst(e)) => 1 + (src(e) + dst(e)) % 3
                                          for e in edges(bn)))
        fw = fit_ergm_count(bn, [SumTerm(), NonzeroTerm()]; weight=:w)
        @test fw.converged && !fw.separated
        @test fw.model.network !== bn
        @test !(:weight in list_edge_attributes(bn))
        @test get_edge_attribute(fw.model.network, :weight) == get_edge_attribute(bn, :w)
        explicit = copy(bn)
        set_edge_attribute!(explicit, :weight, get_edge_attribute(bn, :w))
        @test coef(fit_ergm_count(explicit, [SumTerm(), NonzeroTerm()])) == coef(fw)
        @test coef(fit_ergm_count(explicit, [SumTerm(), NonzeroTerm()]; weight=:weight)) == coef(fw)
        @test fit_ergm_count(explicit, [SumTerm()]; weight=:weight).model.network === explicit
        err = try; fit_ergm_count(bn, [SumTerm()]; weight=:counts); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("weight=:counts names no edge attribute", sprint(showerror, err))
        @test occursin(":w", sprint(showerror, err))
        # A PARTIALLY weighted network — one edge without `:weight` among
        # weighted ones (a weight column with NAs, a merge that dropped rows) —
        # used to count that edge as 1 silently (sum 12 → 13 on the grader's
        # 6-node network, coef [0.179, -2.025] → [0.063, -1.653]): the same
        # `response="w"` mistake on a subset of the edges. Refused, naming
        # the first bare edge and how many there are
        pw = count_net(6, [(1, 2, 2), (2, 3, 1), (3, 4, 3), (4, 5, 1), (5, 1, 2), (2, 4, 1)])
        full = coef(fit_ergm_count(pw, [SumTerm(), NonzeroTerm()]))
        add_edge!(pw, 6, 1)                          # no `:weight`
        err = try; fit_ergm_count(pw, [SumTerm(), NonzeroTerm()]); nothing catch e; e end
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("1 of the network's 7 edges carries no `:weight` (the first is (6,1))", msg)
        @test occursin("not a count of 1", msg)
        @test occursin("set_edge_attribute!(net, :weight, i, j, w)", msg)
        @test occursin("weight=:w", msg)
        @test_throws ArgumentError CountERGMModel([SumTerm()], pw)
        add_edge!(pw, 6, 2)
        err = try; fit_ergm_count(pw, [SumTerm()]); nothing catch e; e end
        @test occursin("2 of the network's 8 edges carry no `:weight`", sprint(showerror, err))
        # ... the explicit-specification simulator reads the seed's counts and
        # refuses the same seed
        @test_throws ArgumentError simulate_count_ergm(pw, [SumTerm()], [-0.5]; n_sim=1,
                                                       burnin=1, interval=1)
        # setting the count on every edge (1 if a bare edge means one event)
        # is the fix, and the fit is then the fit of those counts
        set_edge_attribute!(pw, :weight, 6, 1, 1); set_edge_attribute!(pw, :weight, 6, 2, 1)
        @test fit_ergm_count(pw, [SumTerm(), NonzeroTerm()]).converged
        rem_edge!(pw, 6, 1); rem_edge!(pw, 6, 2)
        @test coef(fit_ergm_count(pw, [SumTerm(), NonzeroTerm()])) == full
        # Non-integer and non-numeric weights: the dyad, the value and the rule
        for (bad, shown) in ((2.5, "2.5 (Float64)"), ("3", "\"3\" (String)"), (1.5f0, "1.5f0 (Float32)"))
            net = count_net(4, [(1, 2, bad), (2, 3, 1)])
            err = try; fit_ergm_count(net, [SumTerm()]); nothing catch e; e end
            @test err isa ArgumentError
            msg = sprint(showerror, err)
            @test occursin("counts must be integers; got $shown at dyad (1,2)", msg)
            @test occursin("round or rescale", msg)
        end
        # ... whereas an integer-valued Float is a count (2.0 → 2)
        @test coef(fit_ergm_count(count_net(4, [(1, 2, 2.0), (2, 3, 1)]), [SumTerm()])) ==
              coef(fit_ergm_count(count_net(4, [(1, 2, 2), (2, 3, 1)]), [SumTerm()]))
        # A negative count is refused with the rule, not with "pass max_val ≥ -1"
        neg = count_net(10, [(10, 9, -1), (1, 2, 2)])
        err = try; fit_ergm_count(neg, [SumTerm()]); nothing catch e; e end
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("Observed count -1 at dyad (10,9): counts must be non-negative integers", msg)
        @test occursin("DiscUnif2Reference(a, b)", msg)
        @test !occursin("max_val", msg)
        # ... the support check itself says the same when reached directly
        err = try; ERGMCount._check_in_support(PoissonReference(), 0:10, -1, 10, 9); nothing
              catch e; e end
        @test err isa ArgumentError
        @test occursin("counts must be non-negative integers", sprint(showerror, err))
        @test !occursin("max_val", sprint(showerror, err))
        # ... and under DiscUnif2 with y < a the wider-reference hint is kept
        err = try; fit_ergm_count(neg, [SumTerm()]; reference=DiscUnif2Reference(0, 3)); nothing
              catch e; e end
        @test occursin("DiscUnif2Reference(-1, 3)", sprint(showerror, err))
        @test fit_ergm_count(neg, [SumTerm()]; reference=DiscUnif2Reference(-1, 3)).converged
    end

    @testset "A binary ERGM.jl term is refused at construction, by name" begin
        # `fit_ergm_count(net, [SumTerm(), ERGM.Edges()])` — the first attempt
        # of a user with `ergm(net ~ edges + mutual, response="w")` and ERGM.jl
        # loaded — passed validation and died in the design build with a
        # MethodError on `change_stat_count(::Edges, …)`. R refuses `edges` on
        # a valued response with a named error; so does the model constructor,
        # naming the count analogue
        net = random_count_net(6; seed=21)
        for (t, analogue) in ((ERGM.Edges(), "`NonzeroTerm()`"),
                              (ERGM.Mutual(), "`CountMutualTerm()`"),
                              (ERGM.Triangle(), "`TransitiveWeightsTerm()`"),
                              (ERGM.IStar(2), "`NodeISumTerm()`"))
            err = try; fit_ergm_count(net, [SumTerm(), t]); nothing catch e; e end
            @test err isa ArgumentError
            msg = sprint(showerror, err)
            @test occursin("$(nameof(typeof(t))) is a binary ERGM.jl term", msg)
            @test occursin("change_stat_count", msg)
            @test occursin(analogue, msg)
            @test !occursin("Closest candidates", msg)
            @test_throws ArgumentError CountERGMModel([t], net)
            @test_throws ArgumentError simulate_count_ergm(network(4; directed=true), [t], [0.1];
                                                           n_sim=1, burnin=1, interval=1)
        end
        @test !ERGMCount._is_count_term(ERGM.Edges())
        @test all(ERGMCount._is_count_term(t) for t in ALL_TERMS)
        # The swapped-argument form of the simulator is a named ArgumentError,
        # as `fit_ergm_count(terms, net)` already was — not a bare MethodError
        for swapped in ([SumTerm()], (SumTerm(),), SumTerm())
            err = try; simulate_count_ergm(swapped, net, [0.1]); nothing catch e; e end
            @test err isa ArgumentError
            @test occursin("simulate_count_ergm(net, terms, coefficients): the network comes first",
                           sprint(showerror, err))
        end
        @test occursin("[SumTerm()]", sprint(showerror,
              try; simulate_count_ergm(SumTerm(), net, [0.1]); catch e; e end))
    end

    @testset "compute is type-stable through the typed :weight snapshot" begin
        # `_get_weights` read the untyped `Dict{Tuple{Int,Int},Any}`, so
        # `compute` inferred `Any` for five terms and boxed every `+` (8.4 KB
        # on a 259-edge network). Every term now infers Float64 on both
        # directednesses, and a term that returns 0.0 early on an undirected
        # network is still Float64
        for directed in (true, false)
            net = random_count_net(7; directed=directed, seed=5)
            for t in ALL_TERMS
                @test Base.return_types(compute, Tuple{typeof(t), typeof(net)}) == [Float64]
            end
        end
        @test ERGMCount._get_weights(random_count_net(5; seed=1)) isa Dict{Tuple{Int,Int},Int}
        # An integer-valued Float weight converts in the typed read (2.0 → 2)
        @test compute(SumTerm(), count_net(3, [(1, 2, 2.0), (2, 3, 1)])) == 3.0
    end

    @testset "Negative counts: the terms R refuses are refused, the rest match R" begin
        # ergm 4.12 refuses `transitiveweights`/`cyclicalweights` on a network
        # with a negative dyad weight and its `mutual(form="geometric")`
        # returns NaN; round 2 returned 3.0/7.0 (the `best = 0` floor of the
        # two-path search defining a statistic R never produces) and threw a
        # bare DomainError from `sqrt`. The fixture freezes R's errors and the
        # values of every other term on a seeded 6-actor network in -2:2.
        g = load_golden(joinpath(@__DIR__, "fixtures", "count_terms.toml"))
        neg = network(Int(g.values["negative_n"]); directed=true)
        for (s, d, w) in zip(g.values["negative_edge_src"], g.values["negative_edge_dst"],
                             g.values["negative_edge_weight"])
            add_edge!(neg, Int(s), Int(d))
            set_edge_attribute!(neg, :weight, Int(s), Int(d), Int(w))
        end
        @test any(v < 0 for v in values(get_edge_attribute(neg, :weight)))
        neg_terms = [SumTerm(), NonzeroTerm(), GreaterthannTerm(0), CountAtleastnTerm(0),
                     SmallerthanTerm(0), EqualToTerm(-1), InIntervalTerm(-2, 1),
                     InIntervalTerm(-2, 1; open=(false, false)), CountMutualTerm(:min),
                     CountMutualTerm(:nabsdiff), CountMutualTerm(:product)]
        @test [name(t, neg) for t in neg_terms] == g.values["negative_summary_names"]
        ns = [compute(t, neg) for t in neg_terms]
        @test check_golden(g, "negative_summary", ns) || error(golden_report(g, "negative_summary", ns))
        for (k, t) in enumerate(neg_terms)
            @test compute(t, neg) ≈ g.values["negative_summary"][k] atol = 1e-9
        end
        # ... and the fit of those terms on negative data is well defined
        ref = DiscUnif2Reference(-2, 2)
        fneg = fit_ergm_count(neg, [SumTerm(), NonzeroTerm(), CountMutualTerm(:product)]; reference=ref)
        @test fneg.converged && fneg.support_control === :bounded
        # The three R refuses: at `compute` (the summary statistic, as in R),
        # at model construction, at simulation over a support with negative
        # values — never a number
        r_sentence = "may not be used with networks with negative dyad weights"
        @test occursin(r_sentence, g.values["r_error_transitiveweights_negative"])
        @test occursin(r_sentence, g.values["r_error_cyclicalweights_negative"])
        @test g.values["r_mutual_geometric_negative"] == "NaN"
        for t in (TransitiveWeightsTerm(), CyclicalWeightsTerm(), CountMutualTerm(:geometric))
            err = try; compute(t, neg); nothing catch e; e end
            @test err isa ArgumentError
            msg = sprint(showerror, err)
            @test occursin("compute: $(nameof(typeof(t))) (`$(name(t))`) $r_sentence", msg)
            @test occursin(t isa CountMutualTerm ? "returns NaN" : "R's `ergm` refuses", msg)
            @test !occursin("DomainError", msg)
            err = try; CountERGMModel([SumTerm(), t], neg, ref); nothing catch e; e end
            @test err isa ArgumentError
            @test occursin("CountERGMModel: $(nameof(typeof(t)))", sprint(showerror, err))
            @test_throws ArgumentError fit_ergm_count(neg, [SumTerm(), t]; reference=ref)
            # non-negative data under a reference whose support reaches below
            # 0: the estimator would enumerate the negative values
            @test_throws ArgumentError CountERGMModel([SumTerm(), t], random_count_net(5; seed=3), ref)
            err = try
                simulate_count_ergm(network(5; directed=true), [SumTerm(), t], [0.0, 0.1];
                                    reference=ref, n_sim=1, burnin=1, interval=1)
                nothing
            catch e; e end
            @test err isa ArgumentError
            @test occursin("simulate_count_ergm: $(nameof(typeof(t)))", sprint(showerror, err))
            @test !ERGMCount._admits_negative(t)
        end
        # ... whereas on non-negative data (and a non-negative support) the
        # same terms are ordinary, and the other `mutual` forms admit negatives
        pos = random_count_net(6; seed=9)
        @test compute(TransitiveWeightsTerm(), pos) >= 0
        @test fit_ergm_count(pos, [SumTerm(), CountMutualTerm(:geometric)];
                             reference=DiscUnif2Reference(0, 4)).converged
        @test all(ERGMCount._admits_negative(t) for t in
                  (SumTerm(), CountMutualTerm(:min), CountMutualTerm(:nabsdiff),
                   CountMutualTerm(:product), TransitiveTiesTerm(), NodeSumTerm()))
    end

    @testset "Adaptive support stops on a rung whose Newton failed" begin
        # `sum + greaterthan(2) + atleast(3)` is exactly collinear on integer
        # counts. Round 2 ran all 8 doublings on zach (max_val 14 → 3584, the
        # design growing 2^8×), each rung's Newton breaking at iteration 1 on
        # the singular Hessian, and then reported `support_control =
        # :unconverged` with the "not normalisable" sentence — a false
        # diagnosis of a collinear design. The doubling now stops at the
        # first rung whose Newton failed (or whose pseudo-Hessian is
        # numerically singular), and the convergence/conditioning warning
        # speaks alone.
        g = load_golden(joinpath(@__DIR__, "fixtures", "zach_poisson.toml"))
        zach = network(Int(g.values["n_actors"]); directed=false)
        for (a, b, w) in zip(Int.(g.values["edge_src"]), Int.(g.values["edge_dst"]),
                             Int.(g.values["edge_weight"]))
            add_edge!(zach, a, b); set_edge_attribute!(zach, :weight, a, b, w)
        end
        terms = [SumTerm(), GreaterthannTerm(2), CountAtleastnTerm(3)]
        rec = Test.collect_test_logs() do
            fit_ergm_count(zach, terms)
        end
        fc = fit_ergm_count(zach, terms; warn=false)
        @test fc.max_val <= 28                       # the first rung, or the first doubling
        @test fc.support_control === :unconverged && !fc.support_stable
        @test isnan(fc.support_delta) && isnan(fc.omitted_tail)
        @test !fc.converged && fc.hessian_cond > ERGMCount._HESSIAN_COND_TOL
        @test fc.collinear == ["greaterthan.2", "atleast.3"]
        warns = [r.message for r in rec[1] if r.level == Base.CoreLogging.Warn]
        @test !any(occursin("normalisable", w) for w in warns)
        @test any(occursin("stopped at max_val = $(fc.max_val)", w) &&
                  occursin("did not converge", w) && occursin("NOT error-controlled", w)
                  for w in warns)
        @test any(occursin("numerically singular", w) && occursin("greaterthan.2, atleast.3", w)
                  for w in warns)
        @test !any(occursin("normalisable", a) for a in approximations(fc))
        @test any(occursin("support doubling stopped at max_val = $(fc.max_val)", a) &&
                  occursin("NOT error-controlled", a) for a in approximations(fc))
        @test any(occursin("numerically singular", a) for a in approximations(fc))
        out = sprint(show, fc)
        @test occursin("adaptive doubling stopped here because this fit did not converge", out)
        @test !occursin("NaN·SE", out)
        @test !is_exact(fc)
        # The rung test itself
        @test !ERGMCount._support_rung_usable((converged=false, hessian_cond=1.0))
        @test !ERGMCount._support_rung_usable((converged=true, hessian_cond=Inf))
        @test !ERGMCount._support_rung_usable((converged=true, hessian_cond=1e9))
        @test ERGMCount._support_rung_usable((converged=true, hessian_cond=1e7))
        # A well-posed model on the same data still doubles to convergence
        # (the fixture model: 14 → 28), so the stop is not a regression of
        # the error-controlled path
        ok = fit_ergm_count(zach, [SumTerm(), NonzeroTerm()])
        @test ok.support_control === :converged && ok.max_val == 28
        # ... and a separated fit — `converged = false` on its asymptote —
        # stops at its first rung too, with the separation caveat alone
        s = count_net(6, [(1,2,1), (2,3,1), (3,4,1), (4,5,1), (5,6,1), (1,3,1), (2,5,1)];
                      directed=false)
        fs = fit_ergm_count(s, [SumTerm(), NonzeroTerm()]; warn=false)
        @test fs.separated && fs.max_val == 10 && fs.support_control === :unconverged
        @test isnan(fs.support_delta)
    end

    @testset "Diagnostic numbers print with three significant digits, no noise" begin
        # `round(x, sigdigits=3)` printed `6.969999999999999e-32` for the
        # docs' seeded example ("these are the numbers you get" — they were
        # not); every diagnostic number goes through `_fmt3` (`%.3g`)
        @test ERGMCount._fmt3(6.969999999999999e-32) == "6.97e-32"
        @test ERGMCount._fmt3(1.6699999999999998e33) == "1.67e+33"
        @test ERGMCount._fmt3(0.218) == "0.218"
        @test ERGMCount._fmt3(47.31) == "47.3"
        @test ERGMCount._fmt3(Inf) == "Inf" && ERGMCount._fmt3(NaN) == "NaN"
        @test ERGMCount._fmt2(100 * BOUNDARY_MASS_TOL) == "0.01"
        noise = r"\d\.\d{5,}e[-+]?\d"       # a mantissa with five or more decimals
        net = random_count_net(12; seed=42, p=0.25)
        for fit in (fit_ergm_count(net, [SumTerm(), NonzeroTerm(), CountMutualTerm()]),
                    fit_ergm_count(net, [SumTerm(), GreaterthannTerm(2), CountAtleastnTerm(3)];
                                   warn=false),
                    fit_ergm_count(net, [SumTerm(), NonzeroTerm()]; max_val=4, warn=false))
            out = sprint(show, fit)
            @test !occursin(noise, out)
            @test !any(occursin(noise, a) for a in approximations(fit))
            @test occursin(r"Boundary mass: \d\.\d\de-\d+|Boundary mass: 0\.\d+", out)
        end
        # ... and the warnings too (the truncation warning of a bound that bites)
        rec = Test.collect_test_logs() do
            fit_ergm_count(net, [SumTerm(), NonzeroTerm()]; max_val=4)
        end
        warns = [r.message for r in rec[1] if r.level == Base.CoreLogging.Warn]
        @test any(occursin("% of the", w) for w in warns)
        @test !any(occursin(noise, w) for w in warns)
    end

    @testset "Every export's docstring carries a runnable example" begin
        # Criterion 5: every export has a docstring with a runnable example.
        # The StatsAPI verbs re-exported from StatsAPI (`coef`, `stderror`,
        # `vcov`, `loglikelihood`, `nobs`, `dof`) carry StatsAPI's docstrings;
        # the ERGMCount-specific ones (`aic`, `bic`, `confint`, `coeftable`)
        # and every other export are documented here, with a ```julia fence
        reexported = (:coef, :stderror, :vcov, :loglikelihood, :nobs, :dof)
        # The raw text of every docstring attached to a binding, looked up in
        # the module that owns it (`aliasof` follows a re-export) — the REPL's
        # `doc` renderer is not a test dependency
        function docstring_text(mod, sym)
            b = Base.Docs.aliasof(Base.Docs.Binding(mod, sym))
            out = String[]
            for m in unique((mod, b.mod))
                d = Base.Docs.meta(m; autoinit=false)
                (d === nothing || !haskey(d, b)) && continue
                for ds in values(d[b].docs)
                    push!(out, join(String[x for x in ds.text if x isa AbstractString]))
                end
            end
            return join(out, "\n")
        end
        for sym in names(ERGMCount)
            sym === :ERGMCount && continue
            doc = docstring_text(ERGMCount, sym)
            @test !isempty(doc) || error("no docstring on $sym")
            sym in reexported && continue
            @test occursin("```julia", doc) || error("no runnable example in the docstring of $sym")
        end
        # ... and every ```julia fence in the source docstrings executes as
        # written (a fresh module per block, warnings silenced)
        src = read(joinpath(dirname(@__DIR__), "src", "ERGMCount.jl"), String)
        # Git may check out source files with CRLF line endings on Windows.
        blocks = [String(m.captures[1]) for m in eachmatch(r"```julia\r?\n(.*?)```"s, src)]
        @test length(blocks) >= 40
        for (k, code) in enumerate(blocks)
            m = Module(Symbol("DocBlock", k))
            ok = try
                Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
                    Core.eval(m, :(using Networks, ERGM, ERGMCount, Random))
                    Core.eval(m, Meta.parseall(code))
                end
                true
            catch e
                @error "docstring block $k failed" code exception = (e, catch_backtrace())
                false
            end
            @test ok
        end
    end

    @testset "show: the model summary and the docstrings describe the current rules" begin
        # `CountERGMModel` prints like `ERGM.ERGMModel`: size, directedness,
        # the R coefficient labels and the reference — never the whole network
        net = random_count_net(6; seed=8)
        model = CountERGMModel([SumTerm(), NonzeroTerm(), CountMutualTerm()], net)
        @test sprint(show, model) ==
              "CountERGMModel{Int64,true}: 6 vertices, $(ne(net)) edges (directed); " *
              "terms: sum + nonzero + mutual.min; reference: PoissonReference(1.0)"
        unet = random_count_net(5; directed=false, seed=8)
        @test sprint(show, CountERGMModel(SumTerm(), unet, BinomialReference(4))) ==
              "CountERGMModel{Int64,false}: 5 vertices, $(ne(unet)) edges (undirected); " *
              "terms: sum; reference: BinomialReference(4)"
        @test !occursin("Network{", sprint(show, model))
        @test is_directed(typeof(model))
        # The user-visible docstrings state the WP2 conjunction (delta AND
        # omitted tail AND boundary mass), not the superseded "or" rule, and
        # the dyad-independent term list is complete
        rdoc = string(@doc CountERGMResult)      # the DocStr repr: raw text, `\n`-escaped
        @test occursin("**and** left at most", rdoc)
        @test occursin("**and** no dyad put more", rdoc)
        @test occursin("`BOUNDARY_MASS_TOL`", rdoc)
        @test !occursin("errors, or left at most", rdoc)
        fdoc = string(@doc fit_ergm_count)
        for t in ("SumTerm", "NonzeroTerm", "GreaterthannTerm", "CountAtleastnTerm",
                  "SmallerthanTerm", "EqualToTerm", "InIntervalTerm")
            @test occursin("`$t`", fdoc)
        end
        @test occursin("weight::Symbol=:weight", fdoc)
        # The result fields the new diagnostics live in
        fit = fit_ergm_count(net, [SumTerm(), NonzeroTerm()])
        @test hasfield(typeof(fit), :separated) && hasfield(typeof(fit), :hessian_cond) &&
              hasfield(typeof(fit), :collinear)
        @test !fit.separated && fit.hessian_cond >= 1 && isempty(fit.collinear)
    end

    # ------------------------------------------------------------------
    # R-parity terms (panel 2026-09 item 31, WP4): the valued `mutual` forms,
    # `transitiveweights`/`cyclicalweights` with R's default (min, max, min)
    # triple, and the dyad-independent `smallerthan`/`equalto`/`ininterval`.
    # Hand values here; the R numbers live in the count_terms fixture below.
    # ------------------------------------------------------------------
    @testset "R-parity terms: mutual forms, weights, thresholds" begin
        # 1→2 (3), 2→1 (2), 2→3 (2), 1→3 (1): one reciprocated pair
        net = count_net(3, [(1, 2, 3), (2, 1, 2), (2, 3, 2), (1, 3, 1)])
        @test compute(CountMutualTerm(), net) == 2.0
        @test compute(CountMutualTerm(:min), net) == 2.0
        @test compute(CountMutualTerm(:nabsdiff), net) == -(1 + 2 + 1)
        @test compute(CountMutualTerm(:geometric), net) ≈ sqrt(6)
        @test compute(CountMutualTerm(:product), net) == 6.0
        @test compute(CountMutualTerm(:threshold; threshold=2), net) == 1.0
        @test compute(CountMutualTerm(:threshold; threshold=3), net) == 0.0
        # threshold <= 0: every pair of an EMPTY network is mutual (R's
        # emptynwstats rule), i.e. the definition is `>=`, not `>`
        @test compute(CountMutualTerm(:threshold; threshold=0), network(4; directed=true)) == 6.0
        @test CountMutualTerm() == CountMutualTerm(:min; threshold=0)
        # R's labels
        @test name(CountMutualTerm()) == "mutual.min"
        @test name(CountMutualTerm(:nabsdiff)) == "mutual.nabsdiff"
        @test name(CountMutualTerm(:geometric)) == "mutual.geom.mean"
        @test name(CountMutualTerm(:product)) == "mutual.product"
        @test name(CountMutualTerm(:threshold; threshold=2)) == "mutual.2"
        @test_throws ArgumentError CountMutualTerm(:sum)
        # Every form is dyad-dependent and directed-only
        @test all(is_dyad_dependent(CountMutualTerm(f)) for f in (:min, :nabsdiff, :geometric, :product, :threshold))
        @test all(requires_directed(CountMutualTerm(f)) for f in (:min, :nabsdiff, :geometric, :product, :threshold))
        @test_throws ArgumentError CountERGMModel([SumTerm(), CountMutualTerm(:product)], network(4; directed=false))

        # transitiveweights(min, max, min): each dyad capped by its strongest
        # two-path. 1→2 (3), 2→3 (2), 1→3 (1): pair (1,3) has the path 1→2→3 of
        # strength min(3,2) = 2, capped by y_13 = 1; no other pair has a path
        tw = count_net(3, [(1, 2, 3), (2, 3, 2), (1, 3, 1)])
        @test compute(TransitiveWeightsTerm(), tw) == 1.0
        @test compute(CyclicalWeightsTerm(), tw) == 0.0         # no 3-cycle
        @test name(TransitiveWeightsTerm()) == "transitiveweights.min.max.min"
        @test name(CyclicalWeightsTerm()) == "cyclicalweights.min.max.min"
        # A 3-cycle 1→2→3→1 with counts 3, 2, 2: each dyad is capped by the
        # two-path closing its cycle (min of the other two), so 2 + 2 + 2
        cy = count_net(3, [(1, 2, 3), (2, 3, 2), (3, 1, 2)])
        @test compute(CyclicalWeightsTerm(), cy) == 6.0
        @test compute(TransitiveWeightsTerm(), cy) == 0.0         # no transitive 2-path
        # Undirected: the cycle and the transitive two-path coincide, so the
        # two statistics are equal there (R allows both terms on an undirected
        # network) — and the sum is over UNORDERED pairs, like R's C code
        un = random_count_net(8; directed=false, seed=21)
        @test compute(CyclicalWeightsTerm(), un) == compute(TransitiveWeightsTerm(), un)
        @test !requires_directed(CyclicalWeightsTerm()) && !requires_directed(TransitiveWeightsTerm())
        @test is_dyad_dependent(TransitiveWeightsTerm()) && is_dyad_dependent(CyclicalWeightsTerm())
        # ... and they are NOT the triple-wise-minimum terms, which have no R
        # counterpart (an undirected triangle counts six times there)
        @test compute(TransitiveTiesTerm(), un) != compute(TransitiveWeightsTerm(), un)
        @test name(TransitiveTiesTerm()) == "transitiveties.count"

        # Dyad-independent threshold terms count the ZERO dyads too
        th = count_net(3, [(1, 2, 1), (2, 3, 2), (3, 1, 3)])   # 6 dyads, 3 empty
        @test compute(SmallerthanTerm(2), th) == 4.0            # 0,0,0,1
        @test compute(SmallerthanTerm(0), th) == 0.0
        @test compute(EqualToTerm(0), th) == 3.0
        @test compute(EqualToTerm(2), th) == 1.0
        @test compute(EqualToTerm(2; tolerance=1), th) == 3.0   # 1, 2, 3
        @test compute(InIntervalTerm(1, 3), th) == 1.0                        # (1,3): only 2
        @test compute(InIntervalTerm(1, 3; open=(false, false)), th) == 3.0   # [1,3]
        @test compute(InIntervalTerm(1, 3; open=(true, false)), th) == 2.0    # (1,3]
        @test compute(InIntervalTerm(1, 3; open=(false, true)), th) == 2.0    # [1,3)
        @test compute(InIntervalTerm(-Inf, Inf), th) == 6.0
        @test compute(InIntervalTerm(0, Inf; open=(false, true)), th) == 6.0
        @test name(SmallerthanTerm(2)) == "smallerthan.2"
        @test name(EqualToTerm(3)) == "equalto.3.pm.0"
        @test name(EqualToTerm(3; tolerance=1)) == "equalto.3.pm.1"
        @test name(InIntervalTerm(1, 3)) == "ininterval(1,3)"
        @test name(InIntervalTerm(1, 3; open=(false, false))) == "ininterval[1,3]"
        @test name(InIntervalTerm(1, 3; open=(true, false))) == "ininterval(1,3]"
        @test name(InIntervalTerm(1.5, Inf)) == "ininterval(1.5,Inf)"
        @test_throws ArgumentError EqualToTerm(3; tolerance=-1)
        @test_throws ArgumentError InIntervalTerm(3, 1)
        for t in (SmallerthanTerm(2), EqualToTerm(3), InIntervalTerm(1, 3))
            @test !is_dyad_dependent(t)
            @test !has_dyad_dependent(CountERGMModel([SumTerm(), t], th))
        end
        # ... so a model made of them is exact under a bounded reference
        fx = fit_ergm_count(random_count_net(7; seed=31), [SumTerm(), InIntervalTerm(1, 3)];
                            reference=BinomialReference(4))
        @test is_exact(fx) && fx.converged

        # A boundary statistic among the new terms, R's drop semantics: on a
        # complete count network every dyad is >= 1, so `smallerthan.1` = 0 is
        # at its smallest attainable value → -Inf; the restricted supports are
        # 1:max_val and `sum` is the zero-truncated Poisson MLE (same closed
        # form as the `nonzero` +Inf case): λ/(1 − e^{-λ}) = ȳ = 9/6
        cn = count_net(3, [(1,2,1), (2,1,2), (1,3,1), (3,1,1), (2,3,3), (3,2,1)])
        fs = @test_logs (:warn, r"smallerthan.1 are at their smallest attainable value") match_mode = :any begin
            fit_ergm_count(cn, [SumTerm(), SmallerthanTerm(1)])
        end
        @test coef(fs)[2] == -Inf && stderror(fs)[2] == 0.0 && fs.p_values[2] == 0.0
        λ = exp(coef(fs)[1])
        @test λ / (1 - exp(-λ)) ≈ 9 / 6 atol = 1e-6
        @test dof(fs) == 1 && !is_exact(fs)
        @test coeftable(fs)["smallerthan.1"].estimate == -Inf
        @test_throws ArgumentError simulate_count_ergm(fs; n_sim=1, burnin=1, interval=1)
    end

    # ------------------------------------------------------------------
    # Golden fixture: `ergm`/`ergm.count` TERM PARITY (count_terms.toml).
    # test/fixtures/r/count_terms.R regenerates it. Every value is a
    # deterministic function of the observed graph, so the comparison is at
    # machine precision (1e-9) and the coefficient NAMES are compared exactly.
    # ------------------------------------------------------------------
    @testset "Golden fixture: ergm term parity on zach and a directed network" begin
        g = load_golden(joinpath(@__DIR__, "fixtures", "count_terms.toml"))
        @test g.provenance["ergm_version"] == "4.12.0"

        function rebuild(prefix; directed)
            net = network(Int(g.values["$(prefix)_n"]); directed=directed)
            for (s, d, w) in zip(g.values["$(prefix)_edge_src"], g.values["$(prefix)_edge_dst"],
                                 g.values["$(prefix)_edge_weight"])
                add_edge!(net, Int(s), Int(d))
                set_edge_attribute!(net, :weight, Int(s), Int(d), Int(w))
            end
            return net
        end

        # Every term of the two formulas has an exact ERGMCount counterpart,
        # in R's order; the names must be R's, the values R's
        zach_terms = [SumTerm(), NonzeroTerm(), GreaterthannTerm(2), GreaterthannTerm(4),
                      CountAtleastnTerm(3), TransitiveWeightsTerm(), CyclicalWeightsTerm(),
                      SmallerthanTerm(2), EqualToTerm(3), InIntervalTerm(1, 3),
                      InIntervalTerm(1, 3; open=(false, false)),
                      # thresholds at or below 0 count the ZERO dyads (all 561)
                      GreaterthannTerm(-1), CountAtleastnTerm(0)]
        zach = rebuild("zach"; directed=false)
        @test !is_directed(zach) && ne(zach) == 78
        @test [name(t, zach) for t in zach_terms] == g.values["zach_summary_names"]
        zs = [compute(t, zach) for t in zach_terms]
        @test check_golden(g, "zach_summary", zs) || error(golden_report(g, "zach_summary", zs))
        # ... row by row too, so a red test names the term
        for (k, t) in enumerate(zach_terms)
            @test compute(t, zach) ≈ g.values["zach_summary"][k] atol = 1e-9
        end

        dir_terms = [SumTerm(), NonzeroTerm(), GreaterthannTerm(2), GreaterthannTerm(4),
                     CountAtleastnTerm(3), TransitiveWeightsTerm(), CyclicalWeightsTerm(),
                     CountMutualTerm(:min), CountMutualTerm(:nabsdiff),
                     CountMutualTerm(:geometric), CountMutualTerm(:product),
                     SmallerthanTerm(2), EqualToTerm(3), InIntervalTerm(1, 3),
                     InIntervalTerm(1, 3; open=(false, false)),
                     InIntervalTerm(1, 3; open=(true, false)),
                     InIntervalTerm(1, 3; open=(false, true))]
        dn = rebuild("directed"; directed=true)
        @test is_directed(dn) && nv(dn) == 8
        @test [name(t, dn) for t in dir_terms] == g.values["directed_summary_names"]
        ds = [compute(t, dn) for t in dir_terms]
        @test check_golden(g, "directed_summary", ds) || error(golden_report(g, "directed_summary", ds))
        for (k, t) in enumerate(dir_terms)
            @test compute(t, dn) ≈ g.values["directed_summary"][k] atol = 1e-9
        end
        # The directed transitive and cyclical weights differ (35 vs 36 in R):
        # ordered-pair sums with different closing two-paths
        @test compute(TransitiveWeightsTerm(), dn) != compute(CyclicalWeightsTerm(), dn)
        # A fitted model labels its coefficients with the same R strings
        fit = fit_ergm_count(dn, [SumTerm(), CountMutualTerm(:nabsdiff), InIntervalTerm(1, 3)];
                             reference=BinomialReference(6))
        @test coeftable(fit).names == ["sum", "mutual.nabsdiff", "ininterval(1,3)"]

        # A tie whose VALUE is 0 is a zero dyad to R's valued terms (nonzero =
        # 1 on a network object holding 2 ties), and it counts with the empty
        # dyads for every threshold term; `compute` agrees — it reads
        # `dyad_value`, never `ne(net)` — so gof's observed statistics and the
        # estimator's view of the data cannot disagree on such a network
        zw = rebuild("zero_weight"; directed=true)
        @test ne(zw) == 2 && get_edge_attribute(zw, :weight, 2, 3) == 0
        zw_terms = [NonzeroTerm(), SumTerm(), GreaterthannTerm(-1), CountAtleastnTerm(0),
                    GreaterthannTerm(0), CountAtleastnTerm(1), SmallerthanTerm(1),
                    EqualToTerm(0), InIntervalTerm(-1, 1)]
        @test [name(t, zw) for t in zw_terms] == g.values["zero_weight_summary_names"]
        zws = [compute(t, zw) for t in zw_terms]
        @test check_golden(g, "zero_weight_summary", zws) ||
              error(golden_report(g, "zero_weight_summary", zws))
        @test compute(NonzeroTerm(), zw) == 1.0 != ne(zw)
        # ... and the estimator's design sees the same dyad values
        @test ERGMCount.dyad_value(zw, get_edge_attribute(zw, :weight, Int), 2, 3) == 0
        zfit = fit_ergm_count(zw, [SumTerm()]; reference=BinomialReference(3))
        @test is_exact(zfit)
        gz = gof(zfit; n_sim=4, burnin=2, interval=1, rng=Xoshiro(1))
        @test gz.statistics[1].observed == [2.0]
        @test gz.statistics[2].observed[1] == 5.0            # five zero dyads, the 0-tie included

        # mutual(form="threshold") is NOT pinned: ergm 4.12.0's own summary
        # fails at C model initialisation, and the fixture records why
        @test occursin("mutual_wt_threshold", g.values["r_error_mutual_threshold"])
        @test haskey(g.provenance, "not_pinned")
        # ... so its value here is the documented definition, by hand: the
        # number of unordered pairs with both counts >= 2
        W = zeros(Int, 8, 8)
        for e in edges(dn)
            W[src(e), dst(e)] = ERGMCount.dyad_value(dn, get_edge_attribute(dn, :weight), src(e), dst(e))
        end
        @test compute(CountMutualTerm(:threshold; threshold=2), dn) ==
              count((W[i, j] >= 2) & (W[j, i] >= 2) for i in 1:8 for j in (i+1):8)
    end

    # ------------------------------------------------------------------
    # Golden fixture: statnet `ergm.count` on Zachary's karate club (issue #8).
    # test/fixtures/r/zach_poisson.R regenerates it.
    #
    # `sum + nonzero` under a Poisson reference is DYAD-INDEPENDENT, so each dyad
    # is an independent draw from a two-parameter law with an EXACT MLE. That is
    # the whole design: ERGMCount.jl's count MPLE enumerates each dyad's full
    # conditional, and for a dyad-independent model the conditional IS the
    # marginal — so the pseudo-likelihood IS the likelihood and ERGMCount.jl is
    # computing the exact MLE. It must therefore agree with R at optimizer
    # precision, not "within Monte-Carlo error".
    # ------------------------------------------------------------------
    @testset "Golden fixture: ergm.count on zach, Poisson reference" begin
        g = load_golden(joinpath(@__DIR__, "fixtures", "zach_poisson.toml"))
        @test g.provenance["ergm_count_version"] == "4.1.3"

        # Rebuild R's zach exactly from the frozen valued edge list.
        n = Int(g.values["n_actors"])
        net = network(n; directed=false)
        s = Int.(g.values["edge_src"])
        d = Int.(g.values["edge_dst"])
        w = Int.(g.values["edge_weight"])
        for k in eachindex(s)
            add_edge!(net, s[k], d[k])
            set_edge_attribute!(net, :weight, s[k], d[k], w[k])
        end
        max_val = Int(g.values["max_val"])

        # Sufficient statistics: deterministic, so machine precision.
        stats = [compute(SumTerm(), net), compute(NonzeroTerm(), net)]
        @test check_golden(g, "summary_statistics", stats) ||
              error(golden_report(g, "summary_statistics", stats))

        fit = fit_ergm_count(net, [SumTerm(), NonzeroTerm()];
                             reference=PoissonReference(), max_val=max_val)
        @test fit.support_control === :fixed

        # --- the DEFAULT path, with no max_val (N6, item 17) ------------------
        # The old default `max(10, 2·7) = 14` landed 1.4e-6/3.9e-6 from the exact
        # MLE — above this tolerance — and only the pinned max_val=30 passed. The
        # error-controlled doubling starts at 14, refits at 28, sees the
        # estimates move by ~2e-5 SE (< support_tol = 1e-3) and stops; the
        # reported fit at 28 is ~2e-12 from the exact MLE.
        dflt = fit_ergm_count(net, [SumTerm(), NonzeroTerm()];
                              reference=PoissonReference())
        @test dflt.support_control === :converged
        @test dflt.support_stable
        @test dflt.max_val >= 20
        @test dflt.max_val == 28               # informative, not a contract
        @test dflt.support_delta <= dflt.support_tol
        @test dflt.omitted_tail <= dflt.support_tol
        @test dflt.converged && is_exact(dflt) == false   # truncated, by design
        # An unmeetable tolerance with one doubling: |Δθ| between 14 and 28 is
        # ~2e-5 SE ≫ 1e-300, so the support is reported unstable, loudly
        unc = @test_logs (:warn, r"error-controlled support did not converge") match_mode = :any begin
            fit_ergm_count(net, [SumTerm(), NonzeroTerm()]; reference=PoissonReference(),
                           support_tol=1e-300, max_doublings=1)
        end
        @test !unc.support_stable && unc.support_control === :unconverged
        @test unc.max_val == 28
        @test unc.support_delta > 1e-300
        @test coef(unc) == coef(dflt)          # same final fit, different verdict
        @test check_golden(g, "exact_coefficients", dflt.coefficients) ||
              error(golden_report(g, "exact_coefficients", dflt.coefficients))
        @test check_golden(g, "exact_std_errors", dflt.std_errors) ||
              error(golden_report(g, "exact_std_errors", dflt.std_errors))
        @test maximum(abs.(dflt.coefficients .- Float64.(g.values["exact_coefficients"]))) < 1e-9
        # ... whereas the OLD default bound, forced, does NOT meet the tolerance:
        # the doubling is what earns it, not luck
        old14 = fit_ergm_count(net, [SumTerm(), NonzeroTerm()];
                               reference=PoissonReference(), max_val=14)
        @test maximum(abs.(old14.coefficients .- Float64.(g.values["exact_coefficients"]))) > 1e-6
        @test check_golden(g, "summary_statistics", stats)

        # --- the assertion with teeth: the EXACT MLE, at 1e-6 ----------------
        # R's value here is not an `optim` output — the score equations were
        # solved analytically, and the fixture freezes the residual score (~1e-14)
        # to prove the golden number is exact enough to police this tolerance.
        # Observed: ERGMCount.jl reproduces it to ~1e-12.
        @test check_golden(g, "exact_coefficients", fit.coefficients) ||
              error(golden_report(g, "exact_coefficients", fit.coefficients))
        @test check_golden(g, "exact_std_errors", fit.std_errors) ||
              error(golden_report(g, "exact_std_errors", fit.std_errors))
        @test fit.loglik ≈ g.values["exact_loglik"] atol = 1e-6

        # --- TRUNCATION: the thing that could silently void the comparison ----
        # The Poisson reference is unbounded and ERGMCount.jl enumerates only
        # 0:max_val, while the exact MLE above truncates nothing. They estimate
        # the same quantity ONLY IF the mass past max_val is negligible. It is:
        # R computes P(y > 30) = 7e-23 under the fitted law, nineteen orders of
        # magnitude below the smallest conditional probability the estimator
        # actually uses. ERGMCount.jl's own reported `boundary_mass` must agree.
        # If a future max_val ever started to bite, this goes red rather than
        # quietly widening the gap above.
        @test fit.truncated
        @test fit.boundary_mass < 1e-15
        @test fit.boundary_mass < 1e4 * Float64(g.values["boundary_tail_mass"])
        @test Float64(g.values["boundary_tail_mass"]) <
              1e-15 * Float64(g.values["min_used_conditional_prob"])

        # --- ergm.count's OWN fit (MCMLE) ------------------------------------
        # statnet has no MPLE for valued ERGMs, so ergm.count reaches this model
        # by MCMC and carries Monte-Carlo error the exact MLE does not. It agrees
        # — but note which way round the error runs: ergm.count's own MCMLE sits
        # ~0.010 from the exact MLE (2-4x its seed-to-seed sd), while
        # ERGMCount.jl sits ~1e-12 from it. The exact value is the reference
        # standard; this is a consistency check on R.
        @test check_golden(g, "mcmle_coefficients", fit.coefficients) ||
              error(golden_report(g, "mcmle_coefficients", fit.coefficients))
        jl_gap = maximum(abs.(fit.coefficients .-
                              Float64.(g.values["exact_coefficients"])))
        @test jl_gap < Float64(g.values["mcmle_vs_exact_max_abs_diff"])
    end
end
