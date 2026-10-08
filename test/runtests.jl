using ERGMCount
using ERGM
using NetworkCore
using Graphs
using Distributions
using Random
using Statistics
using LinearAlgebra: norm, dot, diag
using Test
using Aqua

# Text files read by the tests are compared line by line; a Windows checkout
# (git's core.autocrlf) gives them CRLF endings, so normalise to LF.
_readtext(path) = replace(read(path, String), "\r\n" => "\n")

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

# The workflow testset runs the CI layout step against the sibling checkouts
# that [sources] names. A lone checkout or a registry install (where [sources]
# is not used) has none of them beside it: the step is then not run, with a
# message. Inside the layout every sibling must be present, so a partial
# layout runs the step and fails on the missing one rather than skipping.
function layout_siblings(pkgdir::AbstractString, expected, pkg::AbstractString)
    siblings = sort!([s for s in expected if s != "$pkg.jl"])
    present = filter(s -> isfile(joinpath(dirname(pkgdir), s, "Project.toml")), siblings)
    return (siblings=siblings, present=present, in_layout=!isempty(present))
end

# A count term the impropriety rule knows nothing about (its growth could
# offset another term's), defined at top level for the testset
struct _UnknownCountTerm <: ERGM.AbstractERGMTerm end

const ALL_TERMS = [SumTerm(), NonzeroTerm(), GreaterthannTerm(2),
                   CountAtleastnTerm(2), CountMutualTerm(),
                   TransitiveTiesTerm(), CyclicalTiesTerm(),
                   NodeOSumTerm(), NodeISumTerm(), NodeSumTerm(),
                   # R-parity terms: the (min, max, min) weights, the
                   # other valued `mutual` forms and the threshold terms
                   TransitiveWeightsTerm(), CyclicalWeightsTerm(),
                   CountMutualTerm(:nabsdiff), CountMutualTerm(:geometric),
                   CountMutualTerm(:product), CountMutualTerm(:threshold; threshold=2),
                   SmallerthanTerm(2), EqualToTerm(3), InIntervalTerm(1, 3),
                   InIntervalTerm(1, 3; open=(false, false)),
                   # ergm / ergm.count terms added in 0.2
                   SumTerm(pow=2), SumTerm(pow=0.5), AtmostTerm(2), CMPTerm()]

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
        result = ergm_count(net, [SumTerm(), NonzeroTerm(), CountMutualTerm()]; method=:mple)
        @test result isa CountERGMResult
        @test length(result.coefficients) == 3
        @test isfinite(result.loglik)
        @test all(isfinite, result.coefficients)
    end

    # ------------------------------------------------------------------
    # Allocation regressions on the count MPLE.
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
        # value touches the network, and of NetworkCore.jl's own mutations only
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

    @testset "StatsAPI surface" begin
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
        @test NetworkCore.check_statsapi(result; strict=true) !== nothing
        bin = ergm_count(net, [SumTerm(), NonzeroTerm()]; reference=BinomialReference(4))
        @test NetworkCore.check_statsapi(bin; strict=true) !== nothing
        boot = ergm_count(net, [SumTerm(), NonzeroTerm()]; se=:bootstrap, n_boot=10,
                          rng=Xoshiro(5))
        @test NetworkCore.check_statsapi(boot; strict=true) !== nothing
        @test all(NetworkCore.check_statsapi(result))
        # ... and the optional eleventh verb, `coefnames`: R's labels, in
        # `coef` order, the StatsAPI binding (one `coefnames` under
        # `using ERGM, ERGMCount, StatsAPI`), a fresh vector per call
        for f in (result, bin, boot)
            @test all(values(NetworkCore.check_statsapi(f;
                required=(NetworkCore.STATSAPI_VERBS..., :coefnames), strict=true)))
            @test coefnames(f) == coeftable(f).names == ["sum", "nonzero"]
        end
        @test ERGMCount.coefnames === ERGMCount.StatsAPI.coefnames === NetworkCore.coefnames
        cn = coefnames(result); cn[1] = "changed"
        @test coefnames(result)[1] == "sum"

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
        @test tbl isa NetworkCore.CoefficientTable
        @test tbl.names == [name(t, net) for t in result.model.terms] == ["sum", "nonzero"]
        @test tbl["sum"].estimate == coef(result)[1]
        @test tbl["nonzero"].std_error == stderror(result)[2]
        @test tbl[1].z_value == result.z_values[1]
        @test tbl[2].p_value == result.p_values[2]
        @test tbl.p_values == NetworkCore.z_pvalues(coef(result), stderror(result)).p
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
        @test fit_ergm_count === ergm_count
        # The never-released legacy alias is gone (see the docs' rename table)
        @test !isdefined(ERGMCount, :fit_count_ergm)
        net = random_count_net(5; seed=13)

        # A single term and a Tuple of terms are accepted and give the same fit
        fv = fit_ergm_count(net, [SumTerm(), NonzeroTerm()])
        ft = fit_ergm_count(net, (SumTerm(), NonzeroTerm()))
        @test coef(ft) == coef(fv)
        @test coef(fit_ergm_count(net, SumTerm())) == coef(fit_ergm_count(net, [SumTerm()]))
        @test fv.model.terms isa Tuple

        # An unknown method is refused naming the estimators (ERGM.jl's one
        # resolver); `:mcmle` (ergm.count's estimator) is implemented — see
        # "MCMLE" below
        for m in (:mcmc, :cd, :sa)
            err = try; ergm_count(net, [SumTerm()]; method=m); nothing
                  catch e; e end
            @test err isa ArgumentError
            msg = sprint(showerror, err)
            @test occursin("fit_ergm_count: unknown estimation method", msg)
            @test occursin(":auto", msg) && occursin(":mple", msg) && occursin(":mcmle", msg)
            @test occursin("dyad-independent formula", msg)
        end

        # method=:auto, R's rule: the exact MPLE for a dyad-independent formula,
        # the MCMLE otherwise
        @test fv.method === :mple && fv.mcmc === nothing && !fv.inference_withheld
        @test occursin("Method: mple (maximum pseudo-likelihood, which is the likelihood",
                       sprint(show, fv))
        dnet = count_net(8, [(1,2,2), (2,1,1), (2,3,1), (3,2,2), (3,1,3), (1,4,1),
                             (4,5,2), (5,4,1), (5,6,1), (6,1,2), (2,7,3), (7,8,1),
                             (8,2,2), (8,7,1)])
        dep = fit_ergm_count(dnet, [SumTerm(), CountMutualTerm()]; reference=BinomialReference(4),
                             n_samples=64, bridge_rungs=0, rng=Xoshiro(1), warn=false)
        @test dep.method === :mcmle && dep.mcmc !== nothing
        @test occursin("Method: mcmle (Monte-Carlo maximum likelihood", sprint(show, dep))
        dmp = fit_ergm_count(dnet, [SumTerm(), CountMutualTerm()]; reference=BinomialReference(4),
                             method=:mple)
        @test dmp.method === :mple && dmp.inference_withheld
        @test occursin("the default method=:auto fits the MCMLE here", sprint(show, dmp))
        # A keyword of the other estimator is refused in words, naming the one
        # that takes it — an MPLE keyword on a dyad-dependent formula ...
        err = try
            fit_ergm_count(dnet, [SumTerm(), CountMutualTerm()]; se=:bootstrap, n_boot=10); nothing
        catch e; e end
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("keyword `se`, `n_boot` is not accepted by method=:mcmle", msg)
        @test occursin("method=:auto chose :mcmle because the formula is dyad-dependent", msg)
        @test occursin("pass method=:mple explicitly", msg)
        # ... and an MCMLE keyword on a dyad-independent one
        err = try; fit_ergm_count(net, [SumTerm()]; n_samples=100); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("pass method=:mcmle explicitly", sprint(showerror, err))
        # An explicit method names no `:auto` reason, and a keyword neither
        # estimator takes points at both docstrings
        err = try; fit_ergm_count(net, [SumTerm()]; method=:mcmle, se=:hessian); nothing
              catch e; e end
        @test occursin("pass method=:mple explicitly", sprint(showerror, err)) &&
              !occursin("method=:auto chose", sprint(showerror, err))
        err = try; fit_ergm_count(net, [SumTerm()]; bogus=1); nothing catch e; e end
        @test err isa ArgumentError && occursin("`?count_mcmle`", sprint(showerror, err))

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
        @test fit_ergm_count(unet, [SumTerm(), NodeSumTerm()]; method=:mple) isa CountERGMResult
        @test fit_ergm_count(unet, [SumTerm(), TransitiveTiesTerm()]; method=:mple) isa CountERGMResult
    end

    @testset "Shared contracts are imported, not re-implemented" begin
        # ONE z → p helper, ONE Newton kernel, ONE se validator, ONE dependence
        # predicate, all NetworkCore's or ERGM's
        @test ERGMCount.z_pvalues === NetworkCore.z_pvalues
        @test ERGMCount.newton_fit === NetworkCore.newton_fit
        @test ERGMCount.check_se === NetworkCore.check_se
        @test ERGMCount.has_dyad_dependent === ERGM.has_dyad_dependent
        @test ERGMCount.CoefficientTable === NetworkCore.CoefficientTable
        @test ERGMCount.coeftable === NetworkCore.coeftable
        @test !isdefined(ERGMCount, :_z_pvalues)
        @test !isdefined(ERGMCount, :_has_dyad_dependent)
        @test hasmethod(ERGM.has_dyad_dependent, Tuple{CountERGMModel})
        @test hasmethod(ERGM.requires_directed, Tuple{NodeOSumTerm})
        @test ERGMCount.count_mple isa Function
        # No StatsBase: the Gibbs draw is an inline inverse-CDF loop
        @test !isdefined(ERGMCount, :StatsBase)
        @test !isdefined(ERGMCount, :Weights)

        # No cross-package `Pkg._name` reach-in: ERGM.jl's building blocks
        # come from its extension API, `ERGM.Extension` (semver-covered), and
        # an underscore name of ERGM or NetworkCore is never called
        src = _readtext(joinpath(dirname(@__DIR__), "src", "ERGMCount.jl"))
        code = join(filter(l -> !startswith(strip(l), "#"), split(src, '\n')), '\n')
        @test isempty(collect(eachmatch(r"\b(?:ERGM|NetworkCore)\._\w+", code)))
        ext = sort(unique(String(m.captures[1])
                          for m in eachmatch(r"\bERGM\.Extension\.(\w+)", code)))
        @test ext == ["bridge_integrate", "confidence_test", "ess_sample", "mcmc_defaults",
                      "mcmle_covariance", "mcmle_solve", "n_observed_dyads"]
        @test all(n -> Base.isexported(ERGM.Extension, Symbol(n)), ext)
        # The dyad count is ERGM's one definition, not a private copy
        @test !isdefined(ERGMCount, :_n_dyads)
        # ... so the MCMLE iteration is ERGM's driver, not a local loop
        @test !occursin("_hummel_step", code) && !occursin("cholesky", code)
        @test Base.ispublic(ERGM, :mcmc_convergence) && Base.ispublic(ERGM, :MCMLEConvergence)
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
        # NetworkCore exports a `count_mple` or `dyad_value`
        @test !Base.isexported(ERGM, :count_mple) && !Base.isexported(NetworkCore, :dyad_value)

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
        @test fit.p_values == NetworkCore.z_pvalues(fit.z_values)
        out = sprint(show, fit)
        @test occursin("<1e-16", out)
        @test !occursin("0.0000 ***", out)
    end

    @testset "gof extends the shared NetworkCore.gof generic" begin
        # One generic across the ecosystem: the method is added to
        # NetworkCore.gof, not a package-local function
        @test ERGMCount.gof === NetworkCore.gof

        net = random_count_net(6; seed=19)
        result = fit_ergm_count(net, [SumTerm(), NonzeroTerm()])
        g = ERGMCount.gof(result; n_sim=8, burnin=10, interval=2,
                          rng=Random.Xoshiro(31))
        @test g isa NetworkCore.GOFResult
        @test NetworkCore.n_simulations(g) == 8
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
        @test NetworkCore.missing_policies(fit_ergm_count) == (:error,)
        @test NetworkCore.missing_policies(ergm_count) == (:error,)
        @test NetworkCore.missing_policies(simulate_count_ergm) == (:error,)
        @test NetworkCore.missing_policies(ERGMCount.count_mple) == (:error,)
        @test NetworkCore.supports_missing(fit_ergm_count) == false
        @test NetworkCore.supports_missing(simulate_count_ergm) == false

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
        # observed zero (coef [0.179, -2.025], nobs 30 on a 6-node network)
        # while `fit_ergm_count` on the same network refused
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
        dep = fit_ergm_count(net, [SumTerm(), NodeOSumTerm()]; method=:mple,
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
        # shared `NetworkCore.bootstrap_cov` loop — same API as `ERGM.mple`.
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

        hess = fit_ergm_count(net, terms; method=:mple, reference=ref, se=:hessian)
        boot = fit_ergm_count(net, terms; method=:mple, reference=ref, se=:bootstrap,
                              n_boot=100, rng=MersenneTwister(7))

        # The bootstrap replaces the COVARIANCE, not the point estimate
        @test coef(boot) == coef(hess)
        @test loglikelihood(boot) == loglikelihood(hess)
        @test stderror(boot) != stderror(hess)
        @test vcov(boot) != vcov(hess)
        @test all(isfinite, stderror(boot))

        # Reproducible under a fixed rng
        boot2 = fit_ergm_count(net, terms; method=:mple, reference=ref, se=:bootstrap,
                               n_boot=100, rng=MersenneTwister(7))
        @test stderror(boot2) == stderror(boot)
        @test vcov(boot2) == vcov(boot)
        @test stderror(fit_ergm_count(net, terms; method=:mple, reference=ref, se=:bootstrap,
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
        # ... but the POINT ESTIMATE is a pseudo-likelihood estimate whatever
        # the standard errors are, and `show` keeps saying so under the
        # bootstrap (it used to print no caveat at all there)
        @test occursin("still pseudo-likelihood", out_b)
        @test occursin("not the MLE that R's ergm.count reports", out_b)
        @test occursin("method=:mcmle", out_b) && occursin("method=:mcmle", out_h)

        # The approximations list agrees with the printed prose (one predicate)
        @test any(occursin("anticonservative", a) for a in approximations(hess))
        @test !any(occursin("anticonservative", a) for a in approximations(boot))
        @test any(occursin("parametric bootstrap", a) for a in approximations(boot))
        # The POINT ESTIMATE is a pseudo-likelihood estimate either way, and both
        # fits still say so
        for f in (hess, boot)
            @test any(occursin("not the MLE that R's ergm.count (MCMLE) reports", a)
                      for a in approximations(f))
        end

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
        @test_throws ArgumentError fit_ergm_count(net, terms; method=:mple, reference=ref,
                                                  se=:sandwich)
        @test_throws ArgumentError fit_ergm_count(net, terms; method=:mple, reference=ref,
                                                  se=:bootstrap, n_boot=1)
        # The bootstrap's Gibbs controls default through the shared rule too
        @test_throws ArgumentError fit_ergm_count(net, terms; method=:mple, reference=ref,
                                                  se=:bootstrap, n_boot=4,
                                                  boot_interval=0)
    end

    @testset "Bootstrap is thread-count independent (fresh process)" begin
        # The replicates are drawn serially from the caller's `rng`; only the
        # refits run under `Threads.@threads` (in `NetworkCore.bootstrap_cov`), and
        # each refit is deterministic — so the standard errors are bit-identical
        # whatever the thread count. Prove it in a fresh process with a
        # DIFFERENT thread count, as ERGM.jl does for `mcmle`.
        net = count_net(8, [(1,2,2), (2,1,1), (2,3,1), (3,2,2), (3,1,3), (1,4,1),
                            (4,5,2), (5,4,1), (5,6,1), (6,1,2), (2,7,3), (7,8,1),
                            (8,2,2), (8,7,1)])
        boot = fit_ergm_count(net, [SumTerm(), CountMutualTerm()];
                              reference=BinomialReference(3), method=:mple, se=:bootstrap,
                              n_boot=20, rng=Xoshiro(7))
        # ... and so is the MCMLE: its chain is serial, and every path-sampling
        # grid point runs on its own generator seeded from `rng` up front
        mle = fit_ergm_count(net, [SumTerm(), CountMutualTerm()];
                             reference=BinomialReference(3), method=:mcmle,
                             n_samples=256, rng=Xoshiro(7))
        other_threads = Threads.nthreads() == 1 ? 4 : 1
        script = """
            using ERGMCount, NetworkCore, Graphs, Random
            net = network(8; directed=true)
            for (i, j, w) in [(1,2,2), (2,1,1), (2,3,1), (3,2,2), (3,1,3), (1,4,1),
                              (4,5,2), (5,4,1), (5,6,1), (6,1,2), (2,7,3), (7,8,1),
                              (8,2,2), (8,7,1)]
                add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
            end
            f = fit_ergm_count(net, [SumTerm(), CountMutualTerm()];
                               reference=BinomialReference(3), method=:mple, se=:bootstrap,
                               n_boot=20, rng=Xoshiro(7))
            println(Threads.nthreads()); println(repr(coef(f))); println(repr(stderror(f)))
            m = fit_ergm_count(net, [SumTerm(), CountMutualTerm()];
                               reference=BinomialReference(3), method=:mcmle,
                               n_samples=256, rng=Xoshiro(7))
            println(repr(coef(m))); println(repr(loglikelihood(m)))
            """
        cmd = `$(Base.julia_cmd()) --startup-file=no --threads=$other_threads --project=$(dirname(@__DIR__)) -e $script`
        lines = split(strip(read(pipeline(cmd; stderr=devnull), String)), '\n')
        @test length(lines) == 5
        @test parse(Int, lines[1]) == other_threads
        @test lines[2] == repr(coef(boot))
        @test lines[3] == repr(stderror(boot))
        @test lines[4] == repr(coef(mle))
        @test lines[5] == repr(loglikelihood(mle)) && isfinite(loglikelihood(mle))
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
        # What the exclusion does to the standard errors is said in the same
        # words as every bootstrap of the ERGM family, in all three places
        bias = "The standard errors are conditional on a finite refit: the excluded " *
               "replicates are the extreme ones, so the standard errors are biased downward."
        @test ERGMCount._BOOT_EXCLUSION_BIAS == bias
        rec = Test.collect_test_logs() do
            fit_ergm_count(net, terms; reference=BinomialReference(4),
                           se=:bootstrap, n_boot=40, rng=Xoshiro(1))
        end
        @test any(occursin(bias, r.message) for r in rec[1] if r.level == Base.CoreLogging.Warn)
        @test any(occursin(bias, a) for a in approximations(boot))
        @test occursin(bias, replace(sprint(show, boot), r"\s+" => " "))
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
            fit_ergm_count(net, [SumTerm(), CountMutualTerm()]; method=:mple, reference=BinomialReference(3))
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
        # (the model holds a dyad-dependent term, so the default withholds the
        # naive intervals; `se=:hessian` is the written opt-in)
        @test_throws ArgumentError confint(fit)
        @test isnan(fit.z_values[1]) && fit.inference_withheld
        naive = fit_ergm_count(net, [SumTerm(), CountMutualTerm()]; method=:mple,
                               reference=BinomialReference(3), se=:hessian, warn=false)
        @test confint(naive)[2, :] == [-Inf, -Inf]
        @test isfinite(naive.z_values[1]) && !naive.inference_withheld
        # A limit is not a maximizer: the fit is not exact, whatever the terms
        @test !is_exact(fit)
        # `warn=false` silences R's sentence (the bootstrap refits use it); the
        # result is identical and still records the boundary
        quiet = @test_logs min_level = Base.CoreLogging.Warn begin
            fit_ergm_count(net, [SumTerm(), CountMutualTerm()]; method=:mple,
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
            fit_ergm_count(net, [SumTerm(), CountMutualTerm()]; method=:mple, reference=BinomialReference(3),
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
            fit_ergm_count(net, terms; method=:mple, maxiter=1, max_val=10)
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
            fit_ergm_count(net, terms; method=:mple, maxiter=1, max_val=10)
        end
        warns = [r.message for r in rec[1] if r.level == Base.CoreLogging.Warn]
        @test any(occursin("maxiter = 1", w) && occursin("tol = 1.0e-8", w) &&
                  occursin("pseudo-score norm", w) && occursin("boundary or non-identified", w)
                  for w in warns)
        # `warn=false` silences the warning; the record is unchanged
        quiet = @test_logs min_level = Base.CoreLogging.Warn begin
            fit_ergm_count(net, terms; method=:mple, maxiter=1, max_val=10, warn=false)
        end
        @test !quiet.converged && coef(quiet) == coef(unc)
        # The budget is the caller's: running out of it is not rescued by the
        # coordinate-ascent start (that exists for a Newton that cannot MOVE)
        conv = fit_ergm_count(net, terms; method=:mple, max_val=10)
        @test conv.converged && conv.iterations > 1
        @test conv.gradient_norm < 1e-4
        @test loglikelihood(conv) > loglikelihood(unc)
        @test is_exact(fit_ergm_count(net, [SumTerm()]; reference=BinomialReference(4)))
    end

    @testset "Error-controlled count support" begin
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
        # ERGM.jl (`ERGM.Extension.mcmc_defaults`), converted from toggles to sweeps
        for nd in (30, 100, 870, 5000)
            d = ERGM.Extension.mcmc_defaults(nd)
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
        # The explicit-specification form resolves burnin/interval through
        # the same rule (its support is adaptive; at a mean of 0.6 no
        # conditional comes near the starting bound, so the draw equals the
        # one on a fixed 0:20)
        e_default = simulate_count_ergm(net34, [SumTerm()], [-0.5]; n_sim=1, rng=Xoshiro(8))[1]
        e_spelled = simulate_count_ergm(net34, [SumTerm()], [-0.5]; n_sim=1, burnin=20,
                                        interval=1, max_val=20, rng=Xoshiro(8))[1]
        @test weights_of(e_default) == weights_of(e_spelled)
    end

    @testset "Two-mode (bipartite) networks are refused everywhere" begin
        # The estimator enumerates and the sampler resamples EVERY off-diagonal
        # dyad: on a two-mode network the within-mode dyads are structurally
        # impossible, and they used to be fitted as observed zeros (nobs = 15 on
        # `network(6; bipartite=3)`, wrong pseudo-likelihood, BIC and
        # conditionals) while the README claimed there was no entry point.
        # Now every entry point refuses, as ERGM.jl's `ERGM.Extension.require_supported_network` does.
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
        # The single-column boundary test cannot see a combination; Newton
        # alone "converges" in one step to θ ≈ (−21, 21) with SE 2e4. The
        # shared verdict (NetworkCore.clogit_separation on the compressed
        # design) decides it from the data.
        s = count_net(6, [(1,2,1), (2,3,1), (3,4,1), (4,5,1), (5,6,1), (1,3,1), (2,5,1)];
                      directed=false)
        @test compute(SumTerm(), s) == compute(NonzeroTerm(), s) == 7
        fs = @test_logs (:warn, r"count_mple: the MPLE does not exist \(separation\)") match_mode = :any begin
            fit_ergm_count(s, [SumTerm(), NonzeroTerm()])
        end
        @test !fs.converged
        @test fs.separated
        @test fs.separated_terms == ["sum", "nonzero"]
        @test !is_exact(fs)
        @test coef(fs)[1] < -15 && coef(fs)[2] > 15          # the asymptote Newton stopped on
        # The policy: inference withheld on every coefficient, NaN intervals
        @test all(isnan, fs.z_values) && all(isnan, fs.p_values)
        @test all(isnan, confint(fs))
        @test any(occursin("MPLE does not exist", a) && occursin("`sum`, `nonzero`", a)
                  for a in approximations(fs))
        out = sprint(show, fs)
        @test occursin("Converged: false", out)
        @test occursin("the MPLE does not exist", out) && occursin("`sum`, `nonzero`", out)
        @test occursin("R ergm warns", out)
        # The separation caveat subsumes the conditioning one (one message
        # about meaningless standard errors, not two)
        @test !occursin("numerically singular", out)
        @test !any(occursin("numerically singular", a) for a in approximations(fs))
        rec = Test.collect_test_logs() do
            fit_ergm_count(s, [SumTerm(), NonzeroTerm()])
        end
        warns = [r.message for r in rec[1] if r.level == Base.CoreLogging.Warn]
        @test any(occursin("The MPLE does not exist!", w) && occursin("`sum`, `nonzero`", w)
                  for w in warns)
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
        @test fm.separated_terms == ["sum", "atleast.2"]
        @test coef(fm)[2] > 2 * abs(coef(fm)[1]) * 0.9      # the (−1, 2) direction
        # The bootstrap is refused outright: there is no fitted model to
        # simulate replicates from
        err = try
            fit_ergm_count(s, [SumTerm(), NonzeroTerm()]; se=:bootstrap, n_boot=4,
                           warn=false); nothing
        catch e; e end
        @test err isa ArgumentError && occursin("separation on `sum`, `nonzero`", err.msg)
        # ... and so is the MCMLE, whose start is the MPLE
        dep_s = count_net(4, [(1,2,1), (2,1,1), (2,3,1), (3,2,1), (3,4,1), (4,1,1)])
        fsd = fit_ergm_count(dep_s, [SumTerm(), NonzeroTerm(), CountMutualTerm()];
                             method=:mple, warn=false)
        @test fsd.separated
        err = try
            fit_ergm_count(dep_s, [SumTerm(), NonzeroTerm(), CountMutualTerm()]; warn=false); nothing
        catch e; e end
        @test err isa ArgumentError && occursin("count_mcmle", err.msg) &&
              occursin("separation", err.msg)
        # Well-posed fits are NOT flagged: the fixture model on zach, a wide
        # truncated support whose tilt alone spans > 18 nats (the case that
        # made a Newton-based heuristic need a second signature), a
        # dyad-dependent model, and the bounded references
        g = load_golden(joinpath(@__DIR__, "fixtures", "zach_poisson.toml"))
        zach = network(Int(g.values["n_actors"]); directed=false)
        for (a, b, w) in zip(Int.(g.values["edge_src"]), Int.(g.values["edge_dst"]),
                             Int.(g.values["edge_weight"]))
            add_edge!(zach, a, b); set_edge_attribute!(zach, :weight, a, b, w)
        end
        for kw in ((;), (; max_val=60), (; reference=GeometricReference()),
                   (; reference=BinomialReference(8)))
            f = fit_ergm_count(zach, [SumTerm(), NonzeroTerm()]; kw...)
            @test f.converged && !f.separated && isempty(f.separated_terms)
            @test f.hessian_cond < 1e4
            @test isempty(f.collinear)
        end
        fd = fit_ergm_count(zach, [SumTerm(), NonzeroTerm(), TransitiveWeightsTerm()];
                            max_val=14, method=:mple)
        @test fd.converged && !fd.separated && fd.hessian_cond < 1e4
        # The verdict itself, on the compressed design: it reads the data, not
        # a Newton iterate, so it holds wherever the optimizer stopped
        D = ERGMCount._count_design(fs.model, 0:fs.max_val)
        mask = fill(true, length(D.support), length(D.n_tot))
        v = ERGMCount._count_separation_verdict(D, [1, 2], mask)
        @test v.separated && v.certified && v.terms == [1, 2]
        @test v.direction[2] ≈ -v.direction[1]               # the sum − nonzero direction
        Dz = ERGMCount._count_design(fd.model, 0:14)
        mz = fill(true, 15, length(Dz.n_tot))
        @test !ERGMCount._count_separation_verdict(Dz, [1, 2, 3], mz).separated
        # Where Newton stopped does not matter: cut off after two iterations
        # (far from the asymptote, `converged == false` for the budget alone)
        # the design is still declared separated. A detector that read the
        # Newton iterate needed a converged fit first and called this one
        # merely unconverged.
        cut = fit_ergm_count(s, [SumTerm(), NonzeroTerm()]; maxiter=2, warn=false)
        @test cut.separated && cut.separated_terms == ["sum", "nonzero"]
        @test abs(coef(cut)[1]) < 15                         # nowhere near the asymptote
        # Adding ONE dyad observed at 2 breaks the 0/1 pattern: no longer
        # separated (the observed 2 scores below 0 and 1 along sum − nonzero)
        s2 = copy(s); add_edge!(s2, 4, 6); set_edge_attribute!(s2, :weight, 4, 6, 2)
        f2 = fit_ergm_count(s2, [SumTerm(), NonzeroTerm()])
        @test f2.converged && !f2.separated
        # A separated design with a column fixed by the boundary drop: the
        # verdict runs on the RESTRICTED supports (after `atmost(1)` at its
        # largest value restricts every dyad to 0:1), where `sum − nonzero`
        # is no longer a recession direction (it is 0 on both values) — the
        # verdict sees the restricted design, not the full one
        fb = fit_ergm_count(s, [SumTerm(), AtmostTerm(1)]; warn=false)
        @test isinf(coef(fb)[2]) && !fb.separated
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
        # A 0/1 network fit as counts (`:weight` = 1 on every edge) — `sum ≡
        # nonzero` on the data, which used to be printed as
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
            fit_ergm_count(bn, [SumTerm(), NonzeroTerm(), CountMutualTerm()]; method=:mple)
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
        # used to count that edge as 1 silently (sum 12 → 13 on a 6-node
        # network, coef [0.179, -2.025] → [0.063, -1.653]): the same
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
        # returns NaN; an earlier version returned 3.0/7.0 (the `best = 0` floor of the
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
        fneg = fit_ergm_count(neg, [SumTerm(), NonzeroTerm(), CountMutualTerm(:product)]; method=:mple, reference=ref)
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
        @test fit_ergm_count(pos, [SumTerm(), CountMutualTerm(:geometric)]; method=:mple,
                             reference=DiscUnif2Reference(0, 4)).converged
        @test all(ERGMCount._admits_negative(t) for t in
                  (SumTerm(), CountMutualTerm(:min), CountMutualTerm(:nabsdiff),
                   CountMutualTerm(:product), TransitiveTiesTerm(), NodeSumTerm()))
    end

    @testset "Adaptive support stops on a rung whose Newton failed" begin
        # `sum + greaterthan(2) + atleast(3)` is exactly collinear on integer
        # counts. The doubling used to run all 8 rungs on zach (max_val 14 → 3584, the
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
        for fit in (fit_ergm_count(net, [SumTerm(), NonzeroTerm(), CountMutualTerm()]; method=:mple),
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
        # Every export has a docstring with a runnable example.
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
        src = _readtext(joinpath(dirname(@__DIR__), "src", "ERGMCount.jl"))
        # Git may check out source files with CRLF line endings on Windows.
        blocks = [String(m.captures[1]) for m in eachmatch(r"```julia\r?\n(.*?)```"s, src)]
        @test length(blocks) >= 40
        for (k, code) in enumerate(blocks)
            m = Module(Symbol("DocBlock", k))
            ok = try
                Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
                    Core.eval(m, :(using NetworkCore, ERGM, ERGMCount, Random))
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
        # The user-visible docstrings state the stopping conjunction (delta AND
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
              hasfield(typeof(fit), :collinear) && hasfield(typeof(fit), :separated_terms) &&
              hasfield(typeof(fit), :improper)
        @test !fit.separated && fit.hessian_cond >= 1 && isempty(fit.collinear)
        @test isempty(fit.separated_terms) && !fit.improper
    end

    # ------------------------------------------------------------------
    # R-parity terms: the valued `mutual` forms,
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
        fit = fit_ergm_count(dn, [SumTerm(), CountMutualTerm(:nabsdiff), InIntervalTerm(1, 3)]; method=:mple,
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

        # --- the DEFAULT path, with no max_val ---------------------------------
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

    # ------------------------------------------------------------------
    # The sampler's support: an unbounded reference used to
    # be simulated on a literal 0:20 (explicit form) or on the fit's bound
    # with no diagnostic, and a model whose joint distribution is not
    # normalisable was certified `:converged` and then simulated,
    # bootstrapped and GOF-tested on a chain pinned to the bound.
    # ------------------------------------------------------------------
    @testset "Simulation support is adaptive, tracked, and never silently binding" begin
        draws(sims) = [ERGMCount.dyad_value(s, get_edge_attribute(s, :weight, Int), i, j)
                       for s in sims for i in 1:nv(s) for j in 1:nv(s) if i != j]
        seed10 = network(10; directed=true)

        # A Poisson mean of 25: the default support widens until no
        # conditional reaches its top (the old literal 0:20 gave mean 18.0
        # with 28% of the draws at 20), silently
        sims = @test_logs min_level = Base.CoreLogging.Warn simulate_count_ergm(
            seed10, [SumTerm()], [log(25.0)]; n_sim=50, rng=Xoshiro(1))
        v = draws(sims)
        @test abs(mean(v) - 25) < 5 * 5 / sqrt(length(v))      # 5 MC standard errors
        @test maximum(v) > 30
        # ... and so does a geometric law with mean 99.5 (old: 9.6)
        sims = simulate_count_ergm(seed10, [SumTerm()], [-0.01];
                                   reference=GeometricReference(), n_sim=50, rng=Xoshiro(1))
        v = draws(sims)
        @test abs(mean(v) - 99.5) < 5 * 100 / sqrt(length(v))
        # the draws agree with the reference's own sampler on the mean
        @test abs(mean(draws(simulate_count_ergm(seed10, [SumTerm()], [log(3.0)];
                                                 n_sim=50, rng=Xoshiro(2)))) - 3) < 0.15

        # A FIXED bound that bites is the caller's truncated family: the draws
        # are returned, with the share of draws on the bound in a warning
        trunc = @test_logs (:warn, r"fixed at max_val = 20.*TRUNCATED at 0:20"s) simulate_count_ergm(
            seed10, [SumTerm()], [log(25.0)]; n_sim=50, max_val=20, rng=Xoshiro(1))
        @test maximum(draws(trunc)) == 20 && mean(draws(trunc)) < 19
        # ... and a fixed bound that does not bite is silent
        @test_logs min_level = Base.CoreLogging.Warn simulate_count_ergm(
            seed10, [SumTerm()], [log(2.0)]; n_sim=5, max_val=40, rng=Xoshiro(1))

        # An improper model cannot be simulated on the adaptive support: the
        # chain reaches its cap and the call is refused in words
        err = try
            simulate_count_ergm(seed10, [SumTerm()], [0.1]; reference=GeometricReference(),
                                n_sim=2, rng=Xoshiro(1)); nothing
        catch e; e end
        @test err isa ArgumentError
        @test occursin("NOT normalisable", err.msg) && occursin("max_val=k", err.msg)
        # ... including one whose dyad conditionals are all proper Poisson laws
        # (pair weight e^{θ k²}/k!²: k² outgrows 2 log k!)
        err = try
            simulate_count_ergm(seed10, [SumTerm(), CountMutualTerm(:product)], [0.3, 0.3];
                                n_sim=5, burnin=200, rng=Xoshiro(3)); nothing
        catch e; e end
        @test err isa ArgumentError && occursin("Every dyad conditional can still be proper", err.msg)
        @test occursin("mutual.product", err.msg) && occursin("max_val=k", err.msg)
        @test_throws ArgumentError simulate_count_ergm(seed10, [SumTerm()], [0.0];
                                                       max_doublings=0)

        # Poisson + mutual.product with a positive coefficient. Every
        # conditional at the observed network is a proper Poisson, so the
        # doubling settles — and the fit used to be reported `:converged`,
        # simulate to a mean of 19.6 (bound 48) or 77.6 (bound 200) against
        # 2.05 observed, and bootstrap silently. The analytic rule flags it
        # (`improper`), and the probe finds the mode on the bound as well.
        rng = Xoshiro(5)
        n = 12
        net = network(n; directed=true)
        for i in 1:n, j in 1:n
            i == j && continue
            w = rand(rng, 0:4)
            w > 0 && (add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w))
        end
        for i in 1:2:(n - 1), (a, b) in ((i, i + 1), (i + 1, i))
            has_edge(net, a, b) || add_edge!(net, a, b)
            set_edge_attribute!(net, :weight, a, b, 6)
        end
        terms = [SumTerm(), CountMutualTerm(:product)]
        fit = @test_logs (:warn, r"NOT normalisable") fit_ergm_count(net, terms; method=:mple)
        @test coef(fit)[2] > 0
        @test fit.improper && fit.boundary_mode
        @test fit.support_control === :improper && !fit.support_stable
        @test fit.boundary_mass < 1e-20          # the conditional check alone sees nothing
        @test !is_exact(fit)
        @test any(occursin("NOT normalisable", a) for a in approximations(fit))
        @test any(occursin("NOT error-controlled", a) for a in approximations(fit))
        out = sprint(show, fit)
        @test occursin("WARNING: the fitted model is NOT normalisable", out)
        @test occursin("fitted model is NOT normalisable on the unbounded support;", out)
        # everything downstream refuses
        for f in (() -> simulate_count_ergm(fit; n_sim=2, rng=Xoshiro(1)),
                  () -> gof(fit; n_sim=5, rng=Xoshiro(1)),
                  () -> fit_ergm_count(net, terms; method=:mple, se=:bootstrap, n_boot=5,
                                       rng=Xoshiro(1), warn=false),
                  () -> fit_ergm_count(net, terms; method=:mcmle, rng=Xoshiro(1)))
            err = try; f(); nothing catch e; e end
            @test err isa ArgumentError
            @test occursin("NOT normalisable", err.msg)
        end
        # ... unless the truncated family is asked for in writing, and then
        # the draws carry the warning
        s48 = @test_logs (:warn, r"TRUNCATED at 0:48") simulate_count_ergm(
            fit; n_sim=3, burnin=300, max_val=48, rng=Xoshiro(3))
        @test maximum(draws(s48)) == 48
        # A caller-fixed bound on the same model keeps `:fixed`, records the
        # mode, and the bootstrap and gof refuse the chain pinned to the bound
        fixed = fit_ergm_count(net, terms; max_val=48, method=:mple, warn=false)
        @test fixed.support_control === :fixed && fixed.boundary_mode && fixed.improper
        for f in (() -> gof(fixed; n_sim=5, burnin=300, rng=Xoshiro(1)),
                  () -> fit_ergm_count(net, terms; max_val=48, se=:bootstrap, n_boot=5,
                                       boot_burnin=300, rng=Xoshiro(1), warn=false,
                                       method=:mple))
            err = try; f(); nothing catch e; e end
            @test err isa ArgumentError && occursin("this is refused", err.msg)
        end

        # Proper dyad-dependent models are not flagged and simulate silently
        ok = fit_ergm_count(net, [SumTerm(), CountMutualTerm(:min)]; method=:mple)
        @test !ok.boundary_mode && !ok.improper
        @test ok.support_control === :converged && ok.support_stable
        @test_logs min_level = Base.CoreLogging.Warn gof(ok; n_sim=10, rng=Xoshiro(1))
        neg = fit_ergm_count(net, [SumTerm(), NodeOSumTerm()]; method=:mple)  # squared out-strength
        @test coef(neg)[2] < 0 && !neg.boundary_mode && !neg.improper
        # ... and a positive squared-strength coefficient is: the probe runs
        # on coefficients, not on what the data happened to estimate
        model = CountERGMModel([SumTerm(), NodeOSumTerm()], net, PoissonReference())
        @test ERGMCount._boundary_mode_probe(model, [0.1, 0.01], 0:20)
        @test !ERGMCount._boundary_mode_probe(model, [0.5, -0.01], 0:20)

        # The tracked dyad update draws exactly what the plain one draws, and
        # allocates nothing while the support stands
        tt = (SumTerm(), NonzeroTerm(), CountMutualTerm())
        θ = [0.2, -0.5, 0.3]
        a = copy(net); wa = get_edge_attribute(a, :weight, Int)
        b = copy(net); wb = get_edge_attribute(b, :weight, Int)
        st = ERGMCount._ChainSupport(PoissonReference(), 30; adaptive=false, cap=30)
        support = 0:30
        log_h = [log_reference(PoissonReference(), y) for y in support]
        η = zeros(31); buf = zeros(31)
        ra, rb = Xoshiro(9), Xoshiro(9)
        for _ in 1:3
            ERGMCount._gibbs_sweep!(ra, a, wa, tt, θ, support, log_h, η, buf)
            ERGMCount._gibbs_sweep_tracked!(rb, b, wb, tt, θ, st)
        end
        @test get_edge_attribute(a, :weight) == get_edge_attribute(b, :weight)
        @test st.updates == 3 * n * (n - 1) && st.growths == 0
        tracked(rng, net, w, terms, θ, st) =
            @allocated ERGMCount._gibbs_update_tracked!(rng, net, w, terms, θ, 1, 2, st)
        worst = 0
        for _ in 1:200
            old = ERGMCount.dyad_value(b, wb, 1, 2)
            bytes = tracked(rb, b, wb, tt, θ, st)
            # an edge insertion may grow Graphs' adjacency vector; nothing else allocates
            (old == 0 && ERGMCount.dyad_value(b, wb, 1, 2) != 0) || (worst = max(worst, bytes))
        end
        @test worst == 0
    end

    @testset "Improper models are flagged analytically, whatever the data" begin
        # 25-actor directed counts (each dyad: a shared reciprocity indicator
        # with probability 0.05 plus two Bernoulli(0.15) draws), as in the
        # reproduction that found the hole: on seeds 2 and 10 the Poisson +
        # mutual.product MPLE has a SMALL positive product coefficient (0.061,
        # 0.104). The probe at the fitted bound 0:20 saw nothing, the escape
        # lies at 80-160, and the fits were certified `:converged`, then
        # simulated and bootstrapped silently.
        function countnet(rng, n, lam, mutual)
            g = network(n; directed=true)
            for i in 1:n, j in (i+1):n
                s = rand(rng) < mutual ? 1 : 0
                for (a, b) in ((i, j), (j, i))
                    w = s + Int(rand(rng) < lam) + Int(rand(rng) < lam)
                    w > 0 && (add_edge!(g, a, b); set_edge_attribute!(g, :weight, a, b, w))
                end
            end
            return g
        end
        terms = [SumTerm(), CountMutualTerm(:product)]
        for (seed, θp, escape) in ((2, 0.061, 160), (10, 0.104, 80))
            g = countnet(Xoshiro(seed), 25, 0.15, 0.05)
            fit = @test_logs (:warn, r"NOT normalisable on the unbounded support of PoissonReference, whatever the data") match_mode = :any fit_ergm_count(
                g, terms; method=:mple)
            @test coef(fit)[2] ≈ θp atol = 5e-4
            @test fit.max_val == 20
            @test fit.improper && fit.support_control === :improper && !fit.support_stable
            @test any(occursin("NOT normalisable", a) for a in approximations(fit))
            @test occursin("WARNING: the fitted model is NOT normalisable", sprint(show, fit))
            # The probe on its own: blind at the fitted bound; at the bound
            # where the sampler escaped the conditional modes stick to the
            # top but still weigh less than the data (θ_p·Y² has not yet
            # outgrown 2·log Y! there); at the bound where the adaptive
            # sampler gives up (2^8 × 20) the mode on the bound outweighs them
            @test !ERGMCount._boundary_mode_probe(fit.model, coef(fit), 0:20)
            icm = ERGMCount._icm_from_top(fit.model, coef(fit), 0:escape)
            @test icm.stuck && icm.logweight < icm.observed
            @test !ERGMCount._boundary_mode_probe(fit.model, coef(fit), 0:escape)
            @test ERGMCount._boundary_mode_probe(fit.model, coef(fit), 0:(2^8 * 20))
            @test fit.boundary_mode
            # Everything downstream refuses, naming the term
            for f in (() -> simulate_count_ergm(fit; n_sim=2, rng=Xoshiro(1)),
                      () -> gof(fit; n_sim=5, rng=Xoshiro(1)),
                      () -> fit_ergm_count(g, terms; method=:mple, se=:bootstrap,
                                           n_boot=5, rng=Xoshiro(1), warn=false),
                      () -> fit_ergm_count(g, terms; rng=Xoshiro(1), warn=false))
                err = try; f(); nothing catch e; e end
                @test err isa ArgumentError
                @test occursin("NOT normalisable", err.msg) && occursin("mutual.product", err.msg)
            end
            # The written opt-in to the truncated family still works
            @test length(simulate_count_ergm(fit; n_sim=1, burnin=2, max_val=20,
                                             rng=Xoshiro(1))) == 1
        end

        # The rule itself, on coefficients (no data enters it)
        imp(terms, θ; ref=PoissonReference(), n=10, directed=true) =
            ERGMCount._improper_direction(Tuple(terms), θ, ref, n, directed)
        @test imp(terms, [-1.0, 1e-6]).config === :mutual          # any positive product
        @test imp(terms, [-1.0, -1e-6]) === nothing
        @test imp(terms, [5.0, 0.0]) === nothing                   # linear: Poisson wins
        @test imp([SumTerm(pow=2)], [1e-4]).config === :dyad       # y² beats y log y
        @test imp([SumTerm(pow=1.5)], [1e-4]).order == (1.5, false)
        @test imp([SumTerm(pow=0.5)], [10.0]) === nothing
        @test imp([NodeOSumTerm()], [1e-4]) !== nothing
        @test imp([NodeISumTerm()], [1e-4]) !== nothing
        @test imp([NodeSumTerm()], [1e-4]; directed=false) !== nothing
        # One super-linear term can offset another: Σ out² ≤ Σ (out + in)²
        @test imp([NodeOSumTerm(), NodeSumTerm()], [0.1, -0.1]) === nothing
        @test imp([NodeOSumTerm(), NodeSumTerm()], [0.1, -0.01]) !== nothing
        # ... and a negative squared strength holds down a positive product
        # along the reciprocated pair (2·(2Y)² against Y²), but not when weaker
        @test imp([CountMutualTerm(:product), NodeSumTerm()], [1.0, -0.2]) === nothing
        @test imp([CountMutualTerm(:product), NodeSumTerm()], [1.0, -0.1]) !== nothing
        # CMP: (y!)^θ against the Poisson 1/y! — improper beyond θ = 1; the
        # geometric counting measure has no y! to beat
        @test imp([SumTerm(), CMPTerm()], [-5.0, 0.99]) === nothing
        @test imp([SumTerm(), CMPTerm()], [-5.0, 1.01]).order == (1.0, true)
        @test imp([SumTerm(), CMPTerm()], [-5.0, 0.01]; ref=GeometricReference()) !== nothing
        @test imp([SumTerm(), CMPTerm()], [-5.0, -0.5]; ref=GeometricReference()) === nothing
        # Bounded references, non-finite coefficients and terms of unknown
        # growth: the rule is silent
        @test imp(terms, [-1.0, 1.0]; ref=BinomialReference(5)) === nothing
        @test imp(terms, [-Inf, 1.0]) === nothing
        @test imp([_UnknownCountTerm(), CountMutualTerm(:product)], [0.0, 1.0]) === nothing
        # Every built-in term has a known growth, so a new term cannot be
        # added without deciding it
        for t in ALL_TERMS, d in (true, false), c in ERGMCount._configurations(6, d)
            @test ERGMCount._growth(t, c, 6, d) !== nothing
        end
        # The closed forms of the squared strengths, against `compute` on the
        # configurations themselves (counts Y; Y² coefficient = value / Y²)
        Y = 7
        for (d, configs) in ((true, (:dyad, :mutual, :outstar, :instar, :all)),
                             (false, (:dyad, :star, :all))), c in configs
            n = 5
            net = network(n; directed=d)
            pairs = c === :dyad ? [(1, 2)] : c === :mutual ? [(1, 2), (2, 1)] :
                    c in (:outstar, :star) ? [(1, j) for j in 2:n] :
                    c === :instar ? [(j, 1) for j in 2:n] :
                    [(i, j) for i in 1:n for j in 1:n if (d ? i != j : i < j)]
            for (i, j) in pairs
                add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, Y)
            end
            @test ERGMCount._n_config_dyads(c, n, d) == length(pairs)
            for t in (d ? (NodeOSumTerm(), NodeISumTerm(), NodeSumTerm(),
                           CountMutualTerm(:product), SumTerm(pow=2)) :
                          (NodeSumTerm(), SumTerm(pow=2)))
                @test ERGMCount._growth(t, c, n, d).c ≈ compute(t, net) / Y^2
            end
        end

        # A dyad-INDEPENDENT improper model: `sum(pow=2)` with a small positive
        # coefficient. Its conditional is proper on 0:20 up to a vanishing top
        # mass, so the doubling settles; only the rule sees it.
        # under-dispersed counts (0, 1 and 2 in equal shares: variance 2/3
        # against a mean of 1), so the fitted sum2 < 0
        net = count_net(6, [(i, j, (i + j) % 3) for i in 1:6 for j in (i + 1):6
                            if (i + j) % 3 > 0]; directed=false)
        sq = ERGMCount._improper_direction((SumTerm(), SumTerm(pow=2)), [-2.0, 0.01],
                                           PoissonReference(), 6, false)
        @test sq.config === :dyad && sq.terms == ["sum2"]
        err = try
            simulate_count_ergm(net, [SumTerm(), SumTerm(pow=2)], [-2.0, 0.01]; n_sim=1)
            nothing
        catch e; e end
        @test err isa ArgumentError && occursin("sum2", err.msg) && occursin("max_val=k", err.msg)
        # A proper model is untouched
        good = fit_ergm_count(net, [SumTerm(), SumTerm(pow=2)])
        @test coef(good)[2] < 0 && !good.improper && good.support_control === :converged
    end

    @testset "Geometric reference: linear growth decides propriety, the probe weighs its mode" begin
        # Under the geometric counting measure (h = 1) a model whose statistics
        # grow at most linearly is normalisable exactly when the log-weight
        # θ'g decays along every direction. The rule reads the exact linear
        # coefficient along each configuration: sum + mutual(:min) on a
        # reciprocated pair at (Y, Y) grows like (2θ_sum + θ_mutual)·Y
        imp(terms, θ; n=10, directed=true) =
            ERGMCount._improper_direction(Tuple(terms), θ, GeometricReference(), n, directed)
        sm = [SumTerm(), CountMutualTerm(:min)]
        d = imp(sm, [-1.0, 2.5])
        @test d.config === :mutual && d.order == (1.0, false) && d.coef ≈ 0.5 &&
              d.terms == ["mutual.min"]
        @test occursin("grows like 0.5·Y as Y", ERGMCount._improper_message(d, GeometricReference()))
        @test imp(sm, [-1.0, 1.5]) === nothing             # 2θs + θm = −0.5: proper
        @test imp([SumTerm()], [0.1]).config === :dyad     # a positive sum, one dyad
        @test imp([SumTerm()], [-0.1]) === nothing
        # nabsdiff (−|y_ij − y_ji|) is 0 on a reciprocated pair, so it cannot
        # hold down a positive sum there; with a negative sum, a negative
        # nabsdiff coefficient rewards asymmetry and grows along one dyad
        @test imp([SumTerm(), CountMutualTerm(:nabsdiff)], [0.2, 0.5]).config === :mutual
        @test imp([SumTerm(), CountMutualTerm(:nabsdiff)], [-0.2, 0.5]) === nothing
        @test imp([SumTerm(), CountMutualTerm(:nabsdiff)], [-0.2, -0.5]).config === :dyad
        # the two-path terms grow only when every dyad does: θ_sum·n(n−1) +
        # θ_tw·n(n−1) on 10 directed actors
        @test imp([SumTerm(), TransitiveWeightsTerm()], [-0.5, 0.4]) === nothing
        @test imp([SumTerm(), TransitiveWeightsTerm()], [-0.5, 0.6]).config === :all
        # bounded statistics add nothing; a term of unknown linear growth
        # leaves the configuration undecided
        @test imp([SumTerm(), NonzeroTerm(), GreaterthannTerm(2)], [-0.1, 5.0, 5.0]) === nothing
        @test imp([_UnknownCountTerm(), SumTerm()], [0.0, 0.1]) === nothing
        # Poisson is untouched: −Y log Y beats any linear growth
        @test ERGMCount._improper_direction((SumTerm(), CountMutualTerm(:min)), [-1.0, 2.5],
                                            PoissonReference(), 10, true) === nothing
        # The linear coefficients, against `compute` on the configurations
        # themselves: g(2Y) − g(Y) = a·Y for statistics homogeneous of degree
        # 1, and 0 for bounded ones once Y is past their thresholds
        Y = 50
        for (dd, configs) in ((true, (:dyad, :mutual, :outstar, :instar, :all)),
                              (false, (:dyad, :star, :all))), c in configs
            n = 5
            at(y) = begin
                net = network(n; directed=dd)
                pairs = c === :dyad ? [(1, 2)] : c === :mutual ? [(1, 2), (2, 1)] :
                        c in (:outstar, :star) ? [(1, j) for j in 2:n] :
                        c === :instar ? [(j, 1) for j in 2:n] :
                        [(i, j) for i in 1:n for j in 1:n if (dd ? i != j : i < j)]
                for (i, j) in pairs
                    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, y)
                end
                net
            end
            lo, hi = at(Y), at(2Y)
            for t in ALL_TERMS
                a = ERGMCount._linear_growth(t, c, n, dd)
                (a === nothing || (t isa SumTerm && t.pow < 1)) && continue
                (!dd && ERGM.requires_directed(t)) && continue
                @test (compute(t, hi) - compute(t, lo)) / Y ≈ a atol = 1e-9
            end
        end

        # A reciprocity-heavy construction: 20 actors, reciprocated pairs, the
        # geometric reference. Each MPLE sits in the regime θ_sum + θ_mutual >
        # 0 > 2θ_sum + θ_mutual, where the all-top network is a fixed point of
        # iterated conditional modes; the model is proper, and the fit is
        # neither refused nor flagged
        function recip_net(seed, n, lam, mu)
            rng = Xoshiro(seed)
            g = network(n; directed=true)
            for i in 1:n, j in (i + 1):n
                s = rand(rng) < mu ? 1 : 0
                for (a, b) in ((i, j), (j, i))
                    w = s + Int(rand(rng) < lam) + Int(rand(rng) < lam)
                    w > 0 && (add_edge!(g, a, b); set_edge_attribute!(g, :weight, a, b, w))
                end
            end
            return g
        end
        for (seed, lam, mu) in ((1, 0.05, 0.6), (2, 0.05, 0.8), (3, 0.02, 0.9))
            g = recip_net(seed, 20, lam, mu)
            fm = fit_ergm_count(g, sm; reference=GeometricReference(), method=:mple,
                                warn=false)
            θ = coef(fm)
            @test θ[1] + θ[2] > 0 > 2θ[1] + θ[2]
            @test !fm.boundary_mode && !fm.improper
            @test fm.support_control === :converged && fm.support_stable
            # the iterated conditional modes still end on the bound ...
            top = 2^8 * fm.max_val
            icm = ERGMCount._icm_from_top(fm.model, θ, 0:top)
            @test icm.stuck && icm.logweight < icm.observed
            net = copy(g); w = get_edge_attribute(net, :weight, Int)
            for i in 1:20, j in 1:20
                i == j || ERGMCount._set_dyad!(net, w, i, j, top)
            end
            lh = zeros(top + 1)
            @test ERGMCount._joint_logweight(fm.model.terms, θ, net, w, lh, 0:top,
                                             GeometricReference()) <
                  ERGMCount._joint_logweight(fm.model.terms, θ, g,
                                             get_edge_attribute(g, :weight, Int), lh, 0:top,
                                             GeometricReference())
            # ... but carry less weight than the data, so the probe says no
            @test !ERGMCount._boundary_mode_probe(fm.model, θ, 0:top)
            # simulation and gof run on the adaptive support
            @test length(simulate_count_ergm(fm; n_sim=2, rng=Xoshiro(1))) == 2
        end
        # An improper geometric model is still caught: by the rule, and by
        # the probe (its mode on the bound outweighs the data)
        g = recip_net(1, 20, 0.05, 0.6)
        m = CountERGMModel(sm, g, GeometricReference())
        @test ERGMCount._boundary_mode_probe(m, [-1.0, 2.5], 0:200)
        err = try
            simulate_count_ergm(g, sm, [-1.0, 2.5]; reference=GeometricReference(), n_sim=1)
            nothing
        catch e; e end
        @test err isa ArgumentError && occursin("NOT normalisable", err.msg) &&
              occursin("mutual.min", err.msg)
    end

    @testset "Count MPLE of a dyad-dependent model withholds naive inference unless asked" begin
        net = count_net(8, [(1,2,2), (2,1,1), (2,3,1), (3,2,2), (3,1,3), (1,4,1),
                            (4,5,2), (5,4,1), (5,6,1), (6,1,2), (2,7,3), (7,8,1),
                            (8,2,2), (8,7,1)])
        terms = [SumTerm(), CountMutualTerm()]
        ref = BinomialReference(3)
        dflt = fit_ergm_count(net, terms; method=:mple, reference=ref)
        @test dflt.method === :mple && dflt.inference_withheld
        @test all(isfinite, stderror(dflt)) && se_method(dflt) === :hessian
        @test all(isnan, dflt.z_values) && all(isnan, dflt.p_values)
        @test all(isnan, [r.p_value for r in (coeftable(dflt)[1], coeftable(dflt)[2])])
        out = sprint(show, dflt)
        @test occursin("z values and p-values are not reported (NaN)", out)
        @test occursin("under-cover", out) && occursin("se=:hessian requests", out)
        @test occursin("method=:mcmle", out)
        err = try; confint(dflt); nothing catch e; e end
        @test err isa ArgumentError && occursin("se=:bootstrap", err.msg) &&
              occursin("method=:mcmle", err.msg)
        # The default fit of a dyad-dependent formula is the MCMLE, so the
        # refusal names the fit it concerns: an MPLE with the default `se`
        @test occursin("an MPLE fit (method=:mple)", err.msg) &&
              occursin("with the default se", err.msg) && !occursin("default count MPLE", err.msg)
        @test any(occursin("withheld", a) for a in approximations(dflt))
        @test NetworkCore.check_statsapi(dflt) !== nothing     # the non-strict surface holds

        # Written opt-in: R-style naive Wald table, with the caveat
        naive = fit_ergm_count(net, terms; method=:mple, reference=ref, se=:hessian)
        @test !naive.inference_withheld && coef(naive) == coef(dflt) &&
              stderror(naive) == stderror(dflt)
        @test naive.z_values ≈ coef(naive) ./ stderror(naive)
        @test size(confint(naive)) == (2, 2)
        @test occursin("anticonservative", sprint(show, naive))
        @test !any(occursin("withheld", a) for a in approximations(naive))
        @test NetworkCore.check_statsapi(naive; strict=true) !== nothing

        # The bootstrap is calibrated inference: reported in full
        boot = fit_ergm_count(net, terms; method=:mple, reference=ref, se=:bootstrap, n_boot=40,
                              rng=Xoshiro(2))
        @test !boot.inference_withheld && all(isfinite, boot.z_values)
        @test size(confint(boot)) == (2, 2)

        # A dyad-independent model is untouched: its MPLE is the MLE
        indep = fit_ergm_count(net, [SumTerm(), NonzeroTerm()]; reference=ref)
        @test !indep.inference_withheld && all(isfinite, indep.z_values)
        @test size(confint(indep)) == (2, 2)
        @test !occursin("not reported", sprint(show, indep))

        # Why: the naive standard errors are too small. 60 networks simulated
        # at ergm.count's MLE of the zach model, the MPLE refit on each — the
        # empirical sd of the transitiveweights estimate is about twice the
        # mean naive standard error (0.121 vs 0.058 over 300 replicates, Wald
        # coverage 0.70: the figures the docs quote)
        g = load_golden(joinpath(@__DIR__, "fixtures", "count_mcmle.toml"))
        zach = count_net(Int(g.values["zach_n"]),
                         zip(Int.(g.values["zach_edge_src"]), Int.(g.values["zach_edge_dst"]),
                             Int.(g.values["zach_edge_weight"])); directed=false)
        zterms = [SumTerm(), NonzeroTerm(), TransitiveWeightsTerm()]
        θR = Float64.(g.values["zach_hp_coefficients"])
        sims = simulate_count_ergm(zach, zterms, θR; n_sim=60, burnin=50, interval=10,
                                   rng=Xoshiro(11))
        fits = [fit_ergm_count(s, zterms; method=:mple, max_val=20, se=:hessian, warn=false) for s in sims]
        est = [coef(f)[3] for f in fits]
        @test mean(stderror(f)[3] for f in fits) < 0.75 * std(est)
        @test abs(mean(est) - θR[3]) < 3 * std(est) / sqrt(60) + 0.03   # nearly unbiased
    end

    # ------------------------------------------------------------------
    # MCMLE — ergm.count's estimator, on the exact Gibbs sampler.
    # ------------------------------------------------------------------
    @testset "MCMLE: exact against enumeration on a 3-actor network" begin
        # 6 directed dyads, Binomial(2) reference: 3^6 = 729 states, so the
        # likelihood, its maximizer and the normalizing constant are sums
        ties = [(1,2,2), (2,1,1), (2,3,1), (3,1,1)]
        net = count_net(3, ties)
        terms = [SumTerm(), CountMutualTerm(:min)]
        ref = BinomialReference(2)
        dyads = [(i, j) for i in 1:3 for j in 1:3 if i != j]
        G = Matrix{Float64}(undef, 3^6, 2)
        logh = Vector{Float64}(undef, 3^6)
        for (r, ys) in enumerate(Iterators.product(ntuple(_ -> 0:2, 6)...))
            s = count_net(3, [(i, j, y) for ((i, j), y) in zip(dyads, ys) if y > 0])
            G[r, :] .= (compute(terms[1], s), compute(terms[2], s))
            logh[r] = sum(log_reference(ref, y) for y in ys)
        end
        g_obs = [compute(t, net) for t in terms]
        logh_obs = sum(log_reference(ref, y) for y in
                       (ERGMCount.dyad_value(net, get_edge_attribute(net, :weight, Int), i, j)
                        for (i, j) in dyads))
        function exact(θ)
            η = logh .+ G * θ
            m = maximum(η); w = exp.(η .- m); Z = sum(w); w ./= Z
            μ = G' * w
            Σ = G' * (w .* G) - μ * μ'
            return (dot(θ, g_obs) + logh_obs - (m + log(Z)), g_obs .- μ, -Σ)
        end
        mle = NetworkCore.newton_fit(exact, zeros(2))
        @test mle.converged

        fit = fit_ergm_count(net, terms; reference=ref, method=:mcmle, n_samples=4000,
                             bridge_samples=4000, rng=Xoshiro(3))
        @test fit.method === :mcmle && fit.converged
        mc = fit.mcmc
        @test mc.convergence isa ERGM.MCMLEConvergence
        @test all(abs.(coef(fit) .- mle.θ) .<= 5 .* mc.mc_std_errors)
        @test maximum(abs.(coef(fit) .- mle.θ) ./ mle.se) < 0.3
        # standard errors: inverse Fisher information, plus the MC component.
        # Measured over 30 seeds at these settings (n_samples=4000): the
        # relative error of each SE against the exact one has mean ≤ 0.002 and
        # seed-to-seed sd 0.010-0.011 (largest 0.023), so 0.05 is a 4.5-sd band
        @test stderror(fit) ≈ mle.se rtol = 0.05
        @test all(stderror(fit) .>= sqrt.(diag(inv(cov(mc.samples)))))
        # the log-likelihood, reference measure included
        @test abs(loglikelihood(fit) - mle.loglik) <= 5 * mc.loglik_mc_se + 0.02
        @test mc.loglik_mc_se < 0.1
        @test aic(fit) ≈ -2 * loglikelihood(fit) + 4
        # the MPLE it started from is a different number
        mple = fit_ergm_count(net, terms; method=:mple, reference=ref)
        @test mc.start == coef(mple)
        @test maximum(abs.(coef(mple) .- mle.θ)) > 10 * maximum(mc.mc_std_errors)

        # The record and the protocol
        @test objective(fit) === :likelihood && objective(mple) === :pseudolikelihood
        @test se_method(fit) === :fisher && fit.se_type === :mcmc
        @test !is_exact(fit) && !fit.inference_withheld
        @test all(isfinite, fit.z_values) && size(confint(fit)) == (2, 2)
        @test NetworkCore.check_statsapi(fit; strict=true) !== nothing
        @test all(values(NetworkCore.check_statsapi(fit;
            required=(NetworkCore.STATSAPI_VERBS..., :coefnames), strict=true)))
        @test coefnames(fit) == coeftable(fit).names == coefnames(mple)
        @test any(occursin("MCMLE: the likelihood is approximated", a) for a in approximations(fit))
        out = sprint(show, fit)
        @test occursin("Monte-Carlo maximum likelihood", out)
        @test occursin("Log-likelihood:", out) && !occursin("Pseudo-log-likelihood", out)
        @test occursin("inverse Fisher information + Monte-Carlo error", out)
        @test !occursin("pseudo-likelihood", out)       # no MPLE caveat on an MLE
        # reproducible from `rng` alone
        again = fit_ergm_count(net, terms; reference=ref, method=:mcmle, n_samples=4000,
                               bridge_samples=4000, rng=Xoshiro(3))
        @test coef(again) == coef(fit) && loglikelihood(again) == loglikelihood(fit)
        # gof and simulation run from an MCMLE fit
        gf = gof(fit; n_sim=200, rng=Xoshiro(4))
        @test all(p -> p > 0.05, gf.statistics[1].p_values)   # the MLE matches its statistics
        @test_throws ArgumentError fit_ergm_count(net, terms; reference=ref, method=:mcmle,
                                                  n_samples=8)
        err = try
            fit_ergm_count(net, terms; reference=ref, method=:mcmle, se=:bootstrap); nothing
        catch e; e end
        @test err isa ArgumentError && occursin("pass method=:mple explicitly", err.msg)
        err = try
            count_mcmle(CountERGMModel(terms, net, ref); se=:bootstrap); nothing
        catch e; e end
        @test err isa ArgumentError && occursin("option of the MPLE", err.msg)

        # A dyad-INDEPENDENT model needs no Monte Carlo: the exact MPLE
        di = fit_ergm_count(net, [SumTerm()]; reference=ref, method=:mcmle)
        dm = fit_ergm_count(net, [SumTerm()]; reference=ref)
        @test di.method === :mcmle && di.mcmc === nothing && is_exact(di)
        @test coef(di) == coef(dm) && loglikelihood(di) == loglikelihood(dm)
        @test occursin("exact: the model is dyad-independent", sprint(show, di))

        # A statistic at the boundary of its attainable range: R's drop. On a
        # 3-cycle no pair is reciprocated, so `mutual.min` is at its minimum
        # on every dyad's support; its coefficient is fixed at -Inf and `sum`
        # is the MLE of the model restricted to the networks without a
        # reciprocated pair — computed exactly over the 729 states
        cyc = count_net(3, [(1,2,2), (2,3,1), (3,1,1)])
        keep = G[:, 2] .== 0
        Gk, hk = G[keep, 1:1], logh[keep]
        g_cyc = compute(terms[1], cyc)
        logh_cyc = sum(log_reference(ref, y) for y in
                       (ERGMCount.dyad_value(cyc, get_edge_attribute(cyc, :weight, Int), i, j)
                        for (i, j) in dyads))
        function exact_restricted(θ)
            η = hk .+ Gk * θ
            m = maximum(η); w = exp.(η .- m); Z = sum(w); w ./= Z
            μ = Gk' * w
            Σ = Gk' * (w .* Gk) - μ * μ'
            return (θ[1] * g_cyc + logh_cyc - (m + log(Z)), [g_cyc] .- μ, -Σ)
        end
        rmle = NetworkCore.newton_fit(exact_restricted, zeros(1))
        @test rmle.converged
        dropped = @test_logs (:warn, r"count_mcmle: observed statistic\(s\) mutual.min are at their smallest.*fixed at -Inf.*drop=TRUE.*drop=false"s) fit_ergm_count(
            cyc, terms; reference=ref, method=:mcmle, n_samples=4000, rng=Xoshiro(3))
        @test dropped.converged && coef(dropped)[2] == -Inf
        @test abs(coef(dropped)[1] - rmle.θ[1]) <= 5 * dropped.mcmc.mc_std_errors[1]
        @test abs(coef(dropped)[1] - rmle.θ[1]) < 0.3 * rmle.se[1]
        @test stderror(dropped)[1] ≈ rmle.se[1] rtol = 0.05
        @test stderror(dropped)[2] == 0 && dropped.z_values[2] == -Inf &&
              dropped.p_values[2] == 0
        # the held statistic never moves in the sample, and the start is the
        # restricted MPLE
        @test all(==(0), dropped.mcmc.samples[:, 2])
        @test dropped.mcmc.start == coef(fit_ergm_count(cyc, terms; reference=ref,
                                                         method=:mple, warn=false))
        # one coefficient estimated; the log-likelihood is not, and says why
        @test dof(dropped) == 1 && isnan(loglikelihood(dropped))
        @test any(occursin("maximum-likelihood estimate exists (R ergm fixes it the same " *
                           "way under its default drop=TRUE)", a) for a in approximations(dropped))
        @test occursin("log-likelihood is not estimated", sprint(show, dropped))
        @test coefnames(dropped) == ["sum", "mutual.min"]
        # strict mode: R's control.ergm(drop=FALSE) is refused, by both
        # estimators, before any draw
        for m in (:mcmle, :mple)
            err = try
                fit_ergm_count(cyc, terms; reference=ref, method=m, drop=false); nothing
            catch e; e end
            @test err isa ArgumentError && occursin("drop=false", err.msg) &&
                  occursin("mutual.min", err.msg) && occursin("smallest attainable", err.msg)
        end
        # every statistic at its bound: nothing to estimate
        empty = count_net(3, Tuple{Int,Int,Int}[])
        err = try
            fit_ergm_count(empty, terms; reference=ref, method=:mcmle, warn=false); nothing
        catch e; e end
        @test err isa ArgumentError && occursin("nothing to estimate", err.msg)
        # A fit that cannot converge in its budget says so, quoting the
        # stopping rule it was held to — and only that rule: under
        # `:confidence` the equivalence test and the step length (the t-ratios
        # and the Hotelling test are diagnostics, not the rule) ...
        function nonconv_texts(termination; kw...)
            fit = nothing
            rec = Test.collect_test_logs() do
                fit = fit_ergm_count(net, terms; reference=ref, method=:mcmle,
                                     n_samples=64, maxiter=1, termination=termination,
                                     bridge_rungs=0, rng=Xoshiro(1), kw...)
            end
            warns = [r.message for r in rec[1] if r.level == Base.CoreLogging.Warn &&
                                                  occursin("MCMLE did not converge", r.message)]
            caveats = [a for a in approximations(fit) if occursin("MCMLE did not converge", a)]
            return fit, (only(warns), only(caveats), sprint(show, fit))
        end
        uc, texts = nonconv_texts(:confidence; conv_precision=1e-6)
        @test !uc.converged && isnan(loglikelihood(uc))
        @test uc.mcmc.conv_precision == 1e-6 && uc.mcmc.conv_confidence == 0.99
        for t in texts
            @test occursin("MCMLE did not converge in 1 iteration (", t)
            @test occursin("99% equivalence test p ", t) && occursin("needs < 0.01", t)
            @test occursin("tolerance precision 1e-06, $(uc.mcmc.n_samples) draws", t)
            @test occursin("step length γ", t)
            @test !occursin("t-ratio", t) && !occursin("Hotelling", t)
        end
        # ... and under `:hotelling` its p-value and the largest t-ratio
        uh, texts = nonconv_texts(:hotelling; conv_threshold=1e-6)
        @test !uh.converged
        for t in texts
            @test occursin("Hotelling p ", t) && occursin("max t-ratio", t)
            @test occursin("step length γ", t) && !occursin("equivalence test", t)
        end
        # The diagnostics describe a sample at the returned coefficients: on a
        # fit that did not converge, the fresh sample drawn at the last iterate
        @test size(uc.mcmc.samples, 1) == uc.mcmc.n_samples
    end

    @testset "Golden fixture: ergm.count MCMLE on zach and a directed network" begin
        g = load_golden(joinpath(@__DIR__, "fixtures", "count_mcmle.toml"))
        @test g.provenance["ergm_count_version"] == "4.1.3"
        band_sd = g.tolerance["band_sd"]
        floor_frac = g.tolerance["floor_se_fraction"]
        function check_fit(prefix, net, terms, seeds)
            v(key) = Float64.(g.values["$(prefix)_$key"])
            @test [name(t, net) for t in terms] == g.values["$(prefix)_term_names"]
            @test [compute(t, net) for t in terms] ≈ v("summary_statistics") atol = 1e-9
            n_seeds = length(g.values["$(prefix)_seeds"])
            Rm, Rsd, Sm, Ssd = v("coefficients_mean"), v("coefficients_sd"),
                               v("std_errors_mean"), v("std_errors_sd")
            band = max.(band_sd .* Rsd .* sqrt(1 + 1 / n_seeds), floor_frac .* Sm)
            se_band = max.(band_sd .* Ssd .* sqrt(1 + 1 / n_seeds), floor_frac .* Sm)
            # R reports the log-likelihood relative to the reference measure
            # (logLik at θ = 0 is 0): Σ log h(y_obs) − N·log Σ_y h(y), with
            # Σ_y 1/y! = e for the Poisson reference
            w = get_edge_attribute(net, :weight, Int)
            ys = [ERGMCount.dyad_value(net, w, i, j) for i in 1:nv(net) for j in 1:nv(net)
                  if (is_directed(net) ? i != j : i < j)]
            ll0 = -sum(ERGMCount._logfactorial, ys) - length(ys)
            ll_band = max(band_sd * sqrt(2) * g.values["$(prefix)_loglik_sd"],
                          g.tolerance["loglik_floor"])
            fits = [fit_ergm_count(net, terms; method=:mcmle, rng=Xoshiro(s)) for s in seeds]
            for f in fits
                @test f.converged && f.method === :mcmle
                @test all(abs.(coef(f) .- Rm) .<= band)
                @test all(abs.(stderror(f) .- Sm) .<= se_band)
                @test abs((loglikelihood(f) - ll0) - g.values["$(prefix)_loglik_mean"]) <= ll_band
                # the Monte-Carlo error is at R's precision (R's seed sd is
                # about 4% of an SE; measured here over six seeds per model:
                # 3.4-4.4%, set by the ESS target of n_samples/2 = 512)
                @test maximum(f.mcmc.mc_std_errors ./ stderror(f)) <= 0.055
                @test f.mcmc.termination === :confidence
                # ... and R's long-chain fit lies inside the same band
                @test all(abs.(coef(f) .- v("hp_coefficients")) .<= band)
            end
            return fits
        end
        zach = count_net(Int(g.values["zach_n"]),
                         zip(Int.(g.values["zach_edge_src"]), Int.(g.values["zach_edge_dst"]),
                             Int.(g.values["zach_edge_weight"])); directed=false)
        zterms = [SumTerm(), NonzeroTerm(), TransitiveWeightsTerm()]
        zfit = check_fit("zach", zach, zterms, (1,))[1]
        # The MPLE of the same model sits 1.7-1.8 R standard errors from the
        # MLE on `sum` and `transitiveweights` (the figure the docs quote)
        mple = fit_ergm_count(zach, zterms; method=:mple)
        gap = abs.(coef(mple) .- Float64.(g.values["zach_coefficients_mean"])) ./
              Float64.(g.values["zach_std_errors_mean"])
        @test 1.5 < gap[1] < 2.0 && 1.5 < gap[3] < 2.0 && gap[2] < 1.0
        @test maximum(abs.(coef(zfit) .- Float64.(g.values["zach_coefficients_mean"])) ./
                      Float64.(g.values["zach_std_errors_mean"])) < 0.2

        dn = count_net(Int(g.values["directed_n"]),
                       zip(Int.(g.values["directed_edge_src"]), Int.(g.values["directed_edge_dst"]),
                           Int.(g.values["directed_edge_weight"])))
        check_fit("directed", dn, [SumTerm(), NonzeroTerm(), CountMutualTerm(:min)], (1, 2, 3))
    end

    @testset "Network{Int32}: terms, both estimators, simulation and gof" begin
        # Vertex ids of any integer type: `edges(net)` and the neighbour lists
        # yield the network's own type, the sweeps `Int`; the weight
        # dictionary is keyed in the network's type. Every result equals the
        # `Network{Int}` one
        to32(net) = begin
            h = Network{Int32}(; n=nv(net), directed=is_directed(net))
            for e in edges(net)
                a, b = Int32(src(e)), Int32(dst(e))
                add_edge!(h, a, b)
                set_edge_attribute!(h, :weight, a, b, get_edge_attribute(net, :weight, src(e), dst(e)))
            end
            h
        end
        for directed in (true, false)
            d64 = random_count_net(9; directed=directed, p=0.4, seed=11)
            d32 = to32(d64)
            w64, w32 = get_edge_attribute(d64, :weight, Int), get_edge_attribute(d32, :weight, Int)
            @test keytype(w32) == Tuple{Int32,Int32}
            terms = [t for t in ALL_TERMS if directed || !ERGM.requires_directed(t)]
            @test [compute(t, d32) for t in terms] == [compute(t, d64) for t in terms]
            for t in terms, (i, j) in ((1, 2), (3, 1), (Int32(4), Int32(5)))
                @test change_stat_count(t, d32, w32, i, j, 1, 3) ==
                      change_stat_count(t, d64, w64, Int(i), Int(j), 1, 3)
                @test dyad_value(d32, w32, i, j) == dyad_value(d64, w64, i, j)
            end
        end
        d64 = random_count_net(10; p=0.35, seed=4)
        d32 = to32(d64)
        q(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())
        for terms in ([SumTerm(), NonzeroTerm()],
                      [SumTerm(), NonzeroTerm(), CountMutualTerm(), TransitiveWeightsTerm()])
            a = q(() -> fit_ergm_count(d64, terms; method=:mple))
            b = q(() -> fit_ergm_count(d32, terms; method=:mple))
            @test coef(b) == coef(a) && stderror(b) == stderror(a)
        end
        a = fit_ergm_count(d64, [SumTerm(), CountMutualTerm()]; rng=Xoshiro(5), n_samples=128,
                           bridge_rungs=2)
        b = fit_ergm_count(d32, [SumTerm(), CountMutualTerm()]; rng=Xoshiro(5), n_samples=128,
                           bridge_rungs=2)
        @test coef(b) == coef(a) && loglikelihood(b) == loglikelihood(a)
        s64 = simulate_count_ergm(a; n_sim=2, rng=Xoshiro(1))
        s32 = simulate_count_ergm(b; n_sim=2, rng=Xoshiro(1))
        @test [compute(SumTerm(), x) for x in s32] == [compute(SumTerm(), x) for x in s64]
        @test isequal(gof(b; n_sim=5, rng=Xoshiro(2)).p_overall, gof(a; n_sim=5, rng=Xoshiro(2)).p_overall)
    end

    @testset "Valued covariate terms: definitions, refusals, R's values and labels" begin
        # Brute force: each change statistic is the difference of `compute`
        # on edited copies, at every dyad and several value pairs, and the
        # support profile is the per-value definition
        attrnet(directed) = begin
            net = random_count_net(7; directed=directed, p=0.5, seed=8)
            set_vertex_attribute!(net, :g, ["a", "b", "a", "c", "b", "a", "c"])
            set_vertex_attribute!(net, :x, [1.5, 2.0, 0.5, 3.0, 1.0, 2.5, 0.0])
            net
        end
        W = [Float64((3i + 5j) % 7) / 4 for i in 1:7, j in 1:7]
        cov_terms(directed) = [CountNodeMatchTerm(:g), CountNodeMatchTerm(:g; level="a"),
                               CountNodeMatchTerm(:g; form=:nonzero),
                               CountNodeFactorTerm(:g, "b"),
                               CountNodeFactorTerm(:g, "c"; form=:nonzero),
                               CountAbsDiffTerm(:x), CountAbsDiffTerm(:x; form=:nonzero),
                               CountNodeCovTerm(:x), CountEdgeCovTerm(W; name="w"),
                               CountEdgeCovTerm(W; name="w", form=:nonzero),
                               (directed ? [CountNodeOCovTerm(:x), CountNodeICovTerm(:x)] : [])...]
        for directed in (true, false)
            net = attrnet(directed)
            w = get_edge_attribute(net, :weight, Int)
            for t in cov_terms(directed)
                @test !is_dyad_dependent(t)
                for i in 1:7, j in (directed ? (1:7) : (i + 1:7))
                    i == j && continue
                    for (old, new) in ((0, 3), (2, 0), (1, 4))
                        @test change_stat_count(t, net, w, i, j, old, new) ≈
                              brute_change_count(t, net, i, j, old, new) atol = 1e-12
                    end
                end
                dest = zeros(5)
                ERGMCount.change_stats_support!(dest, t, net, w, 1, 2, 1, 0:4)
                @test dest == [change_stat_count(t, net, w, 1, 2, 1, y) for y in 0:4]
            end
            # a model of them fits and simulates (dyad-independent: exact)
            fit = fit_ergm_count(net, [SumTerm(), CountNodeMatchTerm(:g), CountAbsDiffTerm(:x)])
            @test fit.converged && !has_dyad_dependent(fit.model) && fit.method === :mple
            @test coefnames(fit) == ["sum", "nodematch.sum.g", "absdiff.sum.x"]
            @test length(simulate_count_ergm(fit; n_sim=2, rng=Xoshiro(1))) == 2
        end
        # Refusals, in words
        net = attrnet(true)
        bare = random_count_net(7; seed=8)
        for t in (CountNodeMatchTerm(:g), CountAbsDiffTerm(:x), CountNodeCovTerm(:x))
            err = try; fit_ergm_count(bare, [SumTerm(), t]); nothing catch e; e end
            @test err isa ArgumentError && occursin("has no `:", err.msg)
        end
        set_vertex_attribute!(bare, :x, ["p", "q", "r", "s", "t", "u", "v"])
        err = try; fit_ergm_count(bare, [SumTerm(), CountAbsDiffTerm(:x)]); nothing catch e; e end
        @test err isa ArgumentError && occursin("finite number", err.msg)
        err = try; fit_ergm_count(net, [SumTerm(), CountEdgeCovTerm(ones(3, 3))]); nothing catch e; e end
        @test err isa ArgumentError && occursin("3×3", err.msg)
        @test_throws ArgumentError CountNodeMatchTerm(:g; form=:mean)
        @test_throws ArgumentError CountEdgeCovTerm(ones(2, 3))
        und = attrnet(false)
        err = try; fit_ergm_count(und, [SumTerm(), CountNodeOCovTerm(:x)]); nothing catch e; e end
        @test err isa ArgumentError && occursin("CountNodeCovTerm", err.msg)
        @test compute(CountNodeICovTerm(:x), und) == 0
        err = try; fit_ergm_count(net, [SumTerm(), ERGM.NodeMatch(:g)]); nothing catch e; e end
        @test err isa ArgumentError && occursin("CountNodeMatchTerm", err.msg)

        # R's values and labels (ergm 4.12), and the exact MLE of a model of
        # them against ergm.count's MCMLE
        g = load_golden(joinpath(@__DIR__, "fixtures", "count_covariates.toml"))
        n = Int(g.values["n"])
        rebuild(prefix; directed) = begin
            net = count_net(n, zip(Int.(g.values["$(prefix)_edge_src"]),
                                   Int.(g.values["$(prefix)_edge_dst"]),
                                   Int.(g.values["$(prefix)_edge_weight"])); directed=directed)
            set_vertex_attribute!(net, :g, String.(g.values["g"]))
            set_vertex_attribute!(net, :x, Float64.(g.values["x"]))
            net, permutedims(reshape(Float64.(g.values["$(prefix)_dist"]), n, n))
        end
        dn, Wd = rebuild("directed"; directed=true)
        dterms = [CountNodeMatchTerm(:g),
                  [CountNodeMatchTerm(:g; level=l) for l in ("a", "b", "c")]...,
                  CountNodeMatchTerm(:g; form=:nonzero),
                  CountNodeFactorTerm(:g, "b"), CountNodeFactorTerm(:g, "c"),
                  CountNodeFactorTerm(:g, "b"; form=:nonzero),
                  CountNodeFactorTerm(:g, "c"; form=:nonzero),
                  CountAbsDiffTerm(:x), CountAbsDiffTerm(:x; form=:nonzero),
                  CountNodeCovTerm(:x), CountNodeCovTerm(:x; form=:nonzero),
                  CountNodeOCovTerm(:x), CountNodeICovTerm(:x),
                  CountEdgeCovTerm(Wd; name="dist"),
                  CountEdgeCovTerm(Wd; name="dist", form=:nonzero)]
        @test [name(t, dn) for t in dterms] == g.values["directed_names"]
        @test [compute(t, dn) for t in dterms] ≈ Float64.(g.values["directed_summary"]) atol = 1e-9
        un, Wu = rebuild("undirected"; directed=false)
        uterms = [CountNodeMatchTerm(:g),
                  [CountNodeMatchTerm(:g; level=l) for l in ("a", "b", "c")]...,
                  CountNodeFactorTerm(:g, "b"), CountNodeFactorTerm(:g, "c"),
                  CountAbsDiffTerm(:x), CountNodeCovTerm(:x),
                  CountEdgeCovTerm(Wu; name="dist"),
                  CountEdgeCovTerm(Wu; name="dist", form=:nonzero)]
        @test [name(t, un) for t in uterms] == g.values["undirected_names"]
        @test [compute(t, un) for t in uterms] ≈ Float64.(g.values["undirected_summary"]) atol = 1e-9
        fterms = [SumTerm(), NonzeroTerm(), CountNodeMatchTerm(:g), CountAbsDiffTerm(:x),
                  CountEdgeCovTerm(Wd; name="dist")]
        fit = fit_ergm_count(dn, fterms)
        @test coefnames(fit) == g.values["fit_names"]
        @test fit.method === :mple && !has_dyad_dependent(fit.model) && fit.converged
        k = length(g.values["fit_seeds"])
        Sm = Float64.(g.values["fit_std_errors_mean"])
        band = max.(g.tolerance["band_sd"] .* Float64.(g.values["fit_coefficients_sd"]) .*
                    sqrt(1 + 1 / k), g.tolerance["floor_se_fraction"] .* Sm)
        se_band = max.(g.tolerance["band_sd"] .* Float64.(g.values["fit_std_errors_sd"]) .*
                       sqrt(1 + 1 / k), g.tolerance["floor_se_fraction"] .* Sm)
        @test all(abs.(coef(fit) .- Float64.(g.values["fit_coefficients_mean"])) .<= band)
        @test all(abs.(stderror(fit) .- Sm) .<= se_band)
    end

    @testset "Golden fixture: a proper geometric model and R's drop (ergm.count MCMLE)" begin
        g = load_golden(joinpath(@__DIR__, "fixtures", "count_geometric_drop.toml"))
        @test g.provenance["ergm_count_version"] == "4.1.3"
        band_sd = g.tolerance["band_sd"]
        floor_frac = g.tolerance["floor_se_fraction"]
        v(key) = Float64.(g.values[key])
        function bands(prefix)
            n_seeds = length(g.values["$(prefix)_seeds"])
            Rsd, Sm, Ssd = v("$(prefix)_coefficients_sd"), v("$(prefix)_std_errors_mean"),
                           v("$(prefix)_std_errors_sd")
            return (max.(band_sd .* Rsd .* sqrt(1 + 1 / n_seeds), floor_frac .* Sm),
                    max.(band_sd .* Ssd .* sqrt(1 + 1 / n_seeds), floor_frac .* Sm))
        end

        # (c) sum + mutual(:min) under the geometric reference, in the regime
        # θ_sum + θ_mutual > 0 > 2θ_sum + θ_mutual: R fits it and its
        # simulations reproduce the data. The default fit (MCMLE) used to be
        # refused as "a mode on the truncation bound"
        gn = count_net(Int(g.values["geometric_n"]),
                       zip(Int.(g.values["geometric_edge_src"]),
                           Int.(g.values["geometric_edge_dst"]),
                           Int.(g.values["geometric_edge_weight"])))
        gt = [SumTerm(), CountMutualTerm(:min)]
        @test [name(t, gn) for t in gt] == g.values["geometric_term_names"]
        @test [compute(t, gn) for t in gt] ≈ v("geometric_summary_statistics") atol = 1e-9
        Rm = v("geometric_coefficients_mean")
        @test Rm[1] + Rm[2] > 0 > 2Rm[1] + Rm[2]
        @test abs(g.values["geometric_r_sim_mean"] - g.values["geometric_observed_mean"]) <
              g.tolerance["sim_mean_tolerance"]
        mp = fit_ergm_count(gn, gt; reference=GeometricReference(), method=:mple, warn=false)
        @test !mp.boundary_mode && !mp.improper && mp.support_control === :converged
        fit = fit_ergm_count(gn, gt; reference=GeometricReference(), rng=Xoshiro(1))
        band, se_band = bands("geometric")
        @test fit.method === :mcmle && fit.converged && !fit.boundary_mode && !fit.improper
        @test all(abs.(coef(fit) .- Rm) .<= band)
        @test all(abs.(stderror(fit) .- v("geometric_std_errors_mean")) .<= se_band)
        # its simulations reproduce the observed mean count, as R's do
        sims = simulate_count_ergm(fit; n_sim=100, rng=Xoshiro(2))
        n_arcs = 20 * 19
        @test abs(mean(compute(SumTerm(), x) for x in sims) / n_arcs -
                  g.values["geometric_observed_mean"]) < g.tolerance["sim_mean_tolerance"]

        # (d) transitiveweights on a perfect matching: 0, its smallest
        # attainable value, with a change statistic of 0 on every dyad. R
        # fixes it at -Inf and fits the rest with it held at 0
        dn = count_net(Int(g.values["drop_n"]),
                       zip(Int.(g.values["drop_edge_src"]), Int.(g.values["drop_edge_dst"]),
                           Int.(g.values["drop_edge_weight"])); directed=false)
        dt = [SumTerm(), NonzeroTerm(), TransitiveWeightsTerm()]
        @test [name(t, dn) for t in dt] == g.values["drop_term_names"]
        @test [compute(t, dn) for t in dt] ≈ v("drop_summary_statistics") atol = 1e-9
        dfit = @test_logs (:warn, r"transitiveweights.min.max.min are at their smallest.*-Inf"s) fit_ergm_count(
            dn, dt; rng=Xoshiro(1))
        band, se_band = bands("drop")
        Dm = v("drop_coefficients_mean")
        @test Dm[3] == -Inf && coef(dfit)[3] == -Inf && stderror(dfit)[3] == 0
        @test dfit.converged && all(abs.(coef(dfit)[1:2] .- Dm[1:2]) .<= band[1:2])
        @test all(abs.(stderror(dfit)[1:2] .- v("drop_std_errors_mean")[1:2]) .<= se_band[1:2])
        @test all(==(0), dfit.mcmc.samples[:, 3])
        # strict mode refuses where R's drop=FALSE would keep the term
        @test_throws ArgumentError fit_ergm_count(dn, dt; rng=Xoshiro(1), drop=false)
    end

    @testset "Terms added in 0.2: sum(pow=), atmost(), CMP (ergm.count fixture)" begin
        g = load_golden(joinpath(@__DIR__, "fixtures", "count_terms.toml"))
        rebuild(prefix; directed) =
            count_net(Int(g.values["$(prefix)_n"]),
                      zip(Int.(g.values["$(prefix)_edge_src"]), Int.(g.values["$(prefix)_edge_dst"]),
                          Int.(g.values["$(prefix)_edge_weight"])); directed=directed)
        added = [SumTerm(pow=2), SumTerm(pow=0.5), AtmostTerm(2), AtmostTerm(0), CMPTerm()]
        for (prefix, directed) in (("zach", false), ("directed", true))
            net = rebuild(prefix; directed=directed)
            @test [name(t, net) for t in added] == g.values["$(prefix)_added_summary_names"]
            vals = [compute(t, net) for t in added]
            @test check_golden(g, "$(prefix)_added_summary", vals) ||
                  error(golden_report(g, "$(prefix)_added_summary", vals))
        end
        neg = rebuild("negative"; directed=true)
        neg_terms = [SumTerm(pow=2), SumTerm(pow=3), AtmostTerm(0), AtmostTerm(-1)]
        @test [name(t, neg) for t in neg_terms] == g.values["negative_added_summary_names"]
        nv_ = [compute(t, neg) for t in neg_terms]
        @test check_golden(g, "negative_added_summary", nv_) ||
              error(golden_report(g, "negative_added_summary", nv_))
        # Where R returns Inf (CMP) or NaN (a fractional power) on negative
        # counts, the term is refused in words — at `compute`, at the model
        # and at a simulation over a negative support
        @test g.values["r_cmp_negative"] == "Inf" && g.values["r_sum_pow_half_negative"] == "NaN"
        for t in (CMPTerm(), SumTerm(pow=0.5))
            @test_throws ArgumentError compute(t, neg)
            @test_throws ArgumentError CountERGMModel([t], neg, DiscUnif2Reference(-2, 2))
            @test_throws ArgumentError simulate_count_ergm(network(4; directed=true), [t], [0.1];
                                                           reference=DiscUnif2Reference(-2, 2))
        end
        @test_throws ArgumentError SumTerm(pow=0)
        @test_throws ArgumentError SumTerm(pow=-1)
        @test SumTerm() === SumTerm(pow=1) && name(SumTerm(pow=1 // 3)) == "sum0.333333333333333"
        @test !ERGM.is_dyad_dependent(AtmostTerm(1)) && !ERGM.is_dyad_dependent(CMPTerm()) &&
              !ERGM.is_dyad_dependent(SumTerm(pow=2))

        # CMP is a TERM (ergm.count's), not a reference: with it the Poisson
        # reference becomes the Conway-Maxwell-Poisson family, which nests the
        # Poisson (θ_CMP = 0) and the geometric (θ_CMP = 1) laws — so its
        # maximized likelihood is at least theirs, and it IS the geometric
        # fit's likelihood when the coefficient is pinned there
        zach = rebuild("zach"; directed=false)
        # (zach is overdispersed, so the fitted CMP law has a heavy tail and
        # the fixed bound carries mass: these are fits of the family on 0:60)
        cmp = fit_ergm_count(zach, [SumTerm(), CMPTerm()]; max_val=60, warn=false)
        pois = fit_ergm_count(zach, [SumTerm()]; max_val=60)
        geom = fit_ergm_count(zach, [SumTerm()]; reference=GeometricReference(), max_val=60)
        @test cmp.converged && is_exact(cmp) == false && !cmp.inference_withheld
        @test loglikelihood(cmp) >= loglikelihood(pois) - 1e-8
        @test loglikelihood(cmp) >= loglikelihood(geom) - 1e-8
        model = CountERGMModel([SumTerm(), CMPTerm()], zach, PoissonReference())
        D = ERGMCount._count_design(model, 0:60)
        pl = ERGMCount._count_derivatives(D, [1, 2], fill(true, 61, length(D.n_tot)))
        @test pl([coef(geom)[1], 1.0])[1] ≈ loglikelihood(geom) atol = 1e-8
        @test pl([coef(pois)[1], 0.0])[1] ≈ loglikelihood(pois) atol = 1e-8
    end

    @testset "Workflows reconstruct the ecosystem layout from [sources]" begin
        pkgdir = dirname(@__DIR__)
        PKG = "ERGMCount"
        EXPECTED = Set(["ERGMCount.jl", "ERGM.jl", "NetworkCore.jl"])
        layout = layout_siblings(pkgdir, EXPECTED, PKG)
        layout.in_layout || @info "Workflow layout step not run: none of the sibling " *
            "checkouts $(join(layout.siblings, ", ")) is beside $(dirname(pkgdir)) (a " *
            "lone checkout or a registry install, where [sources] is not used)."
        # The predicate itself: no sibling → skip; any sibling → run (and a
        # missing one then fails the exact-set assertion)
        mktempdir() do root
            fake = joinpath(root, "$PKG.jl")
            mkpath(fake)
            @test !layout_siblings(fake, EXPECTED, PKG).in_layout
            mkpath(joinpath(root, "ERGM.jl"))
            @test !layout_siblings(fake, EXPECTED, PKG).in_layout      # a bare directory is no checkout
            touch(joinpath(root, "ERGM.jl", "Project.toml"))
            l = layout_siblings(fake, EXPECTED, PKG)
            @test l.in_layout && l.present == ["ERGM.jl"] && l.siblings == ["ERGM.jl", "NetworkCore.jl"]
        end
        for wf in ("CI.yml", "Documentation.yml")
            yml = _readtext(joinpath(pkgdir, ".github", "workflows", wf))
            @test !occursin(r"for pkg in", yml)              # no hand-kept clone list
            @test !occursin("checkout_sources.jl", yml)
            @test occursin("path: $PKG.jl\n", yml)
            step = match(r"\n      - name: Reconstruct the ecosystem layout from \[sources\]\n        shell: julia[^\n]*\n        run: \|\n((?:          [^\n]*\n)+)", yml)
            @test step !== nothing
            step === nothing && continue
            @test first(findfirst("setup-julia", yml)) < step.offset
            # Run the workflow's own step without cloning: in the layout this
            # suite runs in, it must find exactly the siblings [sources] names.
            if !layout.in_layout
                @test_skip layout.in_layout
                continue
            end
            script = replace(step.captures[1], r"^          "m => "")
            out = mktemp() do path, io
                write(io, script); close(io)
                withenv("GITHUB_WORKSPACE" => dirname(pkgdir),
                        "GITHUB_REPOSITORY" => "statistical-network-analysis-with-Julia/$PKG.jl",
                        "LAYOUT_CHECK_ONLY" => "true", "GITHUB_STEP_SUMMARY" => nothing,
                        # `Pkg.test` runs this suite with a sandbox load path
                        # that hides the stdlibs the step loads (`TOML`); the
                        # workflow runs it with Julia's default load path
                        "JULIA_LOAD_PATH" => nothing, "JULIA_PROJECT" => nothing) do
                    read(`$(Base.julia_cmd()) --startup-file=no $path`, String)
                end
            end
            @test Set(m.captures[1] for m in eachmatch(r"^\| (\S+\.jl) \|"m, out)) == EXPECTED
        end
        # The threaded cell that the thread-count-independence tests need
        ci = _readtext(joinpath(pkgdir, ".github", "workflows", "CI.yml"))
        @test occursin("JULIA_NUM_THREADS", ci)
    end

    @testset "Not implemented (README, docs index) and Known limitations (CHANGELOG) agree" begin
        root = dirname(@__DIR__)
        section(path, head) = begin
            text = _readtext(joinpath(root, path))
            m = match(Regex("\\n##+ " * head * "[^\\n]*\\n(.*?)(?=\\n## |\\z)", "s"), text)
            m === nothing ? "" : String(m.captures[1])
        end
        lists = (section("README.md", "Not implemented"),
                 section("docs/src/index.md", "Not implemented"),
                 section("CHANGELOG.md", "Known limitations"))
        @test all(!isempty, lists)
        # One item per limitation, each named in all three places
        for key in ("missing-data MLE", "proposals", "PoissonReference()",
                    "not an `ergm.count` estimator", "ormalisability", "StdNormal",
                    "nodecovar", "transitiveweights", "TransitiveTiesTerm",
                    "negative counts", "wo-mode")
            for (where_, text) in zip(("README", "docs index", "CHANGELOG"), lists)
                @test occursin(key, text) || error("\"$key\" is missing from the $where_ list")
            end
        end
        # MCMLE is implemented: no list may still call it missing
        @test !any(occursin("maximum pseudo-likelihood only", t) for t in lists)
        # The CHANGELOG is release notes; the README has one docs link and the
        # install recipe that works outside this workspace
        changelog = _readtext(joinpath(root, "CHANGELOG.md"))
        @test !occursin(r"panel|round \d|item \d"i, changelog)
        readme = _readtext(joinpath(root, "README.md"))
        @test !occursin("docs-stable", readme) && !occursin("root workspace", readme)
        @test occursin("prepare_workspace.jl", readme)
        @test occursin("statistical-network-analysis-with-julia.github.io/citing/", readme)
        docs_project = _readtext(joinpath(root, "docs", "Project.toml"))
        @test occursin("[compat]", docs_project) && occursin("Documenter = \"1\"", docs_project)
    end

    @testset "Aqua.jl quality assurance" begin
        # Ambiguities are checked against ERGM and NetworkCore too (below);
        # everything else is Aqua's full battery
        Aqua.test_all(ERGMCount; ambiguities=false)
        @test isempty(Test.detect_ambiguities(ERGMCount))
        @test isempty(Test.detect_ambiguities(ERGMCount, ERGM, NetworkCore))
    end
end
