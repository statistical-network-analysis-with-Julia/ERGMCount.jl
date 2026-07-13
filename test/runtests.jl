using ERGMCount
using ERGM
using Networks
using Graphs
using Distributions
using Random
using Statistics
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

# Set dyad (i,j) to count value y (0 removes the edge)
function set_dyad!(net, i, j, y)
    if y > 0
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
                   NodeOSumTerm(), NodeISumTerm(), NodeSumTerm()]

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
    # Allocation regression on the count-MPLE derivative loop (review finding 15)
    #
    # The conditional moments E[Δg] and E[Δg Δg'] used to be rebuilt per dyad —
    # `zeros(n_terms)`, `zeros(n_terms, n_terms)` and a fresh `p .* (x * x')`
    # outer product for EVERY support value — i.e. n_dyads × (2 + |support|)
    # allocations on every Newton evaluation: 2.1 MB per evaluation on zach at
    # max_val = 30, and it grows with the support, which is exactly the knob a
    # user turns up to make the truncation harmless. The moments are now filled
    # in place, with the same scalar arithmetic in the same order (the golden
    # fixture is an EXACT MLE — this refactor must not move it, and does not).
    # ------------------------------------------------------------------
    @testset "MPLE derivative evaluations allocate O(p²), not O(dyads · |support| · p²)" begin
        function evaluation_allocs(n, max_val)
            net = random_count_net(n; seed=21)
            model = CountERGMModel(AbstractERGMTerm[SumTerm(), NonzeroTerm()], net,
                                   PoissonReference(), is_directed(net))
            core = ERGMCount._count_mple_fit(model, max_val)
            d = ERGMCount._count_derivatives(core.X, core.log_h, core.support,
                                             [1 for _ in 1:core.n_dyads],
                                             core.n_dyads, 2)
            β = [0.2, -0.1]
            d(β)                    # warm up: @allocated on a first call
            return @allocated d(β)  # would measure compilation
        end

        small = evaluation_allocs(6, 5)      # 15 dyads, 6 support values
        big = evaluation_allocs(20, 30)      # 190 dyads, 31 support values
        # 65x the (dyad × support) work, the same allocations: only the (p)
        # gradient and (p×p) Hessian handed to `newton_fit` are new per call.
        @test small <= 512
        @test big <= 512
        @test big <= small + 64
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

    @testset "StatsAPI accessors" begin
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
        @test_throws ArgumentError ergm_count(net, [SumTerm()]; method=:mcmc)
        # Observed count outside a bounded support errors clearly
        @test_throws ArgumentError ergm_count(net, [SumTerm()];
                                              reference=DiscUnifReference(1))
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

    @testset "Missing dyads are rejected" begin
        # Count MPLE enumerates every dyad as observed; a masked dyad would
        # enter the pseudo-likelihood at its face value.
        net = network(5; directed=true)
        add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 2)

        set_missing_dyad!(net, 3, 4)          # absent face
        @test_throws ArgumentError fit_ergm_count(net, [SumTerm()])

        clear_missing_dyads!(net)
        set_missing_dyad!(net, 1, 2)          # present face
        @test_throws ArgumentError fit_ergm_count(net, [SumTerm()])

        clear_missing_dyads!(net)
        @test fit_ergm_count(net, [SumTerm()]) isa CountERGMResult
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
        net = count_net(8, [(1,2,2), (2,3,1), (3,1,3), (1,4,1), (4,5,2),
                            (5,6,1), (6,1,2), (2,7,3), (7,8,1), (8,2,2)])
        terms = [SumTerm(), CountMutualTerm()]      # dyad-DEPENDENT (mutual)
        ref = BinomialReference(3)

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

        # The dependent model's robust SE for `sum` EXCEEDS the Hessian one, by
        # ~20% here — that gap IS the anticonservatism of the pseudo-likelihood
        # Hessian, and it is the whole point of the option.
        #
        # The comparison is made on `sum` because it is the well-identified
        # coefficient: `mutual.count` is separated in this small network (θ̂ ≈ −21
        # with a Hessian SE of ~1.7e4), and a degenerate Hessian SE is *over*stated,
        # not anticonservative — the bootstrap is rightly far smaller there. The
        # direction of the correction is a property of the fit, not a law.
        @test stderror(boot)[1] > stderror(hess)[1]
        @test stderror(boot)[1] / stderror(hess)[1] > 1.1

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

        # Unknown se symbols are rejected, not silently ignored
        @test_throws ArgumentError fit_ergm_count(net, terms; reference=ref,
                                                  se=:sandwich)
        @test_throws ArgumentError fit_ergm_count(net, terms; reference=ref,
                                                  se=:bootstrap, n_boot=1)
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
