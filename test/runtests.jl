using ERGMCount
using ERGM
using Network
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

    @testset "gof extends the shared Network.gof generic" begin
        # One generic across the ecosystem: the method is added to
        # Network.gof, not a package-local function
        @test ERGMCount.gof === Network.gof

        net = random_count_net(6; seed=19)
        result = fit_ergm_count(net, [SumTerm(), NonzeroTerm()])
        g = ERGMCount.gof(result; n_sim=8, burnin=10, interval=2,
                          rng=Random.Xoshiro(31))
        @test g isa Network.GOFResult
        @test Network.n_simulations(g) == 8
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
end
