# benchmark/fixtures.jl — networks shared by benchmarks.jl and
# regression_tests.jl (included by both; no top-level side effects beyond the
# constants below).

using ERGMCount
using Networks
using Random

"Zachary's karate club with ergm.count's valued edges, rebuilt from the golden fixture."
function zach_network()
    g = load_golden(joinpath(@__DIR__, "..", "test", "fixtures", "zach_poisson.toml"))
    net = network(Int(g.values["n_actors"]); directed=false)
    s = Int.(g.values["edge_src"]); d = Int.(g.values["edge_dst"]); w = Int.(g.values["edge_weight"])
    for k in eachindex(s)
        add_edge!(net, s[k], d[k])
        set_edge_attribute!(net, :weight, s[k], d[k], w[k])
    end
    return net
end

"""
Directed count network with expected mean out-degree `mean_degree` whatever
`n`: counts are Poisson(mean_degree / n) per dyad, so the O(degree) work of the
strength and triadic profiles is the same at every size and only the number of
dyads grows.
"""
function sparse_count_network(rng::AbstractRNG, n::Int; mean_degree::Int=6)
    net = network(n; directed=true)
    λ = mean_degree / n
    for i in 1:n, j in 1:n
        i == j && continue
        # inverse-CDF Poisson draw: no Distributions dependency needed here
        u = rand(rng); p = exp(-λ); c = p; y = 0
        while u > c
            y += 1; p *= λ / y; c += p
        end
        if y > 0
            add_edge!(net, i, j)
            set_edge_attribute!(net, :weight, i, j, y)
        end
    end
    return net
end

const SWEEP_TERMS = (SumTerm(), NonzeroTerm(), CountMutualTerm(), NodeOSumTerm())
const SWEEP_MEAN_DEGREE = 6
const SWEEP_SUPPORT = 0:10

"Coefficients keeping the chain at ≈ `SWEEP_MEAN_DEGREE` mean out-degree at size n."
sweep_theta(n::Int) = [log(SWEEP_MEAN_DEGREE / n), 0.0, 0.3, -0.01]

"""
Chain state for one Gibbs sweep at size `n`, warmed up so that a measured
update or sweep allocates nothing but what the kernel allocates: two sweeps
at a DENSE specification (`sum` coefficient 2, mean count ≈ 7) first, so
every adjacency vector and both weight dictionaries reach the capacity a
full network needs (containers never shrink), then `warm` sweeps at the
target sparse specification to bring the chain back to mean degree
`SWEEP_MEAN_DEGREE`.
"""
function sweep_state(n::Int; warm::Int=3)
    rng = Random.Xoshiro(n)
    cur = sparse_count_network(rng, n; mean_degree=SWEEP_MEAN_DEGREE)
    weights = get_edge_attribute(cur, :weight, Int)
    support = SWEEP_SUPPORT
    log_h = [log_reference(PoissonReference(), y) for y in support]
    η = zeros(length(support)); buf = zeros(length(support))
    dense = [2.0, 0.0, 0.0, 0.0]
    for _ in 1:2
        ERGMCount._gibbs_sweep!(rng, cur, weights, SWEEP_TERMS, dense, support, log_h, η, buf)
    end
    θ = sweep_theta(n)
    for _ in 1:warm
        ERGMCount._gibbs_sweep!(rng, cur, weights, SWEEP_TERMS, θ, support, log_h, η, buf)
    end
    return (rng=rng, cur=cur, weights=weights, θ=θ, support=support, log_h=log_h, η=η, buf=buf)
end

"One full Gibbs sweep on a prepared state (mutates the state's chain)."
function one_sweep!(st)
    return ERGMCount._gibbs_sweep!(st.rng, st.cur, st.weights, SWEEP_TERMS, st.θ,
                                   st.support, st.log_h, st.η, st.buf)
end

n_dyads(n::Int) = n * (n - 1)
