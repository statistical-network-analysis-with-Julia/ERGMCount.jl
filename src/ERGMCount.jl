"""
    ERGMCount.jl - ERGMs for Count-Valued Networks

Extends ERGM to handle networks with integer-valued edge weights,
using Poisson, geometric, or binomial reference measures.

The general form is:
    P(Y=y) ∝ h(y) exp(θ' g(y))

where h(y) is the reference measure determining the baseline distribution
for count-valued edges.

Port of the R ergm.count package from the StatNet collection.
"""
module ERGMCount

using Distributions
using ERGM
using Graphs
using LinearAlgebra
using Network
using Optim
using Random
using Statistics
using StatsBase

# Reference measures
export PoissonReference, GeometricReference, BinomialReference
export DiscUnifReference, DiscUnif2Reference

# Count-specific terms
export SumTerm, NonzeroTerm, GreaterthannTerm
export CountMutualTerm, TransitiveTiesTerm, CyclicalTiesTerm
export NodeOSumTerm, NodeISumTerm, NodeSumTerm
export CountAtleastnTerm

# Estimation
export ergm_count, fit_count_ergm

# Simulation
export simulate_count_ergm

# =============================================================================
# Reference Measures
# =============================================================================

"""
    AbstractReferenceMeasure

Base type for reference measures in count ERGMs.
"""
abstract type AbstractReferenceMeasure end

"""
    PoissonReference

Poisson reference measure for count ERGMs.
h(y_ij) = λ^y_ij / y_ij!

# Fields
- `lambda::Float64`: Rate parameter (default 1.0)
"""
struct PoissonReference <: AbstractReferenceMeasure
    lambda::Float64
    PoissonReference(λ::Float64=1.0) = new(λ)
end

function log_reference(ref::PoissonReference, y::Int)
    return y * log(ref.lambda) - logfactorial(y)
end

function sample_reference(ref::PoissonReference)
    return rand(Poisson(ref.lambda))
end

"""
    GeometricReference

Geometric reference measure for count ERGMs.
h(y_ij) = (1-p)^{y_ij}

# Fields
- `prob::Float64`: Success probability (default 0.5)
"""
struct GeometricReference <: AbstractReferenceMeasure
    prob::Float64

    function GeometricReference(p::Float64=0.5)
        0 < p < 1 || throw(ArgumentError("prob must be in (0, 1)"))
        new(p)
    end
end

function log_reference(ref::GeometricReference, y::Int)
    return y * log(1 - ref.prob)
end

function sample_reference(ref::GeometricReference)
    return rand(Geometric(ref.prob))
end

"""
    BinomialReference

Binomial reference measure for count ERGMs (for bounded counts).
h(y_ij) = C(n, y_ij) p^{y_ij} (1-p)^{n-y_ij}

# Fields
- `n::Int`: Number of trials
- `prob::Float64`: Success probability
"""
struct BinomialReference <: AbstractReferenceMeasure
    n::Int
    prob::Float64

    function BinomialReference(n::Int, p::Float64=0.5)
        n > 0 || throw(ArgumentError("n must be positive"))
        0 <= p <= 1 || throw(ArgumentError("prob must be in [0, 1]"))
        new(n, p)
    end
end

function log_reference(ref::BinomialReference, y::Int)
    return logpdf(Binomial(ref.n, ref.prob), y)
end

function sample_reference(ref::BinomialReference)
    return rand(Binomial(ref.n, ref.prob))
end

"""
    DiscUnifReference

Discrete uniform reference measure on {0, 1, ..., max}.
h(y_ij) = 1/(max+1)

# Fields
- `max::Int`: Maximum value
"""
struct DiscUnifReference <: AbstractReferenceMeasure
    max::Int

    function DiscUnifReference(max::Int)
        max >= 0 || throw(ArgumentError("max must be non-negative"))
        new(max)
    end
end

function log_reference(ref::DiscUnifReference, y::Int)
    return -log(ref.max + 1)
end

function sample_reference(ref::DiscUnifReference)
    return rand(0:ref.max)
end

"""
    DiscUnif2Reference

Discrete uniform reference on {a, a+1, ..., b}.
"""
struct DiscUnif2Reference <: AbstractReferenceMeasure
    a::Int
    b::Int

    function DiscUnif2Reference(a::Int, b::Int)
        a <= b || throw(ArgumentError("a must be <= b"))
        new(a, b)
    end
end

function log_reference(ref::DiscUnif2Reference, y::Int)
    return -log(ref.b - ref.a + 1)
end

function sample_reference(ref::DiscUnif2Reference)
    return rand(ref.a:ref.b)
end

# =============================================================================
# Count-Specific Terms
# =============================================================================

"""
    SumTerm <: AbstractERGMTerm

Sum of edge values: ∑_{i,j} y_{ij}
This is the natural sufficient statistic for Poisson reference.
"""
struct SumTerm <: AbstractERGMTerm end

name(::SumTerm) = "sum"

function compute(::SumTerm, net)
    weights = get_edge_attribute(net, :weight)
    isnothing(weights) && return Float64(ne(net))
    return Float64(sum(values(weights)))
end

function change_stat(::SumTerm, net, i::Int, j::Int, delta::Int=1)
    # Change when edge value changes by delta
    return Float64(delta)
end

"""
    NonzeroTerm <: AbstractERGMTerm

Number of non-zero edges: ∑_{i,j} I(y_{ij} > 0)
"""
struct NonzeroTerm <: AbstractERGMTerm end

name(::NonzeroTerm) = "nonzero"

function compute(::NonzeroTerm, net)
    return Float64(ne(net))
end

function change_stat(::NonzeroTerm, net, i::Int, j::Int, old_val::Int, new_val::Int)
    if old_val == 0 && new_val > 0
        return 1.0
    elseif old_val > 0 && new_val == 0
        return -1.0
    else
        return 0.0
    end
end

"""
    GreaterthannTerm <: AbstractERGMTerm

Number of edges with value > n: ∑_{i,j} I(y_{ij} > n)

# Fields
- `threshold::Int`: Threshold value n
"""
struct GreaterthannTerm <: AbstractERGMTerm
    threshold::Int
end

name(t::GreaterthannTerm) = "greaterthan.$(t.threshold)"

function compute(t::GreaterthannTerm, net)
    weights = get_edge_attribute(net, :weight)
    isnothing(weights) && return 0.0
    return Float64(count(w -> w > t.threshold, values(weights)))
end

"""
    CountAtleastnTerm <: AbstractERGMTerm

Number of edges with value >= n.
"""
struct CountAtleastnTerm <: AbstractERGMTerm
    threshold::Int
end

name(t::CountAtleastnTerm) = "atleast.$(t.threshold)"

function compute(t::CountAtleastnTerm, net)
    weights = get_edge_attribute(net, :weight)
    isnothing(weights) && return 0.0
    return Float64(count(w -> w >= t.threshold, values(weights)))
end

"""
    CountMutualTerm <: AbstractERGMTerm

Sum of minimum edge values in reciprocal pairs:
∑_{i<j} min(y_{ij}, y_{ji})
"""
struct CountMutualTerm <: AbstractERGMTerm end

name(::CountMutualTerm) = "mutual.count"

function compute(::CountMutualTerm, net)
    !is_directed(net) && return 0.0

    weights = get_edge_attribute(net, :weight)
    total = 0.0
    n = nv(net)

    for i in 1:n, j in (i+1):n
        w_ij = if !isnothing(weights)
            get(weights, (i, j), 0)
        else
            has_edge(net, i, j) ? 1 : 0
        end
        w_ji = if !isnothing(weights)
            get(weights, (j, i), 0)
        else
            has_edge(net, j, i) ? 1 : 0
        end
        total += min(w_ij, w_ji)
    end

    return total
end

"""
    TransitiveTiesTerm <: AbstractERGMTerm

Weighted transitivity: ∑_{i,j,k} min(y_{ij}, y_{jk}, y_{ik})
"""
struct TransitiveTiesTerm <: AbstractERGMTerm end

name(::TransitiveTiesTerm) = "transitiveties.count"

function compute(::TransitiveTiesTerm, net)
    weights = get_edge_attribute(net, :weight)
    n = nv(net)
    total = 0.0

    get_weight(i, j) = if !isnothing(weights)
        get(weights, (i, j), 0)
    else
        has_edge(net, i, j) ? 1 : 0
    end

    for i in 1:n, j in 1:n, k in 1:n
        i == j || j == k || i == k || continue
        total += min(get_weight(i, j), get_weight(j, k), get_weight(i, k))
    end

    return total
end

"""
    CyclicalTiesTerm <: AbstractERGMTerm

Weighted cyclicality: ∑_{i,j,k} min(y_{ij}, y_{jk}, y_{ki})
"""
struct CyclicalTiesTerm <: AbstractERGMTerm end

name(::CyclicalTiesTerm) = "cyclicalties.count"

function compute(::CyclicalTiesTerm, net)
    !is_directed(net) && return 0.0

    weights = get_edge_attribute(net, :weight)
    n = nv(net)
    total = 0.0

    get_weight(i, j) = if !isnothing(weights)
        get(weights, (i, j), 0)
    else
        has_edge(net, i, j) ? 1 : 0
    end

    for i in 1:n, j in 1:n, k in 1:n
        i == j || j == k || i == k || continue
        total += min(get_weight(i, j), get_weight(j, k), get_weight(k, i))
    end

    return total / 3  # Each cycle counted 3 times
end

"""
    NodeOSumTerm <: AbstractERGMTerm

Sum of outgoing edge values by node: measures activity heterogeneity.
Coefficient on sum of squared out-strengths.
"""
struct NodeOSumTerm <: AbstractERGMTerm end

name(::NodeOSumTerm) = "nodeOSum"

function compute(::NodeOSumTerm, net)
    weights = get_edge_attribute(net, :weight)
    n = nv(net)
    out_strength = zeros(n)

    for e in edges(net)
        w = isnothing(weights) ? 1 : get(weights, (src(e), dst(e)), 1)
        out_strength[src(e)] += w
    end

    return sum(out_strength .^ 2)
end

"""
    NodeISumTerm <: AbstractERGMTerm

Sum of incoming edge values by node: measures popularity heterogeneity.
"""
struct NodeISumTerm <: AbstractERGMTerm end

name(::NodeISumTerm) = "nodeISum"

function compute(::NodeISumTerm, net)
    weights = get_edge_attribute(net, :weight)
    n = nv(net)
    in_strength = zeros(n)

    for e in edges(net)
        w = isnothing(weights) ? 1 : get(weights, (src(e), dst(e)), 1)
        in_strength[dst(e)] += w
    end

    return sum(in_strength .^ 2)
end

"""
    NodeSumTerm <: AbstractERGMTerm

Sum of total (in + out) edge values by node (for undirected or combined).
"""
struct NodeSumTerm <: AbstractERGMTerm end

name(::NodeSumTerm) = "nodeSum"

function compute(::NodeSumTerm, net)
    weights = get_edge_attribute(net, :weight)
    n = nv(net)
    strength = zeros(n)

    for e in edges(net)
        w = isnothing(weights) ? 1 : get(weights, (src(e), dst(e)), 1)
        strength[src(e)] += w
        strength[dst(e)] += w
    end

    return sum(strength .^ 2)
end

# =============================================================================
# Model and Estimation
# =============================================================================

"""
    CountERGMModel{T}

ERGM model for count-valued networks.
"""
struct CountERGMModel{T}
    terms::Vector{AbstractERGMTerm}
    network::Network{T}
    reference::AbstractReferenceMeasure
    directed::Bool
end

"""
    CountERGMResult

Results from fitting a count ERGM.
"""
struct CountERGMResult{T}
    model::CountERGMModel{T}
    coefficients::Vector{Float64}
    std_errors::Vector{Float64}
    loglik::Float64
    converged::Bool
end

function Base.show(io::IO, result::CountERGMResult)
    println(io, "Count ERGM Results")
    println(io, "==================")
    println(io, "Reference: $(typeof(result.model.reference))")
    println(io, "Log-likelihood: $(round(result.loglik, digits=4))")
    println(io, "Converged: $(result.converged)")
    println(io)
    println(io, "Coefficients:")
    for (i, term) in enumerate(result.model.terms)
        println(io, "  $(rpad(name(term), 20)) $(lpad(round(result.coefficients[i], digits=4), 10)) " *
                    "(SE: $(round(result.std_errors[i], digits=4)))")
    end
end

"""
    ergm_count(net::Network, terms; reference=PoissonReference(), kwargs...)

Fit an ERGM for count-valued networks.

# Arguments
- `net`: Network with edge weights
- `terms`: Vector of ERGM terms
- `reference`: Reference measure (default: Poisson)

# Returns
- `CountERGMResult`: Fitted model results
"""
function ergm_count(net::Network{T}, terms::Vector{<:AbstractERGMTerm};
                    reference::AbstractReferenceMeasure=PoissonReference(),
                    method::Symbol=:mple,
                    maxiter::Int=100) where T

    model = CountERGMModel{T}(terms, net, reference, is_directed(net))

    if method == :mple
        return count_mple(model; maxiter=maxiter)
    else
        throw(ArgumentError("Unknown method: $method"))
    end
end

fit_count_ergm = ergm_count

"""
    count_mple(model::CountERGMModel; kwargs...) -> CountERGMResult

Maximum Pseudo-Likelihood Estimation for count ERGM.
"""
function count_mple(model::CountERGMModel{T}; maxiter::Int=100, tol::Float64=1e-6) where T
    n_terms = length(model.terms)
    n = nv(model.network)
    coef = zeros(n_terms)

    # Compute observed statistics
    obs_stats = [compute(term, model.network) for term in model.terms]

    # Newton-Raphson optimization
    for iter in 1:maxiter
        grad = zeros(n_terms)
        hess = zeros(n_terms, n_terms)

        for i in 1:n
            for j in (model.directed ? (1:n) : (i+1:n))
                i == j && continue

                # Current edge value
                weights = get_edge_attribute(model.network, :weight)
                y_obs = if !isnothing(weights)
                    get(weights, (i, j), 0)
                else
                    has_edge(model.network, i, j) ? 1 : 0
                end

                # Compute expected statistics under current parameters
                # For MPLE, we approximate by conditioning on rest of network
                # This is a simplified version

                delta = [change_stat(term, model.network, i, j) for term in model.terms]
                eta = dot(coef, delta)
                prob = 1.0 / (1.0 + exp(-eta))

                grad .+= delta .* (y_obs > 0 ? 1.0 : 0.0 - prob)
                hess .-= (prob * (1 - prob)) .* (delta * delta')
            end
        end

        # Update
        if det(hess) != 0 && !any(isnan, hess)
            step = -hess \ grad
            coef .+= step

            if maximum(abs.(step)) < tol
                se = sqrt.(abs.(diag(pinv(-hess))))
                return CountERGMResult{T}(model, coef, se, NaN, true)
            end
        end
    end

    se = fill(NaN, n_terms)
    return CountERGMResult{T}(model, coef, se, NaN, false)
end

# =============================================================================
# Simulation
# =============================================================================

"""
    simulate_count_ergm(result::CountERGMResult; n_sim=1, kwargs...) -> Vector{Network}

Simulate networks from a fitted count ERGM.
"""
function simulate_count_ergm(result::CountERGMResult{T};
                             n_sim::Int=1,
                             burnin::Int=1000,
                             interval::Int=100,
                             max_val::Int=20) where T
    networks = Network{T}[]

    current = deepcopy(result.model.network)
    n = nv(current)
    directed = is_directed(current)

    for sim in 1:(burnin + n_sim * interval)
        # Gibbs sampling: update each edge
        for i in 1:n
            for j in (directed ? (1:n) : (i+1:n))
                i == j && continue

                # Sample new value for edge (i,j)
                # Compute conditional distribution
                log_probs = Float64[]
                for y in 0:max_val
                    # Temporarily set edge to value y
                    log_p = log_reference(result.model.reference, y)
                    for (k, term) in enumerate(result.model.terms)
                        # This is approximate - would need proper change stat
                        log_p += result.coefficients[k] * y
                    end
                    push!(log_probs, log_p)
                end

                # Normalize and sample
                log_probs .-= maximum(log_probs)
                probs = exp.(log_probs)
                probs ./= sum(probs)

                new_val = sample(0:max_val, Weights(probs))

                # Update network
                if new_val > 0
                    if !has_edge(current, i, j)
                        add_edge!(current, i, j)
                    end
                    set_edge_attribute!(current, i, j, :weight, new_val)
                else
                    if has_edge(current, i, j)
                        rem_edge!(current, i, j)
                    end
                end
            end
        end

        if sim > burnin && (sim - burnin) % interval == 0
            push!(networks, deepcopy(current))
        end
    end

    return networks
end

end # module
