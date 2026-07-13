"""
    ERGMCount.jl - ERGMs for Count-Valued Networks

Extends ERGM to handle networks with integer-valued edge weights,
using Poisson, geometric, or binomial reference measures.

The general form is:
    P(Y=y) ∝ h(y) exp(θ' g(y))

where h(y) is the reference measure determining the baseline distribution
for count-valued edges (Krivitsky 2012).

Port of the R ergm.count package from the StatNet collection.
"""
module ERGMCount

using Distributions
using ERGM
using Graphs
using LinearAlgebra
using Networks
using Random
using StatsBase

import ERGM: name, compute, is_dyad_dependent
# Shared presentation infrastructure (Networks.jl): the ONE `gof` generic all
# model packages extend, plus the common coefficient-table printer and
# GOF containers
import Networks: gof, print_coeftable, GOFStatistic, GOFResult

# The shared result-metadata protocol (Networks.jl `src/results.jl`): the
# generic accessors that say what a fit actually did. Imported by name because
# ERGMCount adds methods for `CountERGMResult`; `fit_metadata(fit)` collects them.
import Networks: estimand, objective, is_exact, se_method, missing_method,
                 approximations
import StatsAPI
import StatsAPI: coef, stderror, vcov, loglikelihood, nobs, dof

# Reference measures
export PoissonReference, GeometricReference, BinomialReference
export DiscUnifReference, DiscUnif2Reference
export log_reference, sample_reference

# Support truncation: which references the `0:max_val` enumeration truncates,
# and the tolerance past which a fit is warned to be leaning on the bound
export is_truncating, BOUNDARY_MASS_TOL

# Count-specific terms
export SumTerm, NonzeroTerm, GreaterthannTerm
export CountMutualTerm, TransitiveTiesTerm, CyclicalTiesTerm
export NodeOSumTerm, NodeISumTerm, NodeSumTerm
export CountAtleastnTerm
export change_stat_count

# Model / estimation
export CountERGMModel, CountERGMResult
export fit_ergm_count, ergm_count, fit_count_ergm

# Simulation
export simulate_count_ergm

# Goodness of fit (method of the shared Networks.jl `gof` generic)
export gof

# StatsAPI methods (re-exported so `coef(fit)` etc. work with just
# `using ERGMCount`)
export coef, stderror, vcov, loglikelihood, nobs, dof

# =============================================================================
# Reference Measures
# =============================================================================

"""
    AbstractReferenceMeasure

Base type for reference measures in count ERGMs.
"""
abstract type AbstractReferenceMeasure end

# log(y!) without external dependencies; exact for the modest counts used
# in dyad supports
_logfactorial(y::Int) = sum(log, 2:y; init=0.0)

"""
    log_reference(ref::AbstractReferenceMeasure, y::Int) -> Float64

Log of the dyadwise reference measure `h(y)` evaluated at count `y`.
Every concrete reference measure ([`PoissonReference`](@ref ERGMCount.PoissonReference),
[`GeometricReference`](@ref ERGMCount.GeometricReference),
[`BinomialReference`](@ref ERGMCount.BinomialReference),
[`DiscUnifReference`](@ref ERGMCount.DiscUnifReference),
[`DiscUnif2Reference`](@ref ERGMCount.DiscUnif2Reference)) implements a
method; values outside a bounded support return `-Inf`.
"""
function log_reference end

"""
    sample_reference(ref::AbstractReferenceMeasure;
                     rng::Random.AbstractRNG=Random.default_rng()) -> Int

Draw a single random count from the baseline distribution associated with
the reference measure `ref` (e.g. `Poisson(lambda)` for
[`PoissonReference`](@ref ERGMCount.PoissonReference)). All randomness flows through the `rng`
keyword, so the same rng state yields identical draws.
"""
function sample_reference end

"""
    PoissonReference

Poisson reference measure for count ERGMs.
h(y_ij) = λ^y_ij / y_ij!

With this reference and a `SumTerm` coefficient θ, each dyad is
conditionally Poisson(λ·e^θ) (Krivitsky 2012).

# Fields
- `lambda::Float64`: Rate parameter (default 1.0)
"""
struct PoissonReference <: AbstractReferenceMeasure
    lambda::Float64
    PoissonReference(λ::Float64=1.0) = new(λ)
end

function log_reference(ref::PoissonReference, y::Int)
    return y * log(ref.lambda) - _logfactorial(y)
end

function sample_reference(ref::PoissonReference;
                          rng::Random.AbstractRNG=Random.default_rng())
    return rand(rng, Poisson(ref.lambda))
end

"""
    GeometricReference

Geometric reference measure for count ERGMs: the *counting measure*
h(y_ij) = 1 on {0, 1, 2, …}, as in `ergm.count`/Krivitsky (2012). The
geometric shape of the dyad distribution comes from a negative `SumTerm`
coefficient, not from the reference itself, so this measure has no free
parameter.
"""
struct GeometricReference <: AbstractReferenceMeasure end

log_reference(::GeometricReference, y::Int) = 0.0

# The counting measure is improper on its own; sample from the geometric
# shape a unit negative Sum coefficient would induce
sample_reference(::GeometricReference;
                 rng::Random.AbstractRNG=Random.default_rng()) =
    rand(rng, Geometric(1 - exp(-1)))

"""
    BinomialReference

Binomial reference measure for count ERGMs (for bounded counts):
h(y_ij) = C(trials, y_ij), matching `ergm.count`'s `Binomial(trials)`.
The success probability is absorbed into the estimated `SumTerm`
coefficient rather than parameterizing the reference.

# Fields
- `trials::Int`: Number of trials (maximum count value)
"""
struct BinomialReference <: AbstractReferenceMeasure
    trials::Int

    function BinomialReference(trials::Int)
        trials > 0 || throw(ArgumentError("trials must be positive"))
        new(trials)
    end
end

function log_reference(ref::BinomialReference, y::Int)
    (0 <= y <= ref.trials) || return -Inf
    return _logfactorial(ref.trials) - _logfactorial(y) -
           _logfactorial(ref.trials - y)
end

function sample_reference(ref::BinomialReference;
                          rng::Random.AbstractRNG=Random.default_rng())
    return rand(rng, Binomial(ref.trials, 0.5))
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

function sample_reference(ref::DiscUnifReference;
                          rng::Random.AbstractRNG=Random.default_rng())
    return rand(rng, 0:ref.max)
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

function sample_reference(ref::DiscUnif2Reference;
                          rng::Random.AbstractRNG=Random.default_rng())
    return rand(rng, ref.a:ref.b)
end

# Dyad-value support used when enumerating conditional distributions.
# `max_val` truncates unbounded supports (Poisson/geometric).
_support(::PoissonReference, max_val::Int) = 0:max_val
_support(::GeometricReference, max_val::Int) = 0:max_val
_support(ref::BinomialReference, ::Int) = 0:ref.trials
_support(ref::DiscUnifReference, ::Int) = 0:ref.max
_support(ref::DiscUnif2Reference, ::Int) = ref.a:ref.b

"""
    is_truncating(ref::AbstractReferenceMeasure) -> Bool

Whether enumerating this reference's support over `0:max_val` **truncates** it.

`true` for the mathematically unbounded references (`PoissonReference`,
`GeometricReference`): there the finite enumeration is a numerical device, and
the fitted model is only an approximation to the documented unbounded family
insofar as negligible conditional mass sits at the bound — which
[`count_mple`](@ref) measures and reports (see `boundary_mass` on
[`CountERGMResult`](@ref)).

`false` for the genuinely bounded references (`BinomialReference`,
`DiscUnifReference`, `DiscUnif2Reference`), whose support is part of the model.
"""
is_truncating(::AbstractReferenceMeasure) = false
is_truncating(::PoissonReference) = true
is_truncating(::GeometricReference) = true

"""
    BOUNDARY_MASS_TOL

Tolerance for the truncation boundary-mass diagnostic. If a fitted model puts
more than this share of any dyad's conditional mass on the top support value,
`count_mple` warns that the truncation is materially changing the model rather
than merely bounding the arithmetic.
"""
const BOUNDARY_MASS_TOL = 1e-4

# =============================================================================
# Dyad values
# =============================================================================

# Canonical edge-attribute key: (i,j) directed, (min,max) undirected
_wkey(net, i::Int, j::Int) = is_directed(net) ? (i, j) : minmax(i, j)

"""
    dyad_value(net, weights, i, j) -> Int

The count value of dyad (i,j): 0 when the edge is absent, its `:weight`
attribute when present (edges without the attribute count as 1).
"""
function dyad_value(net, weights, i::Int, j::Int)
    has_edge(net, i, j) || return 0
    return Int(get(weights, _wkey(net, i, j), 1))
end

_get_weights(net) = get_edge_attribute(net, :weight)

# =============================================================================
# Count-Specific Terms
# =============================================================================
#
# Each term implements:
#   compute(term, net) — the full statistic
#   change_stat_count(term, net, weights, i, j, old, new) — the change in
#     the statistic when dyad (i,j) moves from value `old` to `new`,
#     holding all other dyads fixed. Implementations must not read the
#     dyad's own current value from the network (old/new are authoritative).

"""
    change_stat_count(term, net, weights, i, j, old, new) -> Float64

Change in `compute(term, net)` when dyad (i,j) moves from count `old` to
count `new`, holding all other dyads fixed. `weights` is the `:weight`
edge-attribute dictionary; hot loops should pass the typed snapshot
`get_edge_attribute(net, :weight, Int)` (a `Dict{Tuple{T,T},Int}`) rather
than the untyped `get_edge_attribute(net, :weight)` Dict.
"""
function change_stat_count end

"""
    SumTerm <: AbstractERGMTerm

Sum of edge values: ∑_{i,j} y_{ij}
This is the natural sufficient statistic for the Poisson reference.
"""
struct SumTerm <: AbstractERGMTerm end

name(::SumTerm) = "sum"

function compute(::SumTerm, net)
    weights = _get_weights(net)
    total = 0.0
    for e in edges(net)
        total += Float64(get(weights, _wkey(net, src(e), dst(e)), 1))
    end
    return total
end

change_stat_count(::SumTerm, net, weights, i::Int, j::Int, old::Int, new::Int) =
    Float64(new - old)

"""
    NonzeroTerm <: AbstractERGMTerm

Number of non-zero edges: ∑_{i,j} I(y_{ij} > 0)
"""
struct NonzeroTerm <: AbstractERGMTerm end

name(::NonzeroTerm) = "nonzero"

compute(::NonzeroTerm, net) = Float64(ne(net))

change_stat_count(::NonzeroTerm, net, weights, i::Int, j::Int, old::Int, new::Int) =
    Float64((new > 0) - (old > 0))

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
    weights = _get_weights(net)
    total = 0.0
    for e in edges(net)
        w = get(weights, _wkey(net, src(e), dst(e)), 1)
        w > t.threshold && (total += 1.0)
    end
    return total
end

change_stat_count(t::GreaterthannTerm, net, weights, i::Int, j::Int, old::Int, new::Int) =
    Float64((new > t.threshold) - (old > t.threshold))

"""
    CountAtleastnTerm <: AbstractERGMTerm

Number of edges with value >= n.
"""
struct CountAtleastnTerm <: AbstractERGMTerm
    threshold::Int
end

name(t::CountAtleastnTerm) = "atleast.$(t.threshold)"

function compute(t::CountAtleastnTerm, net)
    weights = _get_weights(net)
    total = 0.0
    for e in edges(net)
        w = get(weights, _wkey(net, src(e), dst(e)), 1)
        w >= t.threshold && (total += 1.0)
    end
    return total
end

change_stat_count(t::CountAtleastnTerm, net, weights, i::Int, j::Int, old::Int, new::Int) =
    Float64((new >= t.threshold) - (old >= t.threshold))

"""
    CountMutualTerm <: AbstractERGMTerm

Sum of minimum edge values in reciprocal pairs:
∑_{i<j} min(y_{ij}, y_{ji})
"""
struct CountMutualTerm <: AbstractERGMTerm end

name(::CountMutualTerm) = "mutual.count"

function compute(::CountMutualTerm, net)
    !is_directed(net) && return 0.0

    weights = _get_weights(net)
    total = 0.0
    n = nv(net)

    for i in 1:n, j in (i+1):n
        total += min(dyad_value(net, weights, i, j),
                     dyad_value(net, weights, j, i))
    end

    return total
end

function change_stat_count(::CountMutualTerm, net, weights, i::Int, j::Int,
                           old::Int, new::Int)
    is_directed(net) || return 0.0
    y_ji = dyad_value(net, weights, j, i)
    return Float64(min(new, y_ji) - min(old, y_ji))
end

"""
    TransitiveTiesTerm <: AbstractERGMTerm

Triadic-minimum transitivity over ordered distinct triples:
∑_{i≠j≠k} min(y_{ij}, y_{jk}, y_{ik})

Note this is a simple valued transitivity; it is related to but not
identical to `ergm.count`'s `transitiveweights(min, max, min)` statistic.
"""
struct TransitiveTiesTerm <: AbstractERGMTerm end

name(::TransitiveTiesTerm) = "transitiveties.count"

function compute(::TransitiveTiesTerm, net)
    weights = _get_weights(net)
    n = nv(net)
    total = 0.0

    w(i, j) = dyad_value(net, weights, i, j)

    for i in 1:n, j in 1:n, k in 1:n
        (i == j || j == k || i == k) && continue
        total += min(w(i, j), w(j, k), w(i, k))
    end

    return total
end

function change_stat_count(::TransitiveTiesTerm, net, weights, i::Int, j::Int,
                           old::Int, new::Int)
    n = nv(net)
    w(a, b) = dyad_value(net, weights, a, b)

    delta = 0.0
    if is_directed(net)
        for k in 1:n
            (k == i || k == j) && continue
            # Dyad (i,j) in role (a,b): triples (i, j, k) use min(y_ij, y_jk, y_ik)
            delta += min(new, w(j, k), w(i, k)) - min(old, w(j, k), w(i, k))
            # Role (b,c): triples (k, i, j) use min(y_ki, y_ij, y_kj)
            delta += min(w(k, i), new, w(k, j)) - min(w(k, i), old, w(k, j))
            # Role (a,c): triples (i, k, j) use min(y_ik, y_kj, y_ij)
            delta += min(w(i, k), w(k, j), new) - min(w(i, k), w(k, j), old)
        end
    else
        # Undirected: every ordered distinct triple over {i, j, k} uses the
        # same three pair values, so pair {i,j} appears in 6 ordered triples
        # per third vertex k
        for k in 1:n
            (k == i || k == j) && continue
            delta += 6 * (min(new, w(j, k), w(i, k)) - min(old, w(j, k), w(i, k)))
        end
    end
    return delta
end

"""
    CyclicalTiesTerm <: AbstractERGMTerm

Weighted cyclicality: (1/3) ∑_{i≠j≠k} min(y_{ij}, y_{jk}, y_{ki})
(each 3-cycle counted once).
"""
struct CyclicalTiesTerm <: AbstractERGMTerm end

name(::CyclicalTiesTerm) = "cyclicalties.count"

function compute(::CyclicalTiesTerm, net)
    !is_directed(net) && return 0.0

    weights = _get_weights(net)
    n = nv(net)
    total = 0.0

    w(i, j) = dyad_value(net, weights, i, j)

    for i in 1:n, j in 1:n, k in 1:n
        (i == j || j == k || i == k) && continue
        total += min(w(i, j), w(j, k), w(k, i))
    end

    return total / 3  # Each cycle counted 3 times (rotations)
end

function change_stat_count(::CyclicalTiesTerm, net, weights, i::Int, j::Int,
                           old::Int, new::Int)
    is_directed(net) || return 0.0
    n = nv(net)
    w(a, b) = dyad_value(net, weights, a, b)

    # Dyad (i,j) appears once in each of the 3 rotations of a cycle
    # {i→j, j→k, k→i}; the statistic divides by 3, so the net change is
    # one un-rotated sum over k
    delta = 0.0
    for k in 1:n
        (k == i || k == j) && continue
        delta += min(new, w(j, k), w(k, i)) - min(old, w(j, k), w(k, i))
    end
    return delta
end

"""
    NodeOSumTerm <: AbstractERGMTerm

Sum of squared out-strengths: measures activity heterogeneity.
Directed networks only (0 for undirected; use `NodeSumTerm` there).
"""
struct NodeOSumTerm <: AbstractERGMTerm end

name(::NodeOSumTerm) = "nodeOSum"

function compute(::NodeOSumTerm, net)
    is_directed(net) || return 0.0
    weights = _get_weights(net)
    n = nv(net)
    out_strength = zeros(n)

    for e in edges(net)
        w = get(weights, _wkey(net, src(e), dst(e)), 1)
        out_strength[src(e)] += w
    end

    return sum(out_strength .^ 2)
end

function _out_strength(net, weights, v::Int)
    s = 0.0
    for u in outneighbors(net, v)
        s += dyad_value(net, weights, v, u)
    end
    return s
end

function _in_strength(net, weights, v::Int)
    s = 0.0
    for u in inneighbors(net, v)
        s += dyad_value(net, weights, u, v)
    end
    return s
end

function change_stat_count(::NodeOSumTerm, net, weights, i::Int, j::Int,
                           old::Int, new::Int)
    is_directed(net) || return 0.0
    # Out-strength of i excluding the dyad's own contribution
    s = _out_strength(net, weights, i) - dyad_value(net, weights, i, j)
    return (s + new)^2 - (s + old)^2
end

"""
    NodeISumTerm <: AbstractERGMTerm

Sum of squared in-strengths: measures popularity heterogeneity.
Directed networks only (0 for undirected; use `NodeSumTerm` there).
"""
struct NodeISumTerm <: AbstractERGMTerm end

name(::NodeISumTerm) = "nodeISum"

function compute(::NodeISumTerm, net)
    is_directed(net) || return 0.0
    weights = _get_weights(net)
    n = nv(net)
    in_strength = zeros(n)

    for e in edges(net)
        w = get(weights, _wkey(net, src(e), dst(e)), 1)
        in_strength[dst(e)] += w
    end

    return sum(in_strength .^ 2)
end

function change_stat_count(::NodeISumTerm, net, weights, i::Int, j::Int,
                           old::Int, new::Int)
    is_directed(net) || return 0.0
    s = _in_strength(net, weights, j) - dyad_value(net, weights, i, j)
    return (s + new)^2 - (s + old)^2
end

"""
    NodeSumTerm <: AbstractERGMTerm

Sum of squared total (in + out) strengths.
"""
struct NodeSumTerm <: AbstractERGMTerm end

name(::NodeSumTerm) = "nodeSum"

function compute(::NodeSumTerm, net)
    weights = _get_weights(net)
    n = nv(net)
    strength = zeros(n)

    for e in edges(net)
        w = get(weights, _wkey(net, src(e), dst(e)), 1)
        strength[src(e)] += w
        strength[dst(e)] += w
    end

    return sum(strength .^ 2)
end

function change_stat_count(::NodeSumTerm, net, weights, i::Int, j::Int,
                           old::Int, new::Int)
    y_ij = dyad_value(net, weights, i, j)
    # Total strength of a vertex, excluding the (i,j) dyad's contribution.
    # For undirected networks the stored edge contributes to both
    # endpoints once via compute's src/dst accumulation.
    if is_directed(net)
        s_i = _out_strength(net, weights, i) + _in_strength(net, weights, i) - y_ij
        s_j = _out_strength(net, weights, j) + _in_strength(net, weights, j) - y_ij
    else
        # Undirected: strength via unique edges
        s_i = _out_strength(net, weights, i) - y_ij
        s_j = _out_strength(net, weights, j) - y_ij
    end
    d_old, d_new = Float64(old), Float64(new)
    return (s_i + d_new)^2 - (s_i + d_old)^2 + (s_j + d_new)^2 - (s_j + d_old)^2
end

# =============================================================================
# Dependence classification
# =============================================================================
#
# Extends `ERGM.is_dyad_dependent` (whose fallback is the conservative `true`).
# A count term is dyad-independent when its change statistic depends only on the
# old and new value of the dyad being changed — which is exactly what makes the
# dyadwise pseudo-likelihood the likelihood. The strength terms are quadratic in
# the dyads and the mutual/triadic terms read other dyads, so they keep the
# conservative default. This is what `is_exact(::CountERGMResult)` reads.

is_dyad_dependent(::SumTerm) = false
is_dyad_dependent(::NonzeroTerm) = false
is_dyad_dependent(::GreaterthannTerm) = false
is_dyad_dependent(::CountAtleastnTerm) = false

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

`loglik` is the maximized *pseudo*-log-likelihood over the (possibly
truncated) dyad support; `vcov` is the inverse negative Hessian of the
pseudo-log-likelihood at the optimum.

# Support and truncation

Each dyad's conditional is enumerated over `0:max_val`. For the unbounded
references (Poisson, geometric) that enumeration **truncates** the model, so
the estimand is the documented unbounded family only if negligible mass sits at
the bound. Two fields record this, so the approximation is inspectable rather
than implicit:

- `truncated` — whether the reference is unbounded and was therefore truncated
  (see [`is_truncating`](@ref)). `false` for genuinely bounded references, whose
  support is part of the model.
- `boundary_mass` — the largest conditional probability any dyad places on the
  top support value, at the fitted coefficients. A value materially above zero
  means the bound is shaping the fit; `count_mple` warns past
  `BOUNDARY_MASS_TOL`. Always `0.0` when `truncated == false`.

# Standard errors

`se_type` records how `std_errors`/`vcov` were actually obtained — `:hessian`
(the inverse negative pseudo-Hessian, anticonservative under dyadic dependence)
or `:bootstrap` (the parametric bootstrap of `count_mple(model; se=:bootstrap)`).
It is what `Networks.se_method(fit)` reports, and what the `show` method reads
before deciding whether an anticonservatism caveat is still warranted.
"""
struct CountERGMResult{T}
    model::CountERGMModel{T}
    coefficients::Vector{Float64}
    std_errors::Vector{Float64}
    vcov::Matrix{Float64}
    loglik::Float64
    converged::Bool
    max_val::Int
    truncated::Bool
    boundary_mass::Float64
    se_type::Symbol
end

# Backwards-compatible constructor: a result built without an `se_type` reports
# the inverse-Hessian standard errors it in fact had.
CountERGMResult{T}(model, coefficients, std_errors, vcov, loglik, converged,
                   max_val, truncated, boundary_mass) where T =
    CountERGMResult{T}(model, coefficients, std_errors, vcov, loglik, converged,
                       max_val, truncated, boundary_mass, :hessian)

# Two-sided normal p-values via the complementary CDF (the naive
# 2(1 − cdf) form underflows to exactly 0 beyond |z| ≈ 8.3); NaN standard
# errors give NaN p-values, which the shared printer renders as "NaN"
_z_pvalues(z::AbstractVector{Float64}) = 2 .* ccdf.(Normal(), abs.(z))

function Base.show(io::IO, result::CountERGMResult)
    println(io, "Count ERGM Results")
    println(io, "==================")
    println(io, "Reference: $(typeof(result.model.reference))")
    # The support is part of the estimand, not an implementation detail: print
    # it, and mark it as a truncation when the reference is really unbounded.
    if result.truncated
        println(io, "Support:   0:$(result.max_val)  (TRUNCATED — reference is unbounded)")
        println(io, "Boundary mass: $(round(result.boundary_mass, sigdigits=3)) " *
                    "(max over dyads, at the fitted coefficients)")
    else
        println(io, "Support:   0:$(result.max_val)  (bounded reference)")
    end
    println(io, "Pseudo-log-likelihood: $(round(result.loglik, digits=4))")
    println(io, "Converged: $(result.converged)")
    println(io, "Std. errors: ", result.se_type === :bootstrap ?
                "parametric bootstrap" : "inverse pseudo-Hessian")
    println(io)
    z = result.coefficients ./ result.std_errors
    print_coeftable(io, [name(term) for term in result.model.terms],
                    result.coefficients, result.std_errors, _z_pvalues(z);
                    z_values=z)

    # Honest-uncertainty caveat, and the prose twin of what
    # `approximations(result)` reports: the pseudo-likelihood multiplies dyad
    # conditionals as if independent, so the inverse-Hessian standard errors of a
    # dyad-dependent model are expected anticonservative. A dyad-independent model
    # needs no caveat (there the pseudo-likelihood is the likelihood), and neither
    # does a bootstrap fit — those standard errors do NOT assume independence, so
    # claiming they are anticonservative would be a lie.
    if _has_dyad_dependent(result.model) && result.se_type === :hessian
        println(io)
        println(io, "Warning: this model contains dyad-dependent terms and was fit by")
        println(io, "maximum pseudolikelihood. The standard errors are the inverse")
        println(io, "pseudo-Hessian and are expected to be anticonservative; refit with")
        println(io, "`se=:bootstrap` for a parametric-bootstrap covariance.")
    end
end

# ============================================================================
# The shared result-metadata protocol (Networks.jl `src/results.jl`)
# ============================================================================
#
# `fit_metadata(fit)` collects these accessors, so the truncation the `show`
# method prints in prose is also machine-readable — the two are derived from
# the same `truncated`/`boundary_mass` fields and cannot disagree.

estimand(::CountERGMResult) = :count_ergm

objective(::CountERGMResult) = :pseudolikelihood

"""
    _has_dyad_dependent(model::CountERGMModel) -> Bool

Whether any term of the count formula is dyad-dependent (see
`ERGM.is_dyad_dependent`). This is THE predicate that decides whether the
dyadwise pseudo-likelihood is an approximation, and it is defined once so that
`show`'s prose caveat, `is_exact` and `approximations` all read the same answer.
"""
_has_dyad_dependent(model::CountERGMModel) =
    any(is_dyad_dependent(t) for t in model.terms)

"""
    is_exact(result::CountERGMResult) -> Bool

`true` only when **both** conditions hold:

1. every term is dyad-independent (see `ERGM.is_dyad_dependent`), so the dyad
   conditionals that the pseudo-likelihood multiplies are the model's own
   conditionals and their product is the likelihood; and
2. the fit was not truncated — for an unbounded reference (Poisson, geometric)
   the estimator enumerates `0:max_val`, which is the likelihood of a
   *different*, truncated family.

A `SumTerm` model under a bounded reference is therefore exact; the same term
under a Poisson reference is not, and neither is any model containing a
strength, mutual or triadic term.
"""
is_exact(result::CountERGMResult) =
    !result.truncated && !_has_dyad_dependent(result.model)

"""
    se_method(result::CountERGMResult) -> Symbol

What the reported standard errors ACTUALLY are: `:hessian` (the inverse negative
pseudo-Hessian) or `:bootstrap` (the parametric bootstrap of
`count_mple(model; se=:bootstrap)`). Read straight off the fit, so it can never
claim an estimator that was not used.
"""
se_method(result::CountERGMResult) = result.se_type

# `fit_ergm_count` calls `require_observed` with the default `:error` policy:
# the count MPLE enumerates every dyad as observed, so masked data is refused.
missing_method(::CountERGMResult) = :rejected

function approximations(result::CountERGMResult)
    out = String[]
    if result.truncated
        push!(out, "count support truncated at 0:$(result.max_val); max boundary " *
                   "mass $(round(result.boundary_mass, sigdigits=3)) " *
                   "(the reference measure is unbounded, so the enumerated support " *
                   "is an approximation to the model's)")
    end
    if _has_dyad_dependent(result.model)
        # The POINT ESTIMATE is a pseudo-likelihood estimate however the standard
        # errors were computed: the bootstrap replaces the covariance, not θ̂.
        push!(out, "maximum pseudo-likelihood of a dyad-dependent model: the dyad " *
                   "conditionals are multiplied as if independent, so the point " *
                   "estimates are biased in finite samples")
        if result.se_type === :hessian
            push!(out, "inverse-Hessian standard errors of the naive pseudo-likelihood: " *
                       "expected anticonservative under dyadic dependence (refit with " *
                       "`se=:bootstrap` for a parametric-bootstrap covariance)")
        end
    end
    if result.se_type === :bootstrap
        push!(out, "standard errors are a parametric bootstrap of the count MPLE " *
                   "(Gibbs-simulate at θ̂, refit, empirical covariance): they do not " *
                   "assume the dyad conditionals are independent, but they are " *
                   "Monte-Carlo estimates and inherit the simulation's dependence on " *
                   "the fitted model being right")
    end
    return out
end

# StatsAPI interface: methods on the shared statistics generics, so results
# interoperate with StatsBase/GLM-style tooling (`coef(fit)`, `vcov(fit)`, ...)

# Number of free dyads (pseudo-likelihood contributions)
function _n_dyads(model::CountERGMModel)
    n = Int(nv(model.network))
    return model.directed ? n * (n - 1) : n * (n - 1) ÷ 2
end

StatsAPI.coef(result::CountERGMResult) = result.coefficients
StatsAPI.stderror(result::CountERGMResult) = result.std_errors
StatsAPI.vcov(result::CountERGMResult) = result.vcov
StatsAPI.loglikelihood(result::CountERGMResult) = result.loglik
StatsAPI.nobs(result::CountERGMResult) = _n_dyads(result.model)
StatsAPI.dof(result::CountERGMResult) = length(result.coefficients)

"""
    fit_ergm_count(net::Network, terms; reference=PoissonReference(), kwargs...)

Fit an ERGM for count-valued networks by maximum pseudo-likelihood: each
dyad's conditional distribution over the count support
`P(y_ij = y | rest) ∝ h(y)·exp(θ'Δg(y))` is used as an independent
likelihood contribution. The reference measure `h` enters the estimator
directly.

[`ergm_count`](@ref) is the R-faithful alias (matching the `ergm.count`
package); `fit_count_ergm` is a legacy alias.

# Arguments
- `net`: Network with `:weight` edge attribute (unweighted edges count as 1)
- `terms`: Vector of ERGM terms
- `reference`: Reference measure (default: Poisson)
- `max_val::Int`: Truncation of unbounded supports (default: twice the
  maximum observed count, at least 10)
- `se`, `n_boot`, `boot_burnin`, `boot_interval`, `rng`: standard-error controls,
  forwarded to [`count_mple`](@ref). `se=:bootstrap` replaces the
  inverse-pseudo-Hessian covariance with a parametric bootstrap (same API as
  `ERGM.mple`); the point estimates are unchanged.

# Returns
- `CountERGMResult`: Fitted model results
"""
function fit_ergm_count(net::Network{T}, terms::Vector{<:AbstractERGMTerm};
                        reference::AbstractReferenceMeasure=PoissonReference(),
                        method::Symbol=:mple,
                        maxiter::Int=100,
                        tol::Float64=1e-8,
                        max_val::Union{Int,Nothing}=nothing,
                        kwargs...) where T

    # Count MPLE enumerates every dyad as observed, so a masked (unobserved)
    # dyad would enter the pseudo-likelihood at its face value. Reject it.
    require_observed(net; context="fit_ergm_count", face_ok=false)

    model = CountERGMModel{T}(collect(AbstractERGMTerm, terms), net,
                              reference, is_directed(net))

    if method == :mple
        return count_mple(model; maxiter=maxiter, tol=tol, max_val=max_val,
                          kwargs...)
    else
        throw(ArgumentError("Unknown method: $method"))
    end
end

"""
    ergm_count(net::Network, terms; kwargs...)

R-faithful alias for [`fit_ergm_count`](@ref) (the same function), matching
the R `ergm.count` package name.
"""
const ergm_count = fit_ergm_count

"""
    fit_count_ergm(net::Network, terms; kwargs...)

Alias for [`fit_ergm_count`](@ref), kept for backward compatibility.
"""
const fit_count_ergm = fit_ergm_count

# All dyads of the network as (i, j) pairs
function _dyads(net)
    n = nv(net)
    if is_directed(net)
        return [(i, j) for i in 1:n for j in 1:n if i != j]
    else
        return [(i, j) for i in 1:n for j in (i+1):n]
    end
end

function _default_max_val(net, weights)
    m = 0
    for e in edges(net)
        m = max(m, Int(get(weights, _wkey(net, src(e), dst(e)), 1)))
    end
    return max(10, 2 * m)
end

"""
    count_mple(model::CountERGMModel; se=:hessian, n_boot=100, boot_burnin=100,
               boot_interval=10, rng=Random.default_rng(), kwargs...)
        -> CountERGMResult

Maximum pseudo-likelihood estimation for count ERGMs. For each dyad the
full conditional over the count support is enumerated, so the score is
`Σ_dyads [Δg(y_obs) − E_θ(Δg)]` and the Hessian is `−Σ_dyads Var_θ(Δg)`.
The pseudo-log-likelihood is maximized with the shared
`ERGM.newton_fit` Newton–Raphson-with-step-halving optimizer (documented
in ERGM.jl).

# Standard errors

- `se=:hessian` (default) — the inverse negative pseudo-Hessian. **Caution:**
  for a model with dyad-dependent terms (`CountMutualTerm`,
  `TransitiveTiesTerm`, the node-strength terms, ...) the pseudo-likelihood
  multiplies dyad conditionals as if independent, so these standard errors are
  expected to be *anticonservative* (too small) and the p-values too optimistic.
  For a dyad-independent model (e.g. `SumTerm` alone) the pseudo-likelihood is
  the likelihood and they are correct.
- `se=:bootstrap` — parametric bootstrap: Gibbs-simulate `n_boot` count networks
  from the fitted model at θ̂ with [`simulate_count_ergm`](@ref), refit the count
  MPLE on each, and report the empirical covariance of the refits. The point
  estimate is unchanged; only the covariance is replaced. This is the same
  option, with the same keywords and the same semantics, as `ERGM.mple`'s, and
  it runs on the ONE shared `Networks.bootstrap_cov` loop.

# Keyword Arguments
- `se::Symbol=:hessian`: `:hessian` or `:bootstrap` (above)
- `n_boot::Int=100`: number of bootstrap replicates (`se=:bootstrap` only)
- `boot_burnin::Int=100`, `boot_interval::Int=10`: Gibbs controls for the
  bootstrap simulations
- `rng::AbstractRNG=Random.default_rng()`: source of the bootstrap randomness —
  a fixed `rng` reproduces the standard errors exactly
- `maxiter::Int=100`, `tol::Float64=1e-8`, `max_val`: as before; the refits
  reuse the observed fit's support, so every replicate is fit on the same model.
"""
function count_mple(model::CountERGMModel{T}; maxiter::Int=100,
                    tol::Float64=1e-8,
                    max_val::Union{Int,Nothing}=nothing,
                    se::Symbol=:hessian,
                    n_boot::Int=100,
                    boot_burnin::Int=100,
                    boot_interval::Int=10,
                    rng::Random.AbstractRNG=Random.default_rng()) where T
    se in (:hessian, :bootstrap) ||
        throw(ArgumentError("se must be :hessian or :bootstrap, got :$se"))

    net = model.network
    weights = get_edge_attribute(net, :weight, Int)
    mv = isnothing(max_val) ? _default_max_val(net, weights) : max_val

    core = _count_mple_fit(model, mv; maxiter=maxiter, tol=tol)
    fit = core.fit
    support = core.support
    ref = model.reference

    # Boundary-mass diagnostic. For an unbounded reference the enumeration
    # `0:max_val` is a TRUNCATION of the model, not the model itself, so the
    # fit is only an approximation to the documented unbounded family insofar
    # as the fitted conditionals put negligible mass at the top of the support.
    # Measure that directly, at the fitted coefficients, and say so out loud.
    truncated = is_truncating(ref)
    boundary = truncated ?
        _max_boundary_mass(fit.θ, core.X, core.log_h, support, core.n_dyads) : 0.0

    if truncated && boundary > BOUNDARY_MASS_TOL
        @warn """
              Count support was truncated at max_val = $(last(support)), but the \
              fitted model places $(round(100 * boundary, sigdigits=3))% of the \
              conditional mass on the boundary value for at least one dyad \
              (tolerance $(round(100 * BOUNDARY_MASS_TOL, sigdigits=2))%).

              $(typeof(ref)) is mathematically UNBOUNDED. With appreciable mass at \
              the bound, this fit does not approximate that model — it is a \
              different, truncated exponential family, and the estimates are \
              biased toward the bound.

              Refit with a larger `max_val` and check the estimates are stable, \
              e.g. `count_mple(model; max_val = $(2 * last(support)))`.
              """ maxlog = 1
    end

    vcov, std_errors = fit.vcov, fit.se
    if se === :bootstrap
        vcov, std_errors = _count_bootstrap_cov(model, fit.θ, mv;
                                                n_boot=n_boot,
                                                boot_burnin=boot_burnin,
                                                boot_interval=boot_interval,
                                                maxiter=maxiter, tol=tol, rng=rng)
    end

    return CountERGMResult{T}(model, fit.θ, std_errors, vcov, fit.loglik,
                              fit.converged, last(support), truncated, boundary,
                              se)
end

# Parametric-bootstrap covariance of the count MPLE: Gibbs-simulate `n_boot`
# count networks at θ̂, refit the count MPLE on each, take the empirical
# covariance. The loop is the shared `Networks.bootstrap_cov`; this supplies only
# the two callbacks that are ERGMCount's. The refits reuse the observed fit's
# `max_val`, so every replicate is fit on the same (possibly truncated) support —
# and the simulator draws from `0:max_val`, so no replicate can fall outside it.
function _count_bootstrap_cov(model::CountERGMModel{T}, θ̂::Vector{Float64},
                              mv::Int; n_boot::Int, boot_burnin::Int,
                              boot_interval::Int, maxiter::Int, tol::Float64,
                              rng::Random.AbstractRNG) where T
    simulate(rng, B) = _simulate_count(model.network, model.terms,
                                       model.reference, θ̂; n_sim=B,
                                       burnin=boot_burnin, interval=boot_interval,
                                       max_val=mv, rng=rng)

    function refit(sim::Network{T})
        boot_model = CountERGMModel{T}(model.terms, sim, model.reference,
                                       model.directed)
        return _count_mple_fit(boot_model, mv; maxiter=maxiter, tol=tol).fit.θ
    end

    boot = bootstrap_cov(refit, simulate, θ̂; n_boot=n_boot, rng=rng)
    return boot.vcov, boot.se
end

# Core count MPLE: build the per-dyad change-statistic array over the support,
# maximize the pseudo-log-likelihood, and return the Newton fit alongside the
# pieces the boundary-mass diagnostic needs. Shared by `count_mple` and by the
# parametric bootstrap's refits (which need only `.fit.θ`).
function _count_mple_fit(model::CountERGMModel{T}, mv::Int;
                         maxiter::Int=100, tol::Float64=1e-8) where T
    net = model.network
    terms = model.terms
    ref = model.reference
    n_terms = length(terms)
    # Typed snapshot of the :weight attribute (see _simulate_count)
    weights = get_edge_attribute(net, :weight, Int)

    support = _support(ref, mv)
    y0 = first(support)

    dyads = _dyads(net)
    n_dyads = length(dyads)

    # Precompute, per dyad: observed value index and the change-stat matrix
    # Δg(y0 → y) for every support value (terms × support)
    log_h = [log_reference(ref, y) for y in support]
    X = Array{Float64}(undef, n_terms, length(support), n_dyads)
    y_obs_idx = Vector{Int}(undef, n_dyads)

    for (d, (i, j)) in enumerate(dyads)
        y_obs = dyad_value(net, weights, i, j)
        if y_obs > last(support) || y_obs < y0
            throw(ArgumentError(
                "Observed count $y_obs at dyad ($i,$j) outside the support " *
                "$(support); increase max_val"))
        end
        y_obs_idx[d] = y_obs - y0 + 1
        for (s, y) in enumerate(support)
            for (k, term) in enumerate(terms)
                X[k, s, d] = change_stat_count(term, net, weights, i, j, y0, y)
            end
        end
    end

    # Maximize with ERGM.jl's shared Newton-with-step-halving optimizer
    derivatives = _count_derivatives(X, log_h, support, y_obs_idx, n_dyads, n_terms)
    fit = newton_fit(derivatives, zeros(n_terms); maxiter=maxiter, tol=tol)

    return (fit=fit, X=X, log_h=log_h, support=support, n_dyads=n_dyads)
end

# The `(ll, grad, hess)` closure of the count pseudo-log-likelihood: each dyad
# contributes its full conditional P(y_ij = y | rest) ∝ h(y)·exp(θ'Δg(y)) over
# the (truncated) support, so the gradient is Δg(y_obs) − E[Δg] and the Hessian
# is −(E[Δg Δg'] − E[Δg]E[Δg]'), summed over dyads.
#
# The workspaces are allocated ONCE (review finding 15). The conditional moments
# used to be rebuilt per dyad — `zeros(n_terms)`, `zeros(n_terms, n_terms)`, and
# a fresh `x * x'` outer product for every support value — i.e. n_dyads ×
# (2 + |support|) allocations on every Newton evaluation, 2.1 MB per evaluation
# on zach at max_val = 30. They are now filled in place with the SAME scalar
# arithmetic (each entry is still `p * (x[k] * x[l])`, accumulated in the same
# order), so the fit is bit-for-bit what it was — the golden fixture is an exact
# MLE and this refactor must not move it.
function _count_derivatives(X::Array{Float64,3}, log_h::Vector{Float64}, support,
                            y_obs_idx::Vector{Int}, n_dyads::Int, n_terms::Int)
    η = Vector{Float64}(undef, length(support))
    eX = Vector{Float64}(undef, n_terms)
    eXX = Matrix{Float64}(undef, n_terms, n_terms)

    return function (β)
        llv = 0.0
        grad = zeros(n_terms)
        hess = zeros(n_terms, n_terms)
        @inbounds for d in 1:n_dyads
            Xd = @view X[:, :, d]
            for s in eachindex(support)
                η[s] = log_h[s] + dot(β, @view Xd[:, s])
            end
            ηmax = maximum(η)
            Z = 0.0
            for s in eachindex(support)
                Z += exp(η[s] - ηmax)
            end
            logZ = ηmax + log(Z)
            llv += η[y_obs_idx[d]] - logZ

            # E[Δg] and E[Δg Δg'] under the conditional
            fill!(eX, 0.0)
            fill!(eXX, 0.0)
            for s in eachindex(support)
                p = exp(η[s] - logZ)
                for l in 1:n_terms
                    xl = Xd[l, s]
                    eX[l] += p * xl
                    for k in 1:n_terms
                        eXX[k, l] += p * (Xd[k, s] * xl)
                    end
                end
            end
            for l in 1:n_terms
                grad[l] += Xd[l, y_obs_idx[d]] - eX[l]
                for k in 1:n_terms
                    hess[k, l] -= eXX[k, l] - eX[k] * eX[l]
                end
            end
        end
        return llv, grad, hess
    end
end

# Largest conditional probability placed on the top support value by any dyad,
# evaluated at the fitted coefficients. Mirrors the `derivatives` conditional
# exactly (same η, same log-sum-exp) so the diagnostic cannot drift from the
# likelihood it is diagnosing.
function _max_boundary_mass(β::Vector{Float64}, X::Array{Float64,3},
                            log_h::Vector{Float64}, support, n_dyads::Int)
    worst = 0.0
    η = Vector{Float64}(undef, length(support))
    top = length(support)
    for d in 1:n_dyads
        Xd = @view X[:, :, d]
        for s in eachindex(support)
            η[s] = log_h[s] + dot(β, @view Xd[:, s])
        end
        ηmax = maximum(η)
        Z = 0.0
        for s in eachindex(support)
            Z += exp(η[s] - ηmax)
        end
        worst = max(worst, exp(η[top] - ηmax) / Z)
    end
    return worst
end

# =============================================================================
# Simulation
# =============================================================================

"""
    simulate_count_ergm(result::CountERGMResult; n_sim=1, burnin=100,
                        interval=10, max_val=nothing,
                        rng=Random.default_rng()) -> Vector{Network}

Simulate networks from a fitted count ERGM by Gibbs sampling: each dyad is
resampled from its full conditional
`P(y_ij = y | rest) ∝ h(y)·exp(θ'Δg(y))`, using each term's proper change
statistic, so structural terms (mutuality, transitivity, node strength)
influence the draws.

All random draws flow through `rng`; the same rng state yields identical
output.
"""
function simulate_count_ergm(result::CountERGMResult{T};
                             n_sim::Int=1,
                             burnin::Int=100,
                             interval::Int=10,
                             max_val::Union{Int,Nothing}=nothing,
                             rng::Random.AbstractRNG=Random.default_rng()) where T
    model = result.model
    mv = isnothing(max_val) ? result.max_val : max_val
    return _simulate_count(model.network, model.terms, model.reference,
                           result.coefficients; n_sim=n_sim, burnin=burnin,
                           interval=interval, max_val=mv, rng=rng)
end

"""
    simulate_count_ergm(net, terms, coefficients; reference=PoissonReference(),
                        rng=Random.default_rng(), kwargs...)

Simulate count networks from an ERGM specification (without fitting).
`net` provides the size, directedness, and starting state.
"""
function simulate_count_ergm(net::Network{T}, terms::Vector{<:AbstractERGMTerm},
                             coefficients::Vector{Float64};
                             reference::AbstractReferenceMeasure=PoissonReference(),
                             n_sim::Int=1,
                             burnin::Int=100,
                             interval::Int=10,
                             max_val::Int=20,
                             rng::Random.AbstractRNG=Random.default_rng()) where T
    return _simulate_count(net, collect(AbstractERGMTerm, terms), reference,
                           coefficients; n_sim=n_sim, burnin=burnin,
                           interval=interval, max_val=max_val, rng=rng)
end

function _simulate_count(net0::Network{T}, terms::Vector{AbstractERGMTerm},
                         ref::AbstractReferenceMeasure,
                         coefficients::Vector{Float64};
                         n_sim::Int, burnin::Int, interval::Int,
                         max_val::Int,
                         rng::Random.AbstractRNG=Random.default_rng()) where T
    networks = Network{T}[]

    current = deepcopy(net0)
    n = nv(current)
    directed = is_directed(current)
    support = _support(ref, max_val)
    log_h = [log_reference(ref, y) for y in support]
    log_probs = Vector{Float64}(undef, length(support))
    probs = Vector{Float64}(undef, length(support))

    # Typed snapshot of the :weight edge attribute. The attribute storage is
    # an untyped Dict{Tuple{T,T},Any}, which makes every dyad_value read in
    # the innermost conditional loop type-unstable; snapshotting once into a
    # Dict{Tuple{T,T},Int} (Networks.jl's typed accessor) and maintaining it
    # incrementally alongside the network keeps the hot loop concretely
    # typed.
    weights = get_edge_attribute(current, :weight, Int)

    # One Gibbs sweep updates every dyad once
    for sweep in 1:(burnin + n_sim * interval)
        for i in 1:n
            for j in (directed ? (1:n) : (i+1:n))
                i == j && continue

                old = dyad_value(current, weights, i, j)

                for (s, y) in enumerate(support)
                    lp = log_h[s]
                    for (k, term) in enumerate(terms)
                        lp += coefficients[k] *
                              change_stat_count(term, current, weights, i, j, old, y)
                    end
                    log_probs[s] = lp
                end

                log_probs .-= maximum(log_probs)
                probs .= exp.(log_probs)
                probs ./= sum(probs)

                new_val = support[sample(rng, eachindex(support), Weights(probs))]

                # Update the network, its :weight attribute, and the typed
                # snapshot in place
                if new_val > 0
                    if !has_edge(current, i, j)
                        add_edge!(current, i, j)
                    end
                    set_edge_attribute!(current, :weight, i, j, new_val)
                    weights[_wkey(current, i, j)] = new_val
                elseif old > 0
                    rem_edge!(current, i, j)  # also drops the edge's attributes
                    delete!(weights, _wkey(current, i, j))
                end
            end
        end

        if sweep > burnin && (sweep - burnin) % interval == 0
            push!(networks, deepcopy(current))
        end
    end

    return networks
end

# =============================================================================
# Goodness of Fit
# =============================================================================

# Maximum dyad count value in a network (0 for an empty network)
function _max_dyad_value(net)
    weights = get_edge_attribute(net, :weight, Int)
    m = 0
    for e in edges(net)
        m = max(m, Int(get(weights, _wkey(net, src(e), dst(e)), 1)))
    end
    return m
end

# Number of dyads with count value k, for k in 0:K (zeros are the dyads
# without an edge)
function _dyad_value_counts(net, K::Int)
    weights = get_edge_attribute(net, :weight, Int)
    counts = zeros(Float64, K + 1)
    for e in edges(net)
        w = clamp(Int(get(weights, _wkey(net, src(e), dst(e)), 1)), 0, K)
        counts[w + 1] += 1
    end
    n = Int(nv(net))
    n_dyads = is_directed(net) ? n * (n - 1) : n * (n - 1) ÷ 2
    counts[1] += n_dyads - ne(net)
    return counts
end

"""
    gof(result::CountERGMResult; n_sim=100, burnin=100, interval=10,
        max_val=nothing, rng=Random.default_rng()) -> GOFResult

Goodness-of-fit assessment of a fitted count ERGM: networks are simulated
from the fitted model with [`simulate_count_ergm`](@ref) and compared with
the observed network on

- the model statistics (one level per term), and
- the distribution of dyad count values (number of dyads with value
  0, 1, 2, ...).

This is a method of the shared `Networks.gof` generic; it returns the
shared `Networks.GOFResult` (observed value, simulation envelope, and
two-sided Monte-Carlo p-value per level).

# Keyword Arguments
- `n_sim::Int=100`: Number of simulated networks
- `burnin`, `interval`, `max_val`, `rng`: passed to
  [`simulate_count_ergm`](@ref)
"""
function gof(result::CountERGMResult; n_sim::Int=100, burnin::Int=100,
             interval::Int=10, max_val::Union{Int,Nothing}=nothing,
             rng::Random.AbstractRNG=Random.default_rng())
    net = result.model.network
    terms = result.model.terms
    sims = simulate_count_ergm(result; n_sim=n_sim, burnin=burnin,
                               interval=interval, max_val=max_val, rng=rng)

    # Model statistics: observed vs simulated
    obs_stats = [compute(term, net) for term in terms]
    sim_stats = [compute(term, s) for s in sims, term in terms]
    stats = GOFStatistic("model statistics", [name(term) for term in terms],
                         obs_stats, sim_stats)

    # Dyad count-value distribution
    K = maximum(_max_dyad_value(s) for s in sims; init=_max_dyad_value(net))
    obs_counts = _dyad_value_counts(net, K)
    sim_counts = Matrix{Float64}(undef, n_sim, K + 1)
    for (r, s) in enumerate(sims)
        sim_counts[r, :] .= _dyad_value_counts(s, K)
    end
    values = GOFStatistic("dyad count values", string.(0:K),
                          obs_counts, sim_counts)

    return GOFResult([stats, values]; model="Count ERGM")
end

end # module
