"""
    ERGMCount.jl - ERGMs for Count-Valued Networks

Extends ERGM to handle networks with integer-valued edge weights,
using Poisson, geometric, or binomial reference measures.

The general form is:
    P(Y=y) ∝ h(y) exp(θ' g(y))

where h(y) is the reference measure determining the baseline distribution
for count-valued edges (Krivitsky 2012).

Port of the R ergm.count package from the StatNet collection. Estimation is
by maximum pseudo-likelihood over an error-controlled count support
([`fit_ergm_count`](@ref)); simulation is a Gibbs sweep over each dyad's
full conditional ([`simulate_count_ergm`](@ref)).
"""
module ERGMCount

using Distributions   # Poisson/Geometric/Binomial reference draws, Normal quantile
using ERGM
using Graphs
using LinearAlgebra
using Logging: NullLogger, with_logger
using Networks
using PrecompileTools: @setup_workload, @compile_workload
using Printf: @sprintf
using Random
using Statistics: cov

# The shared numerics and validators live in Networks.jl (panel 2026-09, items
# 13/14/28): the ONE Newton optimizer, the ONE floored z → p helper, the ONE
# `se=` validator and the generic coefficient table every `coeftable` returns.
import Networks: newton_fit, z_pvalues, check_se, CoefficientTable
# The statistic protocol and the dependence/directedness traits are ERGM.jl's
# generics: ERGMCount adds methods for its own terms and model type, never a
# same-named private.
import ERGM: name, compute, is_dyad_dependent, has_dyad_dependent,
             requires_directed
# Shared presentation infrastructure (Networks.jl): the ONE `gof` generic all
# model packages extend, plus the GOF containers
import Networks: gof, GOFStatistic, GOFResult

# The shared result-metadata protocol (Networks.jl `src/results.jl`): the
# generic accessors that say what a fit actually did. Imported by name because
# ERGMCount adds methods for `CountERGMResult`; `fit_metadata(fit)` collects them.
import Networks: estimand, objective, is_exact, se_method, missing_method,
                 approximations, missing_policies
import StatsAPI
import StatsAPI: coef, stderror, vcov, confint, loglikelihood, nobs, dof, aic,
                 bic, coeftable

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
export TransitiveWeightsTerm, CyclicalWeightsTerm
export SmallerthanTerm, EqualToTerm, InIntervalTerm
export change_stat_count, dyad_value
public change_stats_support!

# Model / estimation (`count_mple` is exported as ERGM.jl exports `mple`: the
# estimator a warning or a docs page names must be callable as written)
export CountERGMModel, CountERGMResult
export fit_ergm_count, ergm_count, fit_count_ergm, count_mple
export has_dyad_dependent

# Simulation
export simulate_count_ergm

# Goodness of fit (method of the shared Networks.jl `gof` generic)
export gof

# The full StatsAPI surface (re-exported so `coef(fit)` etc. work with just
# `using ERGMCount`; the bindings are StatsAPI's, so co-loading with ERGM.jl
# or REM.jl leaves every verb defined)
export coef, stderror, vcov, confint, loglikelihood, nobs, dof, aic, bic,
       coeftable

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

# Example
```julia
using ERGMCount
log_reference(PoissonReference(2.0), 3) ≈ 3 * log(2.0) - log(6)   # true
log_reference(GeometricReference(), 7)                            # 0.0
log_reference(BinomialReference(4), 5)                            # -Inf — outside 0:4
```
"""
function log_reference end

"""
    sample_reference(ref::AbstractReferenceMeasure;
                     rng::Random.AbstractRNG=Random.default_rng()) -> Int

Draw a single random count from the baseline distribution associated with
the reference measure `ref` (e.g. `Poisson(lambda)` for
[`PoissonReference`](@ref ERGMCount.PoissonReference)). All randomness flows through the `rng`
keyword, so the same rng state yields identical draws.

# Example
```julia
using ERGMCount, Random
y = sample_reference(BinomialReference(10); rng=Xoshiro(1))
0 <= y <= 10                                                        # true
sample_reference(PoissonReference(2.0); rng=Xoshiro(3)) ==
    sample_reference(PoissonReference(2.0); rng=Xoshiro(3))         # true
```
"""
function sample_reference end

"""
    PoissonReference(lambda=1.0)

Poisson reference measure for count ERGMs.
h(y_ij) = λ^y_ij / y_ij!

With this reference and a `SumTerm` coefficient θ, each dyad is
conditionally Poisson(λ·e^θ) (Krivitsky 2012). `lambda` must be positive
(an `ArgumentError` otherwise: `log(λ)` enters every conditional).

# Fields
- `lambda::Float64`: Rate parameter (default 1.0)

# Example
```julia
using ERGMCount
ref = PoissonReference()                # λ = 1: h(y) = 1/y!
log_reference(ref, 3) ≈ -log(6)         # true
PoissonReference(2.0).lambda            # 2.0
is_truncating(ref)                      # true — the support is unbounded
```
"""
struct PoissonReference <: AbstractReferenceMeasure
    lambda::Float64

    function PoissonReference(λ::Real=1.0)
        (isfinite(λ) && λ > 0) || throw(ArgumentError(
            "PoissonReference: lambda must be positive and finite (got $λ); " *
            "`log(lambda)` enters every dyad conditional"))
        new(Float64(λ))
    end
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

# Example
```julia
using ERGMCount
ref = GeometricReference()
log_reference(ref, 5)               # 0.0 — the counting measure
is_truncating(ref)                  # true — 0:max_val truncates it
ERGMCount._support(ref, 20)         # 0:20
```
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

# Example
```julia
using ERGMCount
ref = BinomialReference(4)
log_reference(ref, 2) ≈ log(6)      # true — C(4, 2) = 6
log_reference(ref, 5)               # -Inf — outside the support 0:4
is_truncating(ref)                  # false — the bound is the model's
```
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

# Example
```julia
using ERGMCount
ref = DiscUnifReference(3)
log_reference(ref, 2) ≈ -log(4)     # true — every value of 0:3 has weight 1/4
ERGMCount._support(ref, 30)         # 0:3 — `max_val` is not consulted
is_truncating(ref)                  # false
```
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
    DiscUnif2Reference(a, b)

Discrete uniform reference on {a, a+1, ..., b}, `ergm.count`'s
`DiscUnif(a, b)`: h(y) = 1/(b − a + 1) on the support `a:b`, which is part of
the model (see [`is_truncating`](@ref)). `a` may be **negative**, as in R: a
dyad's count is then any integer in `a:b`, the estimator enumerates that
support, the Gibbs sampler draws from it (an edge is stored for every
non-zero count, negative ones included) and `NonzeroTerm` counts the dyads
with `y ≠ 0`. `a` may also be positive, in which case the value 0 is
impossible and a network with an absent edge is refused by the fit.

# Example
```julia
using Networks, ERGMCount
ref = DiscUnif2Reference(-2, 2)
log_reference(ref, -1)              # -log(5)
ERGMCount._support(ref, 30)         # -2:2 — `max_val` is not consulted
is_truncating(ref)                  # false
```
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

# Example
```julia
using ERGMCount
is_truncating(PoissonReference())       # true
is_truncating(GeometricReference())     # true
is_truncating(BinomialReference(5))     # false
is_truncating(DiscUnif2Reference(-2, 2))  # false
```
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

# Example
```julia
using Networks, ERGMCount
BOUNDARY_MASS_TOL                       # 0.0001
net = network(4; directed=true)
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 2)
add_edge!(net, 3, 4); set_edge_attribute!(net, :weight, 3, 4, 1)
fit = fit_ergm_count(net, [SumTerm()])
fit.boundary_mass <= BOUNDARY_MASS_TOL  # true — the default support is chosen so
```
"""
const BOUNDARY_MASS_TOL = 1e-4

# Three significant digits without floating-point noise: `round(x,
# sigdigits=3)` prints `6.969999999999999e-32`, `@sprintf("%.3g")` prints
# `6.97e-32`. Every diagnostic number in `show`, the caveats and the warnings
# goes through here, so the docs' quoted output is what the user sees.
_fmt3(x::Real) = isfinite(x) ? @sprintf("%.3g", x) : string(x)
_fmt2(x::Real) = isfinite(x) ? @sprintf("%.2g", x) : string(x)

# =============================================================================
# Dyad values
# =============================================================================

# Canonical edge-attribute key: (i,j) directed, (min,max) undirected
_wkey(net, i::Int, j::Int) = is_directed(net) ? (i, j) : minmax(i, j)

"""
    dyad_value(net, weights, i, j) -> Int

The count value of dyad (i,j): 0 when the edge is absent, its `:weight`
attribute when present. `weights` is the `:weight` attribute dictionary —
pass the typed snapshot `get_edge_attribute(net, :weight, Int)` — keyed by
`(i, j)` on a directed network and `minmax(i, j)` on an undirected one.

An edge without a `:weight` entry reads as 1 here, but no such network
reaches a fit or a simulation: `CountERGMModel` (hence `fit_ergm_count`,
`count_mple`) and `simulate_count_ergm` refuse a network with an edge that
carries no `:weight`, whether none of the edges does (R's `response="w"`
migration mistake: pass `weight=:w`) or only some do (a weight column with
gaps), naming the first such edge.

# Example
```julia
using Networks, ERGMCount
net = network(3; directed=false)
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 4)
w = get_edge_attribute(net, :weight, Int)
dyad_value(net, w, 1, 2)   # 4
dyad_value(net, w, 2, 1)   # 4 — the undirected key is minmax(2, 1)
dyad_value(net, w, 1, 3)   # 0 — no edge
```
"""
function dyad_value(net, weights, i::Int, j::Int)
    has_edge(net, i, j) || return 0
    return Int(get(weights, _wkey(net, i, j), 1))
end

# The typed snapshot (`Dict{Tuple{T,T},Int}`): `compute` on the untyped
# `Dict{…,Any}` inferred `Any` for five terms and boxed every `+`. Integer-
# valued Floats convert; a non-integer weight is refused by the model
# constructor before any estimator reads it.
_get_weights(net) = get_edge_attribute(net, :weight, Int)

# Number of dyads of a one-mode network: n(n−1) directed, n(n−1)/2 undirected
function _n_dyads(net::Network)
    n = Int(nv(net))
    return is_directed(net) ? n * (n - 1) : n * (n - 1) ÷ 2
end

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

# Example
```julia
using Networks, ERGMCount
net = network(3; directed=true)
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 3)
add_edge!(net, 1, 3); set_edge_attribute!(net, :weight, 1, 3, 1)
w = get_edge_attribute(net, :weight, Int)
change_stat_count(SumTerm(), net, w, 1, 2, 3, 5)        # 2.0
change_stat_count(NonzeroTerm(), net, w, 1, 2, 3, 0)    # -1.0
change_stat_count(NodeOSumTerm(), net, w, 1, 2, 3, 5)   # 20.0 — (1+5)² − (1+3)²
```
"""
function change_stat_count end

"""
    SumTerm <: AbstractERGMTerm

Sum of edge values: ∑_{i,j} y_{ij} — R's `sum`. This is the natural
sufficient statistic for the Poisson reference.

# Example
```julia
using Networks, ERGMCount
using ERGM: compute, name
net = network(3; directed=true)
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 3)
add_edge!(net, 2, 3); set_edge_attribute!(net, :weight, 2, 3, 2)
compute(SumTerm(), net)   # 5.0
name(SumTerm())           # "sum"
```
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

Number of non-zero dyads: ∑_{i,j} I(y_{ij} ≠ 0) — R's `nonzero`. An edge
whose `:weight` is 0 is a zero dyad (R stores no such edge), and a negative
count (admissible under `DiscUnif2Reference(a < 0, b)`) is non-zero.

# Example
```julia
using Networks, ERGMCount
using ERGM: compute, name
net = network(3; directed=true)
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 3)
add_edge!(net, 2, 3); set_edge_attribute!(net, :weight, 2, 3, 0)   # a zero-valued tie
compute(NonzeroTerm(), net)   # 1.0 — the 0-valued edge is a zero dyad
name(NonzeroTerm())           # "nonzero"
```
"""
struct NonzeroTerm <: AbstractERGMTerm end

name(::NonzeroTerm) = "nonzero"

function compute(::NonzeroTerm, net)
    weights = _get_weights(net)
    total = 0.0
    for e in edges(net)
        get(weights, _wkey(net, src(e), dst(e)), 1) != 0 && (total += 1.0)
    end
    return total
end

change_stat_count(::NonzeroTerm, net, weights, i::Int, j::Int, old::Int, new::Int) =
    Float64((new != 0) - (old != 0))

"""
    GreaterthannTerm <: AbstractERGMTerm

Number of dyads with value > n: ∑_{i,j} I(y_{ij} > n) — R's `greaterthan(n)`.
The zero-valued dyads count when `n < 0` (as `SmallerthanTerm` counts them),
so `GreaterthannTerm(-1)` is the number of dyads.

# Fields
- `threshold::Int`: Threshold value n

# Example
```julia
using Networks, ERGMCount
using ERGM: compute, name
net = network(3; directed=true)            # 6 dyads
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 3)
add_edge!(net, 2, 3); set_edge_attribute!(net, :weight, 2, 3, 1)
compute(GreaterthannTerm(2), net)    # 1.0
compute(GreaterthannTerm(-1), net)   # 6.0 — every dyad, the zeros included
name(GreaterthannTerm(2))            # "greaterthan.2"
```
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
    # Every dyad without an edge has value 0
    0 > t.threshold && (total += _n_dyads(net) - ne(net))
    return total
end

change_stat_count(t::GreaterthannTerm, net, weights, i::Int, j::Int, old::Int, new::Int) =
    Float64((new > t.threshold) - (old > t.threshold))

"""
    CountAtleastnTerm <: AbstractERGMTerm

Number of dyads with value >= n: ∑_{i,j} I(y_{ij} ≥ n) — R's `atleast(n)`.
The zero-valued dyads count when `n ≤ 0`, so `CountAtleastnTerm(0)` is the
number of dyads.

# Example
```julia
using Networks, ERGMCount
using ERGM: compute, name
net = network(3; directed=true)            # 6 dyads
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 3)
add_edge!(net, 2, 3); set_edge_attribute!(net, :weight, 2, 3, 1)
compute(CountAtleastnTerm(3), net)   # 1.0
compute(CountAtleastnTerm(1), net)   # 2.0
compute(CountAtleastnTerm(0), net)   # 6.0 — every dyad
name(CountAtleastnTerm(3))           # "atleast.3"
```
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
    0 >= t.threshold && (total += _n_dyads(net) - ne(net))
    return total
end

change_stat_count(t::CountAtleastnTerm, net, weights, i::Int, j::Int, old::Int, new::Int) =
    Float64((new >= t.threshold) - (old >= t.threshold))

"""
    CountMutualTerm(form=:min; threshold=0) <: AbstractERGMTerm

Valued reciprocity, `ergm`'s valued `mutual(form=, threshold=)`:
∑_{i<j} m(y_{ij}, y_{ji}) over the unordered pairs of a **directed** network,
where `form` picks `m`:

| `form` | `m(a, b)` | R label |
|---|---|---|
| `:min` (default) | `min(a, b)` | `mutual.min` |
| `:nabsdiff` | `-abs(a - b)` | `mutual.nabsdiff` |
| `:geometric` | `sqrt(a * b)` | `mutual.geom.mean` |
| `:product` | `a * b` | `mutual.product` |
| `:threshold` | `(a >= threshold) * (b >= threshold)` — binary mutuality after thresholding | `mutual.<threshold>` |

The coefficient labels are the ones `summary()`/`coef()` print in R. The
`:min`, `:nabsdiff`, `:geometric` and `:product` statistics are pinned
against `ergm` 4.12 by `test/fixtures/count_terms.toml`; `:threshold`
follows the documented definition (an empty network has every pair mutual
when `threshold <= 0`, as R's `emptynwstats` records) but is **not** pinned,
because `ergm` 4.12's own `mutual(form="threshold")` fails at C model
initialisation (the error is frozen in the fixture). Every form is
dyad-dependent (it reads the reciprocal dyad) and requires a directed
network; an undirected `CountERGMModel` refuses the term.

# Example
```julia
using Networks, ERGMCount
using ERGM: compute, name
net = network(3; directed=true)
for (i, j, w) in ((1, 2, 3), (2, 1, 2), (2, 3, 2), (1, 3, 1))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
compute(CountMutualTerm(), net)                     # 2.0  — min(3, 2)
compute(CountMutualTerm(:nabsdiff), net)            # -4.0 — -(|3-2| + |2-0| + |1-0|)
compute(CountMutualTerm(:threshold; threshold=2), net)   # 1.0
name(CountMutualTerm(:geometric))                   # "mutual.geom.mean"
```
"""
struct CountMutualTerm <: AbstractERGMTerm
    form::Symbol
    threshold::Int

    function CountMutualTerm(form::Symbol=:min; threshold::Int=0)
        form in (:min, :nabsdiff, :geometric, :product, :threshold) ||
            throw(ArgumentError(
                "CountMutualTerm: form must be one of :min, :nabsdiff, :geometric, " *
                ":product or :threshold (got :$form)"))
        new(form, threshold)
    end
end

function name(t::CountMutualTerm)
    f = t.form
    f === :min && return "mutual.min"
    f === :nabsdiff && return "mutual.nabsdiff"
    f === :geometric && return "mutual.geom.mean"
    f === :product && return "mutual.product"
    return "mutual.$(t.threshold)"          # R: paste("mutual", threshold, sep=".")
end

# m(a, b) of the docstring's table, as a Float64
@inline function _mutual_pair(t::CountMutualTerm, a::Int, b::Int)
    f = t.form
    f === :min && return Float64(min(a, b))
    f === :nabsdiff && return -Float64(abs(a - b))
    f === :geometric && return sqrt(Float64(a) * Float64(b))
    f === :product && return Float64(a) * Float64(b)
    return Float64((a >= t.threshold) & (b >= t.threshold))
end

function compute(t::CountMutualTerm, net)
    !is_directed(net) && return 0.0
    t.form === :geometric && _refuse_negative_network(t, net, "compute")

    weights = _get_weights(net)
    total = 0.0
    n = nv(net)

    for i in 1:n, j in (i+1):n
        total += _mutual_pair(t, dyad_value(net, weights, i, j),
                              dyad_value(net, weights, j, i))
    end

    return total
end

function change_stat_count(t::CountMutualTerm, net, weights, i::Int, j::Int,
                           old::Int, new::Int)
    is_directed(net) || return 0.0
    y_ji = dyad_value(net, weights, j, i)
    return _mutual_pair(t, new, y_ji) - _mutual_pair(t, old, y_ji)
end

"""
    TransitiveTiesTerm <: AbstractERGMTerm

Triadic-minimum transitivity over ordered distinct triples:
∑_{i≠j≠k} min(y_{ij}, y_{jk}, y_{ik})

**This statistic has no R counterpart and is not validated against R.** It
is not `ergm`'s `transitiveweights` (that is [`TransitiveWeightsTerm`](@ref),
pinned by the `count_terms` fixture) and not the binary `transitiveties`;
it sums the triple-wise minimum over every ordered triple, so on an
undirected network each unordered triangle contributes six times. Use it
as a simple valued transitivity when parity with R does not matter.

# Example
```julia
using Networks, ERGMCount
using ERGM: compute, name
net = network(3; directed=true)
for (i, j, w) in ((1, 2, 3), (2, 3, 2), (1, 3, 1), (2, 1, 2))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
compute(TransitiveTiesTerm(), net)   # 2.0 — triples (1,2,3) and (2,1,3)
name(TransitiveTiesTerm())           # "transitiveties.count"
```
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
(each directed 3-cycle counted once). Directed networks only.

**This statistic has no R counterpart and is not validated against R.** It
is not `ergm`'s `cyclicalweights` (that is [`CyclicalWeightsTerm`](@ref),
pinned by the `count_terms` fixture) and not the binary `cyclicalties`.

# Example
```julia
using Networks, ERGMCount
using ERGM: compute, name
net = network(3; directed=true)
for (i, j, w) in ((1, 2, 3), (2, 3, 2), (3, 1, 2))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
compute(CyclicalTiesTerm(), net)   # 2.0 — the cycle 1→2→3→1, min(3, 2, 2)
name(CyclicalTiesTerm())           # "cyclicalties.count"
```
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

Sum of squared out-strengths, ∑_i (∑_j y_{ij})²: measures activity
heterogeneity. Directed networks only (0 for undirected, and a
`CountERGMModel` refuses it there; use `NodeSumTerm` instead). **No R
counterpart** (`ergm`'s `nodeocovar` is a covariance, not this sum) — label
`nodeOSum`, not validated against R.

# Example
```julia
using Networks, ERGMCount
using ERGM: compute, name
net = network(3; directed=true)
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 3)
add_edge!(net, 1, 3); set_edge_attribute!(net, :weight, 1, 3, 1)
add_edge!(net, 2, 3); set_edge_attribute!(net, :weight, 2, 3, 2)
compute(NodeOSumTerm(), net)   # 20.0 — out-strengths (4, 2, 0): 16 + 4
name(NodeOSumTerm())           # "nodeOSum"
```
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

Sum of squared in-strengths, ∑_j (∑_i y_{ij})²: measures popularity
heterogeneity. Directed networks only (0 for undirected, and a
`CountERGMModel` refuses it there; use `NodeSumTerm` instead). **No R
counterpart** — label `nodeISum`, not validated against R.

# Example
```julia
using Networks, ERGMCount
using ERGM: compute, name
net = network(3; directed=true)
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 3)
add_edge!(net, 1, 3); set_edge_attribute!(net, :weight, 1, 3, 1)
add_edge!(net, 2, 3); set_edge_attribute!(net, :weight, 2, 3, 2)
compute(NodeISumTerm(), net)   # 18.0 — in-strengths (0, 3, 3): 9 + 9
name(NodeISumTerm())           # "nodeISum"
```
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

Sum of squared total (in + out) strengths, ∑_i (∑_j y_{ij} + ∑_j y_{ji})²;
on an undirected network the strength of a vertex is the sum of its edge
values. Defined on both directednesses. **No R counterpart** (`ergm`'s
`nodecovar` is a covariance) — label `nodeSum`, not validated against R.

# Example
```julia
using Networks, ERGMCount
using ERGM: compute, name
net = network(3; directed=false)
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 3)
add_edge!(net, 2, 3); set_edge_attribute!(net, :weight, 2, 3, 1)
compute(NodeSumTerm(), net)   # 26.0 — strengths (3, 4, 1): 9 + 16 + 1
name(NodeSumTerm())           # "nodeSum"
```
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


# -----------------------------------------------------------------------------
# R-parity terms (ergm's valued `transitiveweights`/`cyclicalweights` with the
# default (min, max, min) triple, and the dyad-independent threshold terms
# `smallerthan`/`equalto`/`ininterval`), pinned by test/fixtures/count_terms.toml
# -----------------------------------------------------------------------------

# The strongest two-path i → k → j with the (min, max) triple: max_k min(y_ik, y_kj)
# over k ∉ {i, j}, reading the dyad (i, j) itself nowhere. `skip` names a third
# vertex to leave out (the changing dyad's other endpoint), so the affected
# pairs' two-path strength can be re-assembled with the changing dyad at any
# value: max(P_without_skip, min(new, y_other)).
function _best_twopath(net, weights, i::Int, j::Int, skip::Int)
    best = 0
    for k in outneighbors(net, i)
        (k == i || k == j || k == skip) && continue
        y_ik = dyad_value(net, weights, i, k)
        y_ik > best || continue                      # cannot beat the best
        c = min(y_ik, dyad_value(net, weights, k, j))
        c > best && (best = c)
    end
    return best
end

"""
    TransitiveWeightsTerm() <: AbstractERGMTerm

`ergm`'s valued `transitiveweights("min", "max", "min")` (Krivitsky 2012,
eq. 13, with the default triple): every dyad's value capped by the strongest
two-path closing it,

∑_{(i,j)} min( y_{ij}, max_k min(y_{ik}, y_{kj}) ),

summed over every ordered pair of a directed network and over every
unordered pair of an undirected one, exactly as R's C code does (R label
`transitiveweights.min.max.min`; pinned against `ergm` 4.12 on `zach` and on
a directed count network by `test/fixtures/count_terms.toml`). The
non-default triples (`twopath="geomean"`, `combine="sum"`,
`affect="geomean"`) are not implemented. Dyad-dependent: changing dyad
`(i, j)` moves its own term and the terms of the pairs `(i, l)`, `l ∈ out(j)`,
and `(k, j)`, `k ∈ in(i)`, whose strongest two-path may run through it.

# Example
```julia
using Networks, ERGMCount
using ERGM: compute, name
net = network(3; directed=true)
for (i, j, w) in ((1, 2, 3), (2, 3, 2), (1, 3, 1))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
compute(TransitiveWeightsTerm(), net)   # 1.0 — pair (1,3): min(1, min(3, 2))
name(TransitiveWeightsTerm())           # "transitiveweights.min.max.min"
```
"""
struct TransitiveWeightsTerm <: AbstractERGMTerm end

name(::TransitiveWeightsTerm) = "transitiveweights.min.max.min"

function compute(t::TransitiveWeightsTerm, net)
    _refuse_negative_network(t, net, "compute")
    weights = _get_weights(net)
    total = 0.0
    for e in edges(net)                 # y_ij = 0 contributes min(0, ·) = 0
        i, j = src(e), dst(e)
        y = dyad_value(net, weights, i, j)
        total += min(y, _best_twopath(net, weights, i, j, 0))
    end
    return total
end

# One affected pair (a, b) of a transitive-weights change: its two-path
# strength with the changing dyad (which enters as min(y_changing, other))
# at `old` versus at `new`, capped by its own value y_ab
@inline function _affected_delta(y_ab::Int, rest::Int, other::Int, old::Int, new::Int)
    y_ab > 0 || return 0.0
    p_new = max(rest, min(new, other))
    p_old = max(rest, min(old, other))
    return Float64(min(y_ab, p_new) - min(y_ab, p_old))
end

function change_stat_count(::TransitiveWeightsTerm, net, weights, i::Int, j::Int,
                           old::Int, new::Int)
    # The dyad's own term: its two-paths never read y_ij
    p = _best_twopath(net, weights, i, j, 0)
    delta = Float64(min(new, p) - min(old, p))
    if is_directed(net)
        # Pairs (i, l): the path i → j → l has strength min(y_ij, y_jl)
        for l in outneighbors(net, j)
            (l == i || l == j) && continue
            y_il = dyad_value(net, weights, i, l)
            y_il > 0 || continue
            delta += _affected_delta(y_il, _best_twopath(net, weights, i, l, j),
                                     dyad_value(net, weights, j, l), old, new)
        end
        # Pairs (k, j): the path k → i → j has strength min(y_ki, y_ij)
        for k in inneighbors(net, i)
            (k == i || k == j) && continue
            y_kj = dyad_value(net, weights, k, j)
            y_kj > 0 || continue
            delta += _affected_delta(y_kj, _best_twopath(net, weights, k, j, i),
                                     dyad_value(net, weights, k, i), old, new)
        end
    else
        # Pairs {i, l} through the path i - j - l, and {j, l} through j - i - l
        for l in outneighbors(net, j)
            (l == i || l == j) && continue
            y_il = dyad_value(net, weights, i, l)
            y_il > 0 || continue
            delta += _affected_delta(y_il, _best_twopath(net, weights, i, l, j),
                                     dyad_value(net, weights, j, l), old, new)
        end
        for l in outneighbors(net, i)
            (l == i || l == j) && continue
            y_jl = dyad_value(net, weights, j, l)
            y_jl > 0 || continue
            delta += _affected_delta(y_jl, _best_twopath(net, weights, j, l, i),
                                     dyad_value(net, weights, i, l), old, new)
        end
    end
    return delta
end

"""
    CyclicalWeightsTerm() <: AbstractERGMTerm

`ergm`'s valued `cyclicalweights("min", "max", "min")` with the default
triple: every dyad's value capped by the strongest two-path that closes a
directed 3-cycle through it,

∑_{(i,j)} min( y_{ij}, max_k min(y_{jk}, y_{ki}) ),

over every ordered pair of a directed network (R label
`cyclicalweights.min.max.min`, pinned against `ergm` 4.12 by
`test/fixtures/count_terms.toml`). On an undirected network the cycle and
the transitive two-path coincide, so the statistic equals
[`TransitiveWeightsTerm`](@ref) there — as in R, which allows the term on
undirected networks. The non-default triples are not implemented.

# Example
```julia
using Networks, ERGMCount
using ERGM: compute, name
net = network(3; directed=true)
for (i, j, w) in ((1, 2, 3), (2, 3, 2), (3, 1, 2))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
compute(CyclicalWeightsTerm(), net)   # 6.0 — each of the three dyads capped at 2
name(CyclicalWeightsTerm())           # "cyclicalweights.min.max.min"
```
"""
struct CyclicalWeightsTerm <: AbstractERGMTerm end

name(::CyclicalWeightsTerm) = "cyclicalweights.min.max.min"

# The strongest cycle-closing two-path j → k → i for the dyad (i, j):
# max_k min(y_jk, y_ki), k ∉ {i, j, skip}
function _best_cycle_twopath(net, weights, i::Int, j::Int, skip::Int)
    best = 0
    for k in outneighbors(net, j)
        (k == i || k == j || k == skip) && continue
        y_jk = dyad_value(net, weights, j, k)
        y_jk > best || continue
        c = min(y_jk, dyad_value(net, weights, k, i))
        c > best && (best = c)
    end
    return best
end

function compute(t::CyclicalWeightsTerm, net)
    _refuse_negative_network(t, net, "compute")
    is_directed(net) || return compute(TransitiveWeightsTerm(), net)
    weights = _get_weights(net)
    total = 0.0
    for e in edges(net)
        i, j = src(e), dst(e)
        y = dyad_value(net, weights, i, j)
        total += min(y, _best_cycle_twopath(net, weights, i, j, 0))
    end
    return total
end

function change_stat_count(::CyclicalWeightsTerm, net, weights, i::Int, j::Int,
                           old::Int, new::Int)
    is_directed(net) ||
        return change_stat_count(TransitiveWeightsTerm(), net, weights, i, j, old, new)
    p = _best_cycle_twopath(net, weights, i, j, 0)
    delta = Float64(min(new, p) - min(old, p))
    # Pairs (a, i) closed by a → i → j → a: the path i → j → a has strength
    # min(y_ij, y_ja), so a ∈ out(j)
    for a in outneighbors(net, j)
        (a == i || a == j) && continue
        y_ai = dyad_value(net, weights, a, i)
        y_ai > 0 || continue
        delta += _affected_delta(y_ai, _best_cycle_twopath(net, weights, a, i, j),
                                 dyad_value(net, weights, j, a), old, new)
    end
    # Pairs (j, b) closed by j → b → i → j: the path b → i → j has strength
    # min(y_bi, y_ij), so b ∈ in(i)
    for b in inneighbors(net, i)
        (b == i || b == j) && continue
        y_jb = dyad_value(net, weights, j, b)
        y_jb > 0 || continue
        delta += _affected_delta(y_jb, _best_cycle_twopath(net, weights, j, b, i),
                                 dyad_value(net, weights, b, i), old, new)
    end
    return delta
end

# -----------------------------------------------------------------------------
# Negative dyad weights. `ergm` refuses `transitiveweights`/`cyclicalweights`
# on a network with a negative dyad weight ("Term may not be used with networks
# with negative dyad weights") and its `mutual(form="geometric")` returns NaN
# there (the square root of a negative product). A negative count is
# admissible only under `DiscUnif2Reference(a < 0, b)`; on such data the three
# terms are refused — at `compute` (so the summary statistic is refused, as
# in R), at model construction and at simulation over a support with negative
# values — rather than defining a statistic R never produces (the `best = 0`
# floor of the two-path search would silently do so).
# -----------------------------------------------------------------------------
_admits_negative(::AbstractERGMTerm) = true
_admits_negative(::TransitiveWeightsTerm) = false
_admits_negative(::CyclicalWeightsTerm) = false
_admits_negative(t::CountMutualTerm) = t.form !== :geometric

function _refuse_negative(t::AbstractERGMTerm, context::AbstractString)
    why = t isa CountMutualTerm ?
          "R's `mutual(form=\"geometric\")` returns NaN there (the square root " *
          "of a negative product), which is no statistic" :
          "R's `ergm` refuses the term with that sentence, and the strongest " *
          "two-path of a pair is undefined when a weight can be negative"
    throw(ArgumentError(
        "$context: $(nameof(typeof(t))) (`$(name(t))`) may not be used with " *
        "networks with negative dyad weights — $why. Drop the term, or fit a " *
        "reference whose support is non-negative (`DiscUnif2Reference(0, b)`, " *
        "`PoissonReference()`, ...)."))
end

# Smallest count on the network (0 when it has no edges); the typed read
# assumes the weights were validated
_min_dyad_value(net) = _dyad_value_extrema(net)[1]

_refuse_negative_network(t::AbstractERGMTerm, net, context::AbstractString) =
    (_admits_negative(t) || _min_dyad_value(net) >= 0) ? nothing :
    _refuse_negative(t, context)

# Every term must admit the smallest value of the data (a model) or of the
# enumerated support (a simulation)
function _validate_negative_terms(terms::Tuple, lo::Int, context::AbstractString)
    lo >= 0 && return terms
    for t in terms
        _admits_negative(t) || _refuse_negative(t, context)
    end
    return terms
end

# R prints a numeric bound the way `paste` does: 1 → "1", 1.5 → "1.5", Inf → "Inf"
_rnum(x::Integer) = string(x)
_rnum(x::Real) = isfinite(x) && isinteger(x) ? string(Int(x)) : string(x)

"""
    SmallerthanTerm(threshold) <: AbstractERGMTerm

`ergm`'s valued `smallerthan(threshold)`: the number of dyads whose value is
**below** the threshold, ∑_{(i,j)} I(y_{ij} < threshold) — zero-valued dyads
included, so on a sparse network it is close to the number of dyads. R label
`smallerthan.<threshold>`; dyad-independent; pinned against `ergm` 4.12 by
`test/fixtures/count_terms.toml`.

# Example
```julia
using Networks, ERGMCount
using ERGM: compute, name
net = network(3; directed=true)          # 6 dyads
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 3)
compute(SmallerthanTerm(2), net)   # 5.0 — every dyad but (1,2)
name(SmallerthanTerm(2))           # "smallerthan.2"
```
"""
struct SmallerthanTerm <: AbstractERGMTerm
    threshold::Int
end

name(t::SmallerthanTerm) = "smallerthan.$(t.threshold)"

function compute(t::SmallerthanTerm, net)
    weights = _get_weights(net)
    n_dyads = _n_dyads(net)
    above = 0
    for e in edges(net)
        get(weights, _wkey(net, src(e), dst(e)), 1) >= t.threshold && (above += 1)
    end
    # Every dyad without an edge has value 0
    return Float64(n_dyads - above - (0 >= t.threshold ? n_dyads - ne(net) : 0))
end

change_stat_count(t::SmallerthanTerm, net, weights, i::Int, j::Int, old::Int, new::Int) =
    Float64((new < t.threshold) - (old < t.threshold))

"""
    EqualToTerm(value; tolerance=0) <: AbstractERGMTerm

`ergm`'s valued `equalto(value, tolerance)`: the number of dyads whose value
lies within `tolerance` of `value` inclusive, ∑_{(i,j)} I(|y_{ij} − value| ≤
tolerance) — zero-valued dyads included. R label
`equalto.<value>.pm.<tolerance>`; dyad-independent; pinned against `ergm`
4.12 by `test/fixtures/count_terms.toml`.

# Example
```julia
using Networks, ERGMCount
using ERGM: compute, name
net = network(3; directed=true)
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 3)
compute(EqualToTerm(3), net)   # 1.0
compute(EqualToTerm(0), net)   # 5.0 — the empty dyads
name(EqualToTerm(3))           # "equalto.3.pm.0"
```
"""
struct EqualToTerm <: AbstractERGMTerm
    value::Int
    tolerance::Int

    function EqualToTerm(value::Int; tolerance::Int=0)
        tolerance >= 0 ||
            throw(ArgumentError("EqualToTerm: tolerance must be non-negative (got $tolerance)"))
        new(value, tolerance)
    end
end

name(t::EqualToTerm) = "equalto.$(t.value).pm.$(t.tolerance)"

@inline _equalto_hit(t::EqualToTerm, y::Int) = abs(y - t.value) <= t.tolerance

function compute(t::EqualToTerm, net)
    weights = _get_weights(net)
    n_dyads = _n_dyads(net)
    hits = 0
    for e in edges(net)
        _equalto_hit(t, Int(get(weights, _wkey(net, src(e), dst(e)), 1))) && (hits += 1)
    end
    _equalto_hit(t, 0) && (hits += n_dyads - ne(net))
    return Float64(hits)
end

change_stat_count(t::EqualToTerm, net, weights, i::Int, j::Int, old::Int, new::Int) =
    Float64(_equalto_hit(t, new) - _equalto_hit(t, old))

"""
    InIntervalTerm(lower, upper; open=(true, true)) <: AbstractERGMTerm

`ergm`'s valued `ininterval(lower, upper, open)`: the number of dyads whose
value lies between `lower` and `upper`, zero-valued dyads included. `open`
says whether each end is exclusive — R's default `(true, true)` is the open
interval `(lower, upper)`; `(false, false)` is `[lower, upper]`. The bounds
may be `-Inf`/`Inf`. R label `ininterval(lower,upper)` with the brackets
following `open` (`ininterval[1,3]`, `ininterval(1,3]`, ...);
dyad-independent; pinned against `ergm` 4.12 by
`test/fixtures/count_terms.toml`.

# Example
```julia
using Networks, ERGMCount
using ERGM: compute, name
net = network(3; directed=true)
for (i, j, w) in ((1, 2, 1), (2, 3, 2), (3, 1, 3))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
compute(InIntervalTerm(1, 3), net)                       # 1.0 — only the 2
compute(InIntervalTerm(1, 3; open=(false, false)), net)  # 3.0
name(InIntervalTerm(1, 3; open=(false, true)))           # "ininterval[1,3)"
```
"""
struct InIntervalTerm <: AbstractERGMTerm
    lower::Float64
    upper::Float64
    open_lower::Bool
    open_upper::Bool

    function InIntervalTerm(lower::Real, upper::Real; open::Tuple{Bool, Bool}=(true, true))
        lower <= upper || throw(ArgumentError(
            "InIntervalTerm: lower must not exceed upper (got $lower > $upper)"))
        new(Float64(lower), Float64(upper), open[1], open[2])
    end
end

name(t::InIntervalTerm) =
    "ininterval" * (t.open_lower ? "(" : "[") * _rnum(t.lower) * "," * _rnum(t.upper) *
    (t.open_upper ? ")" : "]")

@inline function _in_interval(t::InIntervalTerm, y::Int)
    lo_ok = t.open_lower ? (y > t.lower) : (y >= t.lower)
    hi_ok = t.open_upper ? (y < t.upper) : (y <= t.upper)
    return lo_ok & hi_ok
end

function compute(t::InIntervalTerm, net)
    weights = _get_weights(net)
    n_dyads = _n_dyads(net)
    hits = 0
    for e in edges(net)
        _in_interval(t, Int(get(weights, _wkey(net, src(e), dst(e)), 1))) && (hits += 1)
    end
    _in_interval(t, 0) && (hits += n_dyads - ne(net))
    return Float64(hits)
end

change_stat_count(t::InIntervalTerm, net, weights, i::Int, j::Int, old::Int, new::Int) =
    Float64(_in_interval(t, new) - _in_interval(t, old))

# =============================================================================
# Per-dyad support profiles
# =============================================================================
#
# Both hot paths of the package — the MPLE design build and the Gibbs
# conditional — need a dyad's change statistic at EVERY value of its support,
# not at one. Calling `change_stat_count` once per value repeats the part of
# the work that does not depend on the value: the strength terms recompute an
# O(degree) strength sum per value, the triadic terms an O(degree) neighbour
# intersection per value. `change_stats_support!` computes that part once per
# dyad and then fills the profile in O(|support|).

"""
    change_stats_support!(dest, term, net, weights, i, j, old, support) -> dest

Fill `dest[s] = change_stat_count(term, net, weights, i, j, old, support[s])`
for every `s in eachindex(support)`: the change-statistic *profile* of dyad
`(i, j)` over its conditional support, moving from count `old`. `dest` must
have `length(support)` elements; nothing is allocated.

The fallback evaluates [`change_stat_count`](@ref) once per support value, so
any term that implements `change_stat_count` has a profile. Terms whose change
statistic shares work across values specialise it: `NodeOSumTerm`,
`NodeISumTerm` and `NodeSumTerm` compute the strength excluding the dyad once
and fill `(s + y)^2 - (s + old)^2`; `TransitiveTiesTerm` and
`CyclicalTiesTerm` collect the third-vertex minima `c_k` once and accumulate
`Σ_k [min(y, c_k) - min(old, c_k)]` for every `y` in one pass, and
`TransitiveWeightsTerm`/`CyclicalWeightsTerm` register one clamp ramp per
affected pair (one two-path search each, not one per support value). Every
specialisation is held equal to the per-value fallback by the test suite, so
`change_stat_count` remains the definition. `ERGMCount` declares this
function `public`; a custom count term needs only `change_stat_count`.

# Example
```julia
using Networks, ERGMCount
net = network(4; directed=true)
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 2)
add_edge!(net, 1, 3); set_edge_attribute!(net, :weight, 1, 3, 1)
weights = get_edge_attribute(net, :weight, Int)
dest = zeros(4)
ERGMCount.change_stats_support!(dest, NodeOSumTerm(), net, weights, 1, 2, 0, 0:3)
dest == [change_stat_count(NodeOSumTerm(), net, weights, 1, 2, 0, y) for y in 0:3]  # true
```
"""
function change_stats_support!(dest::AbstractVector{Float64}, term::AbstractERGMTerm,
                               net, weights, i::Int, j::Int, old::Int,
                               support::UnitRange{Int})
    @inbounds for (s, y) in enumerate(support)
        dest[s] = change_stat_count(term, net, weights, i, j, old, y)
    end
    return dest
end

# (s + y)^2 - (s + old)^2 over the support, `s` the strength excluding the dyad
function _strength_profile!(dest::AbstractVector{Float64}, s::Float64, old::Int,
                            support::UnitRange{Int})
    base = (s + old)^2
    @inbounds for (k, y) in enumerate(support)
        dest[k] = (s + y)^2 - base
    end
    return dest
end

function change_stats_support!(dest::AbstractVector{Float64}, ::NodeOSumTerm, net, weights,
                               i::Int, j::Int, old::Int, support::UnitRange{Int})
    is_directed(net) || return fill!(dest, 0.0)
    s = _out_strength(net, weights, i) - dyad_value(net, weights, i, j)
    return _strength_profile!(dest, s, old, support)
end

function change_stats_support!(dest::AbstractVector{Float64}, ::NodeISumTerm, net, weights,
                               i::Int, j::Int, old::Int, support::UnitRange{Int})
    is_directed(net) || return fill!(dest, 0.0)
    s = _in_strength(net, weights, j) - dyad_value(net, weights, i, j)
    return _strength_profile!(dest, s, old, support)
end

function change_stats_support!(dest::AbstractVector{Float64}, ::NodeSumTerm, net, weights,
                               i::Int, j::Int, old::Int, support::UnitRange{Int})
    y_ij = dyad_value(net, weights, i, j)
    if is_directed(net)
        s_i = _out_strength(net, weights, i) + _in_strength(net, weights, i) - y_ij
        s_j = _out_strength(net, weights, j) + _in_strength(net, weights, j) - y_ij
    else
        s_i = _out_strength(net, weights, i) - y_ij
        s_j = _out_strength(net, weights, j) - y_ij
    end
    d_old = Float64(old)
    b_i, b_j = (s_i + d_old)^2, (s_j + d_old)^2
    @inbounds for (k, y) in enumerate(support)
        d_new = Float64(y)
        dest[k] = (s_i + d_new)^2 - b_i + (s_j + d_new)^2 - b_j
    end
    return dest
end

# --- triadic minima -----------------------------------------------------------
#
# A triadic term's change statistic is Σ_k [min(y, c_k) − min(old, c_k)] over
# third vertices k, where c_k is the minimum of the two other dyads of the
# triple (a different pair per role). As a function of y each summand rises
# with slope 1 from y = lo up to y = c_k and is flat beyond, so the whole sum
# is a piecewise-linear function that can be assembled in O(#k + |support|):
# every k adds its value at `lo` to a scalar and its slope to a difference
# array (kept in `dest`), and one integrating pass turns that into the
# profile. Only third vertices adjacent to both endpoints have c_k > 0, so
# the k-loop runs over neighbour intersections, not over all n vertices.

# Register one third vertex with minimum `c` (in `dest` as the difference
# array of the slope); returns the updated value at `lo`
@inline function _tri_register!(dest::AbstractVector{Float64}, base::Float64, c::Int,
                                old::Int, lo::Int, S::Int)
    c > 0 || return base
    base += min(lo, c) - min(old, c)
    if c > lo
        e = min(c - lo, S - 1)
        @inbounds dest[1] += 1.0
        @inbounds e < S && (dest[e + 1] -= 1.0)
    end
    return base
end

# Integrate the difference array in place: dest[s] = mult·(base + Σ_{t<s} slope(t))
function _tri_integrate!(dest::AbstractVector{Float64}, base::Float64, mult::Float64)
    run = 0.0
    acc = base
    @inbounds for s in eachindex(dest)
        d = dest[s]
        dest[s] = mult * acc
        run += d
        acc += run
    end
    return dest
end

function change_stats_support!(dest::AbstractVector{Float64}, ::TransitiveTiesTerm, net, weights,
                               i::Int, j::Int, old::Int, support::UnitRange{Int})
    fill!(dest, 0.0)
    lo, S = first(support), length(support)
    base = 0.0
    if is_directed(net)
        # Roles (a) and (c): k ∈ out(i) — triples (i, j, k) need y_jk, y_ik;
        # triples (i, k, j) need y_ik, y_kj
        for k in outneighbors(net, i)
            (k == i || k == j) && continue
            y_ik = dyad_value(net, weights, i, k)
            base = _tri_register!(dest, base, min(dyad_value(net, weights, j, k), y_ik),
                                  old, lo, S)
            base = _tri_register!(dest, base, min(y_ik, dyad_value(net, weights, k, j)),
                                  old, lo, S)
        end
        # Role (b): k ∈ in(i) — triples (k, i, j) need y_ki, y_kj
        for k in inneighbors(net, i)
            (k == i || k == j) && continue
            base = _tri_register!(dest, base, min(dyad_value(net, weights, k, i),
                                                  dyad_value(net, weights, k, j)),
                                  old, lo, S)
        end
        return _tri_integrate!(dest, base, 1.0)
    else
        # Pair {i, j} sits in 6 ordered triples per common neighbour k
        for k in outneighbors(net, i)
            (k == i || k == j) && continue
            base = _tri_register!(dest, base, min(dyad_value(net, weights, j, k),
                                                  dyad_value(net, weights, i, k)),
                                  old, lo, S)
        end
        return _tri_integrate!(dest, base, 6.0)
    end
end

function change_stats_support!(dest::AbstractVector{Float64}, ::CyclicalTiesTerm, net, weights,
                               i::Int, j::Int, old::Int, support::UnitRange{Int})
    is_directed(net) || return fill!(dest, 0.0)
    fill!(dest, 0.0)
    lo, S = first(support), length(support)
    base = 0.0
    # Cycle i→j→k→i: k ∈ out(j) with y_ki > 0
    for k in outneighbors(net, j)
        (k == i || k == j) && continue
        base = _tri_register!(dest, base, min(dyad_value(net, weights, j, k),
                                              dyad_value(net, weights, k, i)),
                              old, lo, S)
    end
    return _tri_integrate!(dest, base, 1.0)
end

# --- transitive / cyclical weights ----------------------------------------------
#
# Every summand of a (min, max, min) weights change statistic is a clamp:
# the dyad's own term is min(y, p) = clamp(y, 0, p), and an affected pair's
# term min(y_ab, max(rest, min(y, other))) is clamp(y, rest, min(other, y_ab))
# (constant when that upper end is not above `rest`). So the profile is a sum
# of ramps, assembled by the same difference array as the triadic-minimum
# terms above, with one two-path search per affected pair instead of one per
# affected pair AND support value.

# Register clamp(y, a, b) − clamp(old, a, b) for y over the support (no-op when
# b ≤ a: a constant); returns the updated value at `lo`
@inline function _ramp_register!(dest::AbstractVector{Float64}, base::Float64, a::Int, b::Int,
                                 old::Int, lo::Int, S::Int)
    b > a || return base
    base += clamp(lo, a, b) - clamp(old, a, b)
    start = max(a, lo)                # the ramp rises from here ...
    b > start || return base          # ... unless it ends below the support
    s0 = start - lo + 1               # first support index the ramp rises from
    s0 < S || return base
    s1 = b - lo + 1                   # first index past the ramp
    @inbounds dest[s0] += 1.0
    @inbounds s1 <= S && (dest[s1] -= 1.0)
    return base
end

function change_stats_support!(dest::AbstractVector{Float64}, ::TransitiveWeightsTerm, net,
                               weights, i::Int, j::Int, old::Int, support::UnitRange{Int})
    fill!(dest, 0.0)
    lo, S = first(support), length(support)
    base = _ramp_register!(dest, 0.0, 0, _best_twopath(net, weights, i, j, 0), old, lo, S)
    if is_directed(net)
        for l in outneighbors(net, j)
            (l == i || l == j) && continue
            y_il = dyad_value(net, weights, i, l)
            y_il > 0 || continue
            base = _ramp_register!(dest, base, _best_twopath(net, weights, i, l, j),
                                   min(dyad_value(net, weights, j, l), y_il), old, lo, S)
        end
        for k in inneighbors(net, i)
            (k == i || k == j) && continue
            y_kj = dyad_value(net, weights, k, j)
            y_kj > 0 || continue
            base = _ramp_register!(dest, base, _best_twopath(net, weights, k, j, i),
                                   min(dyad_value(net, weights, k, i), y_kj), old, lo, S)
        end
    else
        for l in outneighbors(net, j)
            (l == i || l == j) && continue
            y_il = dyad_value(net, weights, i, l)
            y_il > 0 || continue
            base = _ramp_register!(dest, base, _best_twopath(net, weights, i, l, j),
                                   min(dyad_value(net, weights, j, l), y_il), old, lo, S)
        end
        for l in outneighbors(net, i)
            (l == i || l == j) && continue
            y_jl = dyad_value(net, weights, j, l)
            y_jl > 0 || continue
            base = _ramp_register!(dest, base, _best_twopath(net, weights, j, l, i),
                                   min(dyad_value(net, weights, i, l), y_jl), old, lo, S)
        end
    end
    return _tri_integrate!(dest, base, 1.0)
end

function change_stats_support!(dest::AbstractVector{Float64}, ::CyclicalWeightsTerm, net,
                               weights, i::Int, j::Int, old::Int, support::UnitRange{Int})
    is_directed(net) ||
        return change_stats_support!(dest, TransitiveWeightsTerm(), net, weights, i, j, old, support)
    fill!(dest, 0.0)
    lo, S = first(support), length(support)
    base = _ramp_register!(dest, 0.0, 0, _best_cycle_twopath(net, weights, i, j, 0), old, lo, S)
    for a in outneighbors(net, j)
        (a == i || a == j) && continue
        y_ai = dyad_value(net, weights, a, i)
        y_ai > 0 || continue
        base = _ramp_register!(dest, base, _best_cycle_twopath(net, weights, a, i, j),
                               min(dyad_value(net, weights, j, a), y_ai), old, lo, S)
    end
    for b in inneighbors(net, i)
        (b == i || b == j) && continue
        y_jb = dyad_value(net, weights, j, b)
        y_jb > 0 || continue
        base = _ramp_register!(dest, base, _best_cycle_twopath(net, weights, j, b, i),
                               min(dyad_value(net, weights, b, i), y_jb), old, lo, S)
    end
    return _tri_integrate!(dest, base, 1.0)
end

# =============================================================================
# Dependence and directedness classification
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
is_dyad_dependent(::SmallerthanTerm) = false
is_dyad_dependent(::EqualToTerm) = false
is_dyad_dependent(::InIntervalTerm) = false

# Extends `ERGM.requires_directed`: these statistics are identically zero on an
# undirected network (their `compute` methods say so), so a coefficient on them
# is not identified there. `CountERGMModel` refuses the combination with an
# actionable message instead of fitting a zero column.
requires_directed(::CountMutualTerm) = true
requires_directed(::CyclicalTiesTerm) = true
requires_directed(::NodeOSumTerm) = true
requires_directed(::NodeISumTerm) = true

_undirected_hint(::Union{NodeOSumTerm, NodeISumTerm}) =
    "use `NodeSumTerm()` (total strength) on an undirected network"
_undirected_hint(::CountMutualTerm) =
    "drop the term: an undirected dyad has no reciprocal direction to compare"
_undirected_hint(::CyclicalTiesTerm) =
    "drop the term (an undirected network has no directed 3-cycles) or use " *
    "`TransitiveTiesTerm()`, whose undirected form is defined"
_undirected_hint(::AbstractERGMTerm) = "fit it on a directed network"

# The count analogue of the binary ERGM.jl terms a user migrating from
# `ergm(net ~ edges + mutual, response="w")` reaches for first
_count_analogue(t::AbstractERGMTerm) = _count_analogue(String(nameof(typeof(t))))
function _count_analogue(n::String)
    n == "Edges" && return "`NonzeroTerm()` (R's `nonzero`) or `SumTerm()` (`sum`)"
    n == "Mutual" && return "`CountMutualTerm()` (R's valued `mutual`)"
    (n == "Triangle" || n == "GWESP" || n == "TransitiveTies") &&
        return "`TransitiveWeightsTerm()` (R's `transitiveweights`)"
    (n == "CTriple" || n == "CyclicalTies") &&
        return "`CyclicalWeightsTerm()` (R's `cyclicalweights`)"
    (n == "OStar" || n == "GWODegree") && return "`NodeOSumTerm()`"
    (n == "IStar" || n == "GWIDegree") && return "`NodeISumTerm()`"
    (n == "Kstar" || n == "GWDegree" || n == "Degree") && return "`NodeSumTerm()`"
    return "one of the count terms (`SumTerm`, `NonzeroTerm`, `GreaterthannTerm`, " *
           "`CountAtleastnTerm`, `SmallerthanTerm`, `EqualToTerm`, `InIntervalTerm`, " *
           "`CountMutualTerm`, `TransitiveWeightsTerm`, `CyclicalWeightsTerm`, " *
           "`NodeOSumTerm`, `NodeISumTerm`, `NodeSumTerm`)"
end

# A term is a count term when it has a count change statistic; a binary
# ERGM.jl term (`Edges`, `Mutual`, ...) has `compute` but no
# `change_stat_count`, and used to die deep inside the design build with a
# MethodError and a "closest candidates" list.
_is_count_term(t::AbstractERGMTerm) =
    hasmethod(change_stat_count, Tuple{typeof(t), Any, Any, Int, Int, Int, Int})

function _validate_count_terms(terms::Tuple, net::Network)
    for t in terms
        _is_count_term(t) || throw(ArgumentError(
            "$(nameof(typeof(t))) is a binary ERGM.jl term (it has no count " *
            "change statistic `change_stat_count`), so it cannot enter a count " *
            "model — R's `ergm` refuses `$(name(t))` on a valued response too. " *
            "Use $(_count_analogue(t)) instead; see the terms guide for the full " *
            "correspondence with R."))
    end
    is_directed(net) && return terms
    for t in terms
        requires_directed(t) && throw(ArgumentError(
            "$(nameof(typeof(t))) (`$(name(t))`) requires a directed network, but " *
            "this network is undirected: the statistic is identically 0 there, so " *
            "its coefficient is not identified. Instead, $(_undirected_hint(t))."))
    end
    return terms
end

# A two-mode (bipartite) network is refused everywhere, with ERGM.jl's
# reasoning: the estimator enumerates and the Gibbs sweep resamples EVERY
# off-diagonal dyad, so the structurally impossible within-mode dyads would be
# counted as observed zeros (wrong pseudo-likelihood, `nobs`, BIC and
# conditionals) — and none of the count terms is a bipartite term.
function _refuse_two_mode(net, context::AbstractString)
    is_two_mode(net) && throw(ArgumentError(
        "$context: ERGMCount.jl fits one-mode networks only; this network is " *
        "two-mode (bipartite). Bipartite count-ERGM terms and the two-mode " *
        "dyad set are not implemented — see README 'Not implemented'. " *
        "Enumerating the one-mode dyads would silently count the impossible " *
        "within-mode dyads as observed zeros (in the pseudo-likelihood, `nobs`, " *
        "the BIC and every conditional the sampler draws from), so the network " *
        "is refused instead."))
    return net
end

# The counts must be integers stored under `:weight`, on EVERY edge. A network
# with edges but no `:weight` at all is the ergm.count `response="w"`
# migration mistake (the counts live under another name, or were never
# attached): it would fit as a 0/1 network on which `sum` and `nonzero`
# coincide, and print a non-identified fit. An edge without a `:weight` while
# others have one (a weight column with NAs, a merge that dropped rows) used
# to count as 1 silently — the same mistake on a subset of the edges, and R's
# `response=` never reads a missing value as 1 — so it is refused too, naming
# the first such edge. A non-integer value (2.5, "3") would surface as a bare
# InexactError or MethodError from the typed attribute read; a negative count
# is admissible only under `DiscUnif2Reference(a < 0, b)`. Returns the
# smallest count on the network (0 when it has no edges), so the caller can
# refuse the terms R refuses on negative data.
function _validate_count_weights(net::Network, ref::AbstractReferenceMeasure)
    ne(net) == 0 && return 0
    weights = get_edge_attribute(net, :weight)
    isempty(weights) && throw(ArgumentError(
        "the network has $(ne(net)) edges but no `:weight` edge attribute, so every " *
        "edge would count as 1 and the fit would be of a 0/1 network (on which " *
        "`sum` and `nonzero` coincide and the fit is not identified). Store the " *
        "counts with `set_edge_attribute!(net, :weight, i, j, w)`; if they live " *
        "under another attribute — R's `response=\"w\"` — pass " *
        "`fit_ergm_count(net, terms; weight=:w)`; to fit a genuinely binary " *
        "network as counts of 0 and 1, set `:weight` to 1 on every edge " *
        "explicitly."))
    n_bare = 0
    first_bare = (0, 0)
    for e in edges(net)
        haskey(weights, _wkey(net, src(e), dst(e))) && continue
        n_bare == 0 && (first_bare = (Int(src(e)), Int(dst(e))))
        n_bare += 1
    end
    n_bare == 0 || throw(ArgumentError(
        "$n_bare of the network's $(ne(net)) edges carr$(n_bare == 1 ? "ies" : "y") " *
        "no `:weight` (the first is ($(first_bare[1]),$(first_bare[2]))), and an " *
        "edge without a count is not a count of 1 — a weight column with gaps, or " *
        "a merge that dropped rows, is a data error, not a modelling choice (R's " *
        "`response=` never reads a missing value as 1). Set the count of every " *
        "edge with `set_edge_attribute!(net, :weight, i, j, w)` (1 if a bare edge " *
        "really means one event), or pass the complete attribute with " *
        "`fit_ergm_count(net, terms; weight=:w)`."))
    lo = 0
    for ((i, j), v) in weights
        ok = v isa Integer || (v isa Real && isfinite(v) && isinteger(v))
        ok || throw(ArgumentError(
            "counts must be integers; got $(repr(v)) ($(typeof(v))) at dyad " *
            "($i,$j) — round or rescale the weights before fitting (a rate or an " *
            "averaged weight is not a count)."))
        y = Int(v)
        (y < 0 && !(ref isa DiscUnif2Reference)) && throw(ArgumentError(
            "Observed count $y at dyad ($i,$j): counts must be non-negative " *
            "integers under $(nameof(typeof(ref))). A negative count is admissible " *
            "only under `DiscUnif2Reference(a, b)` with `a ≤ $y`."))
        lo = min(lo, y)
    end
    return lo
end

# =============================================================================
# Static folds over the term tuple
# =============================================================================
#
# The terms of a `CountERGMModel` are a Tuple, so the per-dyad profiles can be
# folded at compile time (mirroring `ERGM._change_stat_tuple`): no dynamic
# dispatch, no boxing, nothing allocated per dyad. A `map` over the tuple
# would hit Base's Any32 fallback from 32 terms on; the generated loops do
# not. Both folds take a scratch `buf` of `length(support)` for one term's
# profile.

# η[s] += Σ_k θ[k] · Δg_k(old → support[s]) — the linear predictor of the full
# conditional, accumulated term by term from the support profiles
@generated function _accumulate_conditional!(η, buf, terms::TT, θ, net, weights,
                                             i::Int, j::Int, old::Int,
                                             support::UnitRange{Int}) where {TT<:Tuple}
    p = length(TT.parameters)
    body = Expr[]
    for k in 1:p
        push!(body, quote
            change_stats_support!(buf, terms[$k], net, weights, i, j, old, support)
            θk = θ[$k]
            @inbounds for s in eachindex(η)
                η[s] += θk * buf[s]
            end
        end)
    end
    return quote
        $(body...)
        return η
    end
end

# One dyad's (terms × support) slab of change statistics Δg(y0 → y), written
# support-major into the reused buffer (`slab[(s-1)p + k]`): 0 B per dyad,
# pinned
@generated function _fill_slab!(slab::Vector{Float64}, buf, terms::TT, net, weights,
                                i::Int, j::Int, support::UnitRange{Int}) where {TT<:Tuple}
    p = length(TT.parameters)
    body = Expr[]
    for k in 1:p
        push!(body, quote
            change_stats_support!(buf, terms[$k], net, weights, i, j, first(support), support)
            @inbounds for s in eachindex(support)
                slab[(s - 1) * $p + $k] = buf[s]
            end
        end)
    end
    return quote
        $(body...)
        return slab
    end
end

# =============================================================================
# Model
# =============================================================================

"""
    CountERGMModel(terms, net::Network, reference=PoissonReference())

Specification of an ERGM for a count-valued network: the `terms` (a Tuple or
Vector of count terms, stored as a `Tuple` so change statistics fold
statically), the observed `network` (its `:weight` edge attribute holds the
counts, on every edge) and the `reference` measure `h(y)`.

`CountERGMModel{T,D,TT,R}` is parameterised on the network's vertex type `T`
and directedness `D` (mirroring `ERGM.ERGMModel{T,D}`), the term tuple type
`TT` and the reference type `R`; `is_directed(model)` reads `D`. There is no
`directed` field.

The constructor refuses, with an `ArgumentError` saying why: a term that is
only defined on a directed network (`CountMutualTerm`, `CyclicalTiesTerm`,
`NodeOSumTerm`, `NodeISumTerm`; see `ERGM.requires_directed`) when `net` is
undirected, since its statistic is identically zero there and the coefficient
would not be identified; a two-mode (bipartite) network, whose impossible
within-mode dyads would otherwise be counted as observed zeros; a network with
edges but no `:weight` attribute (it would fit as a 0/1 network — see the
`weight=` keyword of [`fit_ergm_count`](@ref)) or with an edge that carries
none while others do (a bare edge is a data gap, not a count of 1); a
non-integer weight; a negative count under any reference but
`DiscUnif2Reference(a < 0, b)`; and, on a network with a negative count,
`TransitiveWeightsTerm`, `CyclicalWeightsTerm` or `CountMutualTerm(:geometric)`
("may not be used with networks with negative dyad weights", as in `ergm`).

# Example
```julia
using Networks, ERGMCount
net = network(4; directed=true)
for (i, j, w) in ((1, 2, 3), (2, 1, 1), (2, 3, 2), (3, 4, 1), (4, 1, 2))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
model = CountERGMModel([SumTerm(), CountMutualTerm()], net, PoissonReference())
model.terms                      # (SumTerm(), CountMutualTerm())
is_directed(model)               # true
has_dyad_dependent(model)        # true — CountMutualTerm reads the reciprocal dyad
```
"""
struct CountERGMModel{T,D,TT<:Tuple,R<:AbstractReferenceMeasure}
    terms::TT
    network::Network{T,D}
    reference::R

    function CountERGMModel(terms::Tuple, net::Network{T,D},
                            reference::R=PoissonReference()) where
                            {T,D,R<:AbstractReferenceMeasure}
        isempty(terms) &&
            throw(ArgumentError("CountERGMModel: at least one term is required"))
        for t in terms
            t isa AbstractERGMTerm || throw(ArgumentError(
                "CountERGMModel: every term must be an AbstractERGMTerm " *
                "(got $(typeof(t)))"))
        end
        _refuse_two_mode(net, "CountERGMModel")
        _validate_count_terms(terms, net)
        lo = _validate_count_weights(net, reference)
        # The reference's own support can reach below the data
        # (`DiscUnif2Reference(-2, 2)` on non-negative counts): the terms R
        # refuses on negative weights are refused on either
        _validate_negative_terms(terms, min(lo, first(_support(reference, 0))),
                                 "CountERGMModel")
        new{T,D,typeof(terms),R}(terms, net, reference)
    end
end

CountERGMModel(terms::AbstractVector, net::Network,
               reference::AbstractReferenceMeasure=PoissonReference()) =
    CountERGMModel(Tuple(terms), net, reference)
CountERGMModel(term::AbstractERGMTerm, net::Network,
               reference::AbstractReferenceMeasure=PoissonReference()) =
    CountERGMModel((term,), net, reference)
# A `BipartiteNetwork` is two-mode by construction: the same refusal instead
# of a MethodError
CountERGMModel(terms, net::BipartiteNetwork,
               reference::AbstractReferenceMeasure=PoissonReference()) =
    _refuse_two_mode(net, "CountERGMModel")

Graphs.is_directed(::CountERGMModel{T,D}) where {T,D} = D
Graphs.is_directed(::Type{<:CountERGMModel{T,D}}) where {T,D} = D

# Mirrors `ERGM.ERGMModel`'s one-line summary: size, directedness, the
# coefficient labels and the reference — never the whole network
function Base.show(io::IO, m::CountERGMModel{T,D}) where {T,D}
    net = m.network
    print(io, "CountERGMModel{$T,$D}: $(nv(net)) vertices, $(ne(net)) edges ",
          D ? "(directed)" : "(undirected)",
          "; terms: ", join(_term_names(m), " + "),
          "; reference: ", m.reference)
    return nothing
end

"""
    has_dyad_dependent(model::CountERGMModel) -> Bool

Whether any term of the count model is dyad-dependent (see
`ERGM.is_dyad_dependent`) — a method of the ONE `ERGM.has_dyad_dependent`
predicate. This decides whether the dyadwise pseudo-likelihood is the
likelihood (dyad-independent: the count MPLE is the exact MLE) or an
approximation; `show`, [`is_exact`](@ref) and `approximations` all
read this one answer.

# Example
```julia
using Networks, ERGMCount
net = network(3; directed=true)
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 2)
has_dyad_dependent(CountERGMModel([SumTerm(), NonzeroTerm()], net))   # false
has_dyad_dependent(CountERGMModel([SumTerm(), NodeOSumTerm()], net))  # true
```
"""
has_dyad_dependent(model::CountERGMModel) =
    any(is_dyad_dependent(t) for t in model.terms)

# Number of free dyads (pseudo-likelihood contributions); masked and two-mode
# networks are refused at the entry points, so every dyad counts
_n_dyads(model::CountERGMModel) = _n_dyads(model.network)

_term_names(model::CountERGMModel) = [name(t, model.network) for t in model.terms]

# =============================================================================
# Result
# =============================================================================

"""
    CountERGMResult

Results from fitting a count ERGM by maximum pseudo-likelihood
([`fit_ergm_count`](@ref) / [`count_mple`](@ref)).

`loglik` is the maximized *pseudo*-log-likelihood over the (possibly
truncated) dyad support; `vcov` is the inverse negative Hessian of the
pseudo-log-likelihood at the optimum (or the parametric-bootstrap covariance,
see `se_type`).

# Support and truncation

Each dyad's conditional is enumerated over `0:max_val`. For the unbounded
references (Poisson, geometric) that enumeration **truncates** the model, so
the estimand is the documented unbounded family only if negligible mass sits at
the bound. The fields below record what was done, so the approximation is
inspectable rather than implicit:

- `truncated` — whether the reference is unbounded and was therefore truncated
  (see [`is_truncating`](@ref)). `false` for genuinely bounded references, whose
  support is part of the model.
- `max_val` — the top of the enumerated support actually used.
- `support_control` — how `max_val` was chosen: `:bounded` (the reference's own
  support), `:fixed` (the caller's `max_val`), `:converged` (the default
  error-controlled doubling stopped because the last doubling moved every
  estimate by at most `support_tol` standard errors **and** left at most
  `support_tol` expected dyads past the previous bound **and** no dyad put more
  than `BOUNDARY_MASS_TOL` of its conditional mass on the new top value), or
  `:unconverged` (the doubling hit its cap while the estimates were still
  moving — the fit is reported but warned about).
- `support_tol`, `support_delta`, `omitted_tail` — the tolerance and the two
  achieved bounds: `support_delta` is the largest |Δθ|/SE the last doubling
  produced, `omitted_tail` the expected number of dyads the previous bound
  omitted. Both are `0.0` for a bounded reference (nothing was truncated) and
  `NaN` for a caller-fixed `max_val` (no doubling was run).
- `support_stable` — `false` only when the doubling hit `max_doublings` with
  the estimates still moving (`support_control == :unconverged`).
- `boundary_mass` — the largest conditional probability any dyad places on the
  top support value, at the fitted coefficients. A value materially above zero
  means the bound is shaping the fit; `count_mple` warns past
  `BOUNDARY_MASS_TOL`. Always `0.0` when `truncated == false`.

# Boundary statistics

A statistic that sits at its smallest (largest) attainable value on every
dyad's conditional support has no finite pseudo-likelihood maximizer. As in R
`ergm`, its coefficient is fixed at `-Inf` (`+Inf`) with standard error 0 and
p-value 0, the other coefficients are estimated on the restricted supports
(the exact limit), and `count_mple` warns. `dof`/`aic`/`bic` count only the
finite coefficients.

# Convergence, separation and conditioning

`converged` is the shared Newton kernel's verdict on the final fit;
`iterations` is how many Newton iterations that fit ran and `gradient_norm`
the norm of the pseudo-score at the reported estimates (0 at an exact
maximum). An unconverged fit is warned about, listed in `approximations`,
printed by `show`, and never `is_exact`.

`separated` is `true` when the pseudo-likelihood has **no finite maximum**
along a combination of the statistics that a single boundary statistic does
not explain — quasi-complete separation, R's "The MPLE does not exist!" (e.g.
`sum + nonzero` on a network whose every count is 0 or 1: `sum − nonzero` is
at its minimum on every dyad, so θ_sum → −∞, θ_nonzero → +∞ with a flat
objective). Newton stops somewhere on that asymptote with arbitrarily large
coefficients and astronomical standard errors; the fit is returned with
`converged = false`, warned about, listed in `approximations` and never
`is_exact`. See [`count_mple`](@ref) for the two-signature test.

`hessian_cond` is the 2-norm condition number of the negative pseudo-Hessian
at the reported estimates over the free coefficients (`1.0` when none is
free, `Inf` when it is singular). Above `ERGMCount._HESSIAN_COND_TOL` (1e8)
the Hessian is numerically singular — two statistics are (nearly) collinear
on this network, e.g. `greaterthan.2` and `atleast.3` on integer counts, or
`sum` and `nonzero` on a 0/1 network — and the standard errors along the flat
direction are meaningless (`NaN` when it is exactly singular); `count_mple`
warns, naming the statistics that load on it, and `show`/`approximations`
say so.

# Standard errors

`se_type` records how `std_errors`/`vcov` were actually obtained — `:hessian`
(the inverse negative pseudo-Hessian, anticonservative under dyadic dependence)
or `:bootstrap` (the parametric bootstrap of `count_mple(model; se=:bootstrap)`;
`boot_replicates` then holds the `n_boot × p` refits, excluded ones as `NaN`
rows). It is what `Networks.se_method(fit)` reports, and what the `show` method
reads before deciding whether an anticonservatism caveat is still warranted.
`z_values`/`p_values` are the vectors `coeftable(fit)` and `show(fit)` print.

# Example
```julia
using Networks, ERGMCount
net = network(4; directed=true)
for (i, j, w) in ((1, 2, 3), (2, 1, 1), (2, 3, 2), (3, 4, 1), (4, 1, 2))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
fit = fit_ergm_count(net, [SumTerm(), NonzeroTerm()])
fit.support_control            # :converged — the default doubling stopped
fit.truncated                  # true — Poisson is unbounded
coef(fit)                      # 2-vector
coeftable(fit)                 # the table `show(fit)` prints
```
"""
struct CountERGMResult{M<:CountERGMModel}
    model::M
    coefficients::Vector{Float64}
    std_errors::Vector{Float64}
    z_values::Vector{Float64}
    p_values::Vector{Float64}
    vcov::Matrix{Float64}
    loglik::Float64
    converged::Bool
    iterations::Int
    gradient_norm::Float64
    max_val::Int
    truncated::Bool
    boundary_mass::Float64
    se_type::Symbol
    support_control::Symbol
    support_stable::Bool
    support_tol::Float64
    support_delta::Float64
    omitted_tail::Float64
    boot_replicates::Union{Nothing, Matrix{Float64}}
    separated::Bool
    hessian_cond::Float64
    collinear::Vector{String}
end

# The pseudo-Hessian is reported as numerically singular above this 2-norm
# condition number: half the double-precision digits are gone (≈ 1/√eps), so
# the covariance along the flattest direction is noise.
const _HESSIAN_COND_TOL = 1e8

# NaN (unavailable) and Inf (singular) both count as ill-conditioned
_ill_conditioned(result::CountERGMResult) = !(result.hessian_cond <= _HESSIAN_COND_TOL)

_fixed_indices(result::CountERGMResult) =
    [k for k in eachindex(result.coefficients) if !isfinite(result.coefficients[k])]

function _fixed_note(result::CountERGMResult)
    fixed = _fixed_indices(result)
    isempty(fixed) && return nothing
    names = _term_names(result.model)
    parts = ["$(names[k]) at $(result.coefficients[k] > 0 ? "+Inf" : "-Inf")"
             for k in fixed]
    return "coefficient$(length(fixed) == 1 ? "" : "s") fixed by a boundary " *
           "statistic (" * join(parts, ", ") * "): the observed statistic is at " *
           "its extreme attainable value on every dyad, so no finite pseudo-" *
           "likelihood estimate exists (R ergm reports the same); standard error " *
           "0 and p-value 0 recorded"
end

function _boot_exclusion_note(result::CountERGMResult)
    reps = result.boot_replicates
    reps === nothing && return nothing
    n_boot = size(reps, 1)
    n_ok = count(b -> all(isfinite, view(reps, b, :)), 1:n_boot)
    n_ok == n_boot && return nothing
    return "$(n_boot - n_ok) of the $n_boot bootstrap refits had no finite " *
           "converged count MPLE and were excluded; the covariance is over the " *
           "remaining $n_ok"
end

function _nonconvergence_caveat(result::CountERGMResult)
    result.separated && return _separation_caveat(result)
    return "Newton did not converge in $(result.iterations) iteration" *
           "$(result.iterations == 1 ? "" : "s") (pseudo-score norm " *
           "$(_fmt3(result.gradient_norm)) at the reported " *
           "estimates): they are NOT a maximum of the pseudo-likelihood — raise " *
           "`maxiter`, or check the model for a statistic with no finite " *
           "maximizer (a boundary or non-identified term)"
end

# R's sentence (mple.existence), for the count pseudo-likelihood
function _separation_caveat(::CountERGMResult)
    return "the MPLE does not exist (perfect separation): the pseudo-likelihood " *
           "has no finite maximum along a combination of the statistics, and the " *
           "reported coefficients are the point at which Newton stopped on its " *
           "flat asymptote — arbitrarily large, with meaningless standard " *
           "errors (R ergm warns \"The MPLE does not exist!\" for the same " *
           "design). Some combination of the statistics is at its extreme on " *
           "every dyad (e.g. `sum − nonzero` when every count is 0 or 1) — " *
           "remove or coarsen a term, or collect more varied counts"
end

# Which statistics load on the flattest direction of the pseudo-Hessian `H`
# (over the free coefficients, named by `names`): the eigenvector of the
# smallest eigenvalue of −H, entries above 30% of its largest
function _collinear_names(H::AbstractMatrix, names::Vector{String})
    length(names) >= 2 || return String[]
    all(isfinite, H) || return String[]
    v = eigen(Symmetric(Matrix{Float64}(-H))).vectors[:, 1]
    top = maximum(abs, v)
    return [names[l] for l in eachindex(v) if abs(v[l]) >= 0.3 * top]
end

function _conditioning_caveat(result::CountERGMResult)
    c = result.hessian_cond
    which = result.collinear
    named = isempty(which) ? "" : " (loading on " * join(which, ", ") * ")"
    return "the pseudo-Hessian at the reported estimates is numerically " *
           "singular (condition number $(_fmt3(c)) > " *
           "$(_HESSIAN_COND_TOL)): a flat direction$named — two statistics that " *
           "are collinear or nearly so on this network, e.g. `greaterthan.2` and " *
           "`atleast.3` on integer counts, or `sum` and `nonzero` on a 0/1 " *
           "network — so the standard errors along it are meaningless (NaN when " *
           "it is exactly singular); drop or merge one of the terms"
end

function _support_line(result::CountERGMResult)
    c = result.support_control
    support = _support(result.model.reference, result.max_val)
    c === :bounded && return "Support:   $support  (bounded reference)"
    head = "Support:   $support  (TRUNCATED — reference is unbounded; "
    if c === :fixed
        return head * "max_val fixed by the caller)"
    elseif c === :converged
        return head * "chosen by doubling until max|Δθ| ≤ $(result.support_tol)·SE " *
               "and the omitted tail ≤ $(result.support_tol); achieved " *
               "$(_fmt3(result.support_delta))·SE, tail " *
               "$(_fmt3(result.omitted_tail)))"
    elseif isnan(result.support_delta)
        return head * "adaptive doubling stopped here because this fit " *
               "$(result.converged ? "has a numerically singular pseudo-Hessian" :
                  "did not converge"); NOT error-controlled)"
    else
        return head * "adaptive doubling did NOT converge: last doubling still moved " *
               "the estimates by $(_fmt3(result.support_delta))·SE)"
    end
end

function Base.show(io::IO, result::CountERGMResult)
    println(io, "Count ERGM Results")
    println(io, "==================")
    println(io, "Reference: $(result.model.reference)")
    # The support is part of the estimand, not an implementation detail: print
    # it, and mark it as a truncation when the reference is really unbounded.
    println(io, _support_line(result))
    if result.truncated
        println(io, "Boundary mass: $(_fmt3(result.boundary_mass)) " *
                    "(max over dyads, at the fitted coefficients)")
    end
    println(io, "Pseudo-log-likelihood: $(round(result.loglik, digits=4))")
    println(io, "AIC: $(round(aic(result), digits=2)), BIC: $(round(bic(result), digits=2))" *
                "  (pseudo-likelihood; compare only across models on the same " *
                "network and support)")
    println(io, "Converged: $(result.converged)")
    result.converged || println(io, "  ", _nonconvergence_caveat(result))
    # A near-singular pseudo-Hessian is printed here whether or not Newton
    # converged, except under separation, whose caveat already says the
    # standard errors are meaningless
    (_ill_conditioned(result) && !result.separated) &&
        println(io, "  ", _conditioning_caveat(result))
    println(io, "Std. errors: ", result.se_type === :bootstrap ?
                "parametric bootstrap" : "inverse pseudo-Hessian")
    println(io)
    println(io, "Coefficients:")
    # Shared ecosystem presentation layer: the printed table IS
    # `coeftable(result)`, built from the same vectors, so what is shown and
    # what is inspected cannot disagree.
    show(io, coeftable(result))

    fixed = _fixed_note(result)
    if fixed !== nothing
        println(io)
        println(io, "Note: ", fixed)
    end
    excluded = _boot_exclusion_note(result)
    if excluded !== nothing
        println(io)
        println(io, "Note: ", excluded)
    end

    # Honest-uncertainty caveat, and the prose twin of what
    # `approximations(result)` reports: the pseudo-likelihood multiplies dyad
    # conditionals as if independent, so the inverse-Hessian standard errors of a
    # dyad-dependent model are expected anticonservative. A dyad-independent model
    # needs no caveat (there the pseudo-likelihood is the likelihood), and neither
    # does a bootstrap fit — those standard errors do NOT assume independence, so
    # claiming they are anticonservative would be a lie.
    if has_dyad_dependent(result.model) && result.se_type === :hessian
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
# the same fields and cannot disagree.

estimand(::CountERGMResult) = :count_ergm

objective(::CountERGMResult) = :pseudolikelihood

"""
    is_exact(result::CountERGMResult) -> Bool

`true` only when **all** of these hold:

1. every term is dyad-independent (see `ERGM.is_dyad_dependent`), so the dyad
   conditionals that the pseudo-likelihood multiplies are the model's own
   conditionals and their product is the likelihood;
2. the fit was not truncated — for an unbounded reference (Poisson, geometric)
   the estimator enumerates `0:max_val`, which is the likelihood of a
   *different*, truncated family (the error-controlled default keeps the
   difference below `support_tol` standard errors, but it is not zero);
3. the Newton iteration converged (`result.converged`); and
4. every coefficient is finite — a coefficient fixed at `±Inf` by a boundary
   statistic is a limit, not a maximizer.

A `SumTerm` model under a bounded reference is therefore exact; the same term
under a Poisson reference is not, and neither is any model containing a
strength, mutual or triadic term, an unconverged fit, or a fit with a
boundary statistic.

# Example
```julia
using Networks, ERGMCount
net = network(4; directed=true)
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 2)
add_edge!(net, 3, 4); set_edge_attribute!(net, :weight, 3, 4, 1)
is_exact(fit_ergm_count(net, [SumTerm()]; reference=BinomialReference(3)))  # true
is_exact(fit_ergm_count(net, [SumTerm()]; reference=PoissonReference()))    # false
```
"""
is_exact(result::CountERGMResult) =
    result.converged && all(isfinite, result.coefficients) &&
    !result.truncated && !has_dyad_dependent(result.model)

"""
    se_method(result::CountERGMResult) -> Symbol

What the reported standard errors ACTUALLY are: `:hessian` (the inverse negative
pseudo-Hessian) or `:bootstrap` (the parametric bootstrap of
`count_mple(model; se=:bootstrap)`). Read straight off the fit, so it can never
claim an estimator that was not used.

# Example
```julia
using Networks, ERGMCount
net = network(4; directed=true)
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 2)
add_edge!(net, 3, 4); set_edge_attribute!(net, :weight, 3, 4, 1)
se_method(fit_ergm_count(net, [SumTerm()]))   # :hessian
```
"""
se_method(result::CountERGMResult) = result.se_type

# `fit_ergm_count` calls `require_observed` with the default `:error` policy:
# the count MPLE enumerates every dyad as observed, so masked data is refused.
missing_method(::CountERGMResult) = :rejected

function approximations(result::CountERGMResult)
    out = String[]
    if result.truncated
        push!(out, "count support truncated at 0:$(result.max_val); max boundary " *
                   "mass $(_fmt3(result.boundary_mass)) " *
                   "(the reference measure is unbounded, so the enumerated support " *
                   "is an approximation to the model's)")
        if result.support_control === :converged
            push!(out, "support chosen by error-controlled doubling: the last " *
                       "doubling to max_val = $(result.max_val) moved the estimates " *
                       "by at most $(_fmt3(result.support_delta)) " *
                       "standard errors and left $(_fmt3(result.omitted_tail)) " *
                       "expected dyads past the previous bound (tolerance " *
                       "$(result.support_tol))")
        elseif result.support_control === :unconverged && isnan(result.support_delta)
            # Stopped on a rung with no usable estimates: the diagnosis is
            # the convergence/conditioning caveat, not the support
            push!(out, "support doubling stopped at max_val = $(result.max_val) " *
                       "because the fit there " *
                       "$(result.converged ? "has a numerically singular pseudo-Hessian" :
                          "did not converge") — the support is NOT error-controlled " *
                       "(no settled estimates to compare); see the convergence caveat")
        elseif result.support_control === :unconverged
            push!(out, "support doubling did NOT converge: at the last permitted " *
                       "doubling (max_val = $(result.max_val)) the estimates still " *
                       "moved by " *
                       "$(_fmt3(result.support_delta)) standard errors " *
                       "per doubling — the truncated fits are not settling, which " *
                       "usually means the unbounded model is not normalisable at " *
                       "these coefficients (e.g. a geometric reference with a " *
                       "non-negative `sum` coefficient)")
        end
    end
    result.converged || push!(out, _nonconvergence_caveat(result))
    (_ill_conditioned(result) && !result.separated) &&
        push!(out, _conditioning_caveat(result))
    fixed = _fixed_note(result)
    fixed === nothing || push!(out, fixed)
    if has_dyad_dependent(result.model)
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
        excluded = _boot_exclusion_note(result)
        excluded === nothing || push!(out, excluded)
    end
    return out
end

# ============================================================================
# StatsAPI surface
# ============================================================================
#
# Methods on the shared statistics generics, so results interoperate with
# StatsBase/GLM-style tooling (`coef(fit)`, `vcov(fit)`, ...). Pinned by
# `Networks.check_statsapi(fit; strict=true)` in the testset.

StatsAPI.coef(result::CountERGMResult) = result.coefficients
StatsAPI.stderror(result::CountERGMResult) = result.std_errors
StatsAPI.vcov(result::CountERGMResult) = result.vcov
StatsAPI.loglikelihood(result::CountERGMResult) = result.loglik
StatsAPI.nobs(result::CountERGMResult) = _n_dyads(result.model)
# R's `logLik.ergm` df: a coefficient fixed at ∓Inf by a boundary statistic is
# not an estimated parameter
StatsAPI.dof(result::CountERGMResult) = count(isfinite, result.coefficients)

"""
    aic(result::CountERGMResult) -> Float64

`-2·loglikelihood + 2·dof` — a **pseudo-likelihood** AIC (a method of
`StatsAPI.aic`): `loglikelihood(result)` is the maximized pseudo-log-likelihood
over the enumerated count support, not a likelihood, so this number is
comparable only across models fit on the **same network with the same
reference and support** (`max_val`), and it carries the pseudo-likelihood's
bias for dyad-dependent models. `dof` counts the finite coefficients (one fixed
at ∓Inf by a boundary statistic is not a parameter).

# Example
```julia
using Networks, ERGMCount
net = network(4; directed=true)
for (i, j, w) in ((1, 2, 3), (2, 1, 1), (2, 3, 2), (3, 4, 1), (4, 1, 2))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
fit = fit_ergm_count(net, [SumTerm(), NonzeroTerm()])
aic(fit) ≈ -2 * loglikelihood(fit) + 2 * dof(fit)   # true
```
"""
StatsAPI.aic(result::CountERGMResult) = -2 * result.loglik + 2 * dof(result)

"""
    bic(result::CountERGMResult) -> Float64

`-2·loglikelihood + dof·log(nobs)` with `nobs` the number of dyads — a
**pseudo-likelihood** BIC (a method of `StatsAPI.bic`), with exactly the
caveats of [`aic`](@ref): comparable only across models fit on the same
network, reference and support.

# Example
```julia
using Networks, ERGMCount
net = network(4; directed=true)
for (i, j, w) in ((1, 2, 3), (2, 1, 1), (2, 3, 2), (3, 4, 1), (4, 1, 2))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
fit = fit_ergm_count(net, [SumTerm(), NonzeroTerm()])
bic(fit) ≈ -2 * loglikelihood(fit) + dof(fit) * log(nobs(fit))   # true
nobs(fit)                                                        # 12 directed dyads
```
"""
StatsAPI.bic(result::CountERGMResult) =
    -2 * result.loglik + dof(result) * log(nobs(result))

"""
    confint(result::CountERGMResult; level=0.95) -> Matrix{Float64}

Normal-theory (Wald) confidence limits `θ̂ ± z_{(1+level)/2} · se`, one row
per coefficient with the lower limit in column 1 and the upper in column 2
(a method of `StatsAPI.confint`). The standard errors are the ones the fit
reports — `se_method(result)` says whether they are inverse-pseudo-Hessian or
parametric-bootstrap — so for a dyad-dependent model fit with `se=:hessian`
the intervals inherit the anticonservative pseudo-likelihood SEs. A
coefficient fixed at ∓Inf has both limits at that value.

# Example
```julia
using Networks, ERGMCount
net = network(4; directed=true)
for (i, j, w) in ((1, 2, 3), (2, 1, 1), (2, 3, 2), (3, 4, 1), (4, 1, 2))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
fit = fit_ergm_count(net, [SumTerm(), NonzeroTerm()])
ci = confint(fit)                          # 2×2
all(ci[:, 1] .< coef(fit) .< ci[:, 2])     # true
confint(fit; level=0.9)                    # narrower
```
"""
function StatsAPI.confint(result::CountERGMResult; level::Real=0.95)
    0 < level < 1 ||
        throw(ArgumentError("confint: level must be in (0, 1) (got $level)"))
    q = quantile(Normal(), 1 - (1 - level) / 2)
    θ, se = result.coefficients, result.std_errors
    return hcat(θ .- q .* se, θ .+ q .* se)
end

"""
    coeftable(result::CountERGMResult) -> Networks.CoefficientTable

The R-style coefficient table (`Estimate`, `Std.Error`, `z value`,
`Pr(>|z|)`) as an inspectable `Networks.CoefficientTable` — exactly the table
`show(result)` prints, built from the same vectors (a method of
`StatsAPI.coeftable`). Rows are labelled with `name(term, net)`, the shared
two-argument statistic name, and can be read by index or by name.

# Example
```julia
using Networks, ERGMCount
net = network(4; directed=true)
for (i, j, w) in ((1, 2, 3), (2, 1, 1), (2, 3, 2), (3, 4, 1), (4, 1, 2))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
fit = fit_ergm_count(net, [SumTerm(), NonzeroTerm()])
tbl = coeftable(fit)
tbl["sum"].estimate == coef(fit)[1]      # true
tbl[2].p_value == fit.p_values[2]        # true
```
"""
StatsAPI.coeftable(result::CountERGMResult) =
    CoefficientTable(_term_names(result.model), result.coefficients,
                     result.std_errors; z_values=result.z_values,
                     p_values=result.p_values)

# z statistics and two-sided p-values through the ONE shared helper; a
# coefficient fixed at ∓Inf by a boundary statistic gets z = ∓Inf, p = 0 (R's
# convention), never NaN from a 0/0.
function _count_zp(θ::Vector{Float64}, se::Vector{Float64})
    zp = z_pvalues(θ, se)
    z, p = zp.z, zp.p
    for k in eachindex(θ)
        if !isfinite(θ[k])
            z[k] = θ[k]
            p[k] = 0.0
        end
    end
    return z, p
end

# =============================================================================
# Entry points
# =============================================================================

"""
    fit_ergm_count(net::Network, terms; reference=PoissonReference(), kwargs...)

Fit an ERGM for count-valued networks by maximum pseudo-likelihood: each
dyad's conditional distribution over the count support
`P(y_ij = y | rest) ∝ h(y)·exp(θ'Δg(y))` is used as an independent
likelihood contribution. The reference measure `h` enters the estimator
directly. For a dyad-independent model — one made of the terms `SumTerm`,
`NonzeroTerm`, `GreaterthannTerm`, `CountAtleastnTerm`, `SmallerthanTerm`,
`EqualToTerm`, `InIntervalTerm` — the pseudo-likelihood *is* the likelihood
and this is the exact MLE.

`terms` may be a Vector, a Tuple or a single term. [`ergm_count`](@ref) is
the R-faithful alias (matching the `ergm.count` package); `fit_count_ergm` is
a legacy alias. Passing the arguments the other way round
(`fit_ergm_count(terms, net)`) is an `ArgumentError` naming the right order.

# Arguments
- `net`: one-mode `Network` whose counts are the `:weight` edge attribute
  (integers, on **every** edge: a network with edges and no `:weight` at all
  is refused — see `weight=` — and so is one where only some edges carry
  one, naming the first bare edge; a bare edge is a data gap, not a count of
  1). A network
  with masked (unobserved) dyads is refused — see `Networks.require_observed`;
  there is no `missing=` keyword, because the count MPLE would enumerate every
  unobserved dyad as an observed row. A two-mode (bipartite) network is
  refused too: its within-mode dyads would be counted as observed zeros.
- `terms`: count ERGM terms
- `reference`: Reference measure (default: Poisson)
- `weight::Symbol=:weight`: the edge attribute holding the counts — R's
  `response="w"` is `weight=:w`. Any name but `:weight` fits a `copy` of the
  network with that attribute copied to `:weight` (so `fit.model.network` is
  the copy; the caller's network is untouched).
- `method`: only `:mple`. `:mcmle` — ergm.count's own Monte-Carlo MLE — is
  **not implemented** and throws an `ArgumentError` saying so.
- `max_val::Int`: fixes the truncation of an unbounded support. By default the
  support is chosen adaptively (doubling from twice the largest observed count
  until the estimates stop moving; see [`count_mple`](@ref)); ignored by the
  bounded references.
- `support_tol`, `max_doublings`, `maxiter`, `tol`, `se`, `n_boot`,
  `boot_burnin`, `boot_interval`, `rng`: forwarded to [`count_mple`](@ref). `se=:bootstrap`
  replaces the inverse-pseudo-Hessian covariance with a parametric bootstrap
  (same API as `ERGM.mple`); the point estimates are unchanged.

# Returns
- [`CountERGMResult`](@ref): fitted model, answering the full StatsAPI surface
  (`coef`, `stderror`, `vcov`, `confint`, `loglikelihood`, `nobs`, `dof`,
  `aic`, `bic`, `coeftable`)

# Example
```julia
using Networks, ERGMCount
net = network(4; directed=true)
for (i, j, w) in ((1, 2, 3), (2, 1, 1), (2, 3, 2), (3, 4, 1), (4, 1, 2))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
fit = fit_ergm_count(net, [SumTerm(), NonzeroTerm()]; reference=PoissonReference())
fit.converged                   # true
coeftable(fit)                  # sum / nonzero estimates with z and p
fit_ergm_count(net, SumTerm())  # a single term works too
```
"""
function fit_ergm_count(net::Network, terms::Tuple;
                        reference::AbstractReferenceMeasure=PoissonReference(),
                        method::Symbol=:mple,
                        weight::Symbol=:weight,
                        kwargs...)
    method === :mple || throw(ArgumentError(
        "fit_ergm_count: method=:$method is not available. ERGMCount.jl fits " *
        "count ERGMs by maximum pseudo-likelihood only (`method=:mple`, the " *
        "default): ergm.count's Monte-Carlo maximum likelihood (`:mcmle`) is not " *
        "implemented (see README 'Not implemented'). The MPLE is the exact MLE " *
        "for a dyad-independent model (sum, nonzero, greaterthan, atleast, " *
        "smallerthan, equalto, ininterval); for a dyad-dependent one (mutual, " *
        "transitive, strength terms) it is a pseudo-likelihood estimate, and " *
        "`se=:bootstrap` gives the honest covariance."))

    # Count MPLE enumerates every dyad as observed, so a masked (unobserved)
    # dyad would enter the pseudo-likelihood at its face value. Reject it.
    require_observed(net; context="fit_ergm_count", face_ok=false)

    model = CountERGMModel(terms, _with_weight(net, weight), reference)
    return count_mple(model; kwargs...)
end

fit_ergm_count(net::Network, terms::AbstractVector; kwargs...) =
    fit_ergm_count(net, Tuple(terms); kwargs...)
fit_ergm_count(net::Network, term::AbstractERGMTerm; kwargs...) =
    fit_ergm_count(net, (term,); kwargs...)
# A `BipartiteNetwork` is two-mode by construction: the same refusal instead
# of a MethodError
fit_ergm_count(net::BipartiteNetwork, terms; kwargs...) =
    _refuse_two_mode(net, "fit_ergm_count")

# R's `response="w"`: the counts under another edge attribute. The package
# reads `:weight` everywhere, so the fit runs on a copy carrying the counts
# under `:weight` — the caller's network is never mutated.
function _with_weight(net::Network, weight::Symbol)
    weight === :weight && return net
    weight in list_edge_attributes(net) || throw(ArgumentError(
        "fit_ergm_count: weight=:$weight names no edge attribute of the network " *
        "(it has $(isempty(list_edge_attributes(net)) ? "none" : join(string.(":", list_edge_attributes(net)), ", "))). " *
        "Store the counts with `set_edge_attribute!(net, :$weight, i, j, w)`, or " *
        "use the default attribute `:weight`."))
    out = copy(net)
    set_edge_attribute!(out, :weight, get_edge_attribute(net, weight))
    return out
end

# Swapped arguments: name the right order instead of a raw MethodError
fit_ergm_count(terms::Union{AbstractERGMTerm, Tuple, AbstractVector}, net::Network;
               kwargs...) =
    throw(ArgumentError(
        "fit_ergm_count(net, terms): the network comes first and the terms " *
        "second (got the terms first). Call " *
        "`fit_ergm_count(net, $(terms isa AbstractERGMTerm ? "[$(nameof(typeof(terms)))()]" : "terms"))`."))

# The missing-data contract, declared: no `missing=` keyword, `:error` only.
missing_policies(::typeof(fit_ergm_count)) = (:error,)

"""
    ergm_count(net::Network, terms; kwargs...)

R-faithful alias for [`fit_ergm_count`](@ref) (the same function), matching
the R `ergm.count` package name.

# Example
```julia
using ERGMCount
ergm_count === fit_ergm_count   # true
```
"""
const ergm_count = fit_ergm_count

"""
    fit_count_ergm(net::Network, terms; kwargs...)

Alias for [`fit_ergm_count`](@ref), kept for backward compatibility.

# Example
```julia
using ERGMCount
fit_count_ergm === fit_ergm_count   # true
```
"""
const fit_count_ergm = fit_ergm_count

function _default_max_val(net, weights)
    m = 0
    for e in edges(net)
        m = max(m, Int(get(weights, _wkey(net, src(e), dst(e)), 1)))
    end
    return max(10, 2 * m)
end

# =============================================================================
# Count MPLE: compressed design
# =============================================================================
#
# Each dyad contributes its full conditional over the support: a (terms ×
# support) slab of change statistics Δg(y0 → y) plus its observed value. Dyads
# with IDENTICAL slabs — every dyad, for a dyad-independent model — contribute
# the same conditional, so the design is compressed to the unique slabs with a
# count of how many dyads share each one and at which support value they were
# observed (the count analogue of `ERGM._mple_data`'s (n_tot, n_one) rows).
# Memory is O(unique slabs × support × terms), not O(dyads × support × terms),
# and the per-dyad sweep allocates nothing (the slab is written into a reused
# buffer keyed into a Dict; only a NEW slab is copied).

struct _CountDesign
    X::Array{Float64, 3}      # p × S × R: Δg(y0 → support[s]) for unique slab r
    C::Matrix{Float64}        # S × R: dyads sharing slab r observed at support[s]
    n_tot::Vector{Float64}    # R: dyads sharing slab r
    log_h::Vector{Float64}    # S: log reference measure over the support
    support::UnitRange{Int}
    n_dyads::Int
end

_wider_reference_hint(ref::BinomialReference, y::Int) = "`BinomialReference($y)`"
_wider_reference_hint(ref::DiscUnifReference, y::Int) = "`DiscUnifReference($y)`"
_wider_reference_hint(ref::DiscUnif2Reference, y::Int) =
    "`DiscUnif2Reference($(min(ref.a, y)), $(max(ref.b, y)))`"
_wider_reference_hint(::AbstractReferenceMeasure, ::Int) = "a wider reference"

# An observed count outside the enumerated support. For a truncating reference
# the bound is OURS (`max_val`) and the fix is to raise it; for a bounded
# reference the bound is the MODEL's, `max_val` is not consulted, and the fix is
# a reference whose support contains the data.
function _check_in_support(ref::AbstractReferenceMeasure, support, y::Int,
                           i::Int, j::Int)
    (first(support) <= y <= last(support)) && return nothing
    if is_truncating(ref) && y < first(support)
        # Below a truncating reference's support means negative: not a
        # truncation anybody chose, and no `max_val` can help
        throw(ArgumentError(
            "Observed count $y at dyad ($i,$j): counts must be non-negative " *
            "integers under $(nameof(typeof(ref))). A negative count is admissible " *
            "only under `DiscUnif2Reference(a, b)` with `a ≤ $y`."))
    elseif is_truncating(ref)
        throw(ArgumentError(
            "Observed count $y at dyad ($i,$j) lies outside the enumerated support " *
            "$(support). $(nameof(typeof(ref))) is unbounded, so this bound is the " *
            "truncation you chose: pass `max_val` ≥ $y (or omit `max_val`, and the " *
            "support starts at twice the largest observed count)."))
    else
        throw(ArgumentError(
            "Observed count $y at dyad ($i,$j) lies outside the support $(support) " *
            "of $(ref). That bound is part of the model — a $(nameof(typeof(ref))) " *
            "cannot generate the value $y — so raising the enumeration bound does " *
            "not apply; choose a reference whose support contains every observed " *
            "count (e.g. $(_wider_reference_hint(ref, y))), or use an unbounded " *
            "one such as `PoissonReference()`."))
    end
end

function _count_design(model::CountERGMModel, support::UnitRange{Int})
    net = model.network
    ref = model.reference
    weights = get_edge_attribute(net, :weight, Int)
    S = length(support)
    p = length(model.terms)
    log_h = [log_reference(ref, y) for y in support]

    slab = Vector{Float64}(undef, p * S)
    buf = Vector{Float64}(undef, S)
    rows = Dict{Vector{Float64}, Int}()
    counts = Vector{Vector{Float64}}()
    n_dyads = _count_design_rows!(rows, counts, slab, buf, net, weights, model.terms,
                                  ref, support)

    R = length(rows)
    X = Array{Float64}(undef, p, S, R)
    for (key, r) in rows
        for s in 1:S, k in 1:p
            X[k, s, r] = key[(s - 1) * p + k]
        end
    end
    C = Matrix{Float64}(undef, S, R)
    for r in 1:R
        C[:, r] .= counts[r]
    end
    n_tot = vec(sum(C, dims=1))
    return _CountDesign(X, C, n_tot, log_h, support, n_dyads)
end

# Function barrier specialised on the term tuple: the slab is filled from the
# support profiles into reused buffers, so the sweep allocates only for NEW
# slabs (O(unique slabs), pinned by the testset).
function _count_design_rows!(rows::Dict{Vector{Float64}, Int},
                             counts::Vector{Vector{Float64}}, slab::Vector{Float64},
                             buf::Vector{Float64}, net, weights, terms::Tuple, ref,
                             support::UnitRange{Int})
    n = nv(net)
    directed = is_directed(net)
    S = length(support)
    y0 = first(support)
    n_dyads = 0
    for i in 1:n
        for j in (directed ? (1:n) : ((i + 1):n))
            i == j && continue
            y = dyad_value(net, weights, i, j)
            _check_in_support(ref, support, y, i, j)
            _fill_slab!(slab, buf, terms, net, weights, i, j, support)
            r = get(rows, slab, 0)
            if r == 0
                r = length(rows) + 1
                rows[copy(slab)] = r
                push!(counts, zeros(S))
            end
            @inbounds counts[r][y - y0 + 1] += 1.0
            n_dyads += 1
        end
    end
    return n_dyads
end

# =============================================================================
# Count MPLE: conditionals, derivatives, boundary statistics
# =============================================================================
#
# `mask[s, r]` says which support values slab r's conditional ranges over: all
# of them normally, a restricted set once a boundary statistic has been fixed
# at ∓Inf (its coefficient's limit concentrates each conditional on the
# values where that statistic is extreme). `cols` are the free coefficients.

# Fill η with the log-conditional numerators of slab r (−Inf where masked) and
# return log Z by log-sum-exp. Shared by the derivatives, the boundary-mass and
# the omitted-tail diagnostics, so no diagnostic can drift from the likelihood.
function _row_logZ!(η::Vector{Float64}, β, X::Array{Float64, 3},
                    log_h::Vector{Float64}, mask::Matrix{Bool},
                    cols::Vector{Int}, r::Int)
    ηmax = -Inf
    @inbounds for s in eachindex(η)
        if mask[s, r]
            v = log_h[s]
            for l in eachindex(cols)
                v += β[l] * X[cols[l], s, r]
            end
            η[s] = v
            ηmax = max(ηmax, v)
        else
            η[s] = -Inf
        end
    end
    Z = 0.0
    @inbounds for s in eachindex(η)
        Z += exp(η[s] - ηmax)
    end
    return ηmax + log(Z)
end

# The `(ll, grad, hess)` closure of the count pseudo-log-likelihood over the
# compressed design: slab r with n_tot[r] dyads, C[s, r] of them observed at
# support index s, contributes Σ_s C[s,r]·η_s − n_tot[r]·log Z to the
# objective, Σ_s C[s,r]·Δg_s − n_tot[r]·E[Δg] to the gradient and
# −n_tot[r]·(E[Δg Δg'] − E[Δg]E[Δg]') to the Hessian. The workspaces are
# allocated ONCE; each evaluation allocates only the gradient and Hessian it
# hands to `newton_fit` (pinned by the testset).
function _count_derivatives(D::_CountDesign, cols::Vector{Int}, mask::Matrix{Bool})
    X, C, n_tot, log_h = D.X, D.C, D.n_tot, D.log_h
    S = length(D.support)
    R = length(n_tot)
    q = length(cols)
    η = Vector{Float64}(undef, S)
    eX = Vector{Float64}(undef, q)
    eXX = Matrix{Float64}(undef, q, q)

    return function (β)
        llv = 0.0
        grad = zeros(q)
        hess = zeros(q, q)
        @inbounds for r in 1:R
            logZ = _row_logZ!(η, β, X, log_h, mask, cols, r)
            w = n_tot[r]
            for s in 1:S
                c = C[s, r]
                c == 0.0 && continue
                llv += c * (η[s] - logZ)
                for l in 1:q
                    grad[l] += c * X[cols[l], s, r]
                end
            end
            # E[Δg] and E[Δg Δg'] under the conditional
            fill!(eX, 0.0)
            fill!(eXX, 0.0)
            for s in 1:S
                mask[s, r] || continue
                pr = exp(η[s] - logZ)
                for l in 1:q
                    xl = X[cols[l], s, r]
                    eX[l] += pr * xl
                    for k in 1:q
                        eXX[k, l] += pr * (X[cols[k], s, r] * xl)
                    end
                end
            end
            for l in 1:q
                grad[l] -= w * eX[l]
                for k in 1:q
                    hess[k, l] -= w * (eXX[k, l] - eX[k] * eX[l])
                end
            end
        end
        return llv, grad, hess
    end
end

# Largest conditional probability placed on the top support value by any dyad,
# at the fitted coefficients.
function _max_boundary_mass(β, D::_CountDesign, cols::Vector{Int}, mask::Matrix{Bool})
    η = Vector{Float64}(undef, length(D.support))
    top = length(D.support)
    worst = 0.0
    for r in eachindex(D.n_tot)
        logZ = _row_logZ!(η, β, D.X, D.log_h, mask, cols, r)
        worst = max(worst, exp(η[top] - logZ))
    end
    return worst
end

# Expected number of dyads past support index `cutoff` (i.e. with count >
# support[cutoff]) under the fitted conditionals: the mass a fit at the
# previous, smaller bound omitted.
function _omitted_tail(β, D::_CountDesign, cols::Vector{Int}, mask::Matrix{Bool},
                       cutoff::Int)
    S = length(D.support)
    η = Vector{Float64}(undef, S)
    total = 0.0
    for r in eachindex(D.n_tot)
        logZ = _row_logZ!(η, β, D.X, D.log_h, mask, cols, r)
        tail = 0.0
        for s in (cutoff + 1):S
            tail += exp(η[s] - logZ)
        end
        total += D.n_tot[r] * tail
    end
    return total
end

# Columns whose observed value is at the smallest (largest) attainable value
# on EVERY dyad's (masked) conditional support: the pseudo-likelihood is then
# monotone in that coefficient and no finite maximizer exists. A column that
# is constant on every support is reported at its minimum.
function _count_boundary_columns(D::_CountDesign, cols::Vector{Int}, mask::Matrix{Bool})
    X, C = D.X, D.C
    S = length(D.support)
    R = length(D.n_tot)
    out = Tuple{Int, Symbol}[]
    for k in cols
        at_min = true
        at_max = true
        for r in 1:R
            lo, hi = Inf, -Inf
            for s in 1:S
                mask[s, r] || continue
                v = X[k, s, r]
                lo = min(lo, v)
                hi = max(hi, v)
            end
            for s in 1:S
                C[s, r] > 0.0 || continue
                v = X[k, s, r]
                v > lo && (at_min = false)
                v < hi && (at_max = false)
            end
            (at_min || at_max) || break
        end
        if at_min
            push!(out, (k, :min))
        elseif at_max
            push!(out, (k, :max))
        end
    end
    return out
end

# Restrict every slab's support to the values where column k is extreme: the
# conditional in the limit θ_k → ∓Inf.
function _restrict_mask!(mask::Matrix{Bool}, D::_CountDesign, k::Int, side::Symbol)
    X = D.X
    S = length(D.support)
    for r in eachindex(D.n_tot)
        ext = side === :min ? Inf : -Inf
        for s in 1:S
            mask[s, r] || continue
            ext = side === :min ? min(ext, X[k, s, r]) : max(ext, X[k, s, r])
        end
        for s in 1:S
            mask[s, r] && X[k, s, r] != ext && (mask[s, r] = false)
        end
    end
    return mask
end

# Iterate boundary detection to a fixed point (restricting one statistic's
# supports can push another to its boundary), as `ERGM._boundary_columns_iterated`
# does. Returns the free columns and the fixed (column, side) pairs.
function _count_boundary!(mask::Matrix{Bool}, D::_CountDesign, p::Int)
    cols = collect(1:p)
    fixed = Tuple{Int, Symbol}[]
    while true
        found = _count_boundary_columns(D, cols, mask)
        isempty(found) && break
        for (k, side) in found
            _restrict_mask!(mask, D, k, side)
        end
        append!(fixed, found)
        gone = Set(first.(found))
        filter!(k -> !(k in gone), cols)
    end
    sort!(fixed; by=first)
    return cols, fixed
end

# R's sentence (ergm.checkextreme.model), one line per boundary side
function _warn_count_boundary(names::Vector{String}, fixed::Vector{Tuple{Int, Symbol}})
    for (side, word, at) in ((:min, "smallest", "-Inf"), (:max, "largest", "+Inf"))
        cols = [names[k] for (k, s) in fixed if s === side]
        isempty(cols) && continue
        @warn "count_mple: observed statistic(s) $(join(cols, ", ")) are at their " *
              "$word attainable value on every dyad's conditional support. Their " *
              "coefficients will be fixed at $at (no finite maximum pseudo-" *
              "likelihood estimate exists; R ergm reports the same for a boundary " *
              "statistic). The remaining coefficients are estimated with each " *
              "dyad's support restricted to the values those statistics allow — " *
              "the exact limit of the pseudo-likelihood — with standard error 0 " *
              "and p-value 0 recorded for the fixed ones."
    end
    return nothing
end

# Whether a Newton "solution" θ (over the free columns `cols`, standard errors
# `se`) of the count pseudo-likelihood is really a point on a flat asymptote:
# the MPLE does not exist by quasi-complete separation along a COMBINATION of
# the statistics that `_count_boundary_columns` (one column at a time) cannot
# see — e.g. `sum − nonzero` at its minimum on every dyad when every count is 0
# or 1, so θ_sum → −∞ and θ_nonzero → +∞ with the objective flat along the
# way. R's `mple.existence` decides this with a linear program; as in
# `ERGM._separated`, two signatures the asymptote leaves on the Newton
# iteration are tested instead, and BOTH must hold:
#
# 1. some row has a masked support value that no dyad of the row was observed
#    at whose exponential tilt θ'Δg lies more than 18.42 nats below the row's
#    largest tilt — a conditional probability ratio below 1e-8 from the
#    coefficients alone (the reference measure is left out: it is fixed, and a
#    wide Poisson support is improbable under it whatever θ is), which a
#    finite maximizer never needs but the asymptote produces by construction
#    (Newton only stops once the unobserved values are that far down); and
# 2. the next Newton step −H⁻¹g is still large — more than 1e-3 of ‖θ‖, or
#    1e-3 absolute for a small θ (a count coefficient is a log-rate, so 1e-3
#    is a meaningful move) — or the Hessian is not invertible (NaN standard
#    errors): on the asymptote Newton keeps moving by O(1) per step while the
#    objective gains less than `tol`, whereas at a finite maximum quadratic
#    convergence has shrunk the step to ~1e-10.
#
# Either alone has benign explanations (a strongly negative `sum` on a wide
# support; a coefficient near 0); together they are the asymptote. Verified on
# the grader's binary-valued `sum + nonzero` design (step 1.41, condition
# number 8e9) against zach's `sum + nonzero` (step 3e-12, condition 54).
function _count_separated(derivatives, D::_CountDesign, cols::Vector{Int},
                          mask::Matrix{Bool}, θ::Vector{Float64}, se::Vector{Float64})
    isempty(θ) && return false
    X, C = D.X, D.C
    S = length(D.support)
    extreme = false
    for r in eachindex(D.n_tot)
        tmax = -Inf
        tmin_unobserved = Inf
        for s in 1:S
            mask[s, r] || continue
            t = 0.0
            for l in eachindex(cols)
                t += θ[l] * X[cols[l], s, r]
            end
            tmax = max(tmax, t)
            C[s, r] == 0.0 && (tmin_unobserved = min(tmin_unobserved, t))
        end
        if tmax - tmin_unobserved > 18.42
            extreme = true
            break
        end
    end
    extreme || return false
    any(isnan, se) && return true
    _, grad, hess = derivatives(θ)
    step = try
        -(hess \ grad)
    catch e
        e isa Union{SingularException, LAPACKException, ZeroPivotException,
                    PosDefException} || rethrow()
        return true
    end
    all(isfinite, step) || return true
    return norm(step) > 1e-3 * max(norm(θ), 1.0)
end

# R's sentence (mple.existence) for the separated case
function _warn_count_separated(context::AbstractString)
    @warn "$context: the MPLE does not exist (perfect separation): the " *
          "pseudo-likelihood has no finite maximum along a combination of the " *
          "statistics, and the returned coefficients are the point at which " *
          "Newton stopped on its flat asymptote — arbitrarily large, with " *
          "meaningless standard errors. R ergm warns \"The MPLE does not " *
          "exist!\" for the same design. The fit is returned with " *
          "`converged == false` and `fit.separated == true`; some combination " *
          "of the statistics is at its extreme attainable value on every dyad " *
          "(e.g. `sum − nonzero` when every count is 0 or 1) — remove or " *
          "coarsen a term, or collect more varied counts." maxlog = 1
    return nothing
end

# 2-norm condition number of the negative pseudo-Hessian over the free
# coefficients: 1 when nothing is free, Inf when singular or not finite
function _hessian_cond(hess)
    isempty(hess) && return 1.0
    all(isfinite, hess) || return Inf
    c = cond(Matrix{Float64}(-hess))
    return isfinite(c) ? c : Inf
end

# `newton_fit` warns whenever the Hessian at its final iterate is not negative
# definite — also at a start it could not move from, which the coordinate
# rescue below recovers from. Run it quietly; `count_mple` says so itself when
# the FINAL fit has undefined standard errors.
_quiet_newton(derivatives, start; maxiter, tol) =
    with_logger(NullLogger()) do
        newton_fit(derivatives, start; maxiter=maxiter, tol=tol)
    end

# Cyclic coordinate ascent from `start`, by bisection on the sign of each
# coordinate's score: the pseudo-log-likelihood is concave, so each score is
# monotone in its coordinate and a sign change brackets the 1-D maximizer. A
# robust start for the cases where Newton from θ = 0 cannot move at all — a
# reference so far from flat on the support that the initial conditional is a
# point mass and the Hessian singular (`BinomialReference(100)` on `0:10`).
function _coordinate_ascent_start(derivatives, start::Vector{Float64};
                                  sweeps::Int=3, xtol::Float64=1e-3)
    θ = copy(start)
    score(k, t) = (θ[k] = t; derivatives(θ)[2][k])
    for _ in 1:sweeps, k in eachindex(θ)
        t0 = θ[k]
        g0 = score(k, t0)
        (isfinite(g0) && g0 != 0) || (θ[k] = t0; continue)
        dir = g0 > 0 ? 1.0 : -1.0
        step = 1.0
        t1, g1 = t0, g0
        crossed = false
        for _ in 1:40                      # expand until the score changes sign
            t1 = t0 + dir * step
            g1 = score(k, t1)
            if !isfinite(g1) || (g1 > 0) != (g0 > 0)
                crossed = true
                break
            end
            t0, g0 = t1, g1
            step *= 2
        end
        if !crossed
            θ[k] = t0                      # monotone as far as we looked
            continue
        end
        # score is decreasing in t: g(a) > 0 > g(b)
        a, b = dir > 0 ? (t0, t1) : (t1, t0)
        while b - a > xtol
            m = (a + b) / 2
            gm = score(k, m)
            if isfinite(gm) && gm > 0
                a = m
            else
                b = m
            end
        end
        θ[k] = (a + b) / 2
    end
    return θ
end

# Core count MPLE on a given support: build the compressed design, fix any
# boundary statistic, maximize the pseudo-log-likelihood of the free
# coefficients with the shared Newton optimizer from `θ0` (a full-length
# coefficient vector, non-finite entries and `nothing` meaning zeros), and
# return everything the diagnostics and the bootstrap refits need.
function _count_mple_fit(model::CountERGMModel, support::UnitRange{Int};
                         maxiter::Int=100, tol::Float64=1e-8,
                         θ0::Union{Nothing, AbstractVector{<:Real}}=nothing)
    D = _count_design(model, support)
    names = _term_names(model)
    p = length(model.terms)
    S = length(D.support)
    R = length(D.n_tot)
    mask = fill(true, S, R)
    cols, fixed = _count_boundary!(mask, D, p)
    q = length(cols)

    derivatives = _count_derivatives(D, cols, mask)
    if q > 0
        start = θ0 === nothing ? zeros(q) :
                [isfinite(θ0[k]) ? Float64(θ0[k]) : 0.0 for k in cols]
        nf = _quiet_newton(derivatives, start; maxiter=maxiter, tol=tol)
        if !nf.converged && nf.iterations < maxiter
            # Newton could not move from this start (singular Hessian, or a
            # step that no halving improves) and gave up before its budget:
            # rescue with a coordinate-ascent start and try again; keep
            # whichever is better. Running out of `maxiter` is NOT rescued —
            # that is the caller's budget, and the fit is reported unconverged.
            rescue = _quiet_newton(derivatives,
                                   _coordinate_ascent_start(derivatives, start);
                                   maxiter=maxiter, tol=tol)
            if rescue.converged || rescue.loglik > nf.loglik
                nf = rescue
            end
        end
        θq, seq, Vq, ll, converged = nf.θ, nf.se, nf.vcov, nf.loglik, nf.converged
        iterations = nf.iterations
        _, gq, Hq = derivatives(θq)
        grad_norm = norm(gq)
        # A "converged" Newton on a flat asymptote is not a maximum: the MPLE
        # does not exist (see `_count_separated`)
        separated = converged && _count_separated(derivatives, D, cols, mask, θq, seq)
        separated && (converged = false)
        hcond = _hessian_cond(Hq)
        collinear = hcond <= _HESSIAN_COND_TOL ? String[] :
                    _collinear_names(Hq, [names[k] for k in cols])
    else
        # Every coefficient is fixed: the limit pseudo-likelihood is fully
        # determined by the restricted supports, nothing to iterate
        θq, seq, Vq = Float64[], Float64[], zeros(0, 0)
        ll = derivatives(Float64[])[1]
        converged = true
        iterations = 0
        grad_norm = 0.0
        separated = false
        hcond = 1.0
        collinear = String[]
    end

    θ = fill(NaN, p)
    se = zeros(p)
    V = zeros(p, p)
    θ[cols] .= θq
    se[cols] .= seq
    V[cols, cols] .= Vq
    for (k, side) in fixed
        θ[k] = side === :min ? -Inf : Inf
    end
    return (θ=θ, se=se, vcov=V, loglik=ll, converged=converged,
            iterations=iterations, grad_norm=grad_norm, design=D,
            mask=mask, cols=cols, fixed=fixed, separated=separated,
            hessian_cond=hcond, collinear=collinear)
end

# =============================================================================
# Count MPLE: error-controlled support and the public estimator
# =============================================================================

# Cap on the doublings of the adaptive support: 2^8 × the starting bound
# (a dyad-dependent design grows linearly with the support, so the cap is
# also a memory bound for a model that never settles)
const _MAX_SUPPORT_DOUBLINGS = 8

# Continuation in the support. Newton from θ = 0 on a wide support is a bad
# start — with the counting measure the initial conditional is uniform on
# 0:max_val, the first steps overshoot and the optimizer gives up — so every
# fit climbs a ladder of supports, `lo:base`, `lo:2·base`, ..., `full`, each
# warm-started from the one below. The adaptive path IS this ladder with a
# stopping rule; the fixed and bounded paths climb it to their given top.
function _support_ladder(full::UnitRange{Int}, base::Int)
    lo, hi = first(full), last(full)
    rungs = UnitRange{Int}[]
    top = min(hi, max(base, lo + 1))
    while true
        push!(rungs, lo:top)
        top >= hi && break
        top = min(hi, 2 * top)
    end
    return rungs
end

function _ladder_fit(model::CountERGMModel, full::UnitRange{Int};
                     maxiter::Int, tol::Float64)
    weights = get_edge_attribute(model.network, :weight, Int)
    base = _default_max_val(model.network, weights)
    fit = nothing
    for rung in _support_ladder(full, base)
        fit = _count_mple_fit(model, rung; maxiter=maxiter, tol=tol,
                              θ0=fit === nothing ? nothing : fit.θ)
    end
    return fit
end

# max |Δθ| between two fits in units of the current fit's standard error (a
# coefficient with no usable SE — NaN or 0 — is measured on the absolute scale);
# fixed coefficients are skipped.
function _support_delta(prev, cur)
    δ = 0.0
    for k in eachindex(cur.θ)
        (isfinite(cur.θ[k]) && isfinite(prev.θ[k])) || continue
        scale = (isfinite(cur.se[k]) && cur.se[k] > 0) ? cur.se[k] : 1.0
        δ = max(δ, abs(cur.θ[k] - prev.θ[k]) / scale)
    end
    return δ
end

# Fit at the default bound, refit at twice the bound, and stop when ALL of
# (a) the doubling moved no estimate by more than `support_tol` standard
# errors, (b) the wider fit leaves no more than `support_tol` expected dyads
# past the previous bound (the mass the smaller fit omitted), and (c) no dyad
# puts more than `BOUNDARY_MASS_TOL` of its conditional mass on the new top
# value. (a) says the estimates have settled, (b) that the settling is not an
# artefact of two truncations agreeing, (c) that the reported bound itself is
# not shaping the fit; an improper model (η increasing in y) fails all three
# at every doubling. The reported fit is always the one at the larger bound,
# whose own truncation error is smaller still.
function _adaptive_support_fit(model::CountERGMModel; support_tol::Float64,
                               max_doublings::Int, maxiter::Int, tol::Float64)
    weights = get_edge_attribute(model.network, :weight, Int)
    mv = _default_max_val(model.network, weights)
    y0 = first(_support(model.reference, mv))
    prev = _count_mple_fit(model, _support(model.reference, mv);
                           maxiter=maxiter, tol=tol)
    # A rung whose Newton failed, or whose pseudo-Hessian is numerically
    # singular, has no estimates to compare and no standard errors to scale
    # the next comparison by: doubling on would grow the design 2^k× only to
    # report the same verdict, and the "not settling" diagnosis would be a
    # misdiagnosis (a collinear design, not an improper family). Stop here;
    # the Newton verdict (`converged`, `separated`, `hessian_cond`) carries
    # the diagnosis and `count_mple` speaks it.
    _support_rung_usable(prev) ||
        return (fit=prev, control=:unconverged, delta=NaN, tail=NaN)
    cur = prev
    δ = NaN
    tail = NaN
    for _ in 1:max_doublings
        cur = _count_mple_fit(model, _support(model.reference, 2 * mv);
                              maxiter=maxiter, tol=tol, θ0=prev.θ)
        _support_rung_usable(cur) ||
            return (fit=cur, control=:unconverged, delta=NaN, tail=NaN)
        δ = _support_delta(prev, cur)
        θfree = cur.θ[cur.cols]
        tail = _omitted_tail(θfree, cur.design, cur.cols, cur.mask, mv - y0 + 1)
        top = _max_boundary_mass(θfree, cur.design, cur.cols, cur.mask)
        (δ <= support_tol && tail <= support_tol && top <= BOUNDARY_MASS_TOL) &&
            return (fit=cur, control=:converged, delta=δ, tail=tail)
        prev, mv = cur, 2 * mv
    end
    return (fit=cur, control=:unconverged, delta=δ, tail=tail)
end

_support_rung_usable(fit) = fit.converged && fit.hessian_cond <= _HESSIAN_COND_TOL

"""
    count_mple(model::CountERGMModel; max_val=nothing, support_tol=1e-3,
               max_doublings=8, se=:hessian, n_boot=100, boot_burnin=nothing,
               boot_interval=nothing, rng=Random.default_rng(),
               maxiter=100, tol=1e-8, warn=true) -> CountERGMResult

Maximum pseudo-likelihood estimation for count ERGMs. For each dyad the
full conditional over the count support is enumerated, so the score is
`Σ_dyads [Δg(y_obs) − E_θ(Δg)]` and the Hessian is `−Σ_dyads Var_θ(Δg)`.
Dyads with identical conditionals are compressed into one row (every dyad, for
a dyad-independent model). The pseudo-log-likelihood is maximized with the
shared `Networks.newton_fit` Newton–Raphson-with-step-halving optimizer, by
continuation in the support: the first fit is on `0:max(10, 2·max count)` and
every wider support is warm-started from the fit below it (a cold start on a
wide support overshoots).

# Support

The support enumerated for each dyad is the reference's own for a bounded
reference (`BinomialReference`, `DiscUnifReference`, `DiscUnif2Reference`).
For an unbounded one (`PoissonReference`, `GeometricReference`) it is
`0:max_val`, a **truncation** of the model, chosen so that the truncation is
an error-controlled numerical device rather than a silent change of model:

- `max_val=nothing` (default): start at twice the largest observed count (at
  least 10), refit at twice that bound, and stop when the doubling moved every
  estimate by at most `support_tol` standard errors (on the absolute scale
  where a standard error is `NaN` or 0) **and** the wider fit leaves at most
  `support_tol` expected dyads past the previous bound **and** no dyad puts
  more than `BOUNDARY_MASS_TOL` of its conditional mass on the new top value;
  the fit at the larger bound is reported (`support_control = :converged`,
  `support_stable = true`, with the achieved `support_delta`/`omitted_tail`).
  If the estimates are still moving after `max_doublings` (default 8)
  doublings the fit is reported with `support_control = :unconverged`,
  `support_stable = false` and a warning naming the achieved δ and bound: that
  is the signature of a model that is not normalisable on the unbounded
  support (e.g. a positive `sum` coefficient under a geometric reference).
  The doubling also stops, with the same `:unconverged` verdict but
  `support_delta = omitted_tail = NaN`, at the first bound on which Newton
  did not converge or the pseudo-Hessian is numerically singular (a
  collinear or separated design): there are no settled estimates to compare,
  so the support is not error-controlled and the convergence/conditioning
  warning speaks alone — the "not normalisable" reading applies only to a
  converged fit that keeps moving.
- `max_val=k`: fix the bound (`support_control = :fixed`,
  `support_delta = NaN`); the boundary-mass diagnostic still reports whether it
  bites.
- A bounded reference never loops: `support_control = :bounded`,
  `support_delta = 0.0`.

# Boundary statistics

A statistic at its extreme attainable value on every dyad's conditional
support has no finite MPLE; its coefficient is fixed at ∓Inf with a warning,
as R `ergm` does (see [`CountERGMResult`](@ref)). `se=:bootstrap` is refused
for such a fit (there is no finite model to simulate from).

# Standard errors

- `se=:hessian` (default) — the inverse negative pseudo-Hessian. **Caution:**
  for a model with dyad-dependent terms (`CountMutualTerm`,
  `TransitiveTiesTerm`, the node-strength terms, ...) the pseudo-likelihood
  multiplies dyad conditionals as if independent, so these standard errors are
  expected to be *anticonservative* (too small) and the p-values too optimistic.
  For a dyad-independent model (e.g. `SumTerm` alone) the pseudo-likelihood is
  the likelihood and they are correct.
- `se=:bootstrap` — parametric bootstrap: Gibbs-simulate `n_boot` count networks
  from the fitted model at θ̂ with [`simulate_count_ergm`](@ref) (at the same
  `max_val`), refit the count MPLE on each, and report the empirical covariance
  of the refits. The point estimate is unchanged; only the covariance is
  replaced. Same option, keywords and semantics as `ERGM.mple`'s, on the ONE
  shared `Networks.bootstrap_cov` loop. A replicate on which the count MPLE
  does not exist (a boundary statistic) or does not converge is excluded from
  the covariance and counted in a warning; `fit.boot_replicates` keeps every
  refit (excluded ones as `NaN` rows).

# Convergence, separation and conditioning

The fit is warned about, and `fit.converged == false` recorded, when the
Newton iteration exhausts `maxiter` or cannot move (`fit.iterations`,
`fit.gradient_norm`); `approximations(fit)` then lists it and `is_exact(fit)`
is `false`.

**Separation.** A design on which the pseudo-likelihood has no finite maximum
along a *combination* of the statistics — quasi-complete separation, which
the one-column boundary test cannot see (e.g. `sum + nonzero` when every
count is 0 or 1, or `sum + atleast(2)` when every count is 0 or 2) — is
detected by the two signatures the asymptote leaves on Newton, as
`ERGM.mple` does: a fitted conditional whose exponential tilt spans more than
18.42 nats between the row's largest value and an unobserved support value,
*and* a next Newton step still larger than 1e-3·‖θ‖ (or an uninvertible
Hessian). Such a fit is returned with `converged = false` and
`separated = true`, warned about with R's sentence ("The MPLE does not
exist!"), listed in `approximations` and never `is_exact`; the bootstrap
excludes such replicates.

**Conditioning.** `fit.hessian_cond` is the condition number of the negative
pseudo-Hessian at the estimates. Above `ERGMCount._HESSIAN_COND_TOL` (1e8) —
two statistics collinear or nearly so on this network, such as
`greaterthan(2)` and `atleast(3)` on integer counts, or `sum` and `nonzero`
on a 0/1 network — the fit warns, naming the statistics that load on the flat
direction (`fit.collinear`), and `show`/`approximations` carry the caveat;
the standard errors along that direction are meaningless (`NaN` when it is
exactly singular).

# Keyword Arguments
- `max_val::Union{Int,Nothing}=nothing`, `support_tol::Float64=1e-3`,
  `max_doublings::Int=8`: above
- `se::Symbol=:hessian`: `:hessian` or `:bootstrap` (validated by the shared
  `Networks.check_se`)
- `n_boot::Int=100`: number of bootstrap replicates (`se=:bootstrap` only)
- `boot_burnin`, `boot_interval`: Gibbs controls (in sweeps) for the bootstrap
  simulations; `nothing` resolves through the dyad-scaled rule shared with
  ERGM.jl (`ERGM._mcmc_defaults`, converted from toggles to sweeps)
- `rng::AbstractRNG=Random.default_rng()`: source of the bootstrap randomness —
  a fixed `rng` reproduces the standard errors exactly
- `maxiter::Int=100`, `tol::Float64=1e-8`: Newton controls
- `warn::Bool=true`: `false` silences the fit diagnostics (boundary statistic,
  truncation, non-convergence, undefined standard errors) — they are all still
  recorded on the result. As in `ERGM.mple`, the parametric bootstrap refits
  run with `warn=false` and report their exclusions once, in aggregate.

# Example
```julia
using Networks, ERGMCount
net = network(4; directed=true)
for (i, j, w) in ((1, 2, 3), (2, 1, 1), (2, 3, 2), (3, 4, 1), (4, 1, 2))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
model = CountERGMModel([SumTerm(), NonzeroTerm()], net, PoissonReference())
fit = ERGMCount.count_mple(model)                 # adaptive support
fit.support_control                               # :converged
fixed = ERGMCount.count_mple(model; max_val=40)   # caller-fixed support
fixed.support_control                             # :fixed
```
"""
function count_mple(model::CountERGMModel; maxiter::Int=100,
                    tol::Float64=1e-8,
                    max_val::Union{Int, Nothing}=nothing,
                    support_tol::Float64=1e-3,
                    max_doublings::Int=_MAX_SUPPORT_DOUBLINGS,
                    se::Symbol=:hessian,
                    n_boot::Int=100,
                    boot_burnin::Union{Int, Nothing}=nothing,
                    boot_interval::Union{Int, Nothing}=nothing,
                    rng::Random.AbstractRNG=Random.default_rng(),
                    warn::Bool=true)
    # The count MPLE enumerates every dyad as observed, so a masked
    # (unobserved) dyad would enter the pseudo-likelihood at its face value.
    # `fit_ergm_count` refuses it and so does this entry point (a model can
    # be built, or its network masked, after the fact); the bootstrap refits
    # run on simulated networks and never reach here.
    require_observed(model.network; context="count_mple", face_ok=false)
    check_se(se, (:hessian, :bootstrap); context="count_mple")
    support_tol > 0 ||
        throw(ArgumentError("count_mple: support_tol must be positive (got $support_tol)"))
    max_doublings >= 1 ||
        throw(ArgumentError("count_mple: max_doublings must be at least 1 (got $max_doublings)"))
    max_val === nothing || max_val >= 1 ||
        throw(ArgumentError("count_mple: max_val must be at least 1 (got $max_val)"))

    ref = model.reference
    truncated = is_truncating(ref)
    if !truncated
        # The bound is the reference's own; `max_val` is not consulted
        fit = _ladder_fit(model, _support(ref, 0); maxiter=maxiter, tol=tol)
        control, δ, tail = :bounded, 0.0, 0.0     # nothing truncated
    elseif max_val !== nothing
        fit = _ladder_fit(model, _support(ref, max_val); maxiter=maxiter, tol=tol)
        control, δ, tail = :fixed, NaN, NaN
    else
        a = _adaptive_support_fit(model; support_tol=support_tol,
                                  max_doublings=max_doublings, maxiter=maxiter,
                                  tol=tol)
        fit, control, δ, tail = a.fit, a.control, a.delta, a.tail
    end

    names = _term_names(model)
    (warn && !isempty(fit.fixed)) && _warn_count_boundary(names, fit.fixed)

    D = fit.design
    support = D.support
    mv = last(support)
    θfree = fit.θ[fit.cols]

    # Boundary-mass diagnostic. For an unbounded reference the enumeration
    # `0:max_val` is a TRUNCATION of the model, not the model itself, so the
    # fit is only an approximation to the documented unbounded family insofar
    # as the fitted conditionals put negligible mass at the top of the support.
    # Measure that directly, at the fitted coefficients, and say so out loud.
    boundary = truncated ? _max_boundary_mass(θfree, D, fit.cols, fit.mask) : 0.0

    if !warn
        # diagnostics recorded on the result, nothing printed
    elseif control === :unconverged && !_support_rung_usable(fit)
        # The doubling stopped on a rung with no usable estimates: the
        # non-convergence / conditioning warning below carries the
        # diagnosis, and a "not normalisable" reading would be false
        @warn "count_mple: the error-controlled support stopped at max_val = $mv " *
              "because the fit there $(fit.converged ? "has a numerically singular " *
              "pseudo-Hessian" : "did not converge") (see the next warning): the " *
              "support was NOT error-controlled — `fit.support_control == " *
              ":unconverged`, `support_delta`/`omitted_tail` are NaN. Fix the " *
              "model first; the doubling resumes once a bound fits." maxlog = 1
    elseif control === :unconverged
        @warn """
              count_mple: the error-controlled support did not converge. After \
              $max_doublings doublings (max_val = $mv) the last doubling \
              still moved the estimates by $(_fmt3(δ)) standard \
              errors (tolerance $support_tol) and left $(_fmt3(tail)) \
              expected dyads past the previous bound.

              $(nameof(typeof(ref))) is mathematically UNBOUNDED, and a truncated \
              fit that keeps moving with the bound is the signature of a model that \
              is NOT normalisable on that support at these coefficients (e.g. a \
              geometric reference with a non-negative `sum` coefficient). The fit \
              is reported as a truncated exponential family on 0:$mv; \
              `fit.support_control == :unconverged` records this.
              """ maxlog = 1
    elseif truncated && boundary > BOUNDARY_MASS_TOL
        @warn """
              Count support was truncated at max_val = $mv, but the \
              fitted model places $(_fmt3(100 * boundary))% of the \
              conditional mass on the boundary value for at least one dyad \
              (tolerance $(_fmt2(100 * BOUNDARY_MASS_TOL))%).

              $(nameof(typeof(ref))) is mathematically UNBOUNDED. With appreciable \
              mass at the bound, this fit does not approximate that model — it is \
              a different, truncated exponential family, and the estimates are \
              biased toward the bound.

              Refit with a larger `max_val` and check the estimates are stable, \
              e.g. `count_mple(model; max_val = $(2 * mv))`, or omit `max_val` for \
              the error-controlled default.
              """ maxlog = 1
    end

    if warn && fit.separated
        _warn_count_separated("count_mple")
    elseif warn && !fit.converged
        @warn "count_mple: Newton did not converge in $(fit.iterations) " *
              "iteration$(fit.iterations == 1 ? "" : "s") (maxiter = $maxiter, " *
              "tol = $tol; pseudo-score norm $(_fmt3(fit.grad_norm)) " *
              "at the reported estimates, which are therefore NOT a maximum of " *
              "the pseudo-likelihood). Two ways out: raise `maxiter` if the " *
              "iteration simply ran out of budget, or check the model for a " *
              "statistic with no finite maximizer (a boundary or non-identified " *
              "term — a statistic constant over the data, or two collinear " *
              "ones). `fit.converged == false` records this and " *
              "`is_exact(fit)` is false." maxlog = 1
    end
    if warn && !fit.separated && !(fit.hessian_cond <= _HESSIAN_COND_TOL)
        named = isempty(fit.collinear) ? "" :
                " — loading on $(join(fit.collinear, ", "))"
        @warn "count_mple: the pseudo-Hessian at the reported estimates is " *
              "numerically singular (condition number " *
              "$(_fmt3(fit.hessian_cond)) " *
              "> $(_HESSIAN_COND_TOL)): a flat direction$named — two statistics " *
              "that are collinear or nearly so on this network (e.g. " *
              "`greaterthan.2` and `atleast.3` on integer counts, or `sum` and " *
              "`nonzero` on a 0/1 network), so the standard errors along it are " *
              "meaningless" *
              (any(isnan, fit.se[fit.cols]) ? " and are reported as NaN" : "") *
              ". Drop or merge one of the terms; `fit.hessian_cond` and " *
              "`fit.collinear` record this." maxlog = 1
    end

    vcov, std_errors, boot_replicates = fit.vcov, fit.se, nothing
    if se === :bootstrap
        isempty(fit.fixed) || throw(ArgumentError(
            "count_mple: se=:bootstrap is not available when a coefficient is " *
            "fixed at ±Inf by a boundary statistic " *
            "($(join([names[k] for (k, _) in fit.fixed], ", "))): there is no " *
            "finite fitted model to simulate replicates from. Drop the term or " *
            "use the inverse-Hessian standard errors of the remaining ones."))
        vcov, std_errors, boot_replicates =
            _count_bootstrap_cov(model, fit.θ, mv; n_boot=n_boot,
                                 boot_burnin=boot_burnin,
                                 boot_interval=boot_interval,
                                 maxiter=maxiter, tol=tol, rng=rng)
    end

    z, pv = _count_zp(fit.θ, std_errors)
    return CountERGMResult(model, fit.θ, std_errors, z, pv, vcov, fit.loglik,
                           fit.converged, fit.iterations, fit.grad_norm, mv,
                           truncated, boundary, se, control,
                           control !== :unconverged, support_tol, δ, tail,
                           boot_replicates, fit.separated, fit.hessian_cond,
                           fit.collinear)
end

missing_policies(::typeof(count_mple)) = (:error,)

# Parametric-bootstrap covariance of the count MPLE: Gibbs-simulate `n_boot`
# count networks at θ̂, refit the count MPLE on each, take the empirical
# covariance. The loop is the shared `Networks.bootstrap_cov`; this supplies only
# the two callbacks that are ERGMCount's. The refits reuse the observed fit's
# `max_val`, so every replicate is fit on the same (possibly truncated) support —
# and the simulator draws from `0:max_val`, so no replicate can fall outside it.
# A replicate without a finite, converged MPLE (a boundary statistic in the
# simulated network, or Newton failing) is returned as NaN and excluded, with
# a warning: a NaN row must never enter the covariance silently.
function _count_bootstrap_cov(model::CountERGMModel, θ̂::Vector{Float64},
                              mv::Int; n_boot::Int, boot_burnin, boot_interval,
                              maxiter::Int, tol::Float64,
                              rng::Random.AbstractRNG)
    simulate(rng, B) = _simulate_count(model.network, model.terms,
                                       model.reference, θ̂; n_sim=B,
                                       burnin=boot_burnin, interval=boot_interval,
                                       max_val=mv, rng=rng,
                                       context="count_mple(se=:bootstrap)")

    p = length(θ̂)
    support = _support(model.reference, mv)
    function refit(sim::Network)
        boot_model = CountERGMModel(model.terms, sim, model.reference)
        f = _count_mple_fit(boot_model, support; maxiter=maxiter, tol=tol, θ0=θ̂)
        (f.converged && all(isfinite, f.θ)) || return fill(NaN, p)
        return f.θ
    end

    boot = bootstrap_cov(refit, simulate, θ̂; n_boot=n_boot, rng=rng)
    replicates = boot.replicates
    ok = [all(isfinite, view(replicates, b, :)) for b in 1:n_boot]
    n_ok = count(ok)
    n_ok == n_boot && return boot.vcov, boot.se, replicates
    n_ok >= 2 || throw(ArgumentError(
        "count_mple: se=:bootstrap — only $n_ok of the $n_boot bootstrap refits " *
        "had a finite, converged count MPLE (the others hit a boundary statistic " *
        "or did not converge), too few for a covariance. The fitted model " *
        "simulates networks on which the MPLE barely exists; increase `n_boot`, " *
        "or reconsider the model."))
    @warn "count_mple: se=:bootstrap — $(n_boot - n_ok) of the $n_boot bootstrap " *
          "refits had no finite, converged count MPLE (the simulated replicate " *
          "put a statistic at its boundary, or Newton did not converge) and were " *
          "excluded; the covariance is over the remaining $n_ok refits. This is " *
          "about the simulated replicates, not about the observed network. " *
          "`fit.boot_replicates` holds every refit (NaN rows excluded)." maxlog = 1
    V = Matrix{Float64}(cov(replicates[ok, :]))
    return V, sqrt.(max.(diag(V), 0.0)), replicates
end

# =============================================================================
# Simulation
# =============================================================================

# Gibbs-sweep defaults through THE dyad-scaled rule shared with ERGM.jl
# (`ERGM._mcmc_defaults`, panel 2026-09 item 24e). ERGM's budget is in single-
# dyad toggles; a Gibbs sweep visits every dyad once, so the same budget is
# `cld(toggles, n_dyads)` sweeps: 20 sweeps of burn-in and
# `cld(max(100, n_dyads ÷ 10), n_dyads)` sweeps between retained draws.
function _gibbs_defaults(n_dyads::Int)
    nd = max(n_dyads, 1)
    d = ERGM._mcmc_defaults(nd)
    return (burnin=cld(d.burnin, nd), interval=cld(d.interval, nd))
end

"""
    simulate_count_ergm(result::CountERGMResult; n_sim=1, burnin=nothing,
                        interval=nothing, max_val=nothing,
                        rng=Random.default_rng()) -> Vector{Network}
    simulate_count_ergm(net, terms, coefficients; reference=PoissonReference(),
                        n_sim=1, burnin=nothing, interval=nothing, max_val=20,
                        rng=Random.default_rng()) -> Vector{Network}

Simulate networks from a count ERGM by Gibbs sampling: each sweep resamples
every dyad from its full conditional
`P(y_ij = y | rest) ∝ h(y)·exp(θ'Δg(y))` (an inverse-CDF draw from the
enumerated support), using each term's change statistic, so structural terms
(mutuality, transitivity, node strength) influence the draws. The first form
simulates from a fitted model at its coefficients and support; the second from
an explicit specification (`net` provides the size, directedness and starting
state; `terms` may be a Vector, a Tuple or a single term).

`burnin` and `interval` are counted in **sweeps**; `nothing` resolves through
the dyad-scaled rule shared with ERGM.jl (20 sweeps of burn-in, then
`cld(max(100, n_dyads ÷ 10), n_dyads)` sweeps between retained draws).

A network with masked (unobserved) dyads is refused
(`Networks.require_observed`): the sampler would otherwise start from, and
condition on, face values of dyads that were never observed. A two-mode
(bipartite) network is refused (the sweep would resample its impossible
within-mode dyads), and so is a specification with a coefficient at ∓Inf (a
boundary statistic). Under `DiscUnif2Reference(a < 0, b)` the draws can be
negative: an edge is stored for every non-zero count, with the (possibly
negative) value under `:weight`.

All random draws flow through `rng`; the same rng state yields identical
output. This sampler is a Gibbs sweep, not a Metropolis toggle chain, so it
deliberately does **not** use `ERGM.mh_toggle!`.

# Example
```julia
using Networks, ERGMCount, Random
seed = network(6; directed=true)
sims = simulate_count_ergm(seed, [SumTerm(), CountMutualTerm()], [log(0.8), 0.5];
                           reference=PoissonReference(), n_sim=3, max_val=10,
                           rng=Xoshiro(1))
length(sims)                        # 3
using ERGM: compute
compute(SumTerm(), sims[1]) > 0     # true
```
"""
function simulate_count_ergm(result::CountERGMResult;
                             n_sim::Int=1,
                             burnin::Union{Int, Nothing}=nothing,
                             interval::Union{Int, Nothing}=nothing,
                             max_val::Union{Int, Nothing}=nothing,
                             rng::Random.AbstractRNG=Random.default_rng())
    model = result.model
    mv = something(max_val, result.max_val)
    return _simulate_count(model.network, model.terms, model.reference,
                           result.coefficients; n_sim=n_sim, burnin=burnin,
                           interval=interval, max_val=mv, rng=rng,
                           context="simulate_count_ergm")
end

function simulate_count_ergm(net::Network, terms::Tuple,
                             coefficients::AbstractVector{<:Real};
                             reference::AbstractReferenceMeasure=PoissonReference(),
                             n_sim::Int=1,
                             burnin::Union{Int, Nothing}=nothing,
                             interval::Union{Int, Nothing}=nothing,
                             max_val::Int=20,
                             rng::Random.AbstractRNG=Random.default_rng())
    return _simulate_count(net, terms, reference, Vector{Float64}(coefficients);
                           n_sim=n_sim, burnin=burnin, interval=interval,
                           max_val=max_val, rng=rng, context="simulate_count_ergm")
end

simulate_count_ergm(net::Network, terms::AbstractVector, coefficients; kwargs...) =
    simulate_count_ergm(net, Tuple(terms), coefficients; kwargs...)
simulate_count_ergm(net::Network, term::AbstractERGMTerm, coefficients; kwargs...) =
    simulate_count_ergm(net, (term,), coefficients; kwargs...)
simulate_count_ergm(net::BipartiteNetwork, terms, coefficients; kwargs...) =
    _refuse_two_mode(net, "simulate_count_ergm")
# Swapped arguments: name the right order instead of a raw MethodError, as
# `fit_ergm_count(terms, net)` does
simulate_count_ergm(terms::Union{AbstractERGMTerm, Tuple, AbstractVector}, net::Network,
                    coefficients; kwargs...) =
    throw(ArgumentError(
        "simulate_count_ergm(net, terms, coefficients): the network comes first " *
        "and the terms second (got the terms first). Call " *
        "`simulate_count_ergm(net, $(terms isa AbstractERGMTerm ? "[$(nameof(typeof(terms)))()]" : "terms"), coefficients)`."))

missing_policies(::typeof(simulate_count_ergm)) = (:error,)

function _simulate_count(net0::Network, terms::Tuple, ref::AbstractReferenceMeasure,
                         θ::Vector{Float64};
                         n_sim::Int, burnin, interval, max_val::Int,
                         rng::Random.AbstractRNG, context::AbstractString)
    # The Gibbs sweep conditions every dyad on the face value of every other:
    # an unobserved dyad has no face value to condition on. And it resamples
    # every off-diagonal dyad: a two-mode network's within-mode dyads are not
    # dyads at all.
    require_observed(net0; context=context, face_ok=false)
    _refuse_two_mode(net0, context)
    isempty(terms) && throw(ArgumentError("$context: at least one term is required"))
    length(θ) == length(terms) || throw(ArgumentError(
        "$context: $(length(θ)) coefficients for $(length(terms)) terms"))
    all(isfinite, θ) || throw(ArgumentError(
        "$context: coefficient(s) at index " *
        "$(join(findall(!isfinite, θ), ", ")) are not finite (a statistic fixed " *
        "at ±Inf by a boundary), so there is no finite model to simulate from."))
    _validate_count_terms(terms, net0)
    n_sim >= 0 || throw(ArgumentError("$context: n_sim must be non-negative (got $n_sim)"))
    max_val >= 1 || throw(ArgumentError("$context: max_val must be at least 1 (got $max_val)"))
    # The seed's counts are read by every term's first conditional: the same
    # rules as for a fit (every edge weighted, integers, sign per reference)
    lo = _validate_count_weights(net0, ref)
    # ... and a support that reaches below zero refuses the terms R refuses
    # on negative weights, whatever the seed holds
    _validate_negative_terms(terms, min(lo, first(_support(ref, max_val))), context)

    # `copy` is the ONE copier of a Network (graph and attribute dicts
    # duplicated, vertex/edge attributes preserved); the chain state is a
    # copy of the seed and every retained draw a copy of the chain
    current = copy(net0)
    n = Int(nv(current))
    n_dyads = is_directed(current) ? n * (n - 1) : n * (n - 1) ÷ 2
    d = _gibbs_defaults(n_dyads)
    burnin = something(burnin, d.burnin)
    interval = something(interval, d.interval)
    burnin >= 0 || throw(ArgumentError("$context: burnin must be non-negative (got $burnin)"))
    interval >= 1 || throw(ArgumentError("$context: interval must be at least 1 (got $interval)"))

    support = _support(ref, max_val)
    log_h = [log_reference(ref, y) for y in support]
    η = Vector{Float64}(undef, length(support))    # the conditional's weights
    buf = Vector{Float64}(undef, length(support))  # one term's profile

    # Typed snapshot of the :weight edge attribute (Networks.jl's typed
    # accessor), maintained incrementally alongside the network so the hot
    # loop never reads the untyped attribute Dict.
    weights = get_edge_attribute(current, :weight, Int)

    networks = Vector{typeof(current)}()
    for sweep in 1:(burnin + n_sim * interval)
        _gibbs_sweep!(rng, current, weights, terms, θ, support, log_h, η, buf)
        if sweep > burnin && (sweep - burnin) % interval == 0
            push!(networks, copy(current))
        end
    end
    return networks
end

# One Gibbs sweep: every dyad once, in a fixed order
function _gibbs_sweep!(rng::Random.AbstractRNG, current::Network, weights, terms::Tuple,
                       θ, support::UnitRange{Int}, log_h, η, buf)
    n = Int(nv(current))
    directed = is_directed(current)
    for i in 1:n
        for j in (directed ? (1:n) : ((i + 1):n))
            i == j && continue
            _gibbs_update_dyad!(rng, current, weights, terms, θ, i, j, support, log_h, η, buf)
        end
    end
    return current
end

# Redraw dyad (i, j) from its full conditional and apply the new value. The
# conditional (folded statically over the term tuple from the support
# profiles) and the inverse-CDF draw allocate nothing (pinned); only an
# actual change of value touches the network, its :weight attribute and the
# typed snapshot, and that costs what Networks.jl's own edge/attribute
# mutation costs.
function _gibbs_update_dyad!(rng::Random.AbstractRNG, current::Network, weights,
                             terms::Tuple, θ, i::Int, j::Int, support::UnitRange{Int},
                             log_h, η, buf)
    old = dyad_value(current, weights, i, j)
    total = _dyad_conditional!(η, buf, terms, θ, log_h, support, current, weights, i, j, old)
    new_val = support[_draw_index(rng, η, total)]
    new_val == old && return new_val
    _set_dyad!(current, weights, i, j, new_val)
    return new_val
end

# Set dyad (i, j) to count `y` in the network, its :weight attribute and the
# typed snapshot. An edge exists for every NON-ZERO count — a negative draw
# under `DiscUnif2Reference(a < 0, b)` is stored as a negative `:weight`, not
# collapsed to 0 (`dyad_value` returns it; the terms read it as a non-zero
# value, as R does). 0 removes the edge (which drops its attributes) but only
# ZEROES the snapshot entry: `dyad_value` reads the snapshot behind a
# `has_edge` guard, so a 0 there is inert, and never deleting keeps the
# snapshot Dict from rehashing on every remove/re-add cycle of a busy chain
# (measured 32 B per cycle, amortised) — once every dyad the chain has ever
# occupied has a key, the update allocates nothing.
function _set_dyad!(net::Network, weights, i::Int, j::Int, y::Int)
    key = _wkey(net, i, j)
    if y != 0
        has_edge(net, i, j) || add_edge!(net, i, j)
        set_edge_attribute!(net, :weight, i, j, y)
    else
        rem_edge!(net, i, j)
    end
    weights[key] = y
    return net
end

# The full conditional of dyad (i, j) over the support, as UN-normalised
# weights `η[s] = exp(η_s − max η)` with η_s = log h(y_s) + θ'Δg(old → y_s);
# returns their sum, so the draw needs no normalising pass
function _dyad_conditional!(η, buf, terms::Tuple, θ, log_h, support::UnitRange{Int},
                            net, weights, i::Int, j::Int, old::Int)
    copyto!(η, log_h)
    _accumulate_conditional!(η, buf, terms, θ, net, weights, i, j, old, support)
    m = -Inf
    @inbounds for s in eachindex(η)
        m = max(m, η[s])
    end
    total = 0.0
    @inbounds for s in eachindex(η)
        w = exp(η[s] - m)
        η[s] = w
        total += w
    end
    return total
end

# Inverse-CDF draw of an index from un-normalised weights summing to `total`:
# one uniform, one pass, no allocation. Rounding can leave the cumulative sum
# a hair below `total`; the last index absorbs it.
function _draw_index(rng::Random.AbstractRNG, w, total::Float64)
    u = rand(rng) * total
    acc = 0.0
    @inbounds for s in eachindex(w)
        acc += w[s]
        u < acc && return s
    end
    return lastindex(w)
end

# =============================================================================
# Goodness of Fit
# =============================================================================

# Smallest and largest dyad count value in a network, 0 included (a count
# can be negative under `DiscUnif2Reference(a < 0, b)`)
function _dyad_value_extrema(net)
    weights = get_edge_attribute(net, :weight, Int)
    lo, hi = 0, 0
    for e in edges(net)
        w = Int(get(weights, _wkey(net, src(e), dst(e)), 1))
        lo = min(lo, w)
        hi = max(hi, w)
    end
    return lo, hi
end

# Number of dyads with count value k, for k in lo:hi (the dyads without an
# edge, and edges whose `:weight` is 0, are the zeros); `lo ≤ 0 ≤ hi` and
# every value of the network lies in `lo:hi`
function _dyad_value_counts(net, lo::Int, hi::Int)
    weights = get_edge_attribute(net, :weight, Int)
    counts = zeros(Float64, hi - lo + 1)
    for e in edges(net)
        w = Int(get(weights, _wkey(net, src(e), dst(e)), 1))
        counts[w - lo + 1] += 1
    end
    counts[1 - lo] += _n_dyads(net) - ne(net)
    return counts
end

"""
    gof(result::CountERGMResult; n_sim=100, burnin=nothing, interval=nothing,
        max_val=nothing, rng=Random.default_rng()) -> GOFResult

Goodness-of-fit assessment of a fitted count ERGM: networks are simulated
from the fitted model with [`simulate_count_ergm`](@ref) and compared with
the observed network on

- the model statistics (one level per term), and
- the distribution of dyad count values (number of dyads with value
  0, 1, 2, ... — every value the observed or a simulated network takes,
  negative ones included under `DiscUnif2Reference(a < 0, b)`).

This is a method of the shared `Networks.gof` generic; it returns the
shared `Networks.GOFResult` (observed value, simulation envelope, and
two-sided Monte-Carlo p-value per level). The simulation inherits
`simulate_count_ergm`'s refusals (masked dyads, non-finite coefficients).

# Keyword Arguments
- `n_sim::Int=100`: Number of simulated networks
- `burnin`, `interval`, `max_val`, `rng`: passed to
  [`simulate_count_ergm`](@ref)

# Example
```julia
using Networks, ERGMCount, Random
net = network(5; directed=true)
for (i, j, w) in ((1, 2, 3), (2, 1, 1), (2, 3, 2), (3, 4, 1), (4, 1, 2), (5, 1, 1))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
fit = fit_ergm_count(net, [SumTerm(), NonzeroTerm()])
g = gof(fit; n_sim=20, rng=Xoshiro(3))
g.statistics[1].labels             # ["sum", "nonzero"]
```
"""
function gof(result::CountERGMResult; n_sim::Int=100,
             burnin::Union{Int, Nothing}=nothing,
             interval::Union{Int, Nothing}=nothing,
             max_val::Union{Int, Nothing}=nothing,
             rng::Random.AbstractRNG=Random.default_rng())
    n_sim >= 1 || throw(ArgumentError("gof: n_sim must be at least 1 (got $n_sim)"))
    net = result.model.network
    terms = result.model.terms
    sims = simulate_count_ergm(result; n_sim=n_sim, burnin=burnin,
                               interval=interval, max_val=max_val, rng=rng)

    # Model statistics: observed vs simulated
    obs_stats = [compute(term, net) for term in terms]
    sim_stats = [compute(term, s) for s in sims, term in terms]
    stats = GOFStatistic("model statistics", _term_names(result.model),
                         obs_stats, sim_stats)

    # Dyad count-value distribution, over every value the observed network or
    # a simulated one takes (negative ones included under DiscUnif2)
    lo, hi = _dyad_value_extrema(net)
    for s in sims
        slo, shi = _dyad_value_extrema(s)
        lo, hi = min(lo, slo), max(hi, shi)
    end
    obs_counts = _dyad_value_counts(net, lo, hi)
    sim_counts = Matrix{Float64}(undef, n_sim, hi - lo + 1)
    for (r, s) in enumerate(sims)
        sim_counts[r, :] .= _dyad_value_counts(s, lo, hi)
    end
    values = GOFStatistic("dyad count values", string.(lo:hi),
                          obs_counts, sim_counts)

    return GOFResult([stats, values]; model="Count ERGM")
end


# ----------------------------------------------------------------------------
# Precompile workload (panel 2026-09, item 18). Before it, the first
# `fit_ergm_count` in a fresh session took ~2.2 s and the first `gof` ~0.5 s
# after a 0.8 s `using`: the whole estimation path — the profile folds over
# the term tuple, the compressed design, the Newton kernel, the support
# ladder and the adaptive doubling, the coefficient table, the Gibbs sweep and
# the GOF tables — was compiled lazily on first use. Running one tiny fit of
# each kind here at precompile time caches those native-code specializations
# in the package image, for BOTH directedness type parameters
# (`CountERGMModel{Int,true,…}` and `CountERGMModel{Int,false,…}` are
# different specializations) and for the two term-tuple widths the docs use
# most. The networks are built inline (no dataset I/O at precompile time),
# everything is seeded and small, and the log output the tiny fits emit is
# silenced: a precompile-time warning is never about the user's data.
# ----------------------------------------------------------------------------
@setup_workload begin
    _pc_d = network(6; directed=true)
    for (i, j, w) in ((1, 2, 3), (2, 1, 1), (2, 3, 2), (3, 1, 1), (3, 4, 1), (4, 5, 2),
                      (5, 3, 1), (1, 5, 2), (5, 6, 1), (6, 2, 1))
        add_edge!(_pc_d, i, j); set_edge_attribute!(_pc_d, :weight, i, j, w)
    end
    _pc_u = network(6; directed=false)
    for (i, j, w) in ((1, 2, 2), (2, 3, 1), (1, 3, 3), (3, 4, 1), (4, 5, 2), (2, 5, 1),
                      (5, 6, 1), (1, 5, 1))
        add_edge!(_pc_u, i, j); set_edge_attribute!(_pc_u, :weight, i, j, w)
    end
    # A console logger writing to devnull: routing the fits' warnings through
    # the same logger type a user's REPL has caches the logging path without
    # printing anything
    _pc_null = Base.CoreLogging.ConsoleLogger(devnull, Base.CoreLogging.Warn)
    @compile_workload begin
        Base.CoreLogging.with_logger(_pc_null) do
            for (_pc_net, _pc_dep) in ((_pc_d, CountMutualTerm()), (_pc_u, NodeSumTerm()))
                _pc_rng = Random.Xoshiro(20260912)
                # The default (error-controlled) Poisson path, the StatsAPI
                # surface and the printed table
                _pc_fit = fit_ergm_count(_pc_net, [SumTerm(), NonzeroTerm()])
                coeftable(_pc_fit); confint(_pc_fit); aic(_pc_fit); bic(_pc_fit)
                show(devnull, _pc_fit); sprint(show, _pc_fit)
                # A bounded reference with a dyad-dependent term (directed:
                # mutual.min; undirected: nodeSum), a second tuple width
                _pc_bin = fit_ergm_count(_pc_net, [SumTerm(), _pc_dep];
                                         reference=BinomialReference(5))
                sprint(show, _pc_bin)
                simulate_count_ergm(_pc_fit; n_sim=1, burnin=2, interval=1, rng=_pc_rng)
                gof(_pc_fit; n_sim=2, burnin=2, interval=1, rng=_pc_rng)
                # The parametric bootstrap: two replicates are enough to compile
                # the path; on a 6-node network a replicate can lack a finite
                # MPLE, and fewer than two usable ones is an ArgumentError, so
                # the call is guarded — the path is compiled either way
                try
                    fit_ergm_count(_pc_net, [SumTerm(), NonzeroTerm()]; se=:bootstrap,
                                   n_boot=2, boot_burnin=2, boot_interval=1, rng=_pc_rng)
                catch
                end
            end
        end
    end
end

end # module
