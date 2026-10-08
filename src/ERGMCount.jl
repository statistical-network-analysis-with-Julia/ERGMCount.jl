"""
    ERGMCount.jl - ERGMs for Count-Valued Networks

Extends ERGM to handle networks with integer-valued edge weights,
using Poisson, geometric, or binomial reference measures.

The general form is:
    P(Y=y) ∝ h(y) exp(θ' g(y))

where h(y) is the reference measure determining the baseline distribution
for count-valued edges (Krivitsky 2012).

Port of the R ergm.count package from the StatNet collection. Estimation
([`fit_ergm_count`](@ref)) follows R's rule (`method=:auto`): maximum
pseudo-likelihood over an error-controlled count support for a
dyad-independent formula (there the exact MLE), Monte-Carlo maximum
likelihood (ergm.count's estimator) otherwise; `method=:mple`/`:mcmle` choose
one. Simulation is a Gibbs sweep over each dyad's full conditional
([`simulate_count_ergm`](@ref)).
"""
module ERGMCount

using Distributions   # Poisson/Geometric/Binomial reference draws, Normal quantile
using ERGM
using Graphs
using LinearAlgebra
using Logging: NullLogger, with_logger
using NetworkCore
using PrecompileTools: @setup_workload, @compile_workload
using Printf: @sprintf
using Random
using Statistics: cov

# The shared numerics and validators live in NetworkCore.jl: the ONE Newton
# optimizer, the ONE floored z → p helper, the ONE `se=` validator and the
# generic coefficient table every `coeftable` returns.
import NetworkCore: newton_fit, z_pvalues, check_se, CoefficientTable
# The ONE separation verdict and its policy helpers (NetworkCore.jl
# `src/separation.jl`): the count pseudo-likelihood is a conditional logit over
# each dyad class's support, decided exactly by `clogit_separation`
import NetworkCore: clogit_separation, warn_separation, separation_caveat,
                    SeparationVerdict
# The statistic protocol and the dependence/directedness traits are ERGM.jl's
# generics: ERGMCount adds methods for its own terms and model type, never a
# same-named private.
import ERGM: name, compute, is_dyad_dependent, has_dyad_dependent,
             requires_directed
# Shared presentation infrastructure (NetworkCore.jl): the ONE `gof` generic all
# model packages extend, plus the GOF containers
import NetworkCore: gof, GOFStatistic, GOFResult

# The shared result-metadata protocol (NetworkCore.jl `src/results.jl`): the
# generic accessors that say what a fit actually did. Imported by name because
# ERGMCount adds methods for `CountERGMResult`; `fit_metadata(fit)` collects them.
import NetworkCore: estimand, objective, is_exact, se_method, missing_method,
                 approximations, missing_policies
import StatsAPI
import StatsAPI: coef, stderror, vcov, confint, loglikelihood, nobs, dof, aic,
                 bic, coeftable, coefnames

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
export CountAtleastnTerm, AtmostTerm, CMPTerm
export TransitiveWeightsTerm, CyclicalWeightsTerm
export SmallerthanTerm, EqualToTerm, InIntervalTerm
# Valued covariate terms (ergm's form="sum" / form="nonzero")
export CountNodeMatchTerm, CountNodeFactorTerm, CountAbsDiffTerm
export CountNodeCovTerm, CountNodeOCovTerm, CountNodeICovTerm, CountEdgeCovTerm
export change_stat_count, dyad_value
public change_stats_support!

# Model / estimation (`count_mple` is exported as ERGM.jl exports `mple`: the
# estimator a warning or a docs page names must be callable as written)
export CountERGMModel, CountERGMResult
export fit_ergm_count, ergm_count, count_mple, count_mcmle
export has_dyad_dependent

# Simulation
export simulate_count_ergm

# Goodness of fit (method of the shared NetworkCore.jl `gof` generic)
export gof

# The full StatsAPI surface (re-exported so `coef(fit)` etc. work with just
# `using ERGMCount`; the bindings are StatsAPI's, so co-loading with ERGM.jl
# or REM.jl leaves every verb defined)
export coef, stderror, vcov, confint, loglikelihood, nobs, dof, aic, bic,
       coeftable, coefnames

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
using NetworkCore, ERGMCount
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
using NetworkCore, ERGMCount
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

# Canonical edge-attribute key: (i,j) directed, (min,max) undirected, in the
# network's own vertex type — the key type of the typed weight dictionary
# (`Dict{Tuple{T,T},Int}`), whatever integer type the caller's indices have
# (`edges(net)` yields `T`, the sweeps `Int`)
function _wkey(net::Network{T}, i::Integer, j::Integer) where {T}
    a, b = T(i), T(j)
    return is_directed(net) ? (a, b) : minmax(a, b)
end

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
using NetworkCore, ERGMCount
net = network(3; directed=false)
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 4)
w = get_edge_attribute(net, :weight, Int)
dyad_value(net, w, 1, 2)   # 4
dyad_value(net, w, 2, 1)   # 4 — the undirected key is minmax(2, 1)
dyad_value(net, w, 1, 3)   # 0 — no edge
```
"""
function dyad_value(net, weights, i::Integer, j::Integer)
    has_edge(net, i, j) || return 0
    return Int(get(weights, _wkey(net, i, j), 1))
end

# The typed snapshot (`Dict{Tuple{T,T},Int}`): `compute` on the untyped
# `Dict{…,Any}` inferred `Any` for five terms and boxed every `+`. Integer-
# valued Floats convert; a non-integer weight is refused by the model
# constructor before any estimator reads it.
_get_weights(net) = get_edge_attribute(net, :weight, Int)

# Dyad counts at face value, for `compute`: every dyad (the observed ones
# `ERGM.Extension.n_observed_dyads` counts plus the masked ones), and those without an edge.
# `compute` reads every count term at face value, so a masked dyad without an
# edge is a zero here as everywhere else in `compute`; the fits and the
# sampler refuse masked networks before any statistic is read.
_n_face_dyads(net::Network) = ERGM.Extension.n_observed_dyads(net) + n_missing_dyads(net)
_n_empty_dyads(net::Network) = _n_face_dyads(net) - ne(net)

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
using NetworkCore, ERGMCount
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
    SumTerm(; pow=1) <: AbstractERGMTerm

Sum of edge values: ∑_{i,j} y_{ij} — R's `sum`. This is the natural
sufficient statistic for the Poisson reference. `SumTerm(pow=p)` is R's
`sum(pow=p)`, ∑_{i,j} y_{ij}^p (labelled `sum<p>` as R prints it: `sum2`,
`sum0.5`); `p` must be positive, and a non-integer power is refused on a
network or support with negative counts (R returns `NaN` there).
Dyad-independent; pinned against `ergm` 4.12 by `test/fixtures/count_terms.toml`.

# Example
```julia
using NetworkCore, ERGMCount
using ERGM: compute, name
net = network(3; directed=true)
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 3)
add_edge!(net, 2, 3); set_edge_attribute!(net, :weight, 2, 3, 2)
compute(SumTerm(), net)        # 5.0
name(SumTerm())                # "sum"
compute(SumTerm(pow=2), net)   # 13.0
name(SumTerm(pow=2))           # "sum2"
```
"""
struct SumTerm <: AbstractERGMTerm
    pow::Float64
    function SumTerm(pow::Real)
        (isfinite(pow) && pow > 0) || throw(ArgumentError(
            "SumTerm: pow must be a positive finite number (got $pow)"))
        return new(Float64(pow))
    end
end
SumTerm(; pow::Real=1) = SumTerm(pow)

name(t::SumTerm) = t.pow == 1 ? "sum" : "sum" * _rnum(t.pow)

@inline _pow(t::SumTerm, y::Int) = t.pow == 1 ? Float64(y) : Float64(y)^t.pow

function compute(t::SumTerm, net)
    isinteger(t.pow) || _refuse_negative_network(t, net, "compute")
    weights = _get_weights(net)
    total = 0.0
    for e in edges(net)
        total += _pow(t, Int(get(weights, _wkey(net, src(e), dst(e)), 1)))
    end
    return total
end

change_stat_count(t::SumTerm, net, weights, i::Integer, j::Integer, old::Int, new::Int) =
    _pow(t, new) - _pow(t, old)

"""
    CMPTerm() <: AbstractERGMTerm

`ergm.count`'s `CMP` term: ∑_{i,j} log(y_{ij}!). Added to a Poisson-reference
model it turns each dyad's law into the Conway–Maxwell–Poisson family (a
coefficient θ multiplies the reference `1/y!` by `(y!)^θ`: under-dispersion
for θ < 0, over-dispersion for 0 < θ < 1, the geometric law at θ = 1). In
`ergm.count` this is a **term**, not a reference measure, and so it is here.
R label `CMP`; dyad-independent; refused on negative counts (R returns
`Inf`); pinned against `ergm.count` 4.1.3 by `test/fixtures/count_terms.toml`.

# Example
```julia
using NetworkCore, ERGMCount
using ERGM: compute, name
net = network(3; directed=true)
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 3)
add_edge!(net, 2, 3); set_edge_attribute!(net, :weight, 2, 3, 2)
compute(CMPTerm(), net) ≈ log(6) + log(2)   # true: log 3! + log 2!
name(CMPTerm())                             # "CMP"
```
"""
struct CMPTerm <: AbstractERGMTerm end

name(::CMPTerm) = "CMP"

function compute(t::CMPTerm, net)
    _refuse_negative_network(t, net, "compute")
    weights = _get_weights(net)
    total = 0.0
    for e in edges(net)
        total += _logfactorial(Int(get(weights, _wkey(net, src(e), dst(e)), 1)))
    end
    return total
end

change_stat_count(::CMPTerm, net, weights, i::Integer, j::Integer, old::Int, new::Int) =
    _logfactorial(new) - _logfactorial(old)

"""
    NonzeroTerm <: AbstractERGMTerm

Number of non-zero dyads: ∑_{i,j} I(y_{ij} ≠ 0) — R's `nonzero`. An edge
whose `:weight` is 0 is a zero dyad (R stores no such edge), and a negative
count (admissible under `DiscUnif2Reference(a < 0, b)`) is non-zero.

# Example
```julia
using NetworkCore, ERGMCount
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

change_stat_count(::NonzeroTerm, net, weights, i::Integer, j::Integer, old::Int, new::Int) =
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
using NetworkCore, ERGMCount
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
    0 > t.threshold && (total += _n_empty_dyads(net))
    return total
end

change_stat_count(t::GreaterthannTerm, net, weights, i::Integer, j::Integer, old::Int, new::Int) =
    Float64((new > t.threshold) - (old > t.threshold))

"""
    CountAtleastnTerm <: AbstractERGMTerm

Number of dyads with value >= n: ∑_{i,j} I(y_{ij} ≥ n) — R's `atleast(n)`.
The zero-valued dyads count when `n ≤ 0`, so `CountAtleastnTerm(0)` is the
number of dyads.

# Example
```julia
using NetworkCore, ERGMCount
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
    0 >= t.threshold && (total += _n_empty_dyads(net))
    return total
end

change_stat_count(t::CountAtleastnTerm, net, weights, i::Integer, j::Integer, old::Int, new::Int) =
    Float64((new >= t.threshold) - (old >= t.threshold))

"""
    AtmostTerm(threshold) <: AbstractERGMTerm

Number of dyads with value <= threshold: ∑_{i,j} I(y_{ij} ≤ threshold) — R's
`atmost(threshold)`. The zero-valued dyads count when `threshold ≥ 0`. R label
`atmost.<threshold>`; dyad-independent; pinned against `ergm` 4.12 by
`test/fixtures/count_terms.toml`.

# Example
```julia
using NetworkCore, ERGMCount
using ERGM: compute, name
net = network(3; directed=true)            # 6 dyads
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 3)
add_edge!(net, 2, 3); set_edge_attribute!(net, :weight, 2, 3, 1)
compute(AtmostTerm(1), net)   # 5.0 — the four empty dyads and (2,3)
compute(AtmostTerm(0), net)   # 4.0
name(AtmostTerm(1))           # "atmost.1"
```
"""
struct AtmostTerm <: AbstractERGMTerm
    threshold::Int
end

name(t::AtmostTerm) = "atmost.$(t.threshold)"

function compute(t::AtmostTerm, net)
    weights = _get_weights(net)
    total = 0.0
    for e in edges(net)
        w = get(weights, _wkey(net, src(e), dst(e)), 1)
        w <= t.threshold && (total += 1.0)
    end
    0 <= t.threshold && (total += _n_empty_dyads(net))
    return total
end

change_stat_count(t::AtmostTerm, net, weights, i::Integer, j::Integer, old::Int, new::Int) =
    Float64((new <= t.threshold) - (old <= t.threshold))

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
using NetworkCore, ERGMCount
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

function change_stat_count(t::CountMutualTerm, net, weights, i::Integer, j::Integer,
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
using NetworkCore, ERGMCount
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

function change_stat_count(::TransitiveTiesTerm, net, weights, i::Integer, j::Integer,
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
using NetworkCore, ERGMCount
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

function change_stat_count(::CyclicalTiesTerm, net, weights, i::Integer, j::Integer,
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
using NetworkCore, ERGMCount
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

function _out_strength(net, weights, v::Integer)
    s = 0.0
    for u in outneighbors(net, v)
        s += dyad_value(net, weights, v, u)
    end
    return s
end

function _in_strength(net, weights, v::Integer)
    s = 0.0
    for u in inneighbors(net, v)
        s += dyad_value(net, weights, u, v)
    end
    return s
end

function change_stat_count(::NodeOSumTerm, net, weights, i::Integer, j::Integer,
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
using NetworkCore, ERGMCount
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

function change_stat_count(::NodeISumTerm, net, weights, i::Integer, j::Integer,
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
using NetworkCore, ERGMCount
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

function change_stat_count(::NodeSumTerm, net, weights, i::Integer, j::Integer,
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
function _best_twopath(net, weights, i::Integer, j::Integer, skip::Integer)
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
using NetworkCore, ERGMCount
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

function change_stat_count(::TransitiveWeightsTerm, net, weights, i::Integer, j::Integer,
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
using NetworkCore, ERGMCount
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
function _best_cycle_twopath(net, weights, i::Integer, j::Integer, skip::Integer)
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

function change_stat_count(::CyclicalWeightsTerm, net, weights, i::Integer, j::Integer,
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
_admits_negative(t::SumTerm) = isinteger(t.pow)
_admits_negative(::CMPTerm) = false

function _refuse_negative(t::AbstractERGMTerm, context::AbstractString)
    why = t isa CountMutualTerm ?
          "R's `mutual(form=\"geometric\")` returns NaN there (the square root " *
          "of a negative product), which is no statistic" :
          t isa SumTerm ?
          "R's `sum(pow=)` with a non-integer power returns NaN there" :
          t isa CMPTerm ?
          "R's `CMP` returns Inf there (log y! of a negative count)" :
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
_rnum(x::Real) = !isfinite(x) ? string(x) :
                 isinteger(x) ? string(Int(x)) : @sprintf("%.15g", x)

"""
    SmallerthanTerm(threshold) <: AbstractERGMTerm

`ergm`'s valued `smallerthan(threshold)`: the number of dyads whose value is
**below** the threshold, ∑_{(i,j)} I(y_{ij} < threshold) — zero-valued dyads
included, so on a sparse network it is close to the number of dyads. R label
`smallerthan.<threshold>`; dyad-independent; pinned against `ergm` 4.12 by
`test/fixtures/count_terms.toml`.

# Example
```julia
using NetworkCore, ERGMCount
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
    n_dyads = _n_face_dyads(net)
    above = 0
    for e in edges(net)
        get(weights, _wkey(net, src(e), dst(e)), 1) >= t.threshold && (above += 1)
    end
    # Every dyad without an edge has value 0
    return Float64(n_dyads - above - (0 >= t.threshold ? n_dyads - ne(net) : 0))
end

change_stat_count(t::SmallerthanTerm, net, weights, i::Integer, j::Integer, old::Int, new::Int) =
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
using NetworkCore, ERGMCount
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
    n_dyads = _n_face_dyads(net)
    hits = 0
    for e in edges(net)
        _equalto_hit(t, Int(get(weights, _wkey(net, src(e), dst(e)), 1))) && (hits += 1)
    end
    _equalto_hit(t, 0) && (hits += n_dyads - ne(net))
    return Float64(hits)
end

change_stat_count(t::EqualToTerm, net, weights, i::Integer, j::Integer, old::Int, new::Int) =
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
using NetworkCore, ERGMCount
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
    n_dyads = _n_face_dyads(net)
    hits = 0
    for e in edges(net)
        _in_interval(t, Int(get(weights, _wkey(net, src(e), dst(e)), 1))) && (hits += 1)
    end
    _in_interval(t, 0) && (hits += n_dyads - ne(net))
    return Float64(hits)
end

change_stat_count(t::InIntervalTerm, net, weights, i::Integer, j::Integer, old::Int, new::Int) =
    Float64(_in_interval(t, new) - _in_interval(t, old))

# =============================================================================
# Valued covariate terms (ergm's `form="sum"` / `form="nonzero"` terms)
# =============================================================================
#
# ergm's valued versions of the dyadic-covariate terms: each is
# Σ_dyads x_ij · f(y_ij), with x_ij a covariate of the dyad built from vertex
# attributes (or a matrix) and f(y) = y (`form="sum"`, the default) or
# f(y) = 1[y ≠ 0] (`form="nonzero"`). They are dyad-independent: the change
# statistic of dyad (i,j) is x_ij·(f(new) − f(old)), read from nothing but the
# covariate. One statistic per term, as every count term here: R's
# multi-statistic forms (`nodefactor` over its levels, `nodematch(diff=TRUE)`)
# are written as one term per level, with R's labels.

abstract type _CountCovariateTerm <: AbstractERGMTerm end

const _COV_FORMS = (:sum, :nonzero)

function _check_form(form::Symbol, what::AbstractString)
    form in _COV_FORMS || throw(ArgumentError(
        "$what: form must be :sum or :nonzero (R's form=\"sum\"/\"nonzero\"), got :$form"))
    return form
end

@inline _cov_f(t::_CountCovariateTerm, y::Int) = t.form === :sum ? Float64(y) : Float64(y != 0)

# The attribute value of vertex v, refused with the attribute's name when
# missing (validated once per model by `_validate_term_data`)
@inline _vattr(net, attr::Symbol, v::Integer) = get_vertex_attribute(net, attr, v)
@inline _vnum(net, attr::Symbol, v::Integer) = Float64(_vattr(net, attr, v))::Float64

function compute(t::_CountCovariateTerm, net)
    weights = _get_weights(net)
    total = 0.0
    for e in edges(net)
        y = Int(get(weights, _wkey(net, src(e), dst(e)), 1))
        total += _dyad_covariate(t, net, src(e), dst(e)) * _cov_f(t, y)
    end
    return total
end

change_stat_count(t::_CountCovariateTerm, net, weights, i::Integer, j::Integer,
                  old::Int, new::Int) =
    _dyad_covariate(t, net, i, j) * (_cov_f(t, new) - _cov_f(t, old))

is_dyad_dependent(::_CountCovariateTerm) = false

# Every vertex must carry the attribute(s) a covariate term reads; refused at
# model construction, naming the attribute and the first vertex without it
_validate_term_data(::AbstractERGMTerm, net) = nothing
function _require_vertex_attr(t, net, attr::Symbol; numeric::Bool=false)
    for v in 1:nv(net)
        x = _vattr(net, attr, v)
        x === nothing && throw(ArgumentError(
            "$(name(t)): vertex $v has no `:$attr` attribute; every vertex needs " *
            "one (`set_vertex_attribute!(net, :$attr, values)`)."))
        if numeric
            (x isa Real && isfinite(x)) || throw(ArgumentError(
                "$(name(t)): `:$attr` must be a finite number on every vertex " *
                "(vertex $v has $(repr(x)))."))
        end
    end
    return nothing
end

"""
    CountNodeMatchTerm(attr; level=nothing, form=:sum) <: AbstractERGMTerm

`ergm`'s valued `nodematch(attr, form=...)`: the sum of the counts (`form=:sum`)
or the number of non-zero dyads (`form=:nonzero`) between two actors with the
same value of the vertex attribute `attr`. With `level=v`, only pairs that
both have value `v` count — one statistic of R's
`nodematch(attr, diff=TRUE)`, which has one per level; write one term per
level. R labels `nodematch.sum.<attr>` and `nodematch.sum.<attr>.<level>`
(`nodematch.nonzero.…` for `form=:nonzero`); dyad-independent; pinned
against `ergm` 4.12 by `test/fixtures/count_covariates.toml`.

# Example
```julia
using NetworkCore, ERGMCount
using ERGM: compute, name
net = network(3; directed=true)
set_vertex_attribute!(net, :g, ["a", "a", "b"])
for (i, j, w) in ((1, 2, 3), (2, 1, 1), (2, 3, 2))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
compute(CountNodeMatchTerm(:g), net)                  # 4.0 — the 1↔2 counts
compute(CountNodeMatchTerm(:g; form=:nonzero), net)   # 2.0
name(CountNodeMatchTerm(:g; level="a"))               # "nodematch.sum.g.a"
```
"""
struct CountNodeMatchTerm <: _CountCovariateTerm
    attr::Symbol
    level::Any
    form::Symbol
    CountNodeMatchTerm(attr::Symbol; level=nothing, form::Symbol=:sum) =
        new(attr, level, _check_form(form, "CountNodeMatchTerm"))
end

name(t::CountNodeMatchTerm) = "nodematch.$(t.form).$(t.attr)" *
    (t.level === nothing ? "" : ".$(t.level)")

@inline function _dyad_covariate(t::CountNodeMatchTerm, net, i::Integer, j::Integer)
    a, b = _vattr(net, t.attr, i), _vattr(net, t.attr, j)
    return Float64(isequal(a, b) && (t.level === nothing || isequal(a, t.level)))
end
_validate_term_data(t::CountNodeMatchTerm, net) = _require_vertex_attr(t, net, t.attr)

"""
    CountNodeFactorTerm(attr, level; form=:sum) <: AbstractERGMTerm

`ergm`'s valued `nodefactor(attr, form=...)`, one level: the sum over dyads of
the count (`form=:sum`) or of 1[count ≠ 0] (`form=:nonzero`), times the number
of the dyad's two actors whose `attr` equals `level`. R's `nodefactor` has one
statistic per level except the first (its default `levels=-1`); write one
term per level. R label `nodefactor.sum.<attr>.<level>`; dyad-independent;
pinned against `ergm` 4.12 by `test/fixtures/count_covariates.toml`.

# Example
```julia
using NetworkCore, ERGMCount
using ERGM: compute, name
net = network(3; directed=true)
set_vertex_attribute!(net, :g, ["a", "b", "b"])
for (i, j, w) in ((1, 2, 3), (2, 3, 2))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
compute(CountNodeFactorTerm(:g, "b"), net)   # 7.0 — 3·1 + 2·2
name(CountNodeFactorTerm(:g, "b"))           # "nodefactor.sum.g.b"
```
"""
struct CountNodeFactorTerm <: _CountCovariateTerm
    attr::Symbol
    level::Any
    form::Symbol
    CountNodeFactorTerm(attr::Symbol, level; form::Symbol=:sum) =
        new(attr, level, _check_form(form, "CountNodeFactorTerm"))
end

name(t::CountNodeFactorTerm) = "nodefactor.$(t.form).$(t.attr).$(t.level)"

@inline _dyad_covariate(t::CountNodeFactorTerm, net, i::Integer, j::Integer) =
    Float64(isequal(_vattr(net, t.attr, i), t.level) + isequal(_vattr(net, t.attr, j), t.level))
_validate_term_data(t::CountNodeFactorTerm, net) = _require_vertex_attr(t, net, t.attr)

"""
    CountAbsDiffTerm(attr; form=:sum) <: AbstractERGMTerm

`ergm`'s valued `absdiff(attr, form=...)`: the sum over dyads of the count
(`form=:sum`) or of 1[count ≠ 0] (`form=:nonzero`), times `|x_i − x_j|` for
the numeric vertex attribute `attr` (R's default `pow=1`). R label
`absdiff.sum.<attr>`; dyad-independent; pinned against `ergm` 4.12 by
`test/fixtures/count_covariates.toml`.

# Example
```julia
using NetworkCore, ERGMCount
using ERGM: compute, name
net = network(3; directed=false)
set_vertex_attribute!(net, :x, [1.0, 3.0, 2.5])
for (i, j, w) in ((1, 2, 2), (2, 3, 4))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
compute(CountAbsDiffTerm(:x), net)   # 6.0 — 2·2 + 4·0.5
name(CountAbsDiffTerm(:x))           # "absdiff.sum.x"
```
"""
struct CountAbsDiffTerm <: _CountCovariateTerm
    attr::Symbol
    form::Symbol
    CountAbsDiffTerm(attr::Symbol; form::Symbol=:sum) =
        new(attr, _check_form(form, "CountAbsDiffTerm"))
end

name(t::CountAbsDiffTerm) = "absdiff.$(t.form).$(t.attr)"

@inline _dyad_covariate(t::CountAbsDiffTerm, net, i::Integer, j::Integer) =
    abs(_vnum(net, t.attr, i) - _vnum(net, t.attr, j))
_validate_term_data(t::CountAbsDiffTerm, net) =
    _require_vertex_attr(t, net, t.attr; numeric=true)

"""
    CountNodeCovTerm(attr; form=:sum) <: AbstractERGMTerm

`ergm`'s valued `nodecov(attr, form=...)`: the sum over dyads of the count
(`form=:sum`) or of 1[count ≠ 0] (`form=:nonzero`), times `x_i + x_j` for the
numeric vertex attribute `attr`. R label `nodecov.sum.<attr>`;
dyad-independent; pinned against `ergm` 4.12 by
`test/fixtures/count_covariates.toml`. [`CountNodeOCovTerm`](@ref) and
[`CountNodeICovTerm`](@ref) are its sender and receiver halves.

# Example
```julia
using NetworkCore, ERGMCount
using ERGM: compute, name
net = network(3; directed=false)
set_vertex_attribute!(net, :x, [1.0, 3.0, 2.5])
for (i, j, w) in ((1, 2, 2), (2, 3, 4))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
compute(CountNodeCovTerm(:x), net)   # 30.0 — 2·4 + 4·5.5
name(CountNodeCovTerm(:x))           # "nodecov.sum.x"
```
"""
struct CountNodeCovTerm <: _CountCovariateTerm
    attr::Symbol
    form::Symbol
    CountNodeCovTerm(attr::Symbol; form::Symbol=:sum) =
        new(attr, _check_form(form, "CountNodeCovTerm"))
end

"""
    CountNodeOCovTerm(attr; form=:sum) <: AbstractERGMTerm

`ergm`'s valued `nodeocov(attr, form=...)`: the sum over dyads i→j of the
count (or 1[count ≠ 0]) times the SENDER's `x_i`. Directed networks only
(refused on an undirected one). R label `nodeocov.sum.<attr>`;
dyad-independent; pinned against `ergm` 4.12 by
`test/fixtures/count_covariates.toml`.

# Example
```julia
using NetworkCore, ERGMCount
using ERGM: compute, name
net = network(3; directed=true)
set_vertex_attribute!(net, :x, [1.0, 3.0, 2.5])
for (i, j, w) in ((1, 2, 2), (2, 3, 4))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
compute(CountNodeOCovTerm(:x), net)   # 14.0 — 2·1 + 4·3
name(CountNodeOCovTerm(:x))           # "nodeocov.sum.x"
```
"""
struct CountNodeOCovTerm <: _CountCovariateTerm
    attr::Symbol
    form::Symbol
    CountNodeOCovTerm(attr::Symbol; form::Symbol=:sum) =
        new(attr, _check_form(form, "CountNodeOCovTerm"))
end

"""
    CountNodeICovTerm(attr; form=:sum) <: AbstractERGMTerm

`ergm`'s valued `nodeicov(attr, form=...)`: the sum over dyads i→j of the
count (or 1[count ≠ 0]) times the RECEIVER's `x_j`. Directed networks only
(refused on an undirected one). R label `nodeicov.sum.<attr>`;
dyad-independent; pinned against `ergm` 4.12 by
`test/fixtures/count_covariates.toml`.

# Example
```julia
using NetworkCore, ERGMCount
using ERGM: compute, name
net = network(3; directed=true)
set_vertex_attribute!(net, :x, [1.0, 3.0, 2.5])
for (i, j, w) in ((1, 2, 2), (2, 3, 4))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
compute(CountNodeICovTerm(:x), net)   # 16.0 — 2·3 + 4·2.5
name(CountNodeICovTerm(:x))           # "nodeicov.sum.x"
```
"""
struct CountNodeICovTerm <: _CountCovariateTerm
    attr::Symbol
    form::Symbol
    CountNodeICovTerm(attr::Symbol; form::Symbol=:sum) =
        new(attr, _check_form(form, "CountNodeICovTerm"))
end

name(t::CountNodeCovTerm) = "nodecov.$(t.form).$(t.attr)"
name(t::CountNodeOCovTerm) = "nodeocov.$(t.form).$(t.attr)"
name(t::CountNodeICovTerm) = "nodeicov.$(t.form).$(t.attr)"
@inline _dyad_covariate(t::CountNodeCovTerm, net, i::Integer, j::Integer) =
    _vnum(net, t.attr, i) + _vnum(net, t.attr, j)
@inline _dyad_covariate(t::CountNodeOCovTerm, net, i::Integer, j::Integer) =
    _vnum(net, t.attr, i)
@inline _dyad_covariate(t::CountNodeICovTerm, net, i::Integer, j::Integer) =
    _vnum(net, t.attr, j)
requires_directed(::Union{CountNodeOCovTerm, CountNodeICovTerm}) = true
function compute(t::Union{CountNodeOCovTerm, CountNodeICovTerm}, net)
    is_directed(net) || return 0.0
    return invoke(compute, Tuple{_CountCovariateTerm, Any}, t, net)
end
_validate_term_data(t::Union{CountNodeCovTerm, CountNodeOCovTerm, CountNodeICovTerm}, net) =
    _require_vertex_attr(t, net, t.attr; numeric=true)

"""
    CountEdgeCovTerm(W; name="W", form=:sum) <: AbstractERGMTerm

`ergm`'s valued `edgecov(x, form=...)`: the sum over dyads of the count
(`form=:sum`) or of 1[count ≠ 0] (`form=:nonzero`), times the dyadic
covariate `W[i, j]` (an n×n numeric matrix; on an undirected network the
entry `W[min(i,j), max(i,j)]`, R's upper triangle). R labels the statistic by
the covariate's name, `edgecov.sum.<name>` — pass `name` as R would print it
(the network attribute's name for `edgecov("w")`). Dyad-independent; pinned
against `ergm` 4.12 by `test/fixtures/count_covariates.toml`.

# Example
```julia
using NetworkCore, ERGMCount
using ERGM: compute, name
net = network(3; directed=true)
for (i, j, w) in ((1, 2, 2), (2, 3, 4))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
W = [0.0 0.5 1.0; 2.0 0.0 0.25; 1.0 1.0 0.0]
compute(CountEdgeCovTerm(W; name="dist"), net)   # 2.0 — 2·0.5 + 4·0.25
name(CountEdgeCovTerm(W; name="dist"))           # "edgecov.sum.dist"
```
"""
struct CountEdgeCovTerm <: _CountCovariateTerm
    W::Matrix{Float64}
    label::String
    form::Symbol
    function CountEdgeCovTerm(W::AbstractMatrix{<:Real}; name::AbstractString="W",
                              form::Symbol=:sum)
        size(W, 1) == size(W, 2) || throw(ArgumentError(
            "CountEdgeCovTerm: the covariate must be a square n×n matrix (got $(size(W)))"))
        all(isfinite, W) || throw(ArgumentError(
            "CountEdgeCovTerm: the covariate must be finite (R's NA is not supported)"))
        new(Matrix{Float64}(W), String(name), _check_form(form, "CountEdgeCovTerm"))
    end
end

name(t::CountEdgeCovTerm) = "edgecov.$(t.form).$(t.label)"
@inline function _dyad_covariate(t::CountEdgeCovTerm, net, i::Integer, j::Integer)
    a, b = is_directed(net) ? (i, j) : minmax(i, j)
    return @inbounds t.W[a, b]
end
function _validate_term_data(t::CountEdgeCovTerm, net)
    size(t.W, 1) == nv(net) || throw(ArgumentError(
        "$(name(t)): the covariate is $(size(t.W, 1))×$(size(t.W, 2)) but the " *
        "network has $(nv(net)) vertices."))
    return nothing
end

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
using NetworkCore, ERGMCount
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
                               net, weights, i::Integer, j::Integer, old::Int,
                               support::UnitRange{Int})
    @inbounds for (s, y) in enumerate(support)
        dest[s] = change_stat_count(term, net, weights, i, j, old, y)
    end
    return dest
end

# log y! − log old! over the support in one pass (the per-value definition
# recomputes each factorial)
function change_stats_support!(dest::AbstractVector{Float64}, ::CMPTerm, net, weights,
                               i::Integer, j::Integer, old::Int, support::UnitRange{Int})
    base = _logfactorial(old)
    lf = _logfactorial(first(support))
    @inbounds for (s, y) in enumerate(support)
        (s > 1 && y >= 2) && (lf += log(y))
        dest[s] = lf - base
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
                               i::Integer, j::Integer, old::Int, support::UnitRange{Int})
    is_directed(net) || return fill!(dest, 0.0)
    s = _out_strength(net, weights, i) - dyad_value(net, weights, i, j)
    return _strength_profile!(dest, s, old, support)
end

function change_stats_support!(dest::AbstractVector{Float64}, ::NodeISumTerm, net, weights,
                               i::Integer, j::Integer, old::Int, support::UnitRange{Int})
    is_directed(net) || return fill!(dest, 0.0)
    s = _in_strength(net, weights, j) - dyad_value(net, weights, i, j)
    return _strength_profile!(dest, s, old, support)
end

function change_stats_support!(dest::AbstractVector{Float64}, ::NodeSumTerm, net, weights,
                               i::Integer, j::Integer, old::Int, support::UnitRange{Int})
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
                               i::Integer, j::Integer, old::Int, support::UnitRange{Int})
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
                               i::Integer, j::Integer, old::Int, support::UnitRange{Int})
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
                               weights, i::Integer, j::Integer, old::Int, support::UnitRange{Int})
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
                               weights, i::Integer, j::Integer, old::Int, support::UnitRange{Int})
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
is_dyad_dependent(::AtmostTerm) = false
is_dyad_dependent(::CMPTerm) = false
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
_undirected_hint(::Union{CountNodeOCovTerm, CountNodeICovTerm}) =
    "use `CountNodeCovTerm(attr)` (both ends) on an undirected network"

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
    n in ("NodeMatch", "MaterializedNodeMatch") &&
        return "`CountNodeMatchTerm(attr)` (R's valued `nodematch`)"
    n in ("NodeFactor", "MaterializedNodeFactor") &&
        return "`CountNodeFactorTerm(attr, level)` (R's valued `nodefactor`)"
    n == "AbsDiff" && return "`CountAbsDiffTerm(attr)` (R's valued `absdiff`)"
    n == "NodeCov" && return "`CountNodeCovTerm(attr)` (R's valued `nodecov`)"
    n == "NodeOCov" && return "`CountNodeOCovTerm(attr)` (R's valued `nodeocov`)"
    n == "NodeICov" && return "`CountNodeICovTerm(attr)` (R's valued `nodeicov`)"
    n == "EdgeCov" && return "`CountEdgeCovTerm(W)` (R's valued `edgecov`)"
    return "one of the count terms (`SumTerm`, `NonzeroTerm`, `GreaterthannTerm`, " *
           "`CountAtleastnTerm`, `AtmostTerm`, `SmallerthanTerm`, `EqualToTerm`, " *
           "`InIntervalTerm`, `CMPTerm`, " *
           "`CountMutualTerm`, `TransitiveWeightsTerm`, `CyclicalWeightsTerm`, " *
           "`NodeOSumTerm`, `NodeISumTerm`, `NodeSumTerm`, `CountNodeMatchTerm`, " *
           "`CountNodeFactorTerm`, `CountAbsDiffTerm`, `CountNodeCovTerm`, " *
           "`CountNodeOCovTerm`, `CountNodeICovTerm`, `CountEdgeCovTerm`)"
end

# A term is a count term when it has a count change statistic; a binary
# ERGM.jl term (`Edges`, `Mutual`, ...) has `compute` but no
# `change_stat_count`, and used to die deep inside the design build with a
# MethodError and a "closest candidates" list.
_is_count_term(t::AbstractERGMTerm) =
    hasmethod(change_stat_count, Tuple{typeof(t), Any, Any, Int, Int, Int, Int})

function _validate_count_terms(terms::Tuple, net::Network)
    for t in terms
        _validate_term_data(t, net)
    end
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
                                             i::Integer, j::Integer, old::Int,
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
                                i::Integer, j::Integer, support::UnitRange{Int}) where {TT<:Tuple}
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
using NetworkCore, ERGMCount
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
using NetworkCore, ERGMCount
net = network(3; directed=true)
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 2)
has_dyad_dependent(CountERGMModel([SumTerm(), NonzeroTerm()], net))   # false
has_dyad_dependent(CountERGMModel([SumTerm(), NodeOSumTerm()], net))  # true
```
"""
has_dyad_dependent(model::CountERGMModel) =
    any(is_dyad_dependent(t) for t in model.terms)


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
  than `BOUNDARY_MASS_TOL` of its conditional mass on the new top value),
  `:unconverged` (the doubling hit its cap while the estimates were still
  moving — the fit is reported but warned about), `:boundary_mode` (the
  doubling settled, but the joint distribution has a mode on the bound — see
  "Joint support" below) or `:improper` (the doubling settled, but the fitted
  coefficients make the unbounded model not normalisable — see `improper`
  below).
- `support_tol`, `support_delta`, `omitted_tail` — the tolerance and the two
  achieved bounds: `support_delta` is the largest |Δθ|/SE the last doubling
  produced, `omitted_tail` the expected number of dyads the previous bound
  omitted. Both are `0.0` for a bounded reference (nothing was truncated) and
  `NaN` for a caller-fixed `max_val` (no doubling was run).
- `support_stable` — `false` when the doubling hit `max_doublings` with
  the estimates still moving (`support_control == :unconverged`), the joint
  distribution has a mode on the bound (`:boundary_mode`) or the fitted model
  is not normalisable (`:improper`).
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
not explain — separation, R's "The MPLE does not exist!" (e.g. `sum +
nonzero` on a network whose every count is 0 or 1: `sum − nonzero` is at its
minimum on every dyad, so θ_sum → −∞, θ_nonzero → +∞ with a flat objective).
It is decided exactly from the data by the shared verdict
(`NetworkCore.clogit_separation`), not from where Newton stopped;
`separated_terms` names the coefficients that carry the direction. Newton
stops somewhere on that asymptote with arbitrarily large coefficients and
astronomical standard errors; the fit is returned with `converged = false`,
warned about, listed in `approximations`, never `is_exact`, and its z values,
p-values and confidence intervals are `NaN` (the ecosystem's separation
policy).

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
(the inverse negative pseudo-Hessian, anticonservative under dyadic dependence),
`:bootstrap` (the parametric bootstrap of `count_mple(model; se=:bootstrap)`;
`boot_replicates` then holds the `n_boot × p` refits, excluded ones as `NaN`
rows) or `:mcmc` (an MCMLE fit: inverse Fisher information estimated from the
final sample plus the Monte-Carlo error of the estimate). It is what
`NetworkCore.se_method(fit)` reports. `z_values`/`p_values` are the vectors
`coeftable(fit)` and `show(fit)` print. `inference_withheld` is `true` for an
MPLE fit (`method=:mple`) of a dyad-dependent model with the default `se`:
its naive pseudo-Hessian standard
errors under-cover, so `z_values`/`p_values` are `NaN` and `confint` refuses
(see [`count_mple`](@ref), "Inference under dyadic dependence").

# Estimator

`method` is `:mple` (maximum pseudo-likelihood, [`count_mple`](@ref)) or
`:mcmle` (Monte-Carlo maximum likelihood, [`count_mcmle`](@ref) — the
estimator of R's `ergm.count`). For an MCMLE fit `loglik` is the
log-likelihood itself, estimated by path sampling, and `mcmc` holds the
Monte-Carlo record: `convergence` (an `ERGM.MCMLEConvergence`),
`mc_std_errors`, the final `samples` of the statistics, `n_samples`,
`burnin`, `interval` (the adapted one), `loglik_mc_se`, the `start` (MPLE)
coefficients, the stopping rule, its p-value and its settings (`termination`,
`termination_p`, `conv_precision`, `conv_confidence`). The verdict is
`termination_p`: on a converged fit it is the test that passed, on the sample
the last step was taken from; on an unconverged one it is the `:confidence`
test re-run on a fresh sample at the returned coefficients (under
`:hotelling`, the last iteration's test). `samples`, `convergence`'s
t-ratios, Hotelling p-value and effective sizes, and the standard errors all
describe a sample drawn **at the returned coefficients**: one further draw at
θ̂ after convergence, or the driver's fresh draw at the last iterate when it
did not converge. They are diagnostics, not the stopping rule; `show`,
`approximations` and the warning quote only the rule (and the t-ratio only
under `:hotelling`). `mcmc` is `nothing` for an MPLE fit.

# Joint support

`boundary_mode` is `true` when the joint distribution on the enumerated
support has a mode on the truncation bound that carries its weight: started
with every dyad at `max_val`, iterated conditional modes leave at least one
dyad there, and the joint log-weight `log h(y) + θ'g(y)` of the configuration
they stop at exceeds the observed network's (a coordinate-wise mode on the
bound that weighs less than the data — the all-top network of a proper
geometric-reference model with strong reciprocity — does not count). The
support check above looks at the dyad conditionals **at the observed network
only**; a model can pass it while its joint distribution is not normalisable
on the unbounded support (a positive `mutual.product` or squared-strength
coefficient under a Poisson reference), or has its mass far beyond the bound.
Such a fit has `support_control = :boundary_mode` (`support_stable = false`)
on the adaptive path, is warned about and listed in `approximations`, and
`simulate_count_ergm`, `gof` and `se=:bootstrap` refuse it. The probe runs on
`0:(2^max_doublings · max_val)`, the bound at which the adaptive sampler
would give up, not on the fitted bound.

`improper` is `true` when the fitted coefficients make the model **not
normalisable whatever the data**: under a Poisson or geometric reference, a
statistic that grows faster than the reference decays — `mutual.product`,
the squared-strength terms `nodeOSum`/`nodeISum`/`nodeSum`, `sum(pow=p)` with
`p > 1`, or `CMP` beyond the reference's own `log y!` — has a positive leading
coefficient along some configuration of the counts (one dyad, a reciprocated
pair, a star, every dyad), so the weights `h(y)·exp(θ'g(y))` grow without
bound along it. This is decided analytically from the terms and `θ` (see
[`count_mple`](@ref), "Improper models"); such a fit has `support_control =
:improper` (`support_stable = false`) on the adaptive path, is warned about,
printed with a `WARNING:` line and listed in `approximations`, and
`simulate_count_ergm`, `gof`, `se=:bootstrap` and `count_mcmle` refuse it
unless a `max_val` is given (the truncated family, chosen in writing).

# Example
```julia
using NetworkCore, ERGMCount
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
    separated_terms::Vector{String}
    hessian_cond::Float64
    collinear::Vector{String}
    method::Symbol
    inference_withheld::Bool
    boundary_mode::Bool
    improper::Bool
    mcmc::Union{Nothing, NamedTuple}
end

# The pseudo-Hessian is reported as numerically singular above this 2-norm
# condition number: half the double-precision digits are gone (≈ 1/√eps), so
# the covariance along the flattest direction is noise.
const _HESSIAN_COND_TOL = 1e8

# NaN (unavailable) and Inf (singular) both count as ill-conditioned
_ill_conditioned(result::CountERGMResult) = !(result.hessian_cond <= _HESSIAN_COND_TOL)

_fixed_indices(result::CountERGMResult) =
    [k for k in eachindex(result.coefficients) if !isfinite(result.coefficients[k])]
# ... as `(column, :min | :max)` pairs, the side read off the sign
_fixed_sides(result::CountERGMResult) =
    Tuple{Int, Symbol}[(k, result.coefficients[k] < 0 ? :min : :max)
                       for k in _fixed_indices(result)]

# R's sentence for the statistics `count_mcmle` fixes at ∓Inf (R's default
# `drop=TRUE`), one warning per side
function _warn_count_drop(names::Vector{String}, fixed)
    for (side, word, at) in ((:min, "smallest", "-Inf"), (:max, "largest", "+Inf"))
        cols = [names[k] for (k, s) in fixed if s === side]
        isempty(cols) && continue
        @warn "count_mcmle: observed statistic(s) $(join(cols, ", ")) are at their " *
              "$word attainable value on every dyad's conditional support. Their " *
              "coefficients will be fixed at $at (no finite maximum-likelihood " *
              "estimate exists; R ergm does the same under its default drop=TRUE). " *
              "The remaining coefficients are estimated with these held fixed: the " *
              "sampler never moves the statistic off its observed bound. Pass " *
              "drop=false to refuse such a model instead."
    end
    return nothing
end

function _fixed_note(result::CountERGMResult)
    fixed = _fixed_indices(result)
    isempty(fixed) && return nothing
    names = _term_names(result.model)
    parts = ["$(names[k]) at $(result.coefficients[k] > 0 ? "+Inf" : "-Inf")"
             for k in fixed]
    head = "coefficient$(length(fixed) == 1 ? "" : "s") fixed by a boundary " *
           "statistic (" * join(parts, ", ") * "): the observed statistic is at " *
           "its extreme attainable value on every dyad, so no finite "
    if result.method === :mcmle && result.mcmc !== nothing
        return head * "maximum-likelihood estimate exists (R ergm fixes it the " *
               "same way under its default drop=TRUE); the other coefficients are " *
               "the MLE with it held at its bound — the sampler never moves the " *
               "statistic off its observed value — standard error 0 and p-value 0 " *
               "recorded, and the log-likelihood is not estimated (NaN: the path " *
               "sampler's dyad-independent reference cannot hold the statistic at " *
               "its bound)"
    end
    return head * "pseudo-likelihood estimate exists (R ergm reports the same); " *
           "standard error 0 and p-value 0 recorded"
end

# `drop=false` (R's `control.ergm(drop=FALSE)`): a statistic at the boundary of
# its attainable range is refused instead of fixed at ∓Inf. R's strict mode
# keeps the term and fits a model whose "MLE is poorly defined"; that is not
# implemented, so the strict mode names the statistics and the ways out.
function _refuse_count_no_drop(names::Vector{String}, fixed; context::AbstractString)
    lo = [names[k] for (k, s) in fixed if s === :min]
    hi = [names[k] for (k, s) in fixed if s === :max]
    parts = String[]
    isempty(lo) || push!(parts, "observed statistic(s) $(join(lo, ", ")) are at their " *
                                "smallest attainable value on every dyad's conditional " *
                                "support (coefficient -Inf)")
    isempty(hi) || push!(parts, "observed statistic(s) $(join(hi, ", ")) are at their " *
                                "largest attainable value on every dyad's conditional " *
                                "support (coefficient +Inf)")
    throw(ArgumentError(
        "$context: " * join(parts, "; ") * ". No finite estimate exists, and " *
        "drop=false asks to keep such a statistic in the model (R's `drop=FALSE`, " *
        "whose \"MLE is poorly defined\"), which is not implemented. Use the " *
        "default drop=true — the coefficient fixed at ±Inf and the rest " *
        "estimated, as R ergm does — or remove the term(s)."))
end

function _boot_exclusion_note(result::CountERGMResult)
    reps = result.boot_replicates
    reps === nothing && return nothing
    n_boot = size(reps, 1)
    n_ok = count(b -> all(isfinite, view(reps, b, :)), 1:n_boot)
    n_ok == n_boot && return nothing
    return "$(n_boot - n_ok) of the $n_boot bootstrap refits had no finite " *
           "converged count MPLE and were excluded; the covariance is over the " *
           "remaining $n_ok. " * _BOOT_EXCLUSION_BIAS
end

# The stopping rule's own verdict, in ERGM.jl's words (`_termination_verdict`):
# under `:confidence` only the equivalence test and the step length — the
# t-ratios and Hotelling test are not the rule and are not quoted; under
# `:hotelling` its p-value and the largest t-ratio
function _count_termination_verdict(rule::Symbol, p, precision, confidence,
                                    n::Int, t_ratios, γ)
    detail = rule === :confidence ?
        "$(round(Int, 100 * confidence))% equivalence test p $(_fmt3(p)) (needs < " *
        "$(_fmt3(1 - confidence)); tolerance precision $(_fmt3(precision)), $n draws)" :
        "Hotelling p $(_fmt3(p)) ($n draws), max t-ratio $(_fmt3(maximum(t_ratios; init=NaN)))"
    return detail * ", step length γ $(_fmt3(γ))"
end

function _nonconvergence_caveat(result::CountERGMResult)
    result.separated && return _separation_caveat(result)
    if result.method === :mcmle && result.mcmc !== nothing
        mc = result.mcmc
        c = mc.convergence
        verdict = _count_termination_verdict(mc.termination, mc.termination_p,
                                             mc.conv_precision, mc.conv_confidence,
                                             mc.n_samples, c.t_ratios, c.step_length)
        return "MCMLE did not converge in $(c.iterations) iteration" *
               "$(c.iterations == 1 ? "" : "s") ($verdict): the estimates are the " *
               "last iterate, NOT a maximum of the likelihood — raise `maxiter` or " *
               "`n_samples`"
    end
    return "Newton did not converge in $(result.iterations) iteration" *
           "$(result.iterations == 1 ? "" : "s") (pseudo-score norm " *
           "$(_fmt3(result.gradient_norm)) at the reported " *
           "estimates): they are NOT a maximum of the pseudo-likelihood — raise " *
           "`maxiter`, or check the model for a statistic with no finite " *
           "maximizer (a boundary or non-identified term)"
end

# What the estimator is, after `Method: <method>` (the wording ERGM.jl's
# `show` uses, so the family reads the same)
function _method_gloss(result::CountERGMResult)
    if result.method === :mcmle
        return result.mcmc === nothing ?
            " (exact: the model is dyad-independent, so its count MPLE is the " *
            "maximum likelihood estimate)" :
            " (Monte-Carlo maximum likelihood, Gibbs sampler)"
    end
    return has_dyad_dependent(result.model) ?
        " (maximum pseudo-likelihood: an approximation under dyadic dependence; " *
        "the default method=:auto fits the MCMLE here, as R's ergm.count does)" :
        " (maximum pseudo-likelihood, which is the likelihood: the formula is " *
        "dyad-independent)"
end

# The ONE sentence every bootstrap caller of the ERGM family uses to disclose
# what excluding failed refits does to the standard errors (in the warning,
# `show` and `approximations`); the same text as ERGM.jl's
const _BOOT_EXCLUSION_BIAS =
    "The standard errors are conditional on a finite refit: the excluded " *
    "replicates are the extreme ones, so the standard errors are biased downward."

# The shared separation caveat (NetworkCore), naming the separated terms, with
# R's sentence for the same design
function _separation_caveat(result::CountERGMResult)
    c = separation_caveat(result.separated_terms)
    c === nothing && (c = "separation: some combination of the coefficients is " *
                          "infinite at the maximum of the pseudo-likelihood")
    return "the MPLE does not exist — " * c * " (R ergm warns \"The MPLE does not " *
           "exist!\" for the same design); remove, merge or coarsen a separating " *
           "term, or collect more varied counts"
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
    elseif c === :boundary_mode
        return head * "the dyad conditionals at the observed network settled, " *
               "but the JOINT distribution has a mode on the bound; NOT " *
               "error-controlled)"
    elseif c === :improper
        return head * "the dyad conditionals at the observed network settled, " *
               "but the fitted model is NOT normalisable on the unbounded " *
               "support; NOT error-controlled)"
    elseif isnan(result.support_delta)
        return head * "adaptive doubling stopped here because this fit " *
               "$(result.converged ? "has a numerically singular pseudo-Hessian" :
                  "did not converge"); NOT error-controlled)"
    else
        return head * "adaptive doubling did NOT converge: last doubling still moved " *
               "the estimates by $(_fmt3(result.support_delta))·SE)"
    end
end

# What the zach comparison measured (test/fixtures/count_mcmle.toml): the count
# MPLE of `sum + nonzero + transitiveweights` sits 1.7–1.8 of R's standard
# errors from ergm.count's MCMLE on two of the three coefficients
const _MPLE_GAP_NOTE = "on ergm.count's `zach` (sum + nonzero + transitiveweights) " *
                       "it sits 1.7-1.8 standard errors from the MLE"

function _improper_caveat(result::CountERGMResult)
    dir = _improper_direction(result.model, result.coefficients)
    dir === nothing && return "the fitted model is NOT normalisable on the " *
                              "unbounded support (see the warning of the fit)"
    return _improper_message(dir, result.model.reference) * "; simulation, " *
           "`gof`, `se=:bootstrap` and `method=:mcmle` are refused unless " *
           "`max_val` fixes the truncated family"
end

function _boundary_mode_caveat(result::CountERGMResult)
    top = result.support_control === :boundary_mode ?
          "2^max_doublings × $(result.max_val), the bound at which the adaptive " *
          "sampler gives up" : "$(result.max_val)"
    return "the joint distribution has a mode on the truncation bound (started " *
           "with every dyad at $top, iterated conditional modes leave dyads " *
           "there, and that configuration outweighs the observed network). The support check " *
           "examines the dyad conditionals at the OBSERVED network only; the " *
           "fitted model is either not normalisable on the unbounded support " *
           "(e.g. a positive `mutual.product` or squared-strength coefficient " *
           "under a Poisson reference — every conditional is proper, the joint " *
           "is not; Krivitsky 2012, sec. 3) or has its mass far beyond the bound."
end

function Base.show(io::IO, result::CountERGMResult)
    mcmle = result.method === :mcmle
    println(io, "Count ERGM Results")
    println(io, "==================")
    println(io, "Reference: $(result.model.reference)")
    println(io, "Method: $(result.method)", _method_gloss(result))
    # The support is part of the estimand, not an implementation detail: print
    # it, and mark it as a truncation when the reference is really unbounded.
    println(io, _support_line(result))
    if result.truncated
        println(io, "Boundary mass: $(_fmt3(result.boundary_mass)) " *
                    "(max over dyads, at the fitted coefficients)")
    end
    if mcmle
        mc = result.mcmc
        println(io, "Log-likelihood: $(round(result.loglik, digits=4))" *
                    (mc === nothing || !isfinite(mc.loglik_mc_se) ? "" :
                     "  (path sampling; MC s.e. $(_fmt3(mc.loglik_mc_se)))"))
        println(io, "AIC: $(round(aic(result), digits=2)), BIC: $(round(bic(result), digits=2))")
    else
        println(io, "Pseudo-log-likelihood: $(round(result.loglik, digits=4))")
        println(io, "AIC: $(round(aic(result), digits=2)), BIC: $(round(bic(result), digits=2))" *
                    "  (pseudo-likelihood; compare only across models on the same " *
                    "network and support)")
    end
    println(io, "Converged: $(result.converged)")
    result.converged || println(io, "  ", _nonconvergence_caveat(result))
    (mcmle && result.mcmc !== nothing) &&
        println(io, "Termination: ", result.mcmc.termination === :confidence ?
                    "confidence rule (R ergm 4), convergence test p-value " :
                    "t-ratio and Hotelling rule, Hotelling p-value ",
                    "$(_fmt3(result.mcmc.termination_p)) ($(result.mcmc.n_samples) " *
                    "draws, $(result.mcmc.interval) sweep" *
                    "$(result.mcmc.interval == 1 ? "" : "s") apart)")
    # A near-singular pseudo-Hessian is printed here whether or not Newton
    # converged, except under separation, whose caveat already says the
    # standard errors are meaningless
    (_ill_conditioned(result) && !result.separated) &&
        println(io, "  ", _conditioning_caveat(result))
    println(io, "Std. errors: ",
            result.se_type === :bootstrap ? "parametric bootstrap" :
            result.se_type === :mcmc ? "inverse Fisher information + Monte-Carlo error" :
            "inverse pseudo-Hessian")
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
    if result.improper
        println(io)
        println(io, "WARNING: ", _improper_caveat(result))
    elseif result.boundary_mode
        println(io)
        println(io, "WARNING: ", _boundary_mode_caveat(result))
    end

    # The estimator caveat, and the prose twin of what `approximations(result)`
    # reports. A dyad-independent model needs none (there the pseudo-likelihood
    # is the likelihood). For a dyad-dependent one the POINT ESTIMATE is a
    # pseudo-likelihood estimate whatever the standard errors are, so the
    # caveat stays under `se=:bootstrap` too.
    if !mcmle && has_dyad_dependent(result.model)
        println(io)
        if result.se_type === :bootstrap
            println(io, "Note: this model contains dyad-dependent terms and was fit by maximum")
            println(io, "pseudo-likelihood. The standard errors are parametric-bootstrap estimates;")
            println(io, "the point estimates are still pseudo-likelihood estimates,")
            println(io, "not the MLE that R's ergm.count reports (on its `zach` example the two")
            println(io, "are 1.7-1.8 standard errors apart). Refit with method=:mcmle for the MLE.")
        elseif result.inference_withheld
            println(io, "Note: z values and p-values are not reported (NaN). This model contains")
            println(io, "dyad-dependent terms and was fit by maximum pseudo-likelihood; the")
            println(io, "standard errors shown are the naive inverse pseudo-Hessian ones, which")
            println(io, "treat dependent dyads as independent and under-cover (95% Wald intervals")
            println(io, "covered 0.70-0.91 in simulation), so no test or interval is built on")
            println(io, "them. The point estimates are not the MLE that R's ergm.count reports.")
            println(io, "For inference refit with method=:mcmle (maximum likelihood) or")
            println(io, "se=:bootstrap; se=:hessian requests the naive Wald table explicitly.")
        else
            println(io, "Warning: this model contains dyad-dependent terms and was fit by")
            println(io, "maximum pseudo-likelihood. The standard errors are the inverse")
            println(io, "pseudo-Hessian and are expected to be anticonservative, and the point")
            println(io, "estimates are not the MLE that R's ergm.count reports; refit with")
            println(io, "method=:mcmle, or with se=:bootstrap for a parametric-bootstrap covariance.")
        end
    end
end

# ============================================================================
# The shared result-metadata protocol (NetworkCore.jl `src/results.jl`)
# ============================================================================
#
# `fit_metadata(fit)` collects these accessors, so the truncation the `show`
# method prints in prose is also machine-readable — the two are derived from
# the same fields and cannot disagree.

estimand(::CountERGMResult) = :count_ergm

objective(result::CountERGMResult) =
    result.method === :mcmle ? :likelihood : :pseudolikelihood

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
using NetworkCore, ERGMCount
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
# (an MCMLE fit of a dyad-dependent model is a Monte-Carlo approximation, and
# `method=:mcmle` on a dyad-independent model returns the exact MPLE: the
# same rule covers both)

"""
    se_method(result::CountERGMResult) -> Symbol

What the reported standard errors ACTUALLY are: `:hessian` (the inverse negative
pseudo-Hessian), `:bootstrap` (the parametric bootstrap of
`count_mple(model; se=:bootstrap)`) or `:fisher` (an MCMLE fit: the inverse
Fisher information of the final sample plus the Monte-Carlo error). Read
straight off the fit, so it can never claim an estimator that was not used.

# Example
```julia
using NetworkCore, ERGMCount
net = network(4; directed=true)
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 2)
add_edge!(net, 3, 4); set_edge_attribute!(net, :weight, 3, 4, 1)
se_method(fit_ergm_count(net, [SumTerm()]))   # :hessian
```
"""
se_method(result::CountERGMResult) =
    result.se_type === :mcmc ? :fisher : result.se_type

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
        elseif result.support_control in (:boundary_mode, :improper)
            push!(out, "support NOT error-controlled: the doubling settled on the " *
                       "dyad conditionals at the observed network (last doubling " *
                       "moved the estimates by $(_fmt3(result.support_delta)) " *
                       "standard errors), but that check cannot see the joint " *
                       "distribution")
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
    result.improper && push!(out, _improper_caveat(result))
    (result.boundary_mode && !result.improper) &&
        push!(out, _boundary_mode_caveat(result))
    if result.method === :mcmle && result.mcmc !== nothing
        mc = result.mcmc
        push!(out, "MCMLE: the likelihood is approximated by a Gibbs sample of " *
                   "$(mc.n_samples) networks, so the estimates carry Monte-Carlo " *
                   "error (largest MC standard error " *
                   "$(_fmt3(maximum(mc.mc_std_errors; init=0.0))), included in the " *
                   "reported standard errors)")
        isfinite(result.loglik) && push!(out,
            "log-likelihood estimated by path sampling from the dyad-independent " *
            "part of the model (MC standard error $(_fmt3(mc.loglik_mc_se)))")
    elseif has_dyad_dependent(result.model)
        # The POINT ESTIMATE is a pseudo-likelihood estimate however the standard
        # errors were computed: the bootstrap replaces the covariance, not θ̂.
        push!(out, "maximum pseudo-likelihood of a dyad-dependent model: the dyad " *
                   "conditionals are multiplied as if independent, so the point " *
                   "estimates are not the MLE that R's ergm.count (MCMLE) reports " *
                   "($(_MPLE_GAP_NOTE)); `method=:mcmle` gives the MLE")
        if result.se_type === :hessian
            push!(out, "inverse-Hessian standard errors of the naive pseudo-likelihood: " *
                       "expected anticonservative under dyadic dependence (refit with " *
                       "`se=:bootstrap` for a parametric-bootstrap covariance)")
        end
        result.inference_withheld &&
            push!(out, "z values, p-values and confidence intervals withheld: " *
                       "the naive pseudo-likelihood standard errors are not " *
                       "calibrated under dyadic dependence (refit with " *
                       "method=:mcmle or se=:bootstrap; se=:hessian opts in)")
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
# `NetworkCore.check_statsapi(fit; strict=true)` in the testset.

StatsAPI.coef(result::CountERGMResult) = result.coefficients
StatsAPI.stderror(result::CountERGMResult) = result.std_errors
StatsAPI.vcov(result::CountERGMResult) = result.vcov
StatsAPI.loglikelihood(result::CountERGMResult) = result.loglik
StatsAPI.nobs(result::CountERGMResult) = ERGM.Extension.n_observed_dyads(result.model.network)
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
using NetworkCore, ERGMCount
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
using NetworkCore, ERGMCount
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
reports — `se_method(result)` says whether they are inverse-pseudo-Hessian,
parametric-bootstrap or (MCMLE) Fisher-information based. For an MPLE fit
(`method=:mple`) of a dyad-dependent model with the default `se`, whose
naive pseudo-likelihood standard errors under-cover, `confint` refuses with an `ArgumentError` (see
[`count_mple`](@ref), "Inference under dyadic dependence"); with an explicit
`se=:hessian` the naive intervals are returned. A coefficient fixed at ∓Inf
has both limits at that value.

# Example
```julia
using NetworkCore, ERGMCount
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
    result.inference_withheld && throw(ArgumentError(
        "confint: no interval is reported for an MPLE fit (method=:mple) of a " *
        "dyad-dependent model with the default se — its naive pseudo-likelihood standard errors " *
        "under-cover (95% Wald intervals covered 0.70-0.91 in simulation). Refit " *
        "with method=:mcmle (maximum likelihood) or se=:bootstrap (parametric " *
        "bootstrap), or pass se=:hessian explicitly to accept the naive Wald " *
        "intervals."))
    # The separation policy: no interval for a fit whose MPLE does not exist
    result.separated && return fill(NaN, length(result.coefficients), 2)
    q = quantile(Normal(), 1 - (1 - level) / 2)
    θ, se = result.coefficients, result.std_errors
    return hcat(θ .- q .* se, θ .+ q .* se)
end

"""
    coeftable(result::CountERGMResult) -> NetworkCore.CoefficientTable

The R-style coefficient table (`Estimate`, `Std.Error`, `z value`,
`Pr(>|z|)`) as an inspectable `NetworkCore.CoefficientTable` — exactly the table
`show(result)` prints, built from the same vectors (a method of
`StatsAPI.coeftable`). Rows are labelled with `name(term, net)`, the shared
two-argument statistic name, and can be read by index or by name.

# Example
```julia
using NetworkCore, ERGMCount
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

"""
    coefnames(result::CountERGMResult) -> Vector{String}

The coefficient labels, in `coef(result)` order — R's `names(coef(fit))`
for the same `ergm.count` formula (`"sum"`, `"nonzero"`,
`"transitiveweights.min.max.min"`, …), and the row labels of
`coeftable(result)` (a method of `StatsAPI.coefnames`). A fresh vector on
every call, so changing it cannot change the fit.

# Example
```julia
using NetworkCore, ERGMCount
net = network(4; directed=true)
for (i, j, w) in ((1, 2, 3), (2, 1, 1), (2, 3, 2), (3, 4, 1), (4, 1, 2))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
fit = fit_ergm_count(net, [SumTerm(), NonzeroTerm()])
coefnames(fit)                              # ["sum", "nonzero"]
coefnames(fit) == coeftable(fit).names      # true
```
"""
StatsAPI.coefnames(result::CountERGMResult) = _term_names(result.model)

# z statistics and two-sided p-values through the ONE shared helper; a
# coefficient fixed at ∓Inf by a boundary statistic gets z = ∓Inf, p = 0 (R's
# convention), never NaN from a 0/0.
function _count_zp(θ::Vector{Float64}, se::Vector{Float64},
                   verdict::SeparationVerdict=_NOT_SEPARATED)
    # A separated fit withholds every z and p (the shared policy)
    verdict.separated && return z_pvalues(θ, se, verdict)
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

Fit an ERGM for count-valued networks, choosing the estimator as R's
`ergm.count` does (`method=:auto`, the default):

- a **dyad-independent** formula — made only of `SumTerm`, `NonzeroTerm`,
  `GreaterthannTerm`, `CountAtleastnTerm`, `AtmostTerm`, `SmallerthanTerm`,
  `EqualToTerm`, `InIntervalTerm`, `CMPTerm` — is fit by the count MPLE
  ([`count_mple`](@ref)): each dyad's conditional distribution over the count
  support `P(y_ij = y | rest) ∝ h(y)·exp(θ'Δg(y))`, with the reference measure
  `h` in it, is a likelihood contribution, and here the pseudo-likelihood *is*
  the likelihood, so this is the exact MLE;
- a **dyad-dependent** formula is fit by Monte-Carlo maximum likelihood
  ([`count_mcmle`](@ref)), the estimator R reports. `method=:mple` asks for
  the count MPLE explicitly (fast, but an approximation whose naive inference
  is withheld; see `count_mple`).

R's `ergm` requires `reference=` for a valued model; here it defaults to
`PoissonReference()` (R's `reference=~Poisson`).

`terms` may be a Vector, a Tuple or a single term. [`ergm_count`](@ref) is
the R-faithful alias (matching the `ergm.count` package). Passing the
arguments the other way round (`fit_ergm_count(terms, net)`) is an
`ArgumentError` naming the right order.

# Arguments
- `net`: one-mode `Network` whose counts are the `:weight` edge attribute
  (integers, on **every** edge: a network with edges and no `:weight` at all
  is refused — see `weight=` — and so is one where only some edges carry
  one, naming the first bare edge; a bare edge is a data gap, not a count of
  1). A network
  with masked (unobserved) dyads is refused — see `NetworkCore.require_observed`;
  there is no `missing=` keyword, because the count MPLE would enumerate every
  unobserved dyad as an observed row. A two-mode (bipartite) network is
  refused too: its within-mode dyads would be counted as observed zeros.
- `terms`: count ERGM terms
- `reference`: Reference measure (default: Poisson)
- `weight::Symbol=:weight`: the edge attribute holding the counts — R's
  `response="w"` is `weight=:w`. Any name but `:weight` fits a `copy` of the
  network with that attribute copied to `:weight` (so `fit.model.network` is
  the copy; the caller's network is untouched).
- `method`: `:auto` (default: `:mple` for a dyad-independent formula,
  `:mcmle` otherwise, through `ERGM.resolve_method`), `:mple`
  ([`count_mple`](@ref)) or `:mcmle` (Monte-Carlo maximum likelihood,
  [`count_mcmle`](@ref) — ergm.count's estimator; its keywords `n_samples`,
  `burnin`, `interval`, `termination`, `effective_size`, `bridge_rungs`, ...
  are forwarded). Anything else is an `ArgumentError`, and so is a keyword the
  chosen estimator does not take (e.g. `se=:bootstrap` on a dyad-dependent
  formula, which `:auto` sends to the MCMLE: the error says to pass
  `method=:mple`).
- `max_val::Int`: fixes the truncation of an unbounded support. By default the
  support is chosen adaptively (doubling from twice the largest observed count
  until the estimates stop moving; see [`count_mple`](@ref)); ignored by the
  bounded references.
- `support_tol`, `max_doublings`, `maxiter`, `tol`, `se`, `n_boot`,
  `boot_burnin`, `boot_interval`, `rng`: forwarded to [`count_mple`](@ref). `se=:bootstrap`
  replaces the inverse-pseudo-Hessian covariance with a parametric bootstrap
  (same API as `ERGM.mple`); the point estimates are unchanged. With the
  default `se=nothing` a dyad-dependent MPLE fit (`method=:mple`) reports no
  z, p or interval (see `count_mple`, "Inference under dyadic dependence").

# Returns
- [`CountERGMResult`](@ref): fitted model, answering the full StatsAPI surface
  (`coef`, `stderror`, `vcov`, `confint`, `loglikelihood`, `nobs`, `dof`,
  `aic`, `bic`, `coeftable`)

# Example
```julia
using NetworkCore, ERGMCount
net = network(4; directed=true)
for (i, j, w) in ((1, 2, 3), (2, 1, 1), (2, 3, 2), (3, 4, 1), (4, 1, 2))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
fit = fit_ergm_count(net, [SumTerm(), NonzeroTerm()]; reference=PoissonReference())
fit.method                      # :mple — dyad-independent, so the MPLE is the MLE
fit.converged                   # true
coeftable(fit)                  # sum / nonzero estimates with z and p
fit_ergm_count(net, SumTerm())  # a single term works too
# A dyad-dependent formula: the count MPLE on request (the default would be
# the MCMLE)
dep = fit_ergm_count(net, [SumTerm(), CountMutualTerm()]; method=:mple)
dep.inference_withheld          # true
```
"""
function fit_ergm_count(net::Network, terms::Tuple;
                        reference::AbstractReferenceMeasure=PoissonReference(),
                        method::Symbol=:auto,
                        weight::Symbol=:weight,
                        kwargs...)
    # Count MPLE enumerates every dyad as observed, so a masked (unobserved)
    # dyad would enter the pseudo-likelihood at its face value. Reject it.
    require_observed(net; context="fit_ergm_count", face_ok=false)

    model = CountERGMModel(terms, _with_weight(net, weight), reference)
    # R's rule (ERGM.jl's one resolver): the exact MPLE for a dyad-independent
    # formula, the MCMLE — ergm.count's estimator — otherwise
    est = ERGM.resolve_method(method, has_dyad_dependent(model); exact=:mple,
                              mcmc=:mcmle, context="fit_ergm_count")
    _check_count_keywords(est, method, keys(kwargs))
    return est === :mcmle ? count_mcmle(model; kwargs...) :
                            count_mple(model; kwargs...)
end

# The keywords each estimator accepts, read from its own signature so the
# list cannot drift from the code. `count_mcmle` takes `se` only to refuse it.
_count_keywords(::Val{:mple}) = Base.kwarg_decl(which(count_mple, Tuple{CountERGMModel}))
_count_keywords(::Val{:mcmle}) =
    filter(!=(:se), Base.kwarg_decl(which(count_mcmle, Tuple{CountERGMModel})))

# A keyword the chosen estimator does not take is refused in words — most
# often an MPLE keyword (`se=:bootstrap`) on a dyad-dependent formula, which
# `method=:auto` sends to the MCMLE
function _check_count_keywords(est::Symbol, method::Symbol, keys)
    accepted = _count_keywords(Val(est))
    bad = [k for k in keys if !(k in accepted)]
    isempty(bad) && return nothing
    other = est === :mple ? :mcmle : :mple
    other_accepted = _count_keywords(Val(other))
    listed = join(("`$k`" for k in bad), ", ")
    why = method === :auto ?
        " method=:auto chose $(repr(est)) because the formula is " *
        (est === :mcmle ? "dyad-dependent (R's ergm.count fits the Monte-Carlo MLE there)." :
                          "dyad-independent (the count MPLE is the exact MLE there).") : ""
    hint = all(in(other_accepted), bad) ?
        " $(length(bad) == 1 ? "It is a keyword" : "They are keywords") of " *
        "method=$(repr(other)); pass method=$(repr(other)) explicitly to use " *
        "$(length(bad) == 1 ? "it" : "them")." :
        " See `?count_mple` and `?count_mcmle` for the keywords of each estimator."
    throw(ArgumentError("fit_ergm_count: keyword $listed is not accepted by " *
                        "method=$(repr(est)).$why$hint"))
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
                           i::Integer, j::Integer)
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
# supports can push another to its boundary), as `ERGM.Extension.boundary_columns`
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

# Whether the count MPLE exists. The count pseudo-likelihood is a conditional
# logit: each compressed row r is a stratum whose alternatives are the support
# values its (restricted) support allows, and the values its dyads were
# observed at are the chosen ones (case weights C[s, r]). It has a finite
# maximizer exactly when no direction d of the free coefficients makes every
# observed value score at least as high as every allowed value, `(x_c −
# x_s)'d ≥ 0`, with a strict inequality somewhere — e.g. `sum − nonzero` when
# every count is 0 or 1 (the observed 0 and 1 tie at 0, every y ≥ 2 scores
# 1 − y < 0), or `sum + atleast(2)` along (−1, 2) when every count is 0 or 2.
# NetworkCore's `clogit_separation` decides this exactly, by a linear programme
# certified in rational arithmetic (R's `mple.existence` is the same
# programme), so the verdict depends on the data and the support, never on
# where Newton stopped. It runs on the design actually fitted: the free
# columns, after `_count_boundary!` has fixed the single boundary columns and
# restricted the supports (run on the full design it would flag a boundary
# column itself). The reference measure h(s) is a constant per alternative and
# cannot change the answer.
function _count_separation_verdict(D::_CountDesign, cols::Vector{Int},
                                   mask::Matrix{Bool})
    isempty(cols) && return _NOT_SEPARATED
    X, C = D.X, D.C
    S = length(D.support)
    n = count(mask)
    Xr = Matrix{Float64}(undef, n, length(cols))
    chosen = Vector{Bool}(undef, n)
    strata = Vector{Int}(undef, n)
    k = 0
    for r in eachindex(D.n_tot), s in 1:S
        mask[s, r] || continue
        k += 1
        for (l, c) in enumerate(cols)
            Xr[k, l] = X[c, s, r]
        end
        chosen[k] = C[s, r] > 0
        strata[k] = r
    end
    return clogit_separation(Xr, chosen, strata)
end

const _NOT_SEPARATED = SeparationVerdict(false, Int[], Float64[], Int[], true, :clogit)

# R's own sentence for the same design, appended to the shared message
const _SEPARATION_NOTE = "R ergm warns \"The MPLE does not exist!\" for the same " *
                         "design. `fit.separated == true` and `fit.separated_terms` " *
                         "record it."

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
        # The MPLE exists exactly when the conditional logit over the
        # restricted supports has no direction of recession — decided from
        # the data alone (see `_count_separation_verdict`), wherever Newton
        # stopped; a separated fit is never converged
        verdict = _count_separation_verdict(D, cols, mask)
        separated = verdict.separated
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
        verdict = _NOT_SEPARATED
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
            separated_terms=String[names[cols[t]] for t in verdict.terms],
            verdict=verdict, hessian_cond=hcond, collinear=collinear)
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
               max_doublings=8, se=nothing, n_boot=100, boot_burnin=nothing,
               boot_interval=nothing, rng=Random.default_rng(),
               maxiter=100, tol=1e-8, warn=true, drop=true) -> CountERGMResult

Maximum pseudo-likelihood estimation for count ERGMs. For each dyad the
full conditional over the count support is enumerated, so the score is
`Σ_dyads [Δg(y_obs) − E_θ(Δg)]` and the Hessian is `−Σ_dyads Var_θ(Δg)`.
Dyads with identical conditionals are compressed into one row (every dyad, for
a dyad-independent model). The pseudo-log-likelihood is maximized with the
shared `NetworkCore.newton_fit` Newton–Raphson-with-step-halving optimizer, by
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

# Improper models

The support check is on the dyad conditionals **at the observed network**. A
model can pass it while its joint distribution is not normalisable on the
unbounded support — a positive `mutual.product` or squared-strength
coefficient under a Poisson reference makes every conditional a proper Poisson
law and the joint improper (Krivitsky 2012, sec. 3) — or has its mass far
beyond the bound. Two checks look for that, after every fit under a Poisson or
geometric reference:

- **An analytic rule.** Under a Poisson reference `log h(y) ≈ −y log y` per
  dyad, under the geometric one `0`. A statistic that grows faster than
  linearly in the counts — `mutual.product`, `nodeOSum`/`nodeISum`/`nodeSum`,
  `sum(pow=p)` with `p > 1`, and `CMP` (which beats the Poisson reference
  beyond a coefficient of 1, the geometric one beyond 0) — makes the weights
  grow without bound along some configuration of the counts when its leading
  coefficient is positive. The rule evaluates the leading order of
  `log h(y) + θ'g(y)` with one dyad, a reciprocated pair, a star and every
  dyad set to a count `Y → ∞`, so a positive coefficient on one super-linear
  term can be offset by a larger negative one on another (a positive
  `nodeOSum` under a larger negative `nodeSum`). A positive leading
  coefficient is a proof that the model is improper **whatever the data**; the
  fit is reported with `fit.improper = true` (`support_control = :improper`,
  `support_stable = false` on the adaptive path), warned about, and refused by
  `se=:bootstrap`, `simulate_count_ergm`, `gof` and `count_mcmle`. A term
  whose growth the rule does not know (a user-defined count term) makes it
  silent.
- **A probe.** Started with every dyad at the bound where the adaptive
  sampler would give up, `2^max_doublings × max_val` (the fitted bound itself
  when the caller fixed `max_val`), conditional modes are iterated; if dyads
  stay at the bound the fit is reported with `fit.boundary_mode = true` (and
  `support_control = :boundary_mode`, `support_stable = false` on the adaptive
  path), warned about, and refused the same way. It covers dyad-dependent
  models whose growth is linear (a positive `transitiveweights` sum under a
  geometric reference) or that the rule cannot read.

On a caller-fixed `max_val` both are recorded but nothing is refused: that is
the truncated family, chosen in writing (`simulate_count_ergm(fit; max_val=k)`
likewise).

# This is not ergm.count's estimator

R's `ergm.count` fits by Monte-Carlo maximum likelihood; statnet has no valued
MPLE. For a dyad-independent model the count MPLE is the exact MLE. For a
dyad-dependent one it is a different, less efficient estimator: on
`ergm.count`'s `zach` example (`sum + nonzero + transitiveweights`) it is
(0.797, −4.924, 0.385) against the MLE (0.597, −4.779, 0.562) — 1.7 and 1.8
standard errors apart on two coefficients. [`count_mcmle`](@ref)
(`method=:mcmle`) is the estimator R reports.

# Standard errors

- `se=:hessian` — the inverse negative pseudo-Hessian. For a dyad-independent
  model (e.g. `SumTerm` alone) the pseudo-likelihood is the likelihood and
  they are correct. For a model with dyad-dependent terms (`CountMutualTerm`,
  `TransitiveWeightsTerm`, the node-strength terms, ...) the pseudo-likelihood
  multiplies dyad conditionals as if independent, so they are *naive*: too
  small (see below).
- `se=:bootstrap` — parametric bootstrap: Gibbs-simulate `n_boot` count networks
  from the fitted model at θ̂ (as [`simulate_count_ergm`](@ref) does), refit
  the count MPLE on each, and report the empirical covariance of the refits.
  The point estimate is unchanged; only the covariance is replaced. Same
  option, keywords and semantics as `ERGM.mple`'s, on the ONE shared
  `NetworkCore.bootstrap_cov` loop. A replicate on which the count MPLE does not
  exist (a boundary statistic) or does not converge is excluded from the
  covariance and counted in a warning; `fit.boot_replicates` keeps every
  refit (excluded ones as `NaN` rows). The simulating chain widens an
  error-controlled support on demand; with a caller-fixed `max_val` it stays
  on `0:max_val`, and the bootstrap is refused if the chain puts more than
  `BOUNDARY_MASS_TOL` of a conditional on that bound.
- `se=nothing` (default) — `:hessian`, with the inference of a dyad-dependent
  model withheld:

# Inference under dyadic dependence

The naive pseudo-Hessian standard errors of a dyad-dependent count model are
not calibrated: simulated at the MLE of the `zach` model above, 95 % Wald
intervals built from them covered the truth 0.91, 0.98 and 0.70 of the time
(300 networks; the mean standard error of `transitiveweights` is 0.058, half
the estimator's sampling standard deviation of 0.121), against 0.96, 0.97 and
0.97 for the parametric bootstrap (100 of those networks, `n_boot=60`). So, **by default, a dyad-dependent MPLE fit reports its
point estimates and naive standard errors but no inference built on them**:
`z_values` and `p_values` are `NaN` (the coefficient table shows `NaN`, with a
note saying why), `confint` refuses with an `ArgumentError`,
`fit.inference_withheld` is `true` and `approximations(fit)` records it. For
calibrated inference use `method=:mcmle` or `se=:bootstrap`. Passing
`se=:hessian` **explicitly** is the written opt-in to the naive Wald table,
printed with its caveat. Dyad-independent models are unaffected.

# Convergence, separation and conditioning

The fit is warned about, and `fit.converged == false` recorded, when the
Newton iteration exhausts `maxiter` or cannot move (`fit.iterations`,
`fit.gradient_norm`); `approximations(fit)` then lists it and `is_exact(fit)`
is `false`.

**Separation.** A design on which the pseudo-likelihood has no finite maximum
along a *combination* of the statistics — separation, which the one-column
boundary test cannot see (e.g. `sum + nonzero` when every count is 0 or 1,
or `sum + atleast(2)` when every count is 0 or 2) — is decided exactly from
the data: the count pseudo-likelihood is a conditional logit over each dyad
class's support, and the ecosystem's shared verdict
(`NetworkCore.clogit_separation`, the linear programme of R's
`mple.existence`) runs on the design actually fitted, after the boundary
columns are fixed. Such a fit is returned with `converged = false`,
`separated = true` and the separating terms in `fit.separated_terms`; it is
warned about (naming them, with R's "The MPLE does not exist!"), its z
values, p-values and confidence intervals are `NaN`, it is listed in
`approximations`, never `is_exact`, `se=:bootstrap` is refused, and the
bootstrap excludes such replicates.

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
- `se=nothing`: `:hessian` or `:bootstrap` (validated by the shared
  `NetworkCore.check_se`); `nothing` is `:hessian` with the inference of a
  dyad-dependent model withheld
- `n_boot::Int=100`: number of bootstrap replicates (`se=:bootstrap` only)
- `boot_burnin`, `boot_interval`: Gibbs controls (in sweeps) for the bootstrap
  simulations; `nothing` resolves through the dyad-scaled rule shared with
  ERGM.jl (`ERGM.Extension.mcmc_defaults`, converted from toggles to sweeps)
- `rng::AbstractRNG=Random.default_rng()`: source of the bootstrap randomness —
  a fixed `rng` reproduces the standard errors exactly
- `maxiter::Int=100`, `tol::Float64=1e-8`: Newton controls
- `warn::Bool=true`: `false` silences the fit diagnostics (boundary statistic,
  truncation, non-convergence, undefined standard errors) — they are all still
  recorded on the result. As in `ERGM.mple`, the parametric bootstrap refits
  run with `warn=false` and report their exclusions once, in aggregate.
- `drop::Bool=true`: R's `control.ergm(drop=)`. `false` refuses a model with a
  statistic at the boundary of its attainable range (an `ArgumentError`
  naming it) instead of fixing its coefficient at ∓Inf

# Example
```julia
using NetworkCore, ERGMCount
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
                    se::Union{Symbol, Nothing}=nothing,
                    n_boot::Int=100,
                    boot_burnin::Union{Int, Nothing}=nothing,
                    boot_interval::Union{Int, Nothing}=nothing,
                    rng::Random.AbstractRNG=Random.default_rng(),
                    warn::Bool=true, drop::Bool=true)
    # The count MPLE enumerates every dyad as observed, so a masked
    # (unobserved) dyad would enter the pseudo-likelihood at its face value.
    # `fit_ergm_count` refuses it and so does this entry point (a model can
    # be built, or its network masked, after the fact); the bootstrap refits
    # run on simulated networks and never reach here.
    require_observed(model.network; context="count_mple", face_ok=false)
    # `se=nothing` (the default) is `:hessian` whose inference is withheld
    # under dyadic dependence; an explicit `se=:hessian` opts in to the naive
    # Wald table (see "Inference under dyadic dependence")
    naive_opt_in = se === :hessian
    se = something(se, :hessian)
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
    (drop || isempty(fit.fixed)) ||
        _refuse_count_no_drop(names, fit.fixed; context="count_mple")
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

    # Joint-support checks. Everything above looks at the dyad conditionals AT
    # THE OBSERVED NETWORK; a model can pass while its joint distribution is
    # not normalisable. First the analytic rule (a super-linear statistic with
    # a positive leading coefficient: improper whatever the data), then the
    # probe (see `_boundary_mode_probe`). On the error-controlled path the
    # probe runs at the bound where the adaptive sampler would give up,
    # 2^max_doublings × max_val, not at the fitted bound: a model whose joint
    # escapes between the two (e.g. at 80 when the fit stopped at 20) is the
    # case the fitted bound cannot see. A caller-fixed bound is probed as given
    # (it is the family the caller chose); an unconverged doubling is already
    # at its cap.
    dir = (truncated && !fit.separated) ? _improper_direction(model, fit.θ) : nothing
    improper = dir !== nothing
    probe_support = control === :converged ?
                    _support(ref, 2^max_doublings * mv) : support
    boundary_mode = truncated && has_dyad_dependent(model) && fit.converged &&
                    all(isfinite, fit.θ) &&
                    _boundary_mode_probe(model, fit.θ, probe_support)
    if improper && control === :converged
        control = :improper
    elseif boundary_mode && control === :converged
        control = :boundary_mode
    end

    if !warn
        # diagnostics recorded on the result, nothing printed
    elseif improper
        @warn "count_mple: " * _improper_message(dir, ref) * ". " *
              "`fit.improper == true`" *
              (control === :improper ? " and `fit.support_control == :improper`" : "") *
              " record this; `simulate_count_ergm`, `gof`, `se=:bootstrap` and " *
              "`method=:mcmle` refuse such a fit unless `max_val` fixes the " *
              "truncated family in writing. Change the model: drop the term, or " *
              "constrain its coefficient (e.g. a bounded reference such as " *
              "`BinomialReference`)." maxlog = 1
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
    elseif boundary_mode
        @warn "count_mple: the joint distribution on " *
              "$(first(probe_support)):$(last(probe_support)) has a mode on the " *
              "truncation bound: started with every dyad at $(last(probe_support)), iterated " *
              "conditional modes leave dyads there. The support check examines " *
              "the dyad conditionals at the observed network only, and they are " *
              "proper; the fitted model is either NOT normalisable on the " *
              "unbounded support (e.g. a positive `mutual.product` or " *
              "squared-strength coefficient under a Poisson reference; Krivitsky " *
              "2012, sec. 3) or has its mass far beyond the bound. " *
              "`fit.boundary_mode == true`" *
              (control === :boundary_mode ?
               " and `fit.support_control == :boundary_mode`" : "") *
              " record this; `simulate_count_ergm`, `gof` and `se=:bootstrap` " *
              "refuse such a fit." maxlog = 1
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
        warn_separation("count_mple", fit.verdict, [names[k] for k in fit.cols];
                        estimate="MPLE", note=_SEPARATION_NOTE)
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
        fit.separated && throw(ArgumentError(
            "count_mple: se=:bootstrap is refused — the MPLE does not exist " *
            "(separation on " * join(("`" * t * "`" for t in fit.separated_terms), ", ") *
            "), so there is no fitted model to simulate replicates from. Remove, " *
            "merge or coarsen the separating term(s)."))
        (improper && control !== :fixed) && throw(ArgumentError(
            "count_mple: se=:bootstrap is refused — " * _improper_message(dir, ref) *
            ". Replicates simulated from it would describe the truncation bound. " *
            "Change the model, or fix `max_val` to bootstrap the truncated family " *
            "deliberately."))
        (boundary_mode && control !== :fixed) && throw(ArgumentError(
            "count_mple: se=:bootstrap is refused — the joint distribution of " *
            "the fitted model on 0:$mv has a mode on the truncation bound (the " *
            "model is not normalisable on the unbounded support, or its mass is " *
            "far beyond the bound), so replicates simulated from it would " *
            "describe the bound, not a model. Change the model (see the " *
            "warning of the fit), or fix `max_val` to bootstrap the truncated " *
            "family deliberately."))
        vcov, std_errors, boot_replicates =
            _count_bootstrap_cov(model, fit.θ, mv, control; n_boot=n_boot,
                                 boot_burnin=boot_burnin,
                                 boot_interval=boot_interval,
                                 max_doublings=max_doublings,
                                 maxiter=maxiter, tol=tol, rng=rng)
    end

    z, pv = _count_zp(fit.θ, std_errors, fit.verdict)
    # Under dyadic dependence the naive inverse-Hessian errors under-cover:
    # unless the caller asked for them in writing, no z, p or interval is
    # built on them (a coefficient fixed at ∓Inf keeps R's z = ∓Inf, p = 0)
    withheld = se === :hessian && !naive_opt_in && has_dyad_dependent(model)
    if withheld
        for k in eachindex(fit.θ)
            isfinite(fit.θ[k]) || continue
            z[k] = NaN
            pv[k] = NaN
        end
    end
    return CountERGMResult(model, fit.θ, std_errors, z, pv, vcov, fit.loglik,
                           fit.converged, fit.iterations, fit.grad_norm, mv,
                           truncated, boundary, se, control,
                           !(control in (:unconverged, :boundary_mode, :improper)),
                           support_tol, δ, tail,
                           boot_replicates, fit.separated, fit.separated_terms,
                           fit.hessian_cond, fit.collinear, :mple, withheld,
                           boundary_mode, improper, nothing)
end

# log h over a whole support. The Poisson log-factorial is accumulated in one
# pass (`log_reference` sums it per value, O(y), which is O(|support|²) over
# the probe's wide support); the values agree to rounding, and the probe only
# takes argmaxes.
_log_reference_vector(ref::AbstractReferenceMeasure, support::UnitRange{Int}) =
    [log_reference(ref, y) for y in support]
function _log_reference_vector(ref::PoissonReference, support::UnitRange{Int})
    out = Vector{Float64}(undef, length(support))
    lf = _logfactorial(first(support))
    ll = log(ref.lambda)
    for (k, y) in enumerate(support)
        k > 1 && (lf += log(y))
        out[k] = y * ll - lf
    end
    return out
end

# Joint-support probe: does the joint distribution on `support` have a mode on
# the truncation bound that carries its mass? Start with EVERY dyad at the top
# value and iterate conditional modes (each dyad set to the argmax of its full
# conditional, in sweep order). If a sweep leaves no dyad at the top, the
# bound is not sticky and the answer is no. If the iteration reaches a fixed
# point (or its sweep budget) with dyads still there, that configuration is a
# coordinate-wise mode on the bound, and its joint log-weight
# log h(y) + θ'g(y) is compared with the observed network's: above it, the
# truncated joint puts its weight on the bound — the unbounded family has no
# normalising constant at θ, or its mass lies far beyond the bound — and the
# answer is yes. At or below it, the mode on the bound is a local one carrying
# less weight than the data themselves, and the answer is no: a proper
# geometric-reference model with a strong reciprocity term (`sum +
# mutual(:min)` with θ_sum + θ_mutual > 0 > 2θ_sum + θ_mutual) has its
# all-top configuration as a fixed point — the conditional slope there is
# θ_sum + θ_mutual > 0 — while its log-weight falls like (2θ_sum +
# θ_mutual)·Y per pair. Deterministic (no draws), O(sweeps · dyads ·
# |support| · degree); a dyad-independent model cannot have such a mode
# without its conditional showing it, so only dyad-dependent fits are probed.
function _boundary_mode_probe(model::CountERGMModel, θ::Vector{Float64},
                              support::UnitRange{Int}; max_sweeps::Int=50)
    r = _icm_from_top(model, θ, support; max_sweeps=max_sweeps)
    return r.stuck && r.logweight > r.observed
end

# Iterated conditional modes from the all-top network on `support`: whether
# dyads are still on the bound when it stops (`stuck`), the joint log-weight
# of the configuration it stopped at, and the observed network's
function _icm_from_top(model::CountERGMModel, θ::Vector{Float64},
                       support::UnitRange{Int}; max_sweeps::Int=50)
    net = copy(model.network)
    weights = get_edge_attribute(net, :weight, Int)
    terms = model.terms
    n = Int(nv(net))
    directed = is_directed(net)
    top = last(support)
    log_h = _log_reference_vector(model.reference, support)
    # The observed network's log-weight, before the dyads are moved
    obs = _joint_logweight(terms, θ, net, weights, log_h, support, model.reference)
    for i in 1:n, j in (directed ? (1:n) : ((i + 1):n))
        i == j && continue
        _set_dyad!(net, weights, i, j, top)
    end
    η = Vector{Float64}(undef, length(support))
    buf = Vector{Float64}(undef, length(support))
    at_top = 0
    for _ in 1:max_sweeps
        changed = false
        at_top = 0
        for i in 1:n, j in (directed ? (1:n) : ((i + 1):n))
            i == j && continue
            old = dyad_value(net, weights, i, j)
            _dyad_conditional!(η, buf, terms, θ, log_h, support, net, weights, i, j, old)
            new = support[argmax(η)]
            if new != old
                _set_dyad!(net, weights, i, j, new)
                changed = true
            end
            at_top += (new == top)
        end
        at_top == 0 && return (stuck=false, logweight=NaN, observed=obs)
        changed || break
    end
    return (stuck=at_top > 0,
            logweight=_joint_logweight(terms, θ, net, weights, log_h, support, model.reference),
            observed=obs)
end

# The joint log-weight log h(y) + θ'g(y) of the network `net` (unnormalised),
# with log h read off `log_h` over `support` where the dyad's value lies in it
function _joint_logweight(terms::Tuple, θ::Vector{Float64}, net, weights,
                          log_h::Vector{Float64}, support::UnitRange{Int},
                          ref::AbstractReferenceMeasure)
    n = Int(nv(net))
    directed = is_directed(net)
    lw = 0.0
    for i in 1:n, j in (directed ? (1:n) : ((i + 1):n))
        i == j && continue
        y = dyad_value(net, weights, i, j)
        lw += y in support ? log_h[y - first(support) + 1] : log_reference(ref, y)
    end
    for (k, t) in enumerate(terms)
        θ[k] == 0 && continue
        lw += θ[k] * compute(t, net)
    end
    return lw
end

missing_policies(::typeof(count_mple)) = (:error,)

# -----------------------------------------------------------------------------
# Improper models: an analytic rule
# -----------------------------------------------------------------------------
#
# Under an unbounded reference the model is normalisable only if the weights
# h(y)·exp(θ'g(y)) decay along every way the counts can grow. log h(y) is
# −log y! ≈ −y log y per dyad for Poisson (plus y log λ, which is linear) and 0
# for the geometric counting measure, so any statistic that grows FASTER than
# linearly in the counts wins against the reference whenever its leading
# coefficient is positive: `mutual.product` (y_ij·y_ji), the squared strengths
# `nodeOSum`/`nodeISum`/`nodeSum`, `sum(pow=p)` for p > 1 (y^p), and `CMP`
# (log y! ≈ y log y, which beats the Poisson reference beyond a coefficient of
# 1 and the geometric one beyond 0). Every other built-in statistic grows at
# most linearly (the minima, the thresholds, `mutual.min`/`nabsdiff`/
# `geom.mean`, the triadic weights) and cannot outrun −y log y.
#
# A super-linear term with a NEGATIVE coefficient can compensate another's
# positive one (a positive `nodeOSum` under a larger negative `nodeSum`), so
# the rule evaluates the leading order of the whole log-weight along fixed
# configurations of the counts scaled by Y → ∞: one dyad, a reciprocated pair
# (directed), an out-star and an in-star (a star, undirected), and every dyad.
# A positive leading coefficient along any of them is a proof that the sum over
# the unbounded support diverges — whatever the observed data. (A configuration
# not on the list could diverge while these do not; then the rule is silent and
# the joint-support probe and the sampler's escape check remain.) A term the
# rule does not know — a user-defined count term — makes it silent, since its
# growth could compensate.

# Leading growth of a statistic along a configuration: coefficient `c` of
# Y^p (lg = false) or of Y·log Y (lg = true). `nothing` marks a term whose
# growth is unknown; `_LINEAR` one that grows at most linearly.
struct _Growth
    p::Float64
    lg::Bool
    c::Float64
end
const _LINEAR = _Growth(1.0, false, 0.0)

# The configurations, with the number of dyads they set to Y
function _configurations(n::Int, directed::Bool)
    n >= 2 || return Symbol[]
    return directed ? (n >= 3 ? [:dyad, :mutual, :outstar, :instar, :all] :
                                [:dyad, :mutual, :all]) :
                      (n >= 3 ? [:dyad, :star, :all] : [:dyad, :all])
end

function _n_config_dyads(c::Symbol, n::Int, directed::Bool)
    c === :dyad && return 1
    c === :mutual && return 2
    c in (:outstar, :instar, :star) && return n - 1
    return directed ? n * (n - 1) : n * (n - 1) ÷ 2
end

const _CONFIG_WORDS = Dict(
    :dyad => "a single dyad", :mutual => "both dyads of one pair (i→j and j→i)",
    :outstar => "every dyad out of one actor", :instar => "every dyad into one actor",
    :star => "every dyad of one actor", :all => "every dyad")

_growth(::AbstractERGMTerm, ::Symbol, ::Int, ::Bool) = nothing
_growth(::Union{NonzeroTerm, GreaterthannTerm, CountAtleastnTerm, AtmostTerm,
                SmallerthanTerm, EqualToTerm, InIntervalTerm, TransitiveTiesTerm,
                CyclicalTiesTerm, TransitiveWeightsTerm, CyclicalWeightsTerm},
        ::Symbol, ::Int, ::Bool) = _LINEAR
_growth(::_CountCovariateTerm, ::Symbol, ::Int, ::Bool) = _LINEAR
_growth(t::SumTerm, c::Symbol, n::Int, d::Bool) =
    t.pow > 1 ? _Growth(t.pow, false, _n_config_dyads(c, n, d)) : _LINEAR
_growth(::CMPTerm, c::Symbol, n::Int, d::Bool) = _Growth(1.0, true, _n_config_dyads(c, n, d))
function _growth(t::CountMutualTerm, c::Symbol, n::Int, d::Bool)
    t.form === :product || return _LINEAR
    d || return _Growth(2.0, false, 0.0)          # identically 0 undirected
    c === :mutual && return _Growth(2.0, false, 1.0)
    c === :all && return _Growth(2.0, false, n * (n - 1) / 2)
    return _Growth(2.0, false, 0.0)
end
# Σ_i s_i² over the strengths the configuration gives (out, in, total)
function _strength_squares(which::Symbol, c::Symbol, n::Int, d::Bool)
    m = n - 1
    if !d                                   # undirected: total strength only
        c === :dyad && return 2.0
        c === :star && return m^2 + m
        return n * m^2
    end
    out, inn = if c === :dyad
        (1.0, 1.0)
    elseif c === :mutual
        (2.0, 2.0)
    elseif c === :outstar
        (m^2, m)
    elseif c === :instar
        (m, m^2)
    else
        (n * m^2, n * m^2)
    end
    which === :out && return out
    which === :in && return inn
    # total strength: (out_i + in_i)²
    c === :dyad && return 2.0
    c === :mutual && return 8.0
    c in (:outstar, :instar) && return m^2 + m
    return n * (2m)^2
end
_growth(::NodeOSumTerm, c::Symbol, n::Int, d::Bool) =
    _Growth(2.0, false, d ? _strength_squares(:out, c, n, d) : 0.0)
_growth(::NodeISumTerm, c::Symbol, n::Int, d::Bool) =
    _Growth(2.0, false, d ? _strength_squares(:in, c, n, d) : 0.0)
_growth(::NodeSumTerm, c::Symbol, n::Int, d::Bool) =
    _Growth(2.0, false, _strength_squares(:total, c, n, d))

# The exact coefficient `a` of the LINEAR growth a·Y of a statistic along a
# configuration, for the statistics that are positively homogeneous of degree
# 1 in the counts (sums and minima: g(Y·x) = Y·g(x)) or bounded (thresholds:
# a = 0); `nothing` for a super-linear statistic (its leading order is read by
# `_growth`) or one whose linear coefficient is unknown. It decides the
# geometric reference, whose counting measure does not decay at all: there a
# positive total linear coefficient along a configuration proves the model
# improper — `sum` ≥ 0 along one dyad, or 2θ_sum + θ_mutual > 0 along a
# reciprocated pair under `mutual(:min)` — whatever the data.
_linear_growth(::AbstractERGMTerm, ::Symbol, ::Int, ::Bool) = nothing
# The covariate terms' linear coefficient depends on WHICH dyads the
# configuration sets (their covariates), not only on how many: unknown here
# (the rule stays silent along a configuration that has one with a non-zero
# coefficient), except for `form=:nonzero`, which is bounded
_linear_growth(t::_CountCovariateTerm, ::Symbol, ::Int, ::Bool) =
    t.form === :nonzero ? 0.0 : nothing
_linear_growth(::Union{NonzeroTerm, GreaterthannTerm, CountAtleastnTerm, AtmostTerm,
                       SmallerthanTerm, EqualToTerm, InIntervalTerm},
               ::Symbol, ::Int, ::Bool) = 0.0
_linear_growth(t::SumTerm, c::Symbol, n::Int, d::Bool) =
    t.pow == 1 ? Float64(_n_config_dyads(c, n, d)) : t.pow < 1 ? 0.0 : nothing
function _linear_growth(t::CountMutualTerm, c::Symbol, n::Int, d::Bool)
    d || return 0.0                                  # identically 0 undirected
    f = t.form
    # y_ij·y_ji: super-linear on a reciprocated pair, identically 0 where the
    # configuration has none
    f === :product && return c in (:mutual, :all) ? nothing : 0.0
    if f === :min || f === :geometric                 # min(Y, Y) = √(Y·Y) = Y
        c === :mutual && return 1.0
        c === :all && return n * (n - 1) / 2
        return 0.0
    elseif f === :nabsdiff                            # −|y_ij − y_ji|
        c === :dyad && return -1.0
        c in (:outstar, :instar) && return -Float64(n - 1)
        return 0.0
    end
    return 0.0                                        # a threshold: bounded
end
# Two-path terms: only the configuration with every dyad at Y has a two-path
_linear_growth(::Union{TransitiveWeightsTerm, CyclicalWeightsTerm}, c::Symbol, n::Int,
               d::Bool) =
    (c === :all && n >= 3) ? (d ? Float64(n * (n - 1)) : n * (n - 1) / 2) : 0.0
_linear_growth(::TransitiveTiesTerm, c::Symbol, n::Int, d::Bool) =
    c === :all ? Float64(n * (n - 1) * (n - 2)) : 0.0
_linear_growth(::CyclicalTiesTerm, c::Symbol, n::Int, d::Bool) =
    (d && c === :all) ? n * (n - 1) * (n - 2) / 3 : 0.0

# The reference's own growth per dyad at Y: −Y log Y (Poisson), none
# (geometric); `nothing` for a truncating reference the rule does not know
_reference_growth(::PoissonReference) = _Growth(1.0, true, -1.0)
_reference_growth(::GeometricReference) = _LINEAR
_reference_growth(::AbstractReferenceMeasure) = nothing

_growth_order(g::_Growth) = (g.p, g.lg)

"""
    _improper_direction(terms, θ, reference, n, directed) -> Union{Nothing, NamedTuple}

The analytic impropriety rule (see the comment above): `nothing` when no
listed configuration makes the log-weight grow without bound (or the rule
cannot decide: a bounded reference, a non-finite coefficient, a term of
unknown growth), otherwise `(config, words, order, coef, terms)` — the
configuration, its description, the leading order `(p, log)`, the positive
leading coefficient, and the names of the terms with a positive share in it.
"""
function _improper_direction(terms, θ::AbstractVector{<:Real},
                             ref::AbstractReferenceMeasure, n::Int, directed::Bool)
    is_truncating(ref) || return nothing
    all(isfinite, θ) || return nothing
    href = _reference_growth(ref)
    href === nothing && return nothing
    for c in _configurations(n, directed)
        m = _n_config_dyads(c, n, directed)
        parts = Tuple{Tuple{Float64, Bool}, Float64, String}[]
        href.c != 0 && push!(parts, (_growth_order(href), href.c * m, ""))
        known = true
        for (k, t) in enumerate(terms)
            g = _growth(t, c, n, directed)
            g === nothing && (known = false; break)
            (g === _LINEAR || g.c == 0 || θ[k] == 0) && continue
            push!(parts, (_growth_order(g), θ[k] * g.c, name(t)))
        end
        known || return nothing
        if isempty(parts)
            # At most linear growth, and the reference does not decay faster
            # (the geometric counting measure; a Poisson reference's −Y log Y
            # is a part, so it never reaches here): the exact linear
            # coefficients decide, when every term's is known
            lin = _linear_total(terms, θ, c, n, directed)
            (lin === nothing || lin.total <= 0) && continue
            return (config=c, words=_CONFIG_WORDS[c], order=(1.0, false),
                    coef=lin.total, terms=lin.terms)
        end
        top = maximum(first, parts)
        top > (1.0, false) || continue
        total = sum(x[2] for x in parts if x[1] == top)
        total > 0 || continue
        return (config=c, words=_CONFIG_WORDS[c], order=top, coef=total,
                terms=unique([x[3] for x in parts if x[1] == top && x[2] > 0 &&
                                                    !isempty(x[3])]))
    end
    return nothing
end

# Σ_k θ_k · a_k along configuration `c` (`_linear_growth`), with the names of
# the terms contributing positively; `nothing` when a coefficient is unknown
function _linear_total(terms, θ, c::Symbol, n::Int, directed::Bool)
    total = 0.0
    pos = String[]
    for (k, t) in enumerate(terms)
        θ[k] == 0 && continue
        a = _linear_growth(t, c, n, directed)
        a === nothing && return nothing
        total += θ[k] * a
        θ[k] * a > 0 && push!(pos, name(t))
    end
    return (total=total, terms=unique(pos))
end

_improper_direction(model::CountERGMModel, θ::AbstractVector{<:Real}) =
    _improper_direction(model.terms, θ, model.reference, Int(nv(model.network)),
                        is_directed(model))

function _improper_message(d, ref)
    order = d.order[2] ? "Y·log Y" : d.order[1] == 1 ? "Y" :
            (d.order[1] == round(d.order[1]) ? "Y^$(Int(d.order[1]))" : "Y^$(d.order[1])")
    return "the fitted model is NOT normalisable on the unbounded support of " *
           "$(nameof(typeof(ref))), whatever the data: with $(d.words) at a count " *
           "Y, the log-weight log h(y) + θ'g(y) grows like $(_fmt3(d.coef))·$order " *
           "as Y → ∞ (from " * join(("`" * t * "`" for t in d.terms), ", ") *
           " with a positive coefficient), faster than the reference decays " *
           "(Krivitsky 2012, sec. 3). Every dyad conditional can still be " *
           "proper; the joint distribution is not. The estimates describe the " *
           "family truncated at the enumerated bound, not a model"
end

# Parametric-bootstrap covariance of the count MPLE: Gibbs-simulate `n_boot`
# count networks at θ̂, refit the count MPLE on each, take the empirical
# covariance. The loop is the shared `NetworkCore.bootstrap_cov`; this supplies only
# the two callbacks that are ERGMCount's. Every replicate is refit on the
# support the simulating chain ended on, so no replicate can fall outside it.
# A replicate without a finite, converged MPLE (a boundary statistic in the
# simulated network, or Newton failing) is returned as NaN and excluded, with
# a warning: a NaN row must never enter the covariance silently.
function _count_bootstrap_cov(model::CountERGMModel, θ̂::Vector{Float64},
                              mv::Int, control::Symbol; n_boot::Int, boot_burnin,
                              boot_interval, max_doublings::Int,
                              maxiter::Int, tol::Float64,
                              rng::Random.AbstractRNG)
    # A caller-fixed `max_val` is the truncated family the caller asked for:
    # simulate it on that support, and refuse if the bound shapes the draws.
    # An error-controlled fit stands for the unbounded family: the chain
    # widens its support on demand (and throws if it never settles), and the
    # refits enumerate whatever support the chain ended on.
    fixed = control === :fixed
    top = Ref(mv)
    function simulate(rng, B)
        chain = _simulate_chain(model.network, model.terms, model.reference, θ̂;
                                n_sim=B, burnin=boot_burnin, interval=boot_interval,
                                max_val=fixed ? mv : nothing, start=mv,
                                max_doublings=max_doublings, rng=rng,
                                context="count_mple(se=:bootstrap)",
                                on_boundary=:error)
        top[] = max(mv, chain.max_val)
        return chain.networks
    end

    p = length(θ̂)
    function refit(sim::Network)
        boot_model = CountERGMModel(model.terms, sim, model.reference)
        f = _count_mple_fit(boot_model, _support(model.reference, top[]);
                            maxiter=maxiter, tol=tol, θ0=θ̂)
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
          "excluded; the covariance is over the remaining $n_ok refits. " *
          _BOOT_EXCLUSION_BIAS * " This is " *
          "about the simulated replicates, not about the observed network. " *
          "`fit.boot_replicates` holds every refit (NaN rows excluded)." maxlog = 1
    V = Matrix{Float64}(cov(replicates[ok, :]))
    return V, sqrt.(max.(diag(V), 0.0)), replicates
end

# =============================================================================
# Monte-Carlo maximum likelihood (ergm.count's estimator)
# =============================================================================

"""
    count_mcmle(model::CountERGMModel; n_samples=1024, burnin=nothing,
                interval=nothing, effective_size=:auto, maxiter=60,
                termination=:confidence, conv_precision=0.1, conv_confidence=0.99,
                conv_threshold=0.1, hotelling_alpha=0.05, gamma0=0.1,
                max_step_norm=5.0, max_n_samples=nothing,
                bridge_rungs=16, bridge_samples=nothing,
                max_val=nothing, max_doublings=8, support_tol=1e-3, init=nothing,
                rng=Random.default_rng(), warn=true, drop=true) -> CountERGMResult

Monte-Carlo maximum likelihood for a count ERGM — the estimator R's
`ergm.count` uses (statnet has no valued MPLE). `fit_ergm_count(net, terms;
method=:mcmle)` calls it.

**The iteration is ERGM.jl's** (`ERGM.Extension.mcmle_solve`, the one `ERGM.mcmle`
runs); this function supplies the sampler. Starting from the count MPLE
([`count_mple`](@ref)), each iteration draws `n_samples` networks from the
model at the current coefficients with the exact Gibbs sampler of
[`simulate_count_ergm`](@ref) (after `burnin` sweeps; the chain continues
across iterations) and takes the Hummel-stepped Monte-Carlo Newton step
toward the observed statistics, `θ ← θ + γ·Σ̂⁻¹(g(y_obs) − ḡ)`. The stopping
rule, tested at full step length, is `termination`:

- `:confidence` (default) — R `ergm` 4's equivalence test: with confidence
  `conv_confidence` the estimating equations at the *updated* coefficients
  lie inside the tolerance region `x'(conv_precision·Σ̂)⁻¹x ≤ 1`
  (`ERGM.Extension.confidence_test`); when it fails near the solution the next sample
  is larger, up to `max_n_samples` (default `16·n_samples`), as R boosts its
  sample;
- `:hotelling` — every statistic's t-ratio `|g_obs − ḡ|/sd` below
  `conv_threshold` and a Hotelling T² test at `hotelling_alpha`
  (`ERGM.mcmc_convergence`).

**Sampling is ESS-adaptive** (`ERGM.Extension.ess_sample`, R's `MCMC.effectiveSize`):
the interval between retained draws starts at `interval` sweeps (default: the
dyad-scaled rule) and doubles until the effective sample size of the
`n_samples` draws reaches `effective_size` (default `:auto` = `n_samples ÷ 2`;
a number sets it, `nothing` keeps the interval fixed). On `ergm.count`'s
`zach` example a sweep-per-draw chain has an effective size of about 180 in
1024; the adapted interval is 4–8 sweeps.

The standard errors are the inverse Fisher information plus the Monte-Carlo
error of the estimate (`ERGM.Extension.mcmle_covariance`, as `ERGM.mcmle` reports
them), computed — with `fit.mcmc.convergence`'s t-ratios and Hotelling
p-value — from a sample drawn **at the returned coefficients**. On a
converged fit that is one further sample drawn after the rule passed: the
sample the last step was taken from can sit a standard error away when the
rule passes after one or two steps from the MPLE, and the Fisher information
there is not the one at the estimate. On a fit that did not converge it is
the fresh sample the iteration draws at its last iterate. The stopping
verdict, `fit.mcmc.termination_p`, is the rule's own test (on a converged fit,
the one that passed, on the pre-step sample); the warning, `show` and
`approximations` quote only that rule and the step length γ — plus the
largest t-ratio under `termination=:hotelling`, whose rule it is. At the defaults the seed-to-seed spread
of the estimates and their Monte-Carlo error (about 4 % of a standard error
on `zach`) are those of `ergm.count`.

`loglikelihood(fit)` is the log-likelihood itself, by path sampling: the
dyad-independent part of the fitted model (the dyad-dependent coefficients set
to zero) has an exact likelihood, and `log Z` is integrated from there to `θ̂`
over `bridge_rungs` Simpson segments with `bridge_samples` draws per grid point
(default `max(64, n_samples ÷ 4)`; `fit.mcmc.loglik_mc_se` is its Monte-Carlo
standard error). It includes the reference measure, so it is comparable
with the (pseudo-)log-likelihood of an MPLE fit of a dyad-independent model on
the same network; R's `logLik` is the same number minus the log-likelihood at
`θ = 0`. `bridge_rungs=0` skips it (`NaN`).

A **dyad-independent** model needs no Monte Carlo: its MPLE is the exact MLE
and is returned as is (`fit.mcmc === nothing`).

# Support

An unbounded reference is sampled on an adaptive support: the chain widens
`0:max_val` whenever a dyad's conditional puts more than 1e-10 of its mass on
the top value, starting from the MPLE's error-controlled bound, and an
`ArgumentError` is thrown if it reaches `2^max_doublings` times that bound —
the signature of a model that is not normalisable. `max_val=k` fixes the
truncated family on `0:k` instead; a chain that puts more than
`BOUNDARY_MASS_TOL` of a conditional on `k` is refused.

# Refusals

`ArgumentError` when the MPLE start is separated or did not converge, or when
its joint-support probe found a mode on the truncation bound
(`fit.boundary_mode`); and when `se=` is passed (an MPLE option — MCMLE
standard errors are always Fisher information plus Monte-Carlo error).

# Boundary statistics

A statistic at the boundary of its attainable range on every dyad's
conditional support (the MPLE fixes its coefficient at ∓Inf) is dropped as R
does under its default `control.ergm(drop=TRUE)`: the coefficient is fixed at
∓Inf (standard error 0, p-value 0) with R's warning, and the others are the
MLE with the statistic held at its observed bound — the sampler gives every
value that would move it zero weight. The log-likelihood is then not
estimated (`NaN`). `drop=false` refuses such a model with an `ArgumentError`
instead (R's `drop=FALSE` keeps the term; that is not implemented); a model
with every statistic at its bound is refused (nothing to estimate).

# Example
```julia
using NetworkCore, ERGMCount, Random
net = network(6; directed=true)
for (i, j, w) in ((1, 2, 3), (2, 1, 2), (2, 3, 2), (3, 2, 1), (3, 4, 1), (4, 5, 2),
                  (5, 4, 1), (1, 5, 2), (5, 6, 1), (6, 2, 1), (6, 5, 2))
    add_edge!(net, i, j); set_edge_attribute!(net, :weight, i, j, w)
end
fit = fit_ergm_count(net, [SumTerm(), CountMutualTerm()]; method=:mcmle,
                     n_samples=400, bridge_samples=100, rng=Xoshiro(1))
fit.method                       # :mcmle
fit.mcmc.convergence.iterations  # MCMLE iterations run
coeftable(fit)                   # z and p from Fisher + Monte-Carlo standard errors
```
"""
function count_mcmle(model::CountERGMModel;
                     n_samples::Int=1024,
                     burnin::Union{Int, Nothing}=nothing,
                     interval::Union{Int, Nothing}=nothing,
                     maxiter::Int=60,
                     termination::Symbol=:confidence,
                     conv_precision::Float64=0.1,
                     conv_confidence::Float64=0.99,
                     conv_threshold::Float64=0.1,
                     hotelling_alpha::Float64=0.05,
                     gamma0::Float64=0.1,
                     max_step_norm::Float64=5.0,
                     max_n_samples::Union{Int, Nothing}=nothing,
                     effective_size::Union{Real, Nothing, Symbol}=:auto,
                     bridge_rungs::Int=16,
                     bridge_samples::Union{Int, Nothing}=nothing,
                     max_val::Union{Int, Nothing}=nothing,
                     max_doublings::Int=_MAX_SUPPORT_DOUBLINGS,
                     support_tol::Float64=1e-3,
                     init::Union{Nothing, AbstractVector{<:Real}}=nothing,
                     rng::Random.AbstractRNG=Random.default_rng(),
                     warn::Bool=true, drop::Bool=true,
                     se::Union{Symbol, Nothing}=nothing)
    require_observed(model.network; context="count_mcmle", face_ok=false)
    se === nothing || throw(ArgumentError(
        "count_mcmle: `se=:$se` is an option of the MPLE (`method=:mple`). An " *
        "MCMLE fit reports the inverse Fisher information of its final sample " *
        "plus the Monte-Carlo error; there is nothing to choose."))
    n_samples >= 16 || throw(ArgumentError(
        "count_mcmle: n_samples must be at least 16 (got $n_samples)"))
    maxiter >= 1 || throw(ArgumentError(
        "count_mcmle: maxiter must be at least 1 (got $maxiter)"))
    bridge_rungs >= 0 || throw(ArgumentError(
        "count_mcmle: bridge_rungs must be non-negative (got $bridge_rungs)"))
    termination in (:confidence, :hotelling) || throw(ArgumentError(
        "count_mcmle: termination must be :confidence or :hotelling (got :$termination)"))
    n_max = something(max_n_samples, 16 * n_samples)
    n_max >= n_samples || throw(ArgumentError(
        "count_mcmle: max_n_samples ($n_max) is below n_samples ($n_samples)"))
    ess_target = effective_size === :auto ? n_samples / 2 :
                 effective_size === nothing ? nothing :
                 effective_size isa Real ? Float64(effective_size) :
                 throw(ArgumentError("count_mcmle: effective_size must be a number, " *
                                     "`nothing` or `:auto` (got :$effective_size)"))
    (ess_target === nothing || ess_target >= 8) || throw(ArgumentError(
        "count_mcmle: effective_size must be at least 8 (got $effective_size)"))

    if !has_dyad_dependent(model)
        # The pseudo-likelihood IS the likelihood: the MPLE is the exact MLE,
        # with exact Fisher-information standard errors and log-likelihood
        exact = count_mple(model; max_val=max_val, support_tol=support_tol,
                           max_doublings=max_doublings, se=:hessian, warn=warn,
                           drop=true)
        drop || isempty(_fixed_indices(exact)) ||
            _refuse_count_no_drop(_term_names(model), _fixed_sides(exact);
                                  context="count_mcmle")
        return _with_method(exact, :mcmle)
    end

    start = count_mple(model; max_val=max_val, support_tol=support_tol,
                       max_doublings=max_doublings, se=:hessian, warn=false)
    ref0 = model.reference
    names = _term_names(model)
    p = length(names)
    # A statistic at the boundary of its attainable range on every dyad's
    # conditional support (the MPLE fixed it at ∓Inf): R's default drop=TRUE.
    # The coefficient is fixed at ∓Inf and the others are estimated with the
    # statistic held at its observed bound — the sampler gives every value
    # that would move it off the bound zero weight (±1e300 in the Gibbs
    # conditional) — so the rest is the MLE of the model restricted to the
    # networks where the statistic keeps its observed value
    fixed = _fixed_sides(start)
    if !isempty(fixed)
        drop || _refuse_count_no_drop(names, fixed; context="count_mcmle")
        warn && _warn_count_drop(names, fixed)
    end
    free = [k for k in 1:p if isfinite(start.coefficients[k])]
    isempty(free) && throw(ArgumentError(
        "count_mcmle: every statistic is at the boundary of its attainable range " *
        "($(join(names, ", "))), so every coefficient is fixed at ∓Inf and there " *
        "is nothing to estimate; `method=:mple` reports the limits."))
    start.separated && throw(ArgumentError(
        "count_mcmle: the MPLE start does not exist (separation on " *
        join(("`" * t * "`" for t in start.separated_terms), ", ") * "), so there " *
        "is nothing to start the Monte-Carlo iteration from; the MLE does not " *
        "exist either. Remove, merge or coarsen the separating term(s)."))
    (start.converged || init !== nothing) || throw(ArgumentError(
        "count_mcmle: the MPLE start did not converge; fix the model (see " *
        "`fit_ergm_count(...; method=:mple)`), or pass `init=`."))
    fixed_support = start.support_control === :fixed
    if start.improper && !fixed_support
        throw(ArgumentError("count_mcmle: at the MPLE start, " *
            _improper_message(_improper_direction(model, start.coefficients), ref0) *
            ". There is no likelihood to maximize. Change the model, or fix " *
            "`max_val` to fit the truncated family deliberately."))
    end
    (start.boundary_mode && !fixed_support) && throw(ArgumentError(
        "count_mcmle: " * _boundary_mode_caveat(start) * " There is no " *
        "likelihood to maximize. Change the model, or fix `max_val` to fit the " *
        "truncated family deliberately."))

    net, terms, ref = model.network, model.terms, model.reference
    init === nothing || length(init) == p || throw(ArgumentError(
        "count_mcmle: init has $(length(init)) coefficients for $p terms"))
    θ0 = init === nothing ? start.coefficients[free] : Vector{Float64}(init)[free]
    all(isfinite, θ0) || throw(ArgumentError(
        "count_mcmle: init must be finite on the estimated coefficients (got $init)"))
    g_obs = Float64[compute(t, net) for t in terms]
    # The coefficients the sampler runs at: the estimated ones, and the fixed
    # ones at ±1e300 — finite, so a value that keeps the statistic at its bound
    # (change 0) contributes 0 and one that moves it gets weight exp(−1e300) = 0
    θ_sampler = [isfinite(c) ? 0.0 : (c > 0 ? 1e300 : -1e300) for c in start.coefficients]

    # The sampler: `n` retained states of the Gibbs chain at θ, as their
    # statistics. The chain (its state and its support) continues across
    # draws; nothing else of the iteration lives here — the Newton step, the
    # stopping rule, the sample boost and the covariance are ERGM.jl's
    # `ERGM.Extension.mcmle_solve`, the iteration `ERGM.mcmle` runs.
    state = copy(net)
    mv = Ref(start.max_val)
    top_mass = Ref(0.0)
    last_interval = Ref(0)
    burn = something(burnin, _gibbs_defaults(ERGM.Extension.n_observed_dyads(model.network)).burnin)
    base_interval = something(interval, _gibbs_defaults(ERGM.Extension.n_observed_dyads(model.network)).interval)
    iv = Ref(base_interval)
    function draw(θ, n)
        θv = copy(θ_sampler)
        θv[free] = θ
        chain_kw = (; max_val=fixed_support ? mv[] : nothing, start=mv[],
                    max_doublings=max_doublings, rng=rng, context="count_mcmle",
                    on_boundary=:error, state=state)
        # `m` further draws of the chain, `interval` sweeps apart, as statistics
        function extend(_, m, interval)
            G = Matrix{Float64}(undef, m, p)
            row = Ref(0)
            function record!(cur)
                r = (row[] += 1)
                for (k, t) in enumerate(terms)
                    G[r, k] = compute(t, cur)
                end
                return nothing
            end
            chain = _simulate_chain(net, terms, ref, θv; n_sim=m, burnin=0,
                                    interval=interval, chain_kw..., stat! = record!)
            mv[] = max(mv[], chain.max_val)
            top_mass[] = max(top_mass[], chain.max_top_mass)
            # the fixed statistics stay at their observed values: only the
            # estimated ones enter the iteration
            return length(free) == p ? G : G[:, free]
        end
        # burn in at the new coefficients, then sample
        burned = _simulate_chain(net, terms, ref, θv; n_sim=0, burnin=burn,
                                 interval=1, chain_kw...)
        mv[] = max(mv[], burned.max_val)
        top_mass[] = 0.0
        ess_target === nothing && return extend(1, n, base_interval)
        # ERGM.jl's ESS-adaptive chooser: the interval doubles (the stored
        # size staying n) until the effective sample size reaches the
        # target; the next draw starts at half the interval reached, as R does
        samples, _, reached = ERGM.Extension.ess_sample(extend, ess_target, n; interval=iv[])
        iv[] = max(base_interval, reached ÷ 2)
        last_interval[] = reached
        return samples
    end
    solve() = ERGM.Extension.mcmle_solve(draw, θ0; labels=names[free], n_samples=n_samples,
                                maxiter=maxiter, termination=termination,
                                conv_precision=conv_precision,
                                conv_confidence=conv_confidence,
                                conv_threshold=conv_threshold,
                                hotelling_alpha=hotelling_alpha, gamma0=gamma0,
                                max_step_norm=max_step_norm, max_n_samples=n_max,
                                target=g_obs[free], context="count_mcmle")
    sol = warn ? solve() : with_logger(solve, NullLogger())
    θf = sol.coef
    G = Matrix{Float64}(sol.final.samples)
    n_final = size(G, 1)
    V, se, mcse, tests = sol.vcov, sol.se, sol.mcmc_se, sol.tests
    if sol.converged
        # The driver's covariance is that of the sample the estimate was
        # stepped FROM. Starting at the MPLE the rule can pass after one or two
        # steps, so that sample sits up to a standard error from θ̂ and its
        # Fisher information is not the one at the estimate (measured on zach:
        # standard errors up to 9 % below R's). One more draw AT θ̂ gives the
        # covariance, the Monte-Carlo error and the diagnostics reported.
        G = draw(θf, n_final)
        _, V, se, mcse = with_logger(NullLogger()) do
            ERGM.Extension.mcmle_covariance(G, [n_final], nothing, nothing, length(free))
        end
        tests = ERGM.mcmc_convergence(G, g_obs[free]; conv_threshold=conv_threshold,
                                      hotelling_alpha=hotelling_alpha)
    end
    # Back to every coefficient: a fixed one is ∓Inf with standard error 0
    # (R's convention), its statistic constant in every draw
    θ = copy(start.coefficients)
    θ[free] = θf
    if length(free) < p
        Vf, sef, mcsef = V, se, mcse
        V = zeros(p, p); V[free, free] = Vf
        se = zeros(p); se[free] = sef
        mcse = zeros(p); mcse[free] = mcsef
        Gf = G
        G = repeat(reshape(g_obs, 1, p), size(Gf, 1), 1)
        G[:, free] = Gf
    end
    singular = !all(isfinite, se)
    convergence = ERGM.MCMLEConvergence((sol.iterations, sol.step_length,
                                         tests.t_ratios, tests.hotelling_p,
                                         tests.n_eff))
    if warn && !sol.converged && !singular
        verdict = _count_termination_verdict(termination, sol.termination_p,
                                             conv_precision, conv_confidence,
                                             n_final, convergence.t_ratios,
                                             sol.step_length)
        @warn "count_mcmle: MCMLE did not converge in $(sol.iterations) iteration" *
              "$(sol.iterations == 1 ? "" : "s") ($verdict): the estimates are " *
              "the last iterate, not a maximum of the likelihood. Raise `maxiter` " *
              "or `n_samples`; `fit.converged == false` records this." maxlog = 1
    end

    loglik, ll_se = NaN, NaN
    # With a statistic held at its bound the path sampler's dyad-independent
    # reference cannot hold it there: the log-likelihood is not estimated
    # (NaN; `show` and `approximations` say why)
    if bridge_rungs > 0 && !singular && all(isfinite, θ)
        try
            loglik, ll_se = _count_bridge_loglik(
                model, θ, g_obs, mv[]; nrungs=bridge_rungs,
                n_samples=something(bridge_samples, max(64, n_samples ÷ 4)),
                burnin=burn,
                interval=interval, fixed_support=fixed_support,
                max_doublings=max_doublings, rng=rng)
        catch e
            e isa ArgumentError || rethrow()
            warn && @warn "count_mcmle: the log-likelihood could not be estimated " *
                          "by path sampling ($(e.msg)) and is reported as NaN; the " *
                          "coefficients and standard errors are unaffected." maxlog = 1
        end
    end

    z, pv = _count_zp(θ, se)
    hcond = singular ? Inf : _hessian_cond(cov(G[:, free]))
    # The analytic rule at the MLE: the iteration can cross into the improper
    # region even from a proper start (only possible on a caller-fixed
    # support, where the chain cannot escape)
    dir = _improper_direction(model, θ)
    improper = dir !== nothing
    (warn && improper) &&
        @warn "count_mcmle: " * _improper_message(dir, ref) * ". `fit.improper " *
              "== true` records this; the fit describes the family truncated at " *
              "0:$(mv[])." maxlog = 1
    mcmc = (convergence=convergence, mc_std_errors=mcse, samples=G,
            n_samples=n_final, burnin=burn,
            interval=last_interval[] == 0 ? base_interval : last_interval[],
            loglik_mc_se=ll_se, start=copy(start.coefficients),
            termination=termination, termination_p=sol.termination_p,
            conv_precision=conv_precision, conv_confidence=conv_confidence)
    diff = g_obs .- vec(sum(G, dims=1)) ./ n_final
    return CountERGMResult(model, θ, se, z, pv, V, loglik,
                           sol.converged, sol.iterations, norm(diff), mv[],
                           start.truncated, top_mass[], :mcmc,
                           start.support_control, start.support_stable,
                           start.support_tol, start.support_delta,
                           start.omitted_tail, nothing, false, String[], hcond,
                           String[], :mcmle, false, start.boundary_mode, improper,
                           mcmc)
end

missing_policies(::typeof(count_mcmle)) = (:error,)

# The same fit under another `method` label (a dyad-independent model's exact
# MPLE returned by `method=:mcmle`)
function _with_method(r::CountERGMResult, method::Symbol)
    return CountERGMResult(r.model, r.coefficients, r.std_errors, r.z_values,
                           r.p_values, r.vcov, r.loglik, r.converged, r.iterations,
                           r.gradient_norm, r.max_val, r.truncated, r.boundary_mass,
                           r.se_type, r.support_control, r.support_stable,
                           r.support_tol, r.support_delta, r.omitted_tail,
                           r.boot_replicates, r.separated, r.separated_terms,
                           r.hessian_cond, r.collinear, method, r.inference_withheld,
                           r.boundary_mode, r.improper, r.mcmc)
end

# Variance of the mean of a (possibly autocorrelated) scalar series by batch
# means, √n batches
function _batch_var_of_mean(x::AbstractVector{Float64})
    n = length(x)
    b = max(2, floor(Int, sqrt(n)))
    len = n ÷ b
    len >= 1 || return NaN
    means = [sum(@view x[((k - 1) * len + 1):(k * len)]) / len for k in 1:b]
    m = sum(means) / b
    return sum(abs2, means .- m) / (b - 1) / b
end

# Path-sampling estimate of the log-likelihood θ'g(y_obs) + Σ log h(y_obs) −
# log Z(θ) and its Monte-Carlo standard error. θ₀ is θ with the dyad-dependent
# coordinates zeroed: a dyad-independent model, whose log-likelihood at y_obs
# is the count pseudo-log-likelihood exactly. Along θ(u) = θ₀ + u(θ − θ₀),
# d/du log Z = (θ − θ₀)'E_{θ(u)}[g], integrated by composite Simpson's rule
# over `nrungs` segments from Gibbs samples at each grid point. Each grid point
# has its own RNG seeded from `rng` up front and its own chain, so the estimate
# is reproducible and thread-count independent.
function _count_bridge_loglik(model::CountERGMModel, θ::Vector{Float64},
                              g_obs::Vector{Float64}, mv::Int; nrungs::Int,
                              n_samples::Int, burnin, interval,
                              fixed_support::Bool, max_doublings::Int,
                              rng::Random.AbstractRNG)
    isodd(nrungs) && (nrungs += 1)
    net, terms, ref = model.network, model.terms, model.reference
    p = length(θ)
    θ0 = copy(θ)
    for (k, t) in enumerate(terms)
        is_dyad_dependent(t) && (θ0[k] = 0.0)
    end
    Δ = θ .- θ0

    seeds = rand(rng, UInt64, nrungs + 1)
    vars = Vector{Float64}(undef, nrungs + 1)
    tops = fill(mv, nrungs + 1)
    # E_{θ(u)}[g] at the k-th grid point, from that point's own seeded chain
    function mean_stats(θu, k)
        X = Matrix{Float64}(undef, n_samples, p)
        idx = Ref(0)
        function record!(cur)
            r = (idx[] += 1)
            for (c, t) in enumerate(terms)
                X[r, c] = Δ[c] == 0.0 ? 0.0 : compute(t, cur)
            end
            return nothing
        end
        chain = _simulate_chain(net, terms, ref, Vector{Float64}(θu); n_sim=n_samples,
                                burnin=burnin, interval=interval,
                                max_val=fixed_support ? mv : nothing, start=mv,
                                max_doublings=max_doublings,
                                rng=Random.Xoshiro(seeds[k]),
                                context="count_mcmle (path sampling, grid point $k)",
                                on_boundary=:error, stat! = record!)
        tops[k] = chain.max_val
        vars[k] = _batch_var_of_mean(X * Δ)
        return vec(sum(X, dims=1)) ./ n_samples
    end
    # The integral is ERGM.jl's (composite Simpson over the grid, the rungs on
    # separate tasks); its Monte-Carlo variance combines the per-rung batch
    # variances with the same Simpson weights (1, 4, 2, ..., 4, 1) / (3·nrungs)
    integral = ERGM.Extension.bridge_integrate(mean_stats, θ0, θ; rungs=nrungs, threaded=true)
    variance = 0.0
    for r in 1:(nrungs + 1)
        w = ((r == 1 || r == nrungs + 1) ? 1.0 : (iseven(r) ? 4.0 : 2.0)) / (3 * nrungs)
        variance += w^2 * vars[r]
    end

    # Exact log-likelihood of the dyad-independent reference model at y_obs,
    # on a support no chain needed to exceed
    support = _support(ref, maximum(tops))
    D = _count_design(model, support)
    ll0 = _count_derivatives(D, collect(1:p), fill(true, length(D.support),
                                                   length(D.n_tot)))(θ0)[1]
    return ll0 + dot(Δ, g_obs) - integral, sqrt(max(variance, 0.0))
end

# =============================================================================
# Simulation
# =============================================================================

# Gibbs-sweep defaults through THE dyad-scaled rule shared with ERGM.jl
# (`ERGM.Extension.mcmc_defaults`), so every sampler of the family burns in alike. ERGM's budget is in single-
# dyad toggles; a Gibbs sweep visits every dyad once, so the same budget is
# `cld(toggles, n_dyads)` sweeps: 20 sweeps of burn-in and
# `cld(max(100, n_dyads ÷ 10), n_dyads)` sweeps between retained draws.
function _gibbs_defaults(n_dyads::Int)
    nd = max(n_dyads, 1)
    d = ERGM.Extension.mcmc_defaults(nd)
    return (burnin=cld(d.burnin, nd), interval=cld(d.interval, nd))
end

"""
    simulate_count_ergm(result::CountERGMResult; n_sim=1, burnin=nothing,
                        interval=nothing, max_val=nothing, max_doublings=8,
                        rng=Random.default_rng()) -> Vector{Network}
    simulate_count_ergm(net, terms, coefficients; reference=PoissonReference(),
                        n_sim=1, burnin=nothing, interval=nothing, max_val=nothing,
                        max_doublings=8, rng=Random.default_rng()) -> Vector{Network}

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

# Support of an unbounded reference

A bounded reference is sampled on its own support. For `PoissonReference` and
`GeometricReference` the enumerated support `0:max_val` is a device, and the
sampler never lets it shape the draws silently:

- `max_val=nothing` (default): the support is **adaptive**. It starts at
  `max(10, 2·largest seed count)` (a fit's own `max_val` in the first form)
  and is doubled whenever a dyad's conditional puts more than 1e-10 of its
  mass on the top value, before that dyad is drawn. If it reaches
  `2^max_doublings` times the start, the call throws an `ArgumentError`: a
  chain that keeps climbing is the signature of a model that is not
  normalisable at these coefficients — which every dyad conditional being
  proper does not rule out (a positive `mutual.product` or squared-strength
  coefficient under a Poisson reference).
- `max_val=k`: the family **truncated** at `0:k`, deliberately. If any
  conditional of the retained sweeps puts more than `BOUNDARY_MASS_TOL` on
  `k`, the draws are returned with a warning giving the share that landed on
  the bound. (`gof` and the parametric bootstrap refuse such a chain.)
- From a fit: a fit whose `max_val` the caller fixed simulates that truncated
  family; an error-controlled fit is simulated adaptively; a fit whose
  joint-support probe found a mode on the bound (`fit.boundary_mode`) is
  refused unless `max_val=k` is passed.

A network with masked (unobserved) dyads is refused
(`NetworkCore.require_observed`): the sampler would otherwise start from, and
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
using NetworkCore, ERGMCount, Random
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
                             max_doublings::Int=_MAX_SUPPORT_DOUBLINGS,
                             rng::Random.AbstractRNG=Random.default_rng())
    model = result.model
    chain = _simulate_chain(model.network, model.terms, model.reference,
                            result.coefficients; n_sim=n_sim, burnin=burnin,
                            interval=interval,
                            _result_support(result, max_val, "simulate_count_ergm")...,
                            max_doublings=max_doublings, rng=rng,
                            context="simulate_count_ergm", on_boundary=:warn)
    return chain.networks
end

# How a chain started from a fit treats the count support (`max_val`, `start`
# keywords of `_simulate_chain`):
# - an explicit `max_val=k` is the caller's truncated family on `0:k`;
# - a fit whose own `max_val` was fixed by the caller simulates that family;
# - an error-controlled fit stands for the UNBOUNDED family, so the chain
#   starts at the fitted bound and widens the support whenever a conditional
#   reaches it — unless the fit's boundary-mode probe found the joint
#   distribution escaping the bound, where there is nothing to draw from.
function _result_support(result::CountERGMResult, max_val, context::AbstractString)
    max_val === nothing || return (max_val=max_val, start=nothing)
    result.support_control === :fixed &&
        return (max_val=result.max_val, start=nothing)
    result.improper && throw(ArgumentError(
        "$context: " * _improper_message(_improper_direction(result.model,
                                                             result.coefficients),
                                         result.model.reference) *
        ". There is no distribution to simulate from. Pass `max_val=k` to draw " *
        "from the family truncated at `0:k` explicitly, or change the model."))
    result.boundary_mode && throw(ArgumentError(
        "$context: " * _boundary_mode_caveat(result) * " There is no settled " *
        "distribution to simulate from. Pass `max_val=k` to draw from the " *
        "family truncated at `0:k` explicitly, or change the model."))
    return (max_val=nothing, start=result.max_val)
end

function simulate_count_ergm(net::Network, terms::Tuple,
                             coefficients::AbstractVector{<:Real};
                             reference::AbstractReferenceMeasure=PoissonReference(),
                             n_sim::Int=1,
                             burnin::Union{Int, Nothing}=nothing,
                             interval::Union{Int, Nothing}=nothing,
                             max_val::Union{Int, Nothing}=nothing,
                             max_doublings::Int=_MAX_SUPPORT_DOUBLINGS,
                             rng::Random.AbstractRNG=Random.default_rng())
    # Coefficients that make the unbounded model improper whatever the data
    # have no distribution to draw from; a chain would sit metastably near
    # its start or climb to the doubling cap. `max_val=k` is the written
    # opt-in to the truncated family.
    if max_val === nothing
        dir = _improper_direction(terms, coefficients, reference, Int(nv(net)),
                                  is_directed(net))
        dir === nothing || throw(ArgumentError(
            "simulate_count_ergm: at these coefficients " *
            replace(_improper_message(dir, reference),
                    "the fitted model is" => "the model is") *
            ". Pass `max_val=k` to draw from the family truncated at `0:k` " *
            "explicitly, or change the coefficients."))
    end
    chain = _simulate_chain(net, terms, reference, Vector{Float64}(coefficients);
                            n_sim=n_sim, burnin=burnin, interval=interval,
                            max_val=max_val, max_doublings=max_doublings, rng=rng,
                            context="simulate_count_ergm", on_boundary=:warn)
    return chain.networks
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

# A conditional that puts more than this on the top support value makes an
# adaptive chain widen its support before drawing: the truncation an unbounded
# reference needs is then invisible at the precision of the draws (a Poisson
# tail past a top value of mass 1e-10 is smaller still).
const _GROW_MASS_TOL = 1e-10

# The enumerated support of a chain over an UNBOUNDED reference, with what the
# chain saw at its top value: the truncation diagnostics of the sampler.
mutable struct _ChainSupport{R<:AbstractReferenceMeasure}
    const ref::R
    support::UnitRange{Int}
    log_h::Vector{Float64}
    η::Vector{Float64}       # the conditional's weights
    buf::Vector{Float64}     # one term's profile
    const adaptive::Bool     # widen on demand (`max_val=nothing`)
    const cap::Int           # an adaptive chain never widens past this
    updates::Int             # dyad draws since the last reset
    top_draws::Int           # ... of which landed on the top value
    max_top_mass::Float64    # largest conditional mass on the top value
    growths::Int
end

function _ChainSupport(ref::AbstractReferenceMeasure, mv::Int; adaptive::Bool, cap::Int)
    support = _support(ref, mv)
    S = length(support)
    return _ChainSupport(ref, support, [log_reference(ref, y) for y in support],
                         Vector{Float64}(undef, S), Vector{Float64}(undef, S),
                         adaptive, cap, 0, 0, 0.0, 0)
end

function _widen!(st::_ChainSupport)
    mv = min(st.cap, 2 * last(st.support))
    st.support = _support(st.ref, mv)
    S = length(st.support)
    st.log_h = [log_reference(st.ref, y) for y in st.support]
    st.η = Vector{Float64}(undef, S)
    st.buf = Vector{Float64}(undef, S)
    st.growths += 1
    return st
end

function _reset_counts!(st::_ChainSupport)
    st.updates = 0
    st.top_draws = 0
    st.max_top_mass = 0.0
    return st
end

struct _ChainEscape <: Exception
    max_val::Int
    mass::Float64
end

# One dyad update of a chain over an unbounded reference: the same draw as
# `_gibbs_update_dyad!` (one uniform, same order), with the conditional mass on
# the top support value measured first. An adaptive chain widens the support
# until that mass is below `_GROW_MASS_TOL`; at its cap it stops (`_ChainEscape`).
function _gibbs_update_tracked!(rng::Random.AbstractRNG, current::Network, weights,
                                terms::Tuple, θ, i::Integer, j::Integer, st::_ChainSupport)
    old = dyad_value(current, weights, i, j)
    total = _dyad_conditional!(st.η, st.buf, terms, θ, st.log_h, st.support,
                               current, weights, i, j, old)
    top = @inbounds st.η[end] / total
    if st.adaptive
        while top > _GROW_MASS_TOL
            last(st.support) >= st.cap && throw(_ChainEscape(last(st.support), top))
            _widen!(st)
            total = _dyad_conditional!(st.η, st.buf, terms, θ, st.log_h, st.support,
                                       current, weights, i, j, old)
            top = @inbounds st.η[end] / total
        end
    end
    idx = _draw_index(rng, st.η, total)
    st.updates += 1
    st.top_draws += (idx == lastindex(st.η))
    top > st.max_top_mass && (st.max_top_mass = top)
    new_val = @inbounds st.support[idx]
    new_val == old && return new_val
    _set_dyad!(current, weights, i, j, new_val)
    return new_val
end

function _gibbs_sweep_tracked!(rng::Random.AbstractRNG, current::Network, weights,
                               terms::Tuple, θ, st::_ChainSupport)
    n = Int(nv(current))
    directed = is_directed(current)
    for i in 1:n
        for j in (directed ? (1:n) : ((i + 1):n))
            i == j && continue
            _gibbs_update_tracked!(rng, current, weights, terms, θ, i, j, st)
        end
    end
    return current
end

function _escape_message(context::AbstractString, ref, e::_ChainEscape, start::Int)
    return "$context: the Gibbs chain reached the cap of its adaptive count " *
           "support (max_val = $(e.max_val), from $start) with conditional mass " *
           "$(_fmt3(e.mass)) still on the top value. $(nameof(typeof(ref))) is " *
           "unbounded, and a chain that keeps climbing is the signature of a " *
           "model that is NOT normalisable at these coefficients — every dyad " *
           "conditional can be proper while the joint is not (e.g. a positive " *
           "`mutual.product` or squared-strength coefficient under a Poisson " *
           "reference: the statistic outgrows log y!), see Krivitsky (2012, " *
           "sec. 3) — or of one whose typical counts are far above the starting " *
           "bound. There is no distribution to draw from in the first case; in " *
           "the second raise `max_doublings`. Pass `max_val=k` to simulate the " *
           "family truncated at `0:k` explicitly."
end

function _boundary_draw_message(context::AbstractString, ref, st::_ChainSupport)
    mv = last(st.support)
    return "$context: the count support was fixed at max_val = $mv, and the " *
           "chain put up to $(_fmt3(100 * st.max_top_mass))% of a dyad's " *
           "conditional mass on that value " *
           "($(_fmt3(100 * st.top_draws / max(st.updates, 1)))% of the retained " *
           "sweeps' dyad draws landed on it; tolerance " *
           "$(_fmt2(100 * BOUNDARY_MASS_TOL))%). $(nameof(typeof(ref))) is " *
           "unbounded: these are draws from the family TRUNCATED at 0:$mv, not " *
           "from the unbounded model."
end

# The chain behind `simulate_count_ergm`, `gof`, the parametric bootstrap and
# the MCMLE. Returns the retained networks and the support the chain ended on
# (`max_val`), with the fraction of retained-sweep dyad draws at the top value
# and the largest conditional mass seen there (both 0 for a bounded reference,
# whose support is the model's own).
#
# `max_val=k` fixes the truncation of an unbounded reference; `max_val=nothing`
# makes it adaptive from `start` (default `max(10, 2·largest seed count)`),
# capped at `2^max_doublings` times that. `on_boundary` says what a FIXED bound
# that shapes the draws does: `:warn` (the caller asked for the truncated
# family) or `:error` (the bootstrap and GOF, which would otherwise report
# numbers from a degenerate chain). `stat!`, when given, is called with every
# retained state instead of copying it (`networks` is then empty).
function _simulate_chain(net0::Network, terms::Tuple, ref::AbstractReferenceMeasure,
                         θ::Vector{Float64};
                         n_sim::Int, burnin, interval,
                         max_val::Union{Int, Nothing},
                         start::Union{Int, Nothing}=nothing,
                         max_doublings::Int=_MAX_SUPPORT_DOUBLINGS,
                         rng::Random.AbstractRNG, context::AbstractString,
                         on_boundary::Symbol=:warn,
                         stat!::F=nothing, state=nothing) where {F}
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
    max_val === nothing || max_val >= 1 ||
        throw(ArgumentError("$context: max_val must be at least 1 (got $max_val)"))
    max_doublings >= 1 || throw(ArgumentError(
        "$context: max_doublings must be at least 1 (got $max_doublings)"))
    # The seed's counts are read by every term's first conditional: the same
    # rules as for a fit (every edge weighted, integers, sign per reference)
    lo = _validate_count_weights(net0, ref)

    # `copy` is the ONE copier of a Network (graph and attribute dicts
    # duplicated, vertex/edge attributes preserved); the chain state is a
    # copy of the seed and every retained draw a copy of the chain
    current = state === nothing ? copy(net0) : state
    # Typed snapshot of the :weight edge attribute (NetworkCore.jl's typed
    # accessor), maintained incrementally alongside the network so the hot
    # loop never reads the untyped attribute Dict.
    weights = get_edge_attribute(current, :weight, Int)

    truncating = is_truncating(ref)
    adaptive = truncating && max_val === nothing
    mv0 = max_val === nothing ?
          max(something(start, 0), _default_max_val(current, weights)) : max_val
    # ... and a support that reaches below zero refuses the terms R refuses
    # on negative weights, whatever the seed holds
    _validate_negative_terms(terms, min(lo, first(_support(ref, mv0))), context)

    n = Int(nv(current))
    n_dyads = is_directed(current) ? n * (n - 1) : n * (n - 1) ÷ 2
    d = _gibbs_defaults(n_dyads)
    burnin = something(burnin, d.burnin)
    interval = something(interval, d.interval)
    burnin >= 0 || throw(ArgumentError("$context: burnin must be non-negative (got $burnin)"))
    interval >= 1 || throw(ArgumentError("$context: interval must be at least 1 (got $interval)"))

    networks = Vector{typeof(current)}()
    if !truncating
        # A bounded reference: the support is the model's own, nothing to track
        support = _support(ref, mv0)
        log_h = [log_reference(ref, y) for y in support]
        η = Vector{Float64}(undef, length(support))
        buf = Vector{Float64}(undef, length(support))
        for sweep in 1:(burnin + n_sim * interval)
            _gibbs_sweep!(rng, current, weights, terms, θ, support, log_h, η, buf)
            if sweep > burnin && (sweep - burnin) % interval == 0
                stat! === nothing ? push!(networks, copy(current)) : stat!(current)
            end
        end
        return (networks=networks, max_val=last(support), top_fraction=0.0,
                max_top_mass=0.0, state=current)
    end

    st = _ChainSupport(ref, mv0; adaptive=adaptive,
                       cap=adaptive ? mv0 * 2^max_doublings : mv0)
    try
        for sweep in 1:(burnin + n_sim * interval)
            # the diagnostics describe the retained part of the chain
            sweep == burnin + 1 && _reset_counts!(st)
            _gibbs_sweep_tracked!(rng, current, weights, terms, θ, st)
            if sweep > burnin && (sweep - burnin) % interval == 0
                stat! === nothing ? push!(networks, copy(current)) : stat!(current)
            end
        end
    catch e
        e isa _ChainEscape || rethrow()
        throw(ArgumentError(_escape_message(context, ref, e, mv0)))
    end
    if !adaptive && st.max_top_mass > BOUNDARY_MASS_TOL
        msg = _boundary_draw_message(context, ref, st)
        on_boundary === :error && throw(ArgumentError(
            msg * " Statistics computed from such a chain describe the bound, " *
            "not the model, so this is refused: refit or simulate with a larger " *
            "`max_val` (or the adaptive default), or — if the bounded family is " *
            "the model — use a bounded reference (`DiscUnifReference`, " *
            "`BinomialReference`)."))
        @warn msg * " Raise `max_val`, or omit it for the adaptive support." maxlog = 1
    end
    return (networks=networks, max_val=last(st.support),
            top_fraction=st.top_draws / max(st.updates, 1),
            max_top_mass=st.max_top_mass, state=current)
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
# typed snapshot, and that costs what NetworkCore.jl's own edge/attribute
# mutation costs.
function _gibbs_update_dyad!(rng::Random.AbstractRNG, current::Network, weights,
                             terms::Tuple, θ, i::Integer, j::Integer, support::UnitRange{Int},
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
function _set_dyad!(net::Network, weights, i::Integer, j::Integer, y::Int)
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
                            net, weights, i::Integer, j::Integer, old::Int)
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
    counts[1 - lo] += _n_empty_dyads(net)
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

This is a method of the shared `NetworkCore.gof` generic; it returns the
shared `NetworkCore.GOFResult` (observed value, simulation envelope, and
two-sided Monte-Carlo p-value per level). The simulation inherits
`simulate_count_ergm`'s refusals (masked dyads, non-finite coefficients, a
model whose chain does not settle on a finite support) and adds one: a
caller-fixed `max_val` on which the chain puts more than `BOUNDARY_MASS_TOL`
of a conditional is an `ArgumentError`, not a warning — the p-values would
describe the bound. For an MCMLE fit the model statistics are matched in
expectation by construction; the count-value panel is the informative one.

# Keyword Arguments
- `n_sim::Int=100`: Number of simulated networks
- `burnin`, `interval`, `max_val`, `rng`: passed to
  [`simulate_count_ergm`](@ref)

# Example
```julia
using NetworkCore, ERGMCount, Random
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
    # A fixed bound that shapes the draws is refused here (`on_boundary=:error`):
    # p-values from a chain pinned to its truncation describe the bound
    sims = _simulate_chain(net, terms, result.model.reference, result.coefficients;
                           n_sim=n_sim, burnin=burnin, interval=interval,
                           _result_support(result, max_val, "gof")...,
                           rng=rng, context="gof", on_boundary=:error).networks

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
# Precompile workload. Before it, the first
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
            # The README's Quick Start, call for call: its default fit of a
            # dyad-dependent formula is the MCMLE (`method=:auto`), so a
            # workload that only ran MPLEs left that path to compile at the
            # user's first call
            _pc_r = Random.Xoshiro(1)
            _pc_q = network(20; directed=true)
            for i in 1:20, j in 1:20
                if i != j && rand(_pc_r) < 0.25
                    add_edge!(_pc_q, i, j)
                    set_edge_attribute!(_pc_q, :weight, i, j, rand(_pc_r, 1:4))
                end
            end
            _pc_terms = [SumTerm(), NonzeroTerm(), CountMutualTerm()]
            _pc_mle = fit_ergm_count(_pc_q, _pc_terms; reference=PoissonReference(),
                                     rng=Random.Xoshiro(2))
            coef(_pc_mle); coeftable(_pc_mle); confint(_pc_mle); loglikelihood(_pc_mle)
            coefnames(_pc_mle)
            sprint(show, MIME"text/plain"(), coeftable(_pc_mle)); sprint(show, _pc_mle)
            _pc_mp = fit_ergm_count(_pc_q, _pc_terms; method=:mple)
            coeftable(_pc_mp); aic(_pc_mp)
            sprint(show, MIME"text/plain"(), coeftable(_pc_mp)); sprint(show, _pc_mp)
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
                                         reference=BinomialReference(5), method=:mple)
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
                # The MCMLE path (Gibbs chain with recorded statistics, Newton
                # step, covariance, path sampling), on a tiny budget
                try
                    _pc_ml = fit_ergm_count(_pc_net, [SumTerm(), _pc_dep];
                                            reference=BinomialReference(5), method=:mcmle,
                                            n_samples=32, maxiter=2, bridge_rungs=2,
                                            bridge_samples=16, rng=_pc_rng)
                    sprint(show, _pc_ml)
                catch
                end
            end
        end
    end
end

end # module
