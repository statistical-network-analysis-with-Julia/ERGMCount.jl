# Estimation

ERGMCount.jl estimates count ERGM parameters by Maximum Pseudo-Likelihood Estimation (MPLE) or by Monte-Carlo Maximum Likelihood (MCMLE — the estimator of R's `ergm.count`). The default, `method=:auto`, chooses as R does: the MPLE for a dyad-independent formula, where it is the exact MLE, and the MCMLE for a dyad-dependent one. This page covers both procedures, configuration, diagnostics, and best practices.

## Overview

The estimation process follows these steps:

1. **Compute observed statistics**: Calculate all term values on the observed network
2. **Construct pseudo-likelihood**: Condition on the rest of the network for each dyad
3. **Optimize**: Find parameters that maximize the pseudo-likelihood via Newton-Raphson

## Maximum Pseudo-Likelihood Estimation

### Why MPLE?

The full likelihood of a count ERGM requires computing a normalizing constant that sums over all possible valued networks -- a combinatorially intractable problem. MPLE avoids this by approximating the joint likelihood with a product of conditional likelihoods:

$$\text{PL}(\theta) = \prod_{(i,j)} P(Y_{ij} = y_{ij} \mid Y_{-ij} = y_{-ij}; \theta)$$

Each conditional probability is tractable because it depends on the reference measure and the change statistics for a single dyad.

### How It Works

For each dyad $(i,j)$, the conditional distribution of $y_{ij}$ given the rest of the network is:

$$P(Y_{ij} = y \mid Y_{-ij}; \theta) \propto h(y) \exp\left(\theta^\top \Delta g_{ij}(y)\right)$$

Where $\Delta g_{ij}(y)$ is the vector of change statistics when $y_{ij}$ takes value $y$.

MPLE maximizes the log-pseudo-likelihood with the ecosystem's shared
`NetworkCore.newton_fit` optimizer (Newton-Raphson with step-halving). Dyads with
identical conditionals are compressed into one row, so a dyad-independent
model costs the same whatever the network size.

For a **dyad-independent** model (`SumTerm`, `NonzeroTerm`,
`GreaterthannTerm`, `CountAtleastnTerm`, `SmallerthanTerm`, `EqualToTerm`,
`InIntervalTerm` only) each dyad's conditional is its marginal, so the pseudo-likelihood *is* the likelihood and the MPLE is the
exact MLE. (`AtmostTerm`, `CMPTerm` and `SumTerm(pow=)` are dyad-independent
too.) For a dyad-dependent model the MPLE is a different estimator from the
one R reports — see [MPLE and MCMLE](@ref) below.

## Fitting a Model

The examples below assume the packages are loaded and a small count-valued
network exists:

```julia
using NetworkCore, ERGMCount
using ERGM: compute, name   # generic term interface shared with ERGM.jl
using Statistics: mean, std
using Random

rng = Xoshiro(1)
net = network(12; directed=true)
for i in 1:12, j in 1:12
    if i != j && rand(rng) < 0.2
        add_edge!(net, i, j)
        set_edge_attribute!(net, :weight, i, j, rand(rng, 1:5))
    end
end
terms = [SumTerm(), NonzeroTerm(), CountMutualTerm()]
```

### Basic Usage

```julia
# `mutual.min` is dyad-dependent: the default fits the MCMLE, as R does
mle = ergm_count(net, terms; reference=PoissonReference(1.0), rng=Xoshiro(2))
mle.method        # :mcmle
# The count MPLE, asked for explicitly (used by most examples on this page)
result = ergm_count(net, terms; reference=PoissonReference(1.0), method=:mple)
```

The `Method:` line of `show` says which estimator ran and why. A keyword of
the estimator `:auto` did not choose is an `ArgumentError` that names the
method taking it (`se=:bootstrap` on a dyad-dependent formula: pass
`method=:mple`).

### Full Options

```julia
result = ergm_count(net, terms;
    reference = PoissonReference(1.0),  # Reference measure
    method = :mple,                     # :auto (default), :mple or :mcmle
    max_val = nothing,                  # nothing: error-controlled support
    support_tol = 1e-3,                 # stopping rule of the support doubling
    se = nothing,                       # :hessian, :bootstrap (n_boot, rng, ...);
                                        # nothing withholds z/p under dependence
    maxiter = 100,                      # Maximum Newton iterations
    warn = true                         # false silences the fit diagnostics
)
```

### Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `reference` | `AbstractReferenceMeasure` | Baseline distribution for edge values | `PoissonReference()` |
| `method` | `Symbol` | `:auto` (R's rule, through `ERGM.resolve_method`: `:mple` for a dyad-independent formula, `:mcmle` otherwise), `:mple` (maximum pseudo-likelihood) or `:mcmle` (Monte-Carlo maximum likelihood; its own keywords are listed under [MPLE and MCMLE](@ref)) | `:auto` |
| `weight` | `Symbol` | The edge attribute holding the counts (R's `response="w"` is `weight=:w`); any name but `:weight` fits a copy of the network | `:weight` |
| `max_val` | `Int` or `nothing` | Fixed top of the enumerated support for an unbounded reference; `nothing` chooses it adaptively | `nothing` |
| `support_tol` | `Float64` | Stop doubling the support when the estimates move by at most this many standard errors *and* at most this many expected dyads sit past the previous bound (and no dyad has appreciable mass at the new bound) | `1e-3` |
| `max_doublings` | `Int` | Cap on the support doublings; hitting it reports `support_stable = false` with a warning | `8` |
| `se` | `Symbol` or `nothing` | MPLE only. `:hessian` (inverse pseudo-Hessian) or `:bootstrap` (parametric bootstrap; `n_boot`, `boot_burnin`, `boot_interval`, `rng`). `nothing` reports the Hessian standard errors and, for a dyad-dependent model fit with `method=:mple`, no z value, p-value or interval | `nothing` |
| `maxiter` | `Int` | Maximum Newton-Raphson iterations | `100` |
| `warn` | `Bool` | `false` silences the boundary, truncation, non-convergence, separation and conditioning warnings (all still recorded on the result) | `true` |

`terms` may be a `Vector`, a `Tuple` or a single term.

### The count support is error-controlled

A Poisson or geometric reference is unbounded, but the estimator enumerates
each dyad's conditional over `0:max_val`. With `max_val` unset the fit starts
at twice the largest observed count (at least 10), refits at twice that bound,
and stops when the doubling moved every estimate by at most `support_tol`
standard errors *and* left at most `support_tol` expected dyads past the
previous bound *and* no dyad puts more than `BOUNDARY_MASS_TOL` of its
conditional mass on the new top value; the fit at the larger bound is what you
get. The result records what happened:

```julia
result.support_control   # :converged (or :fixed, :bounded, :unconverged, :boundary_mode, :improper)
result.support_stable    # false for :unconverged, :boundary_mode and :improper
result.improper          # true when the fitted model is not normalisable, whatever the data
result.boundary_mode     # true when the JOINT distribution has a mode on the bound
result.max_val           # the bound actually used
result.support_delta     # max |Δθ| / SE at the last doubling (NaN for a fixed max_val, 0.0 for a bounded reference)
result.omitted_tail      # expected dyads past the previous bound (same conventions)
```

`:unconverged` — the estimates still moving after `max_doublings` (8)
doublings — is warned about (the warning names the achieved δ and the final
`max_val`), listed in `approximations(result)` and printed by `show`: it is the
signature of a model that is not normalisable on the unbounded support (a
positive `sum` coefficient under a geometric reference, say). The doubling
also stops — `:unconverged`, `support_delta`/`omitted_tail` `NaN`, and the
support line of `show` saying so — at the first bound whose fit did not
converge or has a numerically singular pseudo-Hessian (a collinear design,
separation): there are no settled estimates to compare, so the convergence,
separation and conditioning warnings speak alone and the "not normalisable"
reading is never offered for them. A bounded
reference (`BinomialReference`, `DiscUnifReference`, `DiscUnif2Reference`) has
nothing to control (`:bounded`, `support_delta == 0.0`) and never loops. An
explicit `max_val` keeps a single fit (`:fixed`) with the boundary-mass
diagnostic still computed. On Zachary's karate club the default path doubles
14 → 28 and reproduces the exact MLE to ~1e-12 (the golden fixture pins it at
1e-6); the old fixed default of 14 missed that tolerance.

#### What the support check cannot see: joint normalisability

The doubling compares fits, and every fit looks only at the dyad conditionals
**at the observed network**. A model can pass while its joint distribution
has no normalising constant on the unbounded support. The standard example
(Krivitsky 2012, §3) is a positive `mutual.product` coefficient under a
Poisson reference: given its partner, each dyad is a proper Poisson variable,
but a reciprocated pair at `(k, k)` has weight `exp(θ k²) / k!²`, and `k²`
outgrows `2 log k!` — however small `θ` is. Squared-strength terms
(`nodeOSum`, `nodeISum`, `nodeSum`), `sum(pow=p)` with `p > 1`, and `CMP`
beyond the reference's own `log y!` behave the same way. Two checks run after
every fit under a Poisson or geometric reference.

**The analytic rule.** `log h(y)` is about `−y log y` per dyad for the
Poisson reference and `0` for the geometric one, so a statistic that grows
faster than linearly in the counts outruns it whenever its leading
coefficient is positive. The rule evaluates the leading order of
`log h(y) + θ'g(y)` with one dyad, both dyads of a pair, the dyads of one
actor (a star) and every dyad set to a count `Y → ∞`, summing the
contributions of every super-linear term, so a larger negative
squared-strength term can offset a positive one (`Σ out² ≤ Σ (out + in)²`). A
positive leading coefficient along any of these is a proof that the model is
improper whatever the data, and the fit is reported with

```julia
result.improper          # true
result.support_control   # :improper (a caller-fixed max_val stays :fixed)
result.support_stable    # false
```

Under the geometric reference, which does not decay at all, linear growth
decides too: with no super-linear term the rule sums the exact linear
coefficients along the same configurations, so `sum` ≥ 0 along one dyad, or
`2·θ_sum + θ_mutual > 0` on a reciprocated pair under `mutual(:min)`, is
improper. A user-defined count term, whose growth the rule cannot read,
makes it silent (as does a covariate term's `form=:sum` for linear growth:
its coefficient depends on which dyads grow).

**The probe.** A deterministic probe starts with **every dyad at
`2^max_doublings × max_val`** — the bound at which the adaptive sampler would
give up — and iterates conditional modes. (It used to start at the fitted
`max_val`, where a small positive `mutual.product` coefficient escaped notice:
on a 25-actor network with a fitted coefficient of 0.061 the probe said no at
`0:20` and yes at `0:160`.) If a sweep leaves no dyad on the bound, the bound
is not sticky. If dyads stay there, the configuration the modes stop at is
compared with the observed network by their joint log-weight
`log h(y) + θ'g(y)`: when it weighs more, the truncated joint distribution
puts its weight on its boundary — the model is not normalisable, or its mass
lies far beyond the bound. When it weighs less, it is a local mode the data
outweigh: a proper geometric model with strong reciprocity
(`θ_sum + θ_mutual > 0 > 2·θ_sum + θ_mutual`, which R fits) has the all-top
network as such a mode, and is fitted. The probe covers what the rule cannot
read. A flagged fit is reported with
`result.boundary_mode == true` and `support_control == :boundary_mode`.

Either way the fit warns, `approximations(result)` lists it, `show` prints a
`WARNING:` line, and `se=:bootstrap`, `simulate_count_ergm`, `gof` and
`method=:mcmle` refuse it with an `ArgumentError` — unless `max_val` was
fixed by the caller, which is the truncated family chosen in writing. A
model improper only beyond the probe's bound, and metastable below it, is
not detected.

### Boundary statistics

A statistic at the boundary of its attainable range — `mutual.min` on a
network with no reciprocated pair, `nonzero` on a complete network,
`smallerthan.1` when every dyad carries a count — has no finite
pseudo-likelihood maximizer. As in R `ergm` (its default `drop=TRUE`), the
coefficient is fixed at `-Inf` (`+Inf`) with standard error 0 and p-value 0,
`count_mple` warns with R's sentence ("... are at their smallest attainable
value"), the remaining coefficients are estimated with each dyad's support
restricted to the values the fixed statistic allows (the exact limit of the
pseudo-likelihood), and `dof`/`aic`/`bic` count only the finite
coefficients. Such a fit is never `is_exact`, and `se=:bootstrap`,
`simulate_count_ergm` and `gof` refuse it — there is no finite model to
simulate from. Read a `±Inf` coefficient as a statement about the data and
drop the term. A statistic whose change statistic is 0 on every dyad
(`transitiveweights` on a network where no two ties share an actor) is at
its smallest attainable value too, as in R.

The MCMLE (the default for a dyad-dependent formula) applies the same drop:
the coefficient is fixed at `∓Inf` with R's warning, and the others are the
maximum-likelihood estimates with the statistic held at its bound — the
Gibbs sampler gives every value that would move it zero weight. The
log-likelihood is then not estimated (`NaN`). `drop=false`, R's
`control.ergm(drop=FALSE)`, refuses such a model on either estimator
instead (keeping the term is not implemented).

```julia
acyclic = network(5; directed=true)            # 1→2→3→4→5: nothing reciprocated
for (i, j) in ((1, 2), (2, 3), (3, 4), (4, 5))
    add_edge!(acyclic, i, j); set_edge_attribute!(acyclic, :weight, i, j, 2)
end
b = fit_ergm_count(acyclic, [SumTerm(), CountMutualTerm()]; reference=BinomialReference(3),
                   method=:mple)
coef(b)[2]         # -Inf
dof(b)             # 1
is_exact(b)        # false
```

### MPLE and MCMLE

R's `ergm.count` has no MPLE: it fits every valued model by Monte-Carlo
maximum likelihood. The count MPLE is this package's own estimator. For a
dyad-independent model the two coincide (the MPLE is the exact MLE), which is
why the default `method=:auto` uses the MPLE there and the MCMLE everywhere
else, as R does. For a dyad-dependent model they do not. On `ergm.count`'s `zach` example
(`zach ~ sum + nonzero + transitiveweights("min","max","min")`, Poisson
reference; frozen in `test/fixtures/count_mcmle.toml`):

| | `sum` | `nonzero` | `transitiveweights` |
|:--|--:|--:|--:|
| count MPLE | 0.797 | −4.924 | 0.385 |
| `ergm.count` MCMLE, mean of 9 seeds | 0.597 | −4.779 | 0.562 |
| `ergm.count` seed-to-seed sd | 0.005 | 0.011 | 0.004 |
| `method=:mcmle`, mean of 6 seeds | 0.597 | −4.778 | 0.560 |
| `method=:mcmle` seed-to-seed sd | 0.006 | 0.009 | 0.004 |
| standard error (R) | 0.121 | 0.237 | 0.101 |

The MPLE sits 1.7 and 1.8 standard errors from the MLE on `sum` and
`transitiveweights`. Simulating at the MLE and refitting shows what kind of
difference this is (300 networks drawn with `simulate_count_ergm` at
`ergm.count`'s long-chain estimate, the MPLE refit on each):

| | `sum` | `nonzero` | `transitiveweights` |
|:--|--:|--:|--:|
| mean error of the MPLE | −0.029 | +0.007 | +0.018 |
| sampling standard deviation of the MPLE | 0.137 | 0.238 | 0.121 |
| mean naive (pseudo-Hessian) standard error | 0.113 | 0.287 | 0.058 |
| coverage of naive 95 % Wald intervals | 0.91 | 0.98 | 0.70 |
| coverage with `se=:bootstrap` (100 networks, `n_boot=60`) | 0.96 | 0.97 | 0.97 |

The MPLE is nearly unbiased but less efficient than the MLE (whose standard
error for `transitiveweights` is 0.101), and its naive standard error for the
dependence term is half the true sampling standard deviation.

```julia
mle = ergm_count(net, terms; rng=Xoshiro(2))     # the default for this formula
mle.method                          # :mcmle
mle.converged                       # true
mle.mcmc.convergence                # iterations, step length, t-ratios, Hotelling p, ESS
stderror(mle)                       # inverse Fisher information + Monte-Carlo error
mle.mcmc.mc_std_errors              # the Monte-Carlo part alone
loglikelihood(mle)                  # the log-likelihood, by path sampling
confint(mle)
```

How it works: the iteration is ERGM.jl's (`ERGM.Extension.mcmle_solve`, the one
`ERGM.mcmle` runs); this package supplies the sampler. Starting from the
MPLE, each iteration draws `n_samples` networks at the current coefficients
with the Gibbs sampler (the chain continues across iterations) and moves the
coefficients by the Monte-Carlo Newton step `γ·Σ̂⁻¹(g(y_obs) − ḡ)`, with
Hummel's step length `γ`. At full step length the stopping rule is applied:
by default R `ergm` 4's confidence test — with 99 % confidence the estimating
equations at the updated coefficients lie inside a tolerance region a tenth
the size of the statistics' covariance — which enlarges the sample when it
fails near the solution. The draws are thinned adaptively
(`ERGM.Extension.ess_sample`): the interval between them doubles until their effective
sample size reaches `effective_size`. The standard errors and the convergence
diagnostics are computed from one more sample drawn at the returned
coefficients, so the Fisher information is the one at the estimate.

| Keyword | Description | Default |
|---|---|---|
| `n_samples` | Networks drawn per iteration (doubled on demand) | `1024` |
| `burnin`, `interval` | Gibbs sweeps before the first draw of an iteration and (to start with) between draws | dyad-scaled rule (20, 1) |
| `effective_size` | Effective sample size the thinning must reach (`nothing`: fixed interval) | `n_samples ÷ 2` |
| `maxiter` | Maximum MCMLE iterations | `60` |
| `termination` | `:confidence` (R `ergm` 4's rule; `conv_precision`, `conv_confidence`) or `:hotelling` (`conv_threshold`, `hotelling_alpha`) | `:confidence` |
| `max_n_samples` | Cap on the boosted sample | `16·n_samples` |
| `gamma0`, `max_step_norm` | First step length and the cap on a Newton step | `0.1`, `5.0` |
| `bridge_rungs`, `bridge_samples` | Path-sampling grid for the log-likelihood (`0` skips it) | `16`, `max(64, n_samples ÷ 4)` |
| `max_val`, `max_doublings`, `support_tol` | The count support, as for the MPLE | adaptive |
| `init` | Starting coefficients | the MPLE |

`loglikelihood(mle)` includes the reference measure, so it is comparable with
the log-likelihood of a dyad-independent MPLE fit on the same network. R's
`logLik` for a valued ERGM is relative to the reference (`θ = 0`); subtract
`Σ log h(y_obs) − N·log Σ_y h(y)` to compare. A dyad-independent model under
`method=:mcmle` returns the exact fit with no Monte Carlo. A boundary
statistic is dropped as in R (see "Boundary statistics"); a separated or
unconverged MPLE start, and a start that is improper or has `boundary_mode`,
are refused with an `ArgumentError`.

### Standard errors: Hessian vs parametric bootstrap

`se=:hessian` reports the inverse negative pseudo-Hessian. For a
dyad-independent model the pseudo-likelihood is the likelihood and these are
the usual MLE standard errors. For a dyad-dependent model the
pseudo-likelihood multiplies overlapping conditionals as if independent, so
they are **too small** (see the coverage figures above).

With the default `se=nothing`, a dyad-dependent MPLE fit (`method=:mple`) therefore reports its
estimates and naive standard errors and **withholds the inference built on
them**: `z_values` and `p_values` are `NaN`, `show` explains why, `confint`
throws an `ArgumentError`, `result.inference_withheld` is `true`, and
`NetworkCore.approximations(result)` lists it. Passing `se=:hessian` explicitly
is the written opt-in to the naive Wald table (printed with its caveat). A
dyad-independent fit is unaffected.

`se=:bootstrap` Gibbs-simulates `n_boot` networks at the estimates with
`simulate_count_ergm`, refits the count MPLE on each, and reports the
empirical covariance of the refits on the shared `NetworkCore.bootstrap_cov`
loop (`n_boot`, `boot_burnin`, `boot_interval`, `rng`, exactly as
`ERGM.mple`). The point estimates are unchanged; only the covariance is
replaced. Replicates without a finite, converged MPLE are excluded with one
aggregate warning and kept as `NaN` rows of `result.boot_replicates`; the
warning, `show` and `approximations` then say "The standard errors are
conditional on a finite refit: the excluded replicates are the extreme ones,
so the standard errors are biased downward."

```julia
robust = ergm_count(net, terms; reference=PoissonReference(1.0), method=:mple,
                    se=:bootstrap, n_boot=50, rng=Xoshiro(2))
coef(robust) == coef(result)        # true
stderror(robust)                    # typically larger than the Hessian ones for mutual.min
NetworkCore.se_method(robust)          # :bootstrap
confint(robust)                     # reported: the bootstrap covariance is calibrated
```

The bootstrap replaces the covariance, not the estimate: `show(robust)` still
says that the point estimates are pseudo-likelihood estimates, not the MLE.
The simulating chain widens an error-controlled support on demand; with a
caller-fixed `max_val` the bootstrap is refused if the chain puts more than
`BOUNDARY_MASS_TOL` of a conditional on the bound.

### Model comparison with (pseudo-)AIC/BIC

`aic(result) = -2·loglik + 2·dof` and `bic(result) = -2·loglik +
dof·log(nobs)` are computed from the maximised **pseudo**-log-likelihood
(`nobs` is the number of dyads, `dof` the number of finite coefficients).
They are exact information criteria only for a dyad-independent model;
otherwise they are pseudo-likelihood criteria, comparable across models fit
on the same network with the same reference and the same support, and not
across references or supports (the support is part of the estimand).

### Non-convergence is loud

A Newton iteration that exhausts `maxiter` or cannot move is warned about
(the warning names `maxiter`, `tol`, the pseudo-score norm and the two ways
out), recorded as `result.converged == false` with `result.iterations` and
`result.gradient_norm`, listed in `NetworkCore.approximations(result)`, printed
by `show` ("Converged: false" plus the caveat), and never `is_exact`.
Running out of `maxiter` is reported, not silently topped up. `warn=false`
silences the fit diagnostics without changing what the result records.

### Separation: the MPLE does not exist

A boundary statistic is one column at its extreme. The pseudo-likelihood can
also be flat along a *combination* of the statistics — separation, R's "The
MPLE does not exist!" — which no single-column test sees. The textbook case
is `sum + nonzero` on a network whose every count is 0 or 1: `sum − nonzero`
is at its minimum on every dyad, so the objective keeps rising as
θ_sum → −∞, θ_nonzero → +∞ and Newton stops somewhere on that asymptote with
coefficients around ±20 and standard errors in the tens of thousands.

The count pseudo-likelihood is a conditional logit: each compressed row is a
stratum whose alternatives are the support values, and the values its dyads
were observed at are the chosen ones. `count_mple` decides existence exactly,
from the data, with the ecosystem's shared verdict
(`NetworkCore.clogit_separation`, the linear programme of R's
`mple.existence`, certified in exact arithmetic), run on the design actually
fitted (after the boundary columns are fixed). It does not depend on where
Newton stopped. A separated fit follows the ecosystem's separation policy: a
warning naming the separating terms (with R's sentence), `converged = false`,
`separated = true` and the names in `separated_terms`, `NaN` z values,
p-values and `confint`, an entry in `approximations` and a caveat in `show`;
it is never `is_exact`, `se=:bootstrap` and `method=:mcmle` refuse it, and
the bootstrap excludes such replicates.

```julia
binary = network(6; directed=false)          # seven ties, every count 1
for (i, j) in ((1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (1, 3), (2, 5))
    add_edge!(binary, i, j); set_edge_attribute!(binary, :weight, i, j, 1)
end
sep = fit_ergm_count(binary, [SumTerm(), NonzeroTerm()]; warn=false)
sep.separated        # true
sep.converged        # false
sep.separated_terms  # ["sum", "nonzero"]
confint(sep)         # NaN: no interval without an MPLE
```

### Collinear statistics: a numerically singular pseudo-Hessian

Two statistics that coincide on the data — `greaterthan(2)` and
`atleast(3)` on integer counts (the same statistic), or `sum` and `nonzero`
on a 0/1 network — leave the pseudo-Hessian singular or nearly so.
`result.hessian_cond` records its condition number; above `1e8`
(`ERGMCount._HESSIAN_COND_TOL`) the fit warns, names the statistics loading
on the flat direction in `result.collinear`, and carries the caveat in
`approximations` and under `Converged:` in `show`. The standard errors along
that direction are meaningless (`NaN` when it is exactly singular); drop or
merge one of the terms.

### Alternative Syntax

`fit_ergm_count` is the standardized entry point (the ecosystem's
`fit_<model>` naming); `ergm_count` is the R-faithful alias of the same
function:

```julia
# The same function
result = fit_ergm_count(net, terms; reference=PoissonReference(), method=:mple)
result = ergm_count(net, terms; reference=PoissonReference(), method=:mple)
```

## Understanding Results

The `CountERGMResult` object contains:

| Field | Type | Description |
|-------|------|-------------|
| `model` | `CountERGMModel` | The fitted model specification (`model.terms` is a `Tuple`) |
| `coefficients` | `Vector{Float64}` | Estimated coefficients (`±Inf` for a boundary statistic) |
| `std_errors` | `Vector{Float64}` | Standard errors (`se_type` says which kind) |
| `z_values`, `p_values` | `Vector{Float64}` | What `coeftable`/`show` print |
| `loglik` | `Float64` | Log-pseudo-likelihood at convergence |
| `converged`, `iterations`, `gradient_norm` | | Whether the Newton iteration converged, how many iterations the final fit ran, and the pseudo-score norm at the estimates |
| `separated`, `separated_terms` | `Bool`, `Vector{String}` | The MPLE does not exist along a combination of the statistics (then `converged` is `false` and z, p and intervals are `NaN`), and the terms that carry it |
| `hessian_cond`, `collinear` | `Float64`, `Vector{String}` | Condition number of the negative pseudo-Hessian, and the statistics loading on its flat direction when it exceeds `1e8` |
| `max_val`, `truncated`, `support_control`, `support_stable`, `support_tol`, `support_delta`, `omitted_tail`, `boundary_mass` | | The support used and how it was chosen |
| `se_type`, `boot_replicates` | | `:hessian`/`:bootstrap`/`:mcmc`, and the bootstrap refits |
| `method`, `mcmc` | `Symbol`, `NamedTuple` or `nothing` | the estimator that ran, `:mple` or `:mcmle`; the MCMLE's record (`convergence`, `mc_std_errors`, `samples`, `n_samples`, `burnin`, `interval`, `loglik_mc_se`, `start`) |
| `inference_withheld` | `Bool` | The MPLE of a dyad-dependent model with the default `se`: `z_values`/`p_values` are `NaN`, `confint` refuses |
| `boundary_mode` | `Bool` | The joint distribution, probed at `2^max_doublings × max_val`, has a mode on the bound that outweighs the observed network |
| `improper` | `Bool` | The fitted model is not normalisable on the unbounded support, whatever the data (the analytic rule) |

Prefer the StatsAPI verbs to the fields: `coef`, `stderror`, `vcov`,
`confint`, `loglikelihood`, `nobs`, `dof`, `aic`, `bic`, `coeftable`.

### Displaying Results

```julia
println(result)
```

Output (the network above is seeded, so these are the numbers you get):

```text
Count ERGM Results
==================
Reference: PoissonReference(1.0)
Method: mple (maximum pseudo-likelihood: an approximation under dyadic dependence; the default method=:auto fits the MCMLE here, as R's ergm.count does)
Support:   0:40  (TRUNCATED — reference is unbounded; chosen by doubling until max|Δθ| ≤ 0.001·SE and the omitted tail ≤ 0.001; achieved 3.57e-10·SE, tail 1.61e-10)
Boundary mass: 6.97e-32 (max over dyads, at the fitted coefficients)
Pseudo-log-likelihood: -115.8861
AIC: 237.77, BIC: 246.42  (pseudo-likelihood; compare only across models on the same network and support)
Converged: true
Std. errors: inverse pseudo-Hessian

Coefficients:
            Estimate  Std.Error  z value  Pr(>|z|)
sum           1.0529     0.1221      NaN       NaN
nonzero      -4.2415     0.4265      NaN       NaN
mutual.min    0.2183     0.1729      NaN       NaN
---
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Note: z values and p-values are not reported (NaN). This model contains
dyad-dependent terms and was fit by maximum pseudo-likelihood; the
standard errors shown are the naive inverse pseudo-Hessian ones, which
treat dependent dyads as independent and under-cover (95% Wald intervals
covered 0.70-0.91 in simulation), so no test or interval is built on
them. The point estimates are not the MLE that R's ergm.count reports.
For inference refit with method=:mcmle (maximum likelihood) or
se=:bootstrap; se=:hessian requests the naive Wald table explicitly.
```

(The coefficient table is the shared presentation layer from NetworkCore.jl,
used identically across the ecosystem's model packages; it is exactly
`coeftable(result)`.)

### Accessing Results

```julia
coef(result)             # coefficient vector
stderror(result)         # standard errors
vcov(result)             # covariance matrix
confint(mle)             # normal-theory 95% limits, one row per coefficient
                         # (`confint(result)` is refused: dyad-dependent MPLE)
coeftable(result)["sum"] # one labelled row of the table
aic(result), bic(result) # pseudo-likelihood criteria
nobs(result), dof(result)

# Model specification
result.model.terms      # Tuple of terms
result.model.reference  # Reference measure
result.model.network    # Original network
```

## Interpreting Coefficients

### General Interpretation

Coefficients in count ERGMs describe how each unit change in a sufficient statistic affects the log-odds of edge value configurations:

| Coefficient | Meaning |
|-------------|---------|
| $\theta > 0$ | The corresponding statistic is over-represented relative to the reference |
| $\theta < 0$ | The corresponding statistic is under-represented relative to the reference |
| $\theta = 0$ | The statistic value matches the reference expectation |

### Reference-Specific Interpretation

With **Poisson reference** and `SumTerm`:
- The coefficient $\theta_{\text{sum}}$ shifts the conditional Poisson mean
- Conditional mean for each dyad: $\lambda \exp(\theta_{\text{sum}})$
- Example: $\theta_{\text{sum}} = 0.5$ with $\lambda = 1$ gives conditional mean $\exp(0.5) \approx 1.65$

With **Geometric reference** and `SumTerm`:
- Shifts the geometric distribution parameter
- Changes the probability of high vs. low values

### Example Interpretations

| Term | Coefficient | Interpretation |
|------|-------------|----------------|
| SumTerm | 0.3 | Edge values are 35% higher than baseline ($\exp(0.3) \approx 1.35$) |
| NonzeroTerm | -2.5 | Strong penalty on edge existence; sparse network |
| CountMutualTerm | 0.6 | Reciprocal dyads tend to match intensities |
| NodeOSumTerm | 0.01 | Slight heterogeneity in sending activity |

## Convergence

### Checking Convergence

```julia
if result.converged
    println("Model converged")
else
    println("WARNING: Model did not converge")
end
```

### Common Issues

| Issue | Symptom | Solution |
|-------|---------|----------|
| Non-convergence | `converged = false`; a warning naming `maxiter`, `tol` and the pseudo-score norm; listed in `approximations`; `is_exact` false | Raise `maxiter` if Newton ran out of budget (`iterations == maxiter`), otherwise look for a statistic with no finite maximizer (constant over the data, or collinear with another) |
| Boundary statistic | coefficient `±Inf`, SE 0, a warning naming the statistic (MPLE and MCMLE) | The statistic is at its extreme attainable value on every dyad (e.g. `mutual.min` with no reciprocated pair); no finite estimate exists — drop the term, or `drop=false` to refuse |
| Separation | `separated = true`, `converged = false`; a warning "the MPLE does not exist (separation)" naming the terms (`separated_terms`); coefficients around ±20 with astronomical SEs, z, p and `confint` `NaN` | A combination of the statistics is at its extreme on every dyad (`sum + nonzero` when every count is 0 or 1); remove or coarsen a term, or collect more varied counts |
| Improper model | `improper = true`, `support_control == :improper` (warned); simulation, `gof`, `se=:bootstrap` and `method=:mcmle` refused | A positive leading coefficient on a super-linear term (`mutual.product`, a squared strength, `sum(pow>1)`, `CMP` beyond `log y!`) under a Poisson or geometric reference: no distribution has these coefficients. Drop the term, use a bounded reference, or fix `max_val` to work with the truncated family deliberately |
| Collinear statistics | a warning "numerically singular (condition number …)" naming the terms; `hessian_cond > 1e8`, `collinear` non-empty; SEs `NaN` when exactly singular | Two statistics coincide on the data (`greaterthan(2)` and `atleast(3)`; `sum` and `nonzero` on 0/1 counts); drop or merge one |
| Unconverged support | `support_control == :unconverged` (warned) | The truncated fits keep moving with the bound (`support_delta` finite): the unbounded model is not normalisable at these coefficients; reconsider the reference. Or (`support_delta` `NaN`) the doubling stopped at a bound whose fit did not converge or is singular: fix the model first (see the rows above) |
| Negative counts with `transitiveweights`/`cyclicalweights`/`mutual(:geometric)` | `ArgumentError` "may not be used with networks with negative dyad weights" (from `compute`, the model constructor and the simulator) | As in `ergm`: drop the term, or fit a reference whose support is non-negative |
| Some edges without `:weight` | `ArgumentError` naming the first bare edge and how many there are | A bare edge is a data gap, not a count of 1: set every edge's count (`weight=:w` maps R's `response="w"`) |
| A binary ERGM.jl term (`ERGM.Edges()`, `ERGM.Mutual()`, …) | `ArgumentError` naming the count analogue | Count models use `NonzeroTerm()` for `edges`, `CountMutualTerm()` for `mutual`, `TransitiveWeightsTerm()` for `triangle` |
| No `:weight` attribute, or a non-integer / negative one | `ArgumentError` naming the dyad and the rule | Store integer counts under `:weight` (`weight=:w` maps R's `response="w"`); round rates; a negative count needs `DiscUnif2Reference(a < 0, b)` |
| Two-mode network | `ArgumentError` "fits one-mode networks only" | Bipartite count ERGMs are not implemented; nothing is fit |
| Directed-only term on an undirected network | `ArgumentError` naming the term | Use `NodeSumTerm`, or drop the reciprocity/cycle term |

### Handling Non-Convergence

```julia
# Increase iterations
result = ergm_count(net, terms;
    reference=PoissonReference(),
    method=:mple,
    maxiter=500
)

# Check for problematic terms
for (i, term) in enumerate(terms)
    println("$(name(term)): coef=$(round(coef(result)[i], digits=4)), " *
            "SE=$(round(stderror(result)[i], digits=4))")
end
```

## Model Comparison

### Comparing Models

```julia
# Model 1: Basic
terms1 = [SumTerm(), NonzeroTerm()]
result1 = ergm_count(net, terms1; reference=PoissonReference(), method=:mple)

# Model 2: Add reciprocity (the count MPLE too, so both criteria are on the
# pseudo-likelihood scale)
terms2 = [SumTerm(), NonzeroTerm(), CountMutualTerm()]
result2 = ergm_count(net, terms2; reference=PoissonReference(), method=:mple)

println("Model 1: AIC=", round(aic(result1), digits=2), ", BIC=", round(bic(result1), digits=2))
println("Model 2: AIC=", round(aic(result2), digits=2), ", BIC=", round(bic(result2), digits=2))
```

`aic`/`bic` are pseudo-likelihood criteria: compare them only across models
fit on the same network with the same reference and support.

### Comparing Reference Measures

```julia
terms = [SumTerm(), NonzeroTerm()]

result_pois = ergm_count(net, terms; reference=PoissonReference(1.0))
result_geom = ergm_count(net, terms; reference=GeometricReference())

println("Poisson: ", round.(coef(result_pois), digits=3))
println("Geometric: ", round.(coef(result_geom), digits=3))
```

## Simulation-Based Validation

After fitting, validate by simulating networks from the estimated model and comparing summary statistics.

`simulate_count_ergm` is a **Gibbs sampler**: one sweep visits every dyad
once and redraws it from its full conditional `P(y_ij = y | rest) ∝
h(y)·exp(θ'Δg(y))` over the count support. `burnin` and `interval` are
counted in sweeps and default through the dyad-scaled rule shared with
ERGM.jl (`ERGM.Extension.mcmc_defaults`, a budget in single-dyad moves, divided by
the number of dyads): **20 sweeps of burn-in and `cld(max(100, n_dyads ÷
10), n_dyads)` sweeps between retained draws** — 1 for any network with more
than ten nodes — and `gof` and the parametric bootstrap (`boot_burnin`,
`boot_interval`) resolve their defaults the same way.

For an unbounded reference the sampler's support is **adaptive**. It starts
at the fit's bound (the explicit form `simulate_count_ergm(net, terms, coefs)`
starts at `max(10, 2·largest seed count)`) and doubles whenever a dyad's
conditional puts more than 1e-10 of its mass on the top value, before that
dyad is drawn — so a Poisson mean of 25 is simulated as 25, not as the 18 a
fixed `0:20` gives. A chain that reaches `2^max_doublings` times its starting
bound throws an `ArgumentError`: that is how a model that is not normalisable
shows itself. `max_val=k` asks for the family truncated at `0:k` instead; if
any conditional then puts more than `BOUNDARY_MASS_TOL` on `k`, the draws are
returned with a warning that gives the share of draws on the bound, and `gof`
and the bootstrap refuse. Each
retained draw is an independent `copy` of the chain that keeps the seed
network's vertex attributes. A sweep costs ≈ 0.6 µs per dyad and allocates
nothing once the chain has warmed up (see the CHANGELOG for the measured
numbers); a 1000-node directed network is ≈ 0.6 s per sweep.

```julia
# Simulate from fitted model (pass rng for reproducible draws)
using Random
sim_nets = simulate_count_ergm(result;
    n_sim=100,
    burnin=50,        # Gibbs sweeps (default: 20, the dyad-scaled rule shared with ERGM.jl)
    rng=Xoshiro(42)
)

# Compare observed vs simulated
for term in terms
    obs = compute(term, net)
    sim_vals = [compute(term, sn) for sn in sim_nets]
    sim_mean = mean(sim_vals)
    sim_sd = std(sim_vals)

    z = (obs - sim_mean) / sim_sd
    println("$(name(term)): obs=$(round(obs, digits=2)), " *
            "sim=$(round(sim_mean, digits=2)) +/- $(round(sim_sd, digits=2)), " *
            "z=$(round(z, digits=2))")
end
```

A well-fitting model should have $|z| < 2$ for all terms, meaning the observed statistics fall within the range of simulated values.

## Best Practices

1. **Always include SumTerm and NonzeroTerm**: These control baseline scale and density
2. **Check convergence**: Verify `result.converged == true` before interpreting
   (and `result.support_stable`, `result.improper`, `result.boundary_mode`)
3. **Validate with simulation**: Compare observed and simulated network statistics
   (for an MPLE fit of a dyad-dependent model, a mismatch on the model's own
   statistics is the estimator's inefficiency showing; the MCMLE, the
   default for such a model, matches them by construction)
4. **Start simple**: Begin with basic terms, add complexity gradually
5. **Match reference to data**: Choose a reference that reflects your data's properties
6. **Keep the reference fixed**: the only reference parameter is the Poisson rate $\lambda$ (confounded with the `sum` coefficient — leave it at 1 unless a baseline rate is known) and the binomial `trials`; the geometric shape and the binomial success probability come from the estimated `sum` coefficient
7. **Read a `±Inf` coefficient as a statement about the data**: the statistic is at its boundary and no finite estimate exists; drop the term
8. **Inspect standard errors**: NaN or very large SEs indicate estimation problems
