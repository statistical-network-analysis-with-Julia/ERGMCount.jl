# ERGMCount.jl


[![Network Analysis](https://img.shields.io/badge/Network-Analysis-orange.svg)](https://github.com/statistical-network-analysis-with-Julia/ERGMCount.jl)
[![Build Status](https://github.com/statistical-network-analysis-with-Julia/ERGMCount.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/statistical-network-analysis-with-Julia/ERGMCount.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://statistical-network-analysis-with-Julia.github.io/ERGMCount.jl/dev/)
[![Julia](https://img.shields.io/badge/Julia-1.12+-purple.svg)](https://julialang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<p align="center">
  <img src="docs/src/assets/logo.svg" alt="ERGMCount.jl icon" width="160">
</p>

ERGMs for Count-Valued Networks in Julia.

## Overview

ERGMCount.jl fits exponential-family random graph models to networks whose
edges carry **integer counts** — emails exchanged, co-occurrences, contacts —
rather than a presence/absence bit. A count ERGM is

```text
P(Y = y) ∝ h(y) · exp(θ' g(y)),     h(y) = ∏_ij h(y_ij)
```

where the **reference measure** `h` fixes the baseline distribution of a
dyad's count (Poisson, geometric, binomial, discrete uniform) and the
**terms** `g` are valued sufficient statistics (Krivitsky 2012). It is a port
of the R `ergm.count` package from the statnet collection, with the term
statistics validated against `ergm`/`ergm.count` by provenanced golden
fixtures (see [Validation against R](#validation-against-r)).

Two estimators are offered, and the default chooses between them as R does
(`method=:auto`): a **dyad-independent** formula is fit by maximum
pseudo-likelihood over an error-controlled count support, which is then the
exact MLE; a **dyad-dependent** formula is fit by **Monte-Carlo maximum
likelihood** on an exact Gibbs sampler, the estimator `ergm.count` reports.
`method=:mple` asks for the count MPLE explicitly (see
[Estimators](#estimators-mple-and-mcmle)).

## Installation

Requires Julia 1.12 or newer. The packages are not yet registered.

**Recommended: the ecosystem workspace.** It clones every package side by
side, develops them together in one environment, and adds the packages the
examples also use (CSV, DataFrames, Distributions, Graphs, StatsAPI,
StatsBase):

```bash
mkdir network-analysis && cd network-analysis
git clone https://github.com/statistical-network-analysis-with-Julia/statistical-network-analysis-with-Julia.github.io
julia statistical-network-analysis-with-Julia.github.io/tools/prepare_workspace.jl "$PWD" --clone
julia --project=.snippet-env
```

**Only this package, in your own environment.** Add its dependencies first,
in this order:

```julia
using Pkg
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/NetworkCore.jl")
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/ERGM.jl")
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/ERGMCount.jl")
```

The examples below load only `NetworkCore`, `ERGM`, `ERGMCount` and the `Random`
standard library.

## Features

- **Reference measures**: `PoissonReference`, `GeometricReference`,
  `BinomialReference`, `DiscUnifReference`, `DiscUnif2Reference`
- **Valued terms with R's labels**: `sum` (and `sum(pow=)`), `nonzero`,
  `greaterthan.k`, `atleast.k`, `atmost.k`, `smallerthan.k`,
  `equalto.v.pm.t`, `ininterval(a,b)`, `CMP`, the `mutual` forms (`min`,
  `nabsdiff`, `geometric`, `product`, `threshold`),
  `transitiveweights.min.max.min` / `cyclicalweights.min.max.min`, the
  valued covariate terms `nodematch`, `nodefactor`, `absdiff`, `nodecov`,
  `nodeocov`, `nodeicov` and `edgecov` (`form="sum"`/`"nonzero"`), and the
  node-strength terms `nodeOSum`/`nodeISum`/`nodeSum`
- **Estimation**: R's default rule — maximum pseudo-likelihood over an
  error-controlled count support for a dyad-independent model (the exact
  MLE), Monte-Carlo maximum likelihood (pinned against `ergm.count`) for a
  dyad-dependent one;
  parametric-bootstrap standard errors for the MPLE; R's boundary-statistic
  (`±Inf`) semantics; loud non-convergence, separation ("The MPLE does not
  exist!") and collinearity
- **Full StatsAPI surface**: `coef`, `stderror`, `vcov`, `confint`,
  `loglikelihood`, `nobs`, `dof`, `aic`, `bic`, `coeftable`, `coefnames`
- **Simulation and GOF**: a Gibbs sweep over each dyad's full conditional
  (0 bytes per dyad once warmed up) on an adaptive support that refuses a
  model whose chain does not settle, `gof(fit)` on the shared `NetworkCore.gof`
  generic

## Quick Start

```julia
using NetworkCore, ERGMCount, Random

# A directed count network: the counts live in the :weight edge attribute —
# ERGMCount's fixed name for what R's `response="w"` selects per call. Counts
# stored under another name are passed as `fit_ergm_count(net, terms; weight=:w)`;
# a network with edges but no `:weight` at all is refused (it would fit as 0/1),
# and so is one where only SOME edges carry a `:weight` (a bare edge is not a 1).
rng = Xoshiro(1)
net = network(20; directed=true)
for i in 1:20, j in 1:20
    if i != j && rand(rng) < 0.25
        add_edge!(net, i, j)
        set_edge_attribute!(net, :weight, i, j, rand(rng, 1:4))
    end
end

# Terms: overall volume, density, and reciprocity of intensity
terms = [SumTerm(), NonzeroTerm(), CountMutualTerm()]

# Fit under a Poisson reference (the default; R requires `reference=` to be
# spelled out). `mutual.min` is dyad-dependent, so the default method=:auto
# fits by Monte-Carlo maximum likelihood, as R's ergm.count does
mle = fit_ergm_count(net, terms; reference=PoissonReference(), rng=Xoshiro(2))
mle.method         # :mcmle
coef(mle)          # 3-vector, labelled as in R: sum, nonzero, mutual.min
coeftable(mle)     # estimate, SE (Fisher + Monte-Carlo error), z, p
confint(mle)       # Wald intervals
loglikelihood(mle) # the log-likelihood, by path sampling

# The count MPLE, asked for explicitly: fast, but for this dyad-dependent
# formula an approximation, so z and p are NaN (see "Estimators")
fit = fit_ergm_count(net, terms; method=:mple)
coeftable(fit)     # estimate and naive SE; z and p withheld
aic(fit)           # pseudo-likelihood AIC
fit.max_val        # the top of the enumerated support, chosen by doubling
```

`fit_ergm_count` is the standardized entry point (the ecosystem's
`fit_<model>` naming); `ergm_count` is the R-faithful alias of the same
function (`ergm_count === fit_ergm_count`).
`terms` may be a `Vector`, a `Tuple` or a single term.

## Reference measures

The reference measure `h(y)` is the baseline law of a dyad's count; the
`sum` coefficient shifts it (Krivitsky 2012, `ergm.count`'s conventions):

```julia
PoissonReference(1.0)      # h(y) = λ^y / y!, unbounded; with a `sum` coefficient θ
                           # each dyad is conditionally Poisson(λ·e^θ)
GeometricReference()       # the counting measure h(y) = 1, unbounded; a negative
                           # `sum` coefficient θ gives a geometric dyad law with
                           # success probability 1 − e^θ (no reference parameter)
BinomialReference(10)      # h(y) = C(10, y) on 0:10; the success probability is
                           # absorbed into the `sum` coefficient (logistic(θ))
DiscUnifReference(10)      # h(y) = 1 on 0:10
DiscUnif2Reference(1, 10)  # h(y) = 1 on 1:10 (zero is not a valid value)
```

Poisson and geometric are **unbounded** (`is_truncating(ref) == true`): the
estimator enumerates each dyad's conditional on `0:max_val` and must choose
that bound. The other three are bounded by construction and never truncate.

## Terms

Every term is a subtype of `AbstractERGMTerm` implementing `compute`,
`name` and `change_stat_count`; coefficients are labelled with the string R
prints for the same term.

```julia
using ERGM: compute, name   # the generic term interface shared with ERGM.jl

# Dyad-independent (the model is then an exact MLE)
SumTerm()                                   # sum:            Σ y_ij
SumTerm(pow=2)                              # sum2:           Σ y_ij², R's sum(pow=2)
NonzeroTerm()                               # nonzero:        Σ I(y_ij ≠ 0)
GreaterthannTerm(2)                         # greaterthan.2:  Σ I(y_ij > 2)
CountAtleastnTerm(3)                        # atleast.3:      Σ I(y_ij ≥ 3)
AtmostTerm(3)                               # atmost.3:       Σ I(y_ij ≤ 3), zero dyads included
CMPTerm()                                   # CMP:            Σ log(y_ij!) — ergm.count's
                                            # Conway–Maxwell–Poisson TERM (not a reference)
SmallerthanTerm(2)                          # smallerthan.2:  Σ I(y_ij < 2), zero dyads included
EqualToTerm(3)                              # equalto.3.pm.0: Σ I(y_ij = 3)
InIntervalTerm(1, 3)                        # ininterval(1,3): Σ I(1 < y_ij < 3)
InIntervalTerm(1, 3; open=(false, false))   # ininterval[1,3]: closed ends

# Dyad-independent covariate terms (vertex attributes :g, :x; matrix W);
# form=:nonzero counts the non-zero dyads instead of summing the counts
CountNodeMatchTerm(:g)                      # nodematch.sum.g:   Σ y_ij I(g_i = g_j)
CountNodeMatchTerm(:g; level="a")           # nodematch.sum.g.a: one level of nodematch(diff=TRUE)
CountNodeFactorTerm(:g, "b")                # nodefactor.sum.g.b: Σ y_ij (I(g_i = b) + I(g_j = b))
CountAbsDiffTerm(:x)                        # absdiff.sum.x:     Σ y_ij |x_i − x_j|
CountNodeCovTerm(:x)                        # nodecov.sum.x:     Σ y_ij (x_i + x_j); nodeocov/nodeicov
W = rand(5, 5)                              # any n × n dyadic covariate
CountEdgeCovTerm(W; name="dist")            # edgecov.sum.dist:  Σ y_ij W_ij

# Dyad-dependent, directed networks only
CountMutualTerm()                           # mutual.min:       Σ_{i<j} min(y_ij, y_ji)
CountMutualTerm(:nabsdiff)                  # mutual.nabsdiff:  −Σ |y_ij − y_ji|
CountMutualTerm(:geometric)                 # mutual.geom.mean: Σ sqrt(y_ij y_ji)
CountMutualTerm(:product)                   # mutual.product:   Σ y_ij y_ji
CountMutualTerm(:threshold; threshold=2)    # mutual.2: binary mutuality after thresholding at ≥ 2
NodeOSumTerm(); NodeISumTerm()              # nodeOSum / nodeISum: Σ_i (out-/in-strength)²

# Dyad-dependent, any network
TransitiveWeightsTerm()   # transitiveweights.min.max.min: Σ min(y_ij, max_k min(y_ik, y_kj))
CyclicalWeightsTerm()     # cyclicalweights.min.max.min:   Σ min(y_ij, max_k min(y_jk, y_ki))
NodeSumTerm()             # nodeSum: Σ_i (total strength)²

name(CountMutualTerm(:geometric))   # "mutual.geom.mean", as R prints it
```

`TransitiveTiesTerm()` and `CyclicalTiesTerm()` (a triple-wise minimum
summed over ordered triples) are **not** R terms and are not validated
against R; use `TransitiveWeightsTerm`/`CyclicalWeightsTerm` for parity with
`ergm`'s `transitiveweights`/`cyclicalweights`. A directed-only term on an
undirected network is refused at model construction with a hint
(`NodeSumTerm`, or drop the term) rather than fit as an all-zero column.

## Fitting

### The count support is error-controlled

For an unbounded reference the fit is a truncated exponential family on
`0:max_val`. With `max_val` unset, `fit_ergm_count` starts at
`max(10, 2 · largest count)`, refits at twice the bound, and stops at the
first doubling at which **all** of these hold: every estimate moved by at most
`support_tol` (default `1e-3`) standard errors, at most `support_tol`
expected dyads sit past the previous bound, and no dyad puts more than
`BOUNDARY_MASS_TOL` of its conditional mass on the new top value. The fit at
the larger bound is reported, and the result records what happened:

```julia
fit.support_control   # :converged (or :fixed, :bounded, :unconverged, :boundary_mode, :improper)
fit.support_stable    # false for :unconverged, :boundary_mode and :improper
fit.max_val           # the bound actually used
fit.support_delta     # max |Δθ|/SE at the last doubling (NaN for a fixed max_val, 0.0 bounded)
fit.omitted_tail      # expected dyads past the previous bound (same conventions)
fit.boundary_mass     # largest conditional mass any dyad puts on the top value

fixed = fit_ergm_count(net, terms; method=:mple, max_val=40)   # fix the bound instead
fixed.support_control                            # :fixed
```

A support still moving after `max_doublings` (default 8) is warned about and
reported `:unconverged` — the signature of a model that is not normalisable
on the unbounded support (a non-negative `sum` coefficient under a geometric
reference, say). The doubling also stops, `:unconverged` with
`support_delta`/`omitted_tail` `NaN`, at the first bound whose fit did not
converge or has a numerically singular pseudo-Hessian (a collinear design,
say): there are no settled estimates to compare, and the convergence and
conditioning warnings carry the diagnosis alone — never the "not
normalisable" sentence. On Zachary's karate club the default path doubles 14 → 28
and reproduces the exact MLE to ~1e-12 (the fixture pins it at 1e-6; the old
fixed default of 14 missed that tolerance).

**The check is on the dyad conditionals at the observed network.** A model
can pass it while its *joint* distribution is not normalisable on the
unbounded support — with a positive `mutual.product` or squared-strength
coefficient under a Poisson reference every conditional is a proper Poisson
law and the joint is improper (Krivitsky 2012, §3), however small the
coefficient — or has its mass far beyond the bound. Two checks look for it:

- **An analytic rule.** A statistic that grows faster than linearly in the
  counts (`mutual.product`, `nodeOSum`/`nodeISum`/`nodeSum`, `sum(pow=p)` with
  `p > 1`, and `CMP` beyond the reference's own `log y!`) outruns the Poisson
  reference's `−y log y` (and the geometric reference's 0) when its leading
  coefficient is positive. The rule evaluates the leading order of the
  log-weight with one dyad, a reciprocated pair, a star and every dyad at a
  count `Y → ∞` — so a larger negative squared-strength term can offset a
  positive one — and a positive leading coefficient makes the fit improper
  **whatever the data**: `fit.improper == true` (`support_control ==
  :improper` on the adaptive path).
- **A probe.** Started with every dyad at the bound where the adaptive sampler
  would give up (`2^max_doublings` times the fitted bound), conditional modes
  are iterated; if dyads stay there the fit is reported
  `fit.boundary_mode == true` (`support_control == :boundary_mode`). It covers
  what the rule cannot read: linear growth (a positive `transitiveweights`
  sum under a geometric reference) and user-defined terms.

Either way the fit is warned about, listed in `approximations`, and refused
by `se=:bootstrap`, `simulate_count_ergm`, `gof` and `method=:mcmle` unless
`max_val` fixes the truncated family in writing.

### Estimators: MPLE and MCMLE

The default, `method=:auto`, follows R (through ERGM.jl's
`ERGM.resolve_method`): the count MPLE for a dyad-independent formula and the
MCMLE otherwise. `method=:mple` maximizes the product of the dyad
conditionals. For a **dyad-independent** model that product is the
likelihood, the fit is the exact MLE and its Hessian standard errors are
exact. A keyword of the estimator `:auto` did not choose is refused in words
(`se=:bootstrap` on a dyad-dependent formula names `method=:mple`).

For a **dyad-dependent** model the count MPLE is this package's estimator,
not R's: `ergm.count` has no MPLE and fits by Monte-Carlo maximum likelihood,
which is why `:auto` does too. The two can differ materially. On `ergm.count`'s `zach` example
(`sum + nonzero + transitiveweights`, Poisson reference):

| | `sum` | `nonzero` | `transitiveweights` |
|:--|--:|--:|--:|
| count MPLE (`method=:mple`) | 0.797 | −4.924 | 0.385 |
| `ergm.count` MCMLE (mean of 9 seeds) | 0.597 | −4.779 | 0.562 |
| `method=:mcmle` (mean of 6 seeds) | 0.597 | −4.778 | 0.560 |
| standard error (R) | 0.121 | 0.237 | 0.101 |

The MPLE is 1.7 and 1.8 standard errors from the MLE on two coefficients.
Simulation at the MLE (300 networks drawn at `ergm.count`'s long-chain
estimate, the MPLE refit on each) shows it is nearly unbiased but less
efficient (its sampling standard deviation for `transitiveweights` is 0.121
against the MLE's standard error of 0.101), and that its naive
pseudo-Hessian standard errors are too small (0.058 on average for that
coefficient): 95 % Wald intervals built from them covered 0.91, 0.98 and
0.70, against 0.96, 0.97 and 0.97 for the parametric bootstrap (100 of the
networks, `n_boot=60`).
So **a dyad-dependent MPLE fit (`method=:mple`) reports estimates and naive
standard errors but no z, p or interval**: the table shows `NaN` with a note,
`confint` throws an `ArgumentError`, and `approximations(fit)` records it.
There are three ways to inference:

```julia
# 1. Maximum likelihood — what R reports, and the default here
mle = fit_ergm_count(net, terms; rng=Xoshiro(2))
mle.converged, mle.mcmc.convergence.iterations
stderror(mle)                       # inverse Fisher information + Monte-Carlo error
mle.mcmc.mc_std_errors              # the Monte-Carlo part alone

# 2. The MPLE with a parametric bootstrap (calibrated SEs, same point estimates)
robust = fit_ergm_count(net, terms; method=:mple, se=:bootstrap, n_boot=50,
                        rng=Xoshiro(2))
stderror(robust)                    # the bootstrap standard errors
NetworkCore.se_method(robust)          # :bootstrap
coef(robust) == coef(fit)           # true — only the covariance changed
robust.boot_replicates              # the 50 × 3 refits (excluded ones as NaN rows)

# 3. The naive Wald table, asked for in writing (printed with its caveat)
naive = fit_ergm_count(net, terms; method=:mple, se=:hessian)
naive.z_values
```

`method=:mcmle` runs ERGM.jl's MCMLE iteration (`ERGM.Extension.mcmle_solve`, the one
`ERGM.mcmle` runs) on the Gibbs sampler: from the MPLE it draws `n_samples`
(default 1024) networks per iteration, thinned until their effective sample
size reaches half that number, takes Hummel-stepped Monte-Carlo Newton steps,
and stops on R `ergm` 4's confidence rule (`termination=:hotelling` selects
the t-ratio and Hotelling tests instead). The standard errors come from one
more sample at the estimate. `loglikelihood` is estimated by path sampling from the dyad-independent
part of the model (`bridge_rungs=0` skips it) and includes the reference
measure; R's `logLik` is the same number minus the log-likelihood at `θ = 0`.
On a dyad-independent model `method=:mcmle` returns the exact fit with no
Monte Carlo. The bootstrap's replicates on which the count MPLE does not
exist or does not converge are excluded from the covariance with one
aggregate warning, which `show` and `approximations` repeat: the standard
errors are then conditional on a finite refit, and since the excluded
replicates are the extreme ones they are biased downward.

### Boundary statistics

A statistic at the boundary of its attainable range — `mutual.min` on a
network with no reciprocated pair, `nonzero` on a complete network — has no
finite MPLE. As in R `ergm`, its coefficient is fixed at `±Inf` with standard
error 0 (warned with R's sentence), the remaining coefficients are estimated
on the restricted supports (the exact limit), `dof`/`aic`/`bic` count only
the finite ones, and the fit is never `is_exact`. `se=:bootstrap`,
`simulate_count_ergm` and `gof` refuse such a fit — there is no finite model
to simulate from.

```julia
acyclic = network(5; directed=true)            # 1→2→3→4→5, nothing reciprocated
for (i, j) in ((1, 2), (2, 3), (3, 4), (4, 5))
    add_edge!(acyclic, i, j); set_edge_attribute!(acyclic, :weight, i, j, 2)
end
b = fit_ergm_count(acyclic, [SumTerm(), CountMutualTerm()]; reference=BinomialReference(3),
                   method=:mple)
coef(b)[2]         # -Inf
dof(b)             # 1
is_exact(b)        # false
```

### Non-convergence, separation and collinearity are loud

A Newton iteration that exhausts `maxiter` or cannot move is warned about
(naming `maxiter`, `tol` and the pseudo-score norm), recorded as
`fit.converged == false` with `fit.iterations`/`fit.gradient_norm`, listed in
`NetworkCore.approximations(fit)`, printed by `show`, and never `is_exact`.

A design on which the pseudo-likelihood has no finite maximum along a
**combination** of the statistics — separation, which no single-column
boundary test sees; `sum + nonzero` on a network whose every count is 0 or 1
is the textbook case (`sum − nonzero` is at its minimum on every dyad, so
θ_sum → −∞, θ_nonzero → +∞) — is decided exactly from the data by the
ecosystem's shared verdict (`NetworkCore.clogit_separation`: the count
pseudo-likelihood is a conditional logit over each dyad class's support, and
the verdict is the linear programme of R's `mple.existence`). Such a fit is
returned with `converged = false`, `fit.separated = true` and the separating
terms in `fit.separated_terms`, warned about (naming them, with R's "The
MPLE does not exist!"), and its z values, p-values and confidence intervals
are `NaN`; `se=:bootstrap` and `method=:mcmle` are refused. Two statistics that are collinear on
the data — `greaterthan(2)` and `atleast(3)` on integer counts, `sum` and
`nonzero` on a 0/1 network — make the pseudo-Hessian numerically singular:
`fit.hessian_cond` (its condition number) is recorded, and above `1e8` the
fit warns naming the statistics that load on the flat direction
(`fit.collinear`), with the caveat in `approximations` and under
`Converged:` in `show`.

```julia
binary = network(6; directed=false)                      # every count is 1
for (i, j) in ((1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (1, 3), (2, 5))
    add_edge!(binary, i, j); set_edge_attribute!(binary, :weight, i, j, 1)
end
sep = fit_ergm_count(binary, [SumTerm(), NonzeroTerm()]; warn=false)
sep.separated, sep.converged        # (true, false) — no finite MPLE
sep.separated_terms                 # ["sum", "nonzero"]
confint(sep)                        # NaN: no inference without an MPLE
```

`warn=false` silences every fit diagnostic (boundary, truncation,
non-convergence, separation, conditioning) without changing what the result
records.

### StatsAPI surface

```julia
coef(fit); stderror(fit); vcov(fit)
confint(mle; level=0.95)     # normal-theory limits, one row per coefficient
                             # (refused for the dyad-dependent MPLE `fit`)
loglikelihood(fit)           # the maximised pseudo-log-likelihood (MPLE fit)
nobs(fit), dof(fit)          # dyads, finite coefficients
aic(fit), bic(fit)           # pseudo-likelihood criteria: compare only across
                             # models on the same network, reference and support
coeftable(fit)["sum"]        # one labelled row of the table
coefnames(fit)               # ["sum", "nonzero", "mutual.min"]: R's labels
NetworkCore.fit_metadata(fit)   # estimand, objective, is_exact, se_method, ...
```

## Simulation and GOF

`simulate_count_ergm` is a **Gibbs sweep**: every dyad is redrawn once per
sweep from its full conditional `P(y_ij = y | rest) ∝ h(y)·exp(θ'Δg(y))`,
so structural terms (mutuality, transitivity, strength) shape the draws.
`burnin` and `interval` are counted in **sweeps** and default through the
dyad-scaled rule shared with ERGM.jl (`ERGM.Extension.mcmc_defaults` converted from
single-dyad moves to sweeps): 20 sweeps of burn-in and 1 sweep between
retained draws on any network with more than ten nodes. `gof` and the
parametric bootstrap use the same sampler and the same defaults.

For an unbounded reference the sampler's support is **adaptive**: it starts
at the fit's bound (or `max(10, 2·largest seed count)` for an explicit
specification) and doubles whenever a dyad's conditional puts more than
1e-10 of its mass on the top value, so the truncation never shapes a draw. A
chain that reaches `2^max_doublings` times its starting bound is refused with
an `ArgumentError` — the signature of a model that is not normalisable.
`max_val=k` asks for the family truncated at `0:k` instead; if the chain then
puts more than `BOUNDARY_MASS_TOL` of a conditional on `k`, the draws come
with a warning giving the share that landed on the bound, and `gof` and the
bootstrap refuse.

```julia
sims = simulate_count_ergm(mle; n_sim=100, rng=Xoshiro(3))   # default 20 / 1 sweeps
sims = simulate_count_ergm(mle; n_sim=100, burnin=50, interval=5, rng=Xoshiro(3))
compute(SumTerm(), sims[1])

# Goodness of fit: the model statistics and the dyad count-value distribution
g = gof(mle; n_sim=50, rng=Xoshiro(4))
g.statistics[1].labels     # ["sum", "nonzero", "mutual.min"]
```

Under `DiscUnif2Reference(a, b)` with a negative `a` (R's `DiscUnif(a, b)`
admits one) the draws can be negative: an edge is stored for every non-zero
count, negative ones included, `nonzero` counts the dyads with `y ≠ 0`, and
`gof`'s count-value panel runs over every value the networks take. As in
`ergm`, `TransitiveWeightsTerm`, `CyclicalWeightsTerm` and
`CountMutualTerm(:geometric)` "may not be used with networks with negative
dyad weights": `compute`, `CountERGMModel` and `simulate_count_ergm` refuse
them with that sentence on negative data (or over a support reaching below
0) instead of returning a number R never produces.

A warmed-up sweep costs ≈ 0.6 µs and 0 bytes per dyad (pinned by the test
suite and `benchmark/regression_tests.jl`). Draws are reproducible from
`rng`; the stream changed in 0.2.0, so seeded draws differ from 0.1.x.

## Missing dyads

A network with masked (unobserved) dyads is **refused** by `fit_ergm_count`,
`count_mple`, `simulate_count_ergm` and `gof`: the count MPLE would enumerate
every unobserved dyad as an observed row, and the Gibbs sweep would condition
on a face value nobody observed. There is no `missing=` keyword and no
`:face` policy (`NetworkCore.missing_policies(fit_ergm_count) == (:error,)`);
the error names `clear_missing_dyads!` for the case where the dyads really
are observed zeros.

## Validation against R

Three provenanced golden fixtures, generated by checked-in R scripts and
loaded with `NetworkCore.load_golden` (which refuses a fixture without a
`[provenance]` block):

- **`test/fixtures/count_mcmle.toml`** (`test/fixtures/r/count_mcmle.R`,
  ergm.count 4.1.3): the **MCMLE** of two dyad-dependent models —
  `zach ~ sum + nonzero + transitiveweights` and a seeded 14-actor directed
  count network `~ sum + nonzero + mutual(min)`, Poisson reference — each
  refit by `ergm.count` under nine seeds, plus one long-chain fit. Both sides
  are Monte-Carlo estimates of the same MLE, so `method=:mcmle` is held to
  R's nine-seed mean within four of R's own seed-to-seed standard deviations
  (floored at a tenth of a standard error), coefficients and standard errors
  alike, and its path-sampled log-likelihood to R's bridge estimate.
- **`test/fixtures/zach_poisson.toml`** (`test/fixtures/r/zach_poisson.R`,
  ergm.count 4.1.3): `zach ~ sum + nonzero`, Poisson reference, on
  `ergm.count`'s Zachary karate club. The model is **dyad-independent on
  purpose**: each dyad is then an independent draw from a two-parameter law
  with an *exact* MLE, and the count MPLE's dyad-conditional enumeration *is*
  that likelihood — so ERGMCount.jl is held to the exact MLE at **1e-6**
  (reproduced to ~1e-12), on both the pinned `max_val=30` and the default
  error-controlled path. The golden value is solved analytically (R's
  `optim` was only good to 7e-6); `ergm.count`'s own MCMLE sits 0.0103 from
  it and is frozen as a cross-check, not as the standard. The fixture also
  freezes `P(y > 30) = 6.9e-23` under the fitted law, so the truncation is
  measured, not assumed.
- **`test/fixtures/count_terms.toml`** (`test/fixtures/r/count_terms.R`,
  ergm 4.12.0): term parity. The summary statistics *and R's coefficient
  names* for `sum`, `nonzero`, `greaterthan`, `atleast`, `smallerthan`,
  `equalto`, `ininterval` (all four bracket combinations),
  `transitiveweights`/`cyclicalweights` (`min`, `max`, `min`) and the
  `mutual` forms `min`/`nabsdiff`/`geometric`/`product`, plus
  `greaterthan(-1)`/`atleast(0)` (thresholds at or below 0 count the zero
  dyads) and `sum(pow=2)`, `sum(pow=1/2)`, `atmost` and `CMP`, on `zach`
  (undirected) and on a seeded 8-actor directed count network frozen in the
  fixture, at 1e-9; a 3-actor network with a
  **0-valued tie** (R's `nonzero` is 1 there and the threshold terms count
  the tie with the empty dyads); and a seeded 6-actor directed network with
  values in **-2:2**, on which the terms R accepts are pinned and R's
  refusal of `transitiveweights`/`cyclicalweights` ("may not be used with
  networks with negative dyad weights") and its `NaN` for
  `mutual(form="geometric")` are frozen as the behaviour ERGMCount refuses.
  `mutual(form="threshold")` cannot be pinned — ergm 4.12.0's own
  `summary()` fails at C initialisation, and the fixture records the error
  — so `CountMutualTerm(:threshold)` is tested by hand value.

Every change statistic is also brute-force tested against `compute` on
edited networks, and every support profile against the per-value change
statistic. The MCMLE is additionally checked against exact enumeration (a
3-actor network, 729 states: the maximizer, its standard errors and the
log-likelihood).

## Not implemented

What a user of R `ergm.count` will not find here yet, and what happens
instead. Each is an explicit error or an absent name, never a silently wrong
number:

- **A missing-data MLE for masked count networks, contrastive divergence
  and stochastic approximation.** `method` is `:auto`, `:mple` or `:mcmle`;
  anything else is an `ArgumentError`, and a masked network is refused by
  both estimators.
- **`ergm.count`'s MCMC proposals and `control.ergm` tuning.** The sampler is
  a systematic-scan Gibbs sweep over each dyad's enumerated conditional, not
  R's Poisson/zero-inflated Metropolis proposals; the two were verified to
  draw from the same distribution, but seeds and tuning constants do not
  carry over.
- **A default reference where R requires one.** R's `ergm` stops without
  `reference=` on a valued model; `fit_ergm_count` defaults to
  `PoissonReference()`.
- **The count MPLE is not an `ergm.count` estimator.** The default
  (`method=:auto`) uses it only where it is the exact MLE (a dyad-independent
  formula) and fits the MCMLE otherwise, as R does; `method=:mple` on a
  dyad-dependent formula is this package's own approximation (see
  [Estimators](#estimators-mple-and-mcmle)).
- **Joint normalisability is decided analytically only along five
  configurations.** A positive leading coefficient on `mutual.product`, a
  squared-strength term, `sum(pow>1)` or `CMP`, or — under the geometric
  reference — a positive linear growth (`sum` ≥ 0, or 2·`sum` + `mutual.min`
  > 0 on a reciprocated pair) is refused whatever the data. Growth the rule
  cannot read (a user-defined term, a covariate term's `form=:sum`, a
  configuration other than a dyad, a reciprocated pair, a star or every
  dyad) is caught only by the probe at `2^max_doublings` times the fitted
  bound — whose mode on the bound must outweigh the observed network — and
  by the sampler's adaptive support, so a model improper only beyond that
  bound, and metastable below it, is not detected.
- **An MCMLE with a statistic at its bound has no log-likelihood.** The
  coefficient is fixed at ∓Inf and the rest estimated, as R does; the path
  sampler cannot hold the statistic at its bound, so `loglikelihood(fit)`,
  AIC and BIC are `NaN` (R reports a log-likelihood there).
- **The `StdNormal` and continuous `Unif` references**: no reference type to
  call; the estimator enumerates integer supports only. (`CMP` is a *term* in
  `ergm.count`, and is one here: `CMPTerm()`.)
- **The `nodecovar`/`nodeocovar`/`nodeicovar`/`nodesqrtcovar` family, the
  valued `nodemix`, `absdiff(pow≠1)` and `nodematch(keep=, levels=)`**: no
  term to call. The valued `nodematch`, `nodefactor`, `absdiff`, `nodecov`,
  `nodeocov`, `nodeicov` and `edgecov` (`form="sum"`/`"nonzero"`) are
  `CountNodeMatchTerm` … `CountEdgeCovTerm`, one statistic per term: R's
  `nodefactor` and `nodematch(diff=TRUE)` are one term per level.
- **`transitiveweights`/`cyclicalweights` with a non-default triple**
  (`twopath="geomean"`, `combine="sum"`, `affect="geomean"`):
  `TransitiveWeightsTerm`/`CyclicalWeightsTerm` implement the default
  `(min, max, min)` only and take no arguments.
- **`TransitiveTiesTerm`/`CyclicalTiesTerm` are not R terms** (a triple-wise
  minimum over ordered triples, kept from 0.1). They are not validated
  against R; their names say `.count`, not `transitiveweights`.
- **`sum(pow=)` with a non-integer power, and `CMP`, on negative counts**:
  refused (R returns `NaN` and `Inf`).
- **Curved terms, `constraints=`, offsets, and two-mode networks**: none;
  every term is a fixed-parameter statistic, and a two-mode (bipartite)
  network — `network(n; bipartite=k)` or a `BipartiteNetwork` — is
  **refused** with an `ArgumentError` by `fit_ergm_count`, `CountERGMModel`,
  `simulate_count_ergm` and `gof`: enumerating the one-mode dyads would count
  the impossible within-mode dyads as observed zeros, so nothing is fit.

## Documentation

For more detailed documentation, see:

- [Documentation](https://statistical-network-analysis-with-Julia.github.io/ERGMCount.jl/dev/)

## References

1. Krivitsky, P.N. (2012). Exponential-family random graph models for valued networks. *Electronic Journal of Statistics*, 6, 1100-1128.

2. Desmarais, B.A., Cranmer, S.J. (2012). Statistical mechanics of networks: Estimation and uncertainty. *Physica A*, 391(4), 1865-1876.

3. Hunter, D.R., Handcock, M.S., Butts, C.T., Goodreau, S.M., Morris, M. (2008). ergm: A package to fit, simulate and diagnose exponential-family models for networks. *Journal of Statistical Software*, 24(3), 1-29.

## Citation

If you use ERGMCount.jl in your work, please cite it using the entry in
[`CITATION.bib`](CITATION.bib):

```biblatex
@misc{SNWJERGMCountJL,
  author = {Santoni, Simone},
  title = {ERGMCount.jl: Exponential Random Graph Models for Count-Valued Networks in Julia},
  year = {2026},
  url = {https://github.com/statistical-network-analysis-with-Julia/ERGMCount.jl},
  note = {Homepage: https://statistical-network-analysis-with-Julia.github.io/ERGMCount.jl; GitHub: https://github.com/statistical-network-analysis-with-Julia}
}
```

ERGMCount.jl implements methods developed for, and is validated against, the
R package `ergm.count`. **Please also cite the R package and the methods
paper** — Krivitsky (2012) and the `ergm.count` package (Krivitsky 2025;
`citation("ergm.count")` in R gives the current entry); the per-package list
is at
<https://statistical-network-analysis-with-julia.github.io/citing/>.

## License

MIT License - see [LICENSE](LICENSE) for details.
