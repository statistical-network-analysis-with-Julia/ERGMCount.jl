# ERGMCount.jl


[![Network Analysis](https://img.shields.io/badge/Network-Analysis-orange.svg)](https://github.com/statistical-network-analysis-with-Julia/ERGMCount.jl)
[![Build Status](https://github.com/statistical-network-analysis-with-Julia/ERGMCount.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/statistical-network-analysis-with-Julia/ERGMCount.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://statistical-network-analysis-with-Julia.github.io/ERGMCount.jl/stable/)
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

Estimation is by **maximum pseudo-likelihood** over an error-controlled
count support, which is the exact MLE for a dyad-independent model;
`ergm.count`'s Monte-Carlo MLE is not implemented (see
[Not implemented](#not-implemented)).

## Installation

Requires Julia 1.12+. Until the packages are registered, ERGMCount.jl and
its dependencies [Networks.jl](https://github.com/statistical-network-analysis-with-Julia/Networks.jl)
and [ERGM.jl](https://github.com/statistical-network-analysis-with-Julia/ERGM.jl)
are added from GitHub, **in this order** (each depends on the one before):

```julia
using Pkg
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/Networks.jl")
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/ERGM.jl")
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/ERGMCount.jl")
```

For development, clone the ecosystem repositories side by side and start
Julia with the root workspace project (`julia --project=.` in the clone
root): the `[sources]` path dependencies wire the packages together with no
ordered installs.

## Features

- **Reference measures**: `PoissonReference`, `GeometricReference`,
  `BinomialReference`, `DiscUnifReference`, `DiscUnif2Reference`
- **Valued terms with R's labels**: `sum`, `nonzero`, `greaterthan.k`,
  `atleast.k`, `smallerthan.k`, `equalto.v.pm.t`, `ininterval(a,b)`, the
  `mutual` forms (`min`, `nabsdiff`, `geometric`, `product`, `threshold`),
  `transitiveweights.min.max.min` / `cyclicalweights.min.max.min`, and the
  node-strength terms `nodeOSum`/`nodeISum`/`nodeSum`
- **Estimation**: maximum pseudo-likelihood over an error-controlled count
  support (exact MLE for dyad-independent models); inverse-pseudo-Hessian or
  parametric-bootstrap standard errors; R's boundary-statistic (`±Inf`)
  semantics; loud non-convergence, separation ("The MPLE does not exist!")
  and collinearity
- **Full StatsAPI surface**: `coef`, `stderror`, `vcov`, `confint`,
  `loglikelihood`, `nobs`, `dof`, `aic`, `bic`, `coeftable`
- **Simulation and GOF**: a Gibbs sweep over each dyad's full conditional
  (0 bytes per dyad once warmed up), `gof(fit)` on the shared
  `Networks.gof` generic

## Quick Start

```julia
using Networks, ERGMCount, Random

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

# Fit under a Poisson reference (the default); the count support is chosen
# by error control, see below
fit = fit_ergm_count(net, terms; reference=PoissonReference())

coef(fit)          # 3-vector, labelled as in R: sum, nonzero, mutual.min
coeftable(fit)     # the table `show(fit)` prints: estimate, SE, z, p
aic(fit)           # pseudo-likelihood AIC
fit.max_val        # the top of the enumerated support, chosen by doubling
```

`fit_ergm_count` is the standardized entry point (the ecosystem's
`fit_<model>` naming); `ergm_count` is the R-faithful alias of the same
function and `fit_count_ergm` a legacy alias (`ergm_count === fit_ergm_count`).
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
NonzeroTerm()                               # nonzero:        Σ I(y_ij ≠ 0)
GreaterthannTerm(2)                         # greaterthan.2:  Σ I(y_ij > 2)
CountAtleastnTerm(3)                        # atleast.3:      Σ I(y_ij ≥ 3)
SmallerthanTerm(2)                          # smallerthan.2:  Σ I(y_ij < 2), zero dyads included
EqualToTerm(3)                              # equalto.3.pm.0: Σ I(y_ij = 3)
InIntervalTerm(1, 3)                        # ininterval(1,3): Σ I(1 < y_ij < 3)
InIntervalTerm(1, 3; open=(false, false))   # ininterval[1,3]: closed ends

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
fit.support_control   # :converged (or :fixed, :bounded, :unconverged)
fit.support_stable    # false only for :unconverged
fit.max_val           # the bound actually used
fit.support_delta     # max |Δθ|/SE at the last doubling (NaN for a fixed max_val, 0.0 bounded)
fit.omitted_tail      # expected dyads past the previous bound (same conventions)
fit.boundary_mass     # largest conditional mass any dyad puts on the top value

fixed = fit_ergm_count(net, terms; max_val=40)   # fix the bound instead
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

### Standard errors

`se=:hessian` (default) is the inverse negative pseudo-Hessian — correct for
a dyad-independent model, **anticonservative** for any dyad-dependent one,
which `show(fit)` says out loud. `se=:bootstrap` Gibbs-simulates `n_boot`
networks at the estimates, refits each, and reports the empirical
covariance on the shared `Networks.bootstrap_cov` loop; the point estimates
are unchanged.

```julia
robust = fit_ergm_count(net, terms; se=:bootstrap, n_boot=50, rng=Xoshiro(2))
stderror(robust)                    # the bootstrap standard errors
Networks.se_method(robust)          # :bootstrap
coef(robust) == coef(fit)           # true — only the covariance changed
robust.boot_replicates              # the 50 × 3 refits (excluded ones as NaN rows)
```

Replicates on which the count MPLE does not exist or does not converge are
excluded from the covariance with one aggregate warning.

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
b = fit_ergm_count(acyclic, [SumTerm(), CountMutualTerm()]; reference=BinomialReference(3))
coef(b)[2]         # -Inf
dof(b)             # 1
is_exact(b)        # false
```

### Non-convergence, separation and collinearity are loud

A Newton iteration that exhausts `maxiter` or cannot move is warned about
(naming `maxiter`, `tol` and the pseudo-score norm), recorded as
`fit.converged == false` with `fit.iterations`/`fit.gradient_norm`, listed in
`Networks.approximations(fit)`, printed by `show`, and never `is_exact`.

A design on which the pseudo-likelihood has no finite maximum along a
**combination** of the statistics — quasi-complete separation, which no
single-column boundary test sees; `sum + nonzero` on a network whose every
count is 0 or 1 is the textbook case (`sum − nonzero` is at its minimum on
every dyad, so θ_sum → −∞, θ_nonzero → +∞) — is detected as `ERGM.mple`
detects it (an unobserved support value pushed > 18 nats down by the tilt,
*and* a Newton step still O(1)) and returned with `converged = false`,
`fit.separated = true`, R's sentence "The MPLE does not exist!" in the
warning, `approximations` and `show`. Two statistics that are collinear on
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
sep.collinear                       # ["sum", "nonzero"]
```

`warn=false` silences every fit diagnostic (boundary, truncation,
non-convergence, separation, conditioning) without changing what the result
records.

### StatsAPI surface

```julia
coef(fit); stderror(fit); vcov(fit)
confint(fit; level=0.95)     # normal-theory limits, one row per coefficient
loglikelihood(fit)           # the maximised pseudo-log-likelihood
nobs(fit), dof(fit)          # dyads, finite coefficients
aic(fit), bic(fit)           # pseudo-likelihood criteria: compare only across
                             # models on the same network, reference and support
coeftable(fit)["sum"]        # one labelled row of the table
Networks.fit_metadata(fit)   # estimand, objective, is_exact, se_method, ...
```

## Simulation and GOF

`simulate_count_ergm` is a **Gibbs sweep**: every dyad is redrawn once per
sweep from its full conditional `P(y_ij = y | rest) ∝ h(y)·exp(θ'Δg(y))`,
so structural terms (mutuality, transitivity, strength) shape the draws.
`burnin` and `interval` are counted in **sweeps** and default through the
dyad-scaled rule shared with ERGM.jl (`ERGM._mcmc_defaults` converted from
single-dyad moves to sweeps): 20 sweeps of burn-in and 1 sweep between
retained draws on any network with more than ten nodes. `gof` and the
parametric bootstrap use the same sampler and the same defaults.

```julia
sims = simulate_count_ergm(fit; n_sim=100, rng=Xoshiro(3))   # default 20 / 1 sweeps
sims = simulate_count_ergm(fit; n_sim=100, burnin=50, interval=5, rng=Xoshiro(3))
compute(SumTerm(), sims[1])

# Goodness of fit: the model statistics and the dyad count-value distribution
g = gof(fit; n_sim=50, rng=Xoshiro(4))
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
`:face` policy (`Networks.missing_policies(fit_ergm_count) == (:error,)`);
the error names `clear_missing_dyads!` for the case where the dyads really
are observed zeros.

## Validation against R

Two provenanced golden fixtures, generated by checked-in R scripts and loaded
with `Networks.load_golden` (which refuses a fixture without a
`[provenance]` block):

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
  dyads), on `zach` (undirected, 13 rows) and on a seeded 8-actor directed
  count network frozen in the fixture, at 1e-9; a 3-actor network with a
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
statistic.

## Not implemented

What a user of R `ergm.count` will not find here yet, and what happens
instead. Each is an explicit error or an absent name, never a silently wrong
number:

- **Monte-Carlo maximum likelihood (`method=:mcmle`)**, `ergm.count`'s own
  estimator, and with it a missing-data MLE for masked count networks.
  `fit_ergm_count(net, terms; method=:mcmle)` throws an `ArgumentError`
  explaining that ERGMCount.jl fits by maximum pseudo-likelihood only — the
  exact MLE for a dyad-independent model, a pseudo-likelihood estimate with
  `se=:bootstrap` as the honest covariance for a dyad-dependent one:

  ```julia
  msg = try
      fit_ergm_count(net, terms; method=:mcmle)
  catch err
      sprint(showerror, err)
  end
  occursin("maximum pseudo-likelihood only", msg)   # true
  ```
- **The `CMP`, `StdNormal` and continuous `Unif` references**: no reference
  type to call; the estimator enumerates integer supports only.
- **The `nodecovar`/`nodeocovar`/`nodeicovar`/`nodesqrtcovar` family** and
  the valued forms of `nodematch`/`nodefactor`/`absdiff`/`edgecov`
  (`form=`): no term to call.
- **`transitiveweights`/`cyclicalweights` with a non-default triple**
  (`twopath="geomean"`, `combine="sum"`, `affect="geomean"`):
  `TransitiveWeightsTerm`/`CyclicalWeightsTerm` implement the default
  `(min, max, min)` only and take no arguments.
- **`TransitiveTiesTerm`/`CyclicalTiesTerm` are not R terms** (a triple-wise
  minimum over ordered triples, kept from 0.1). They are not validated
  against R; their names say `.count`, not `transitiveweights`.
- **Curved terms, `constraints=`, offsets, and two-mode networks**: none;
  every term is a fixed-parameter statistic, and a two-mode (bipartite)
  network — `network(n; bipartite=k)` or a `BipartiteNetwork` — is
  **refused** with an `ArgumentError` by `fit_ergm_count`, `CountERGMModel`,
  `simulate_count_ergm` and `gof`: enumerating the one-mode dyads would count
  the impossible within-mode dyads as observed zeros, so nothing is fit.

## Documentation

For more detailed documentation, see:

- [Stable Documentation](https://statistical-network-analysis-with-Julia.github.io/ERGMCount.jl/stable/)
- [Development Documentation](https://statistical-network-analysis-with-Julia.github.io/ERGMCount.jl/dev/)

## References

1. Krivitsky, P.N. (2012). Exponential-family random graph models for valued networks. *Electronic Journal of Statistics*, 6, 1100-1128.

2. Desmarais, B.A., Cranmer, S.J. (2012). Statistical mechanics of networks: Estimation and uncertainty. *Physica A*, 391(4), 1865-1876.

3. Hunter, D.R., Handcock, M.S., Butts, C.T., Goodreau, S.M., Morris, M. (2008). ergm: A package to fit, simulate and diagnose exponential-family models for networks. *Journal of Statistical Software*, 24(3), 1-29.

## Citation

If you use ERGMCount.jl in your work, please cite it using the entry in
[`CITATION.bib`](CITATION.bib):

```biblatex
@misc{SNWJERGMCountJL,
  author = {{Statistical Network Analysis with Julia}},
  title = {ERGMCount.jl: Exponential Random Graph Models for Count-Valued Networks in Julia},
  year = {2026},
  url = {https://github.com/statistical-network-analysis-with-Julia/ERGMCount.jl},
  note = {Homepage: https://statistical-network-analysis-with-Julia.github.io/ERGMCount.jl; GitHub: https://github.com/statistical-network-analysis-with-Julia}
}
```

## License

MIT License - see [LICENSE](LICENSE) for details.
