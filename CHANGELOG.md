# Changelog

All notable changes to ERGMCount.jl are documented in this file. The format
is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
package adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - Unreleased

First public release of ERGMCount.jl, a Julia port of R's `ergm.count`
(statnet): exponential-family random graph models for count-valued networks
(Krivitsky 2012), fitted by maximum pseudo-likelihood or Monte-Carlo maximum
likelihood, with an exact Gibbs sampler, goodness of fit, and term values,
labels and estimates checked against R `ergm` 4.12 / `ergm.count` 4.1.3 by
provenanced fixtures. The changes below are relative to 0.1.0, a development
version, never released.

**Dependency renamed:** the foundation package is now `NetworkCore` (developed as `Networks`); write `using NetworkCore` where code said `using Networks`. Types and functions keep their names.

### Highlights

- **Monte-Carlo maximum likelihood** (`method=:mcmle`, `count_mcmle`), the
  estimator R's `ergm.count` uses, on the exact Gibbs sampler. It matches
  `ergm.count` on two dyad-dependent models within R's own seed-to-seed
  spread, and exact enumeration on a small network.
- **R's default estimator.** `method=:auto` (the default, ERGM.jl's
  `resolve_method`) fits a dyad-independent formula by the count MPLE, which
  is then the exact MLE, and a dyad-dependent one by MCMLE, as `ergm.count`
  does.
- **A proper count MPLE**: each dyad's full conditional is enumerated over
  the count support, so the reference measure enters the estimator. For a
  dyad-independent model it is the exact MLE (reproduced to 1e-12 on
  `ergm.count`'s `zach`).
- **Honest inference under dependence.** A dyad-dependent MPLE fit no longer
  prints naive Wald z and p-values by default; `method=:mcmle` and
  `se=:bootstrap` give calibrated ones.
- **The count support is never silently binding.** The estimator chooses the
  truncation of an unbounded reference by error control, the sampler widens
  it on demand, and a model that is not normalisable — decided analytically
  for super-linear terms, probed at the sampler's cap otherwise — is flagged
  and refused downstream.
- **R's terms and labels**: `sum` (and `sum(pow=)`), `nonzero`, `greaterthan`,
  `atleast`, `atmost`, `smallerthan`, `equalto`, `ininterval`, `CMP`, the
  valued `mutual` forms, `transitiveweights`, `cyclicalweights`.
- **Loud failures instead of silent mis-fits**: boundary statistics,
  separation, collinearity, non-convergence, two-mode networks, masked dyads,
  unweighted or non-integer counts and binary ERGM terms are refused or
  reported as R does, in words.

### Breaking

- **The default estimator follows R.** `fit_ergm_count`/`ergm_count` default
  to `method=:auto`: the count MPLE for a dyad-independent formula, the MCMLE
  for a dyad-dependent one (it was the MPLE for every formula). Pass
  `method=:mple` for the count MPLE. A keyword of the estimator `:auto` did
  not choose (e.g. `se=:bootstrap` on a dyad-dependent formula) is an
  `ArgumentError` naming the method that takes it. `show` prints a `Method:`
  line saying which estimator ran and why.
- **`fit_count_ergm` is removed** (a legacy alias of `fit_ergm_count`, never
  released); the docs list it under "Renamed and removed names".
- **MPLE inference under dyadic dependence is withheld.** With the default
  `se=nothing`, a `method=:mple` fit of a dyad-dependent model reports estimates and
  naive pseudo-Hessian standard errors, but `z_values`/`p_values` are `NaN`,
  `show` explains why, `confint` throws an `ArgumentError`, and
  `fit.inference_withheld` is `true` (simulated at the MLE of the `zach`
  model, 95 % Wald intervals from those standard errors covered 0.91, 0.98
  and 0.70; the bootstrap's 0.96, 0.97 and 0.97). `se=:hessian` passed explicitly returns the naive Wald table;
  `se=:bootstrap` and `method=:mcmle` report full inference.
  Dyad-independent fits are unchanged.
- **`show` keeps the point-estimate caveat under `se=:bootstrap`**: the
  bootstrap replaces the covariance, not the estimate, which is still not the
  MLE R reports.
- **`simulate_count_ergm` has an adaptive support.** The explicit form no
  longer truncates an unbounded reference at a literal `max_val=20` (a Poisson
  mean of 25 simulated a mean of 18.0), and the fitted form no longer stays on
  the fit's bound: with `max_val=nothing` the support doubles whenever a
  dyad's conditional puts more than 1e-10 of its mass on the top value, and a
  chain that reaches `2^max_doublings` times its starting bound throws an
  `ArgumentError` (a model that is not normalisable). `max_val=k` simulates
  the truncated family and warns, with the share of draws on the bound, when
  the bound carries more than `BOUNDARY_MASS_TOL` of a conditional; `gof` and
  `se=:bootstrap` refuse such a chain. Draws for models whose conditionals
  stay far below the bound are unchanged.
- **A joint distribution with a mode on the truncation bound is flagged.** The
  support check examines dyad conditionals at the observed network only; a
  model such as Poisson + `mutual.product` with a positive coefficient passed
  it as `support_control = :converged` and was then simulated, bootstrapped
  and GOF-tested on a chain pinned to the bound. A deterministic probe
  (conditional modes iterated from the all-`max_val` configuration) now sets
  `fit.boundary_mode`, reports `support_control = :boundary_mode` with
  `support_stable = false`, warns, and `simulate_count_ergm`, `gof`,
  `se=:bootstrap` and `method=:mcmle` refuse the fit. The probe runs at the
  bound where the adaptive sampler gives up, `2^max_doublings` times the
  fitted bound, not at the fitted bound, where small positive coefficients
  went unseen.
- **A proper geometric model with strong reciprocity is no longer refused.**
  `sum + mutual(:min)` under `GeometricReference()` with θ_sum + θ_mutual > 0
  > 2θ_sum + θ_mutual (R fits it) was flagged as having "a mode on the
  truncation bound" and refused by the default fit, simulation and `gof`.
  The joint-support probe now counts a mode on the bound only when its joint
  log-weight exceeds the observed network's, and the analytic rule reads
  linear growth under the geometric reference, so an improper geometric
  model (2θ_sum + θ_mutual > 0) is still refused, now with the reason.
- **Improper models are flagged analytically.** Under a Poisson or geometric
  reference a positive leading coefficient on a statistic that grows faster
  than the reference decays — `mutual.product`, `nodeOSum`/`nodeISum`/
  `nodeSum`, `sum(pow>1)`, `CMP` beyond `log y!` — makes the model not
  normalisable whatever the data. The rule evaluates the leading order of the
  log-weight along a dyad, a reciprocated pair, a star and every dyad (so a
  larger negative squared-strength term can offset a positive one); a fit it
  flags has `fit.improper = true` and `support_control = :improper`, warns,
  and is refused by `simulate_count_ergm`, `gof`, `se=:bootstrap` and
  `method=:mcmle` unless `max_val` fixes the truncated family. The explicit
  `simulate_count_ergm(net, terms, θ)` refuses such coefficients too. A
  Poisson + `mutual.product` MPLE of 0.06 used to be certified `:converged`
  and simulated silently.
- **An MCMLE that does not converge quotes its own stopping rule.** The
  warning, `show` and `approximations` give the equivalence test p-value
  (with its threshold, tolerance precision and sample size) and the step
  length γ under `termination=:confidence`, and the Hotelling p-value and the
  largest t-ratio only under `:hotelling`; `fit.mcmc` records
  `conv_precision`/`conv_confidence`. The docstrings say which sample each
  diagnostic describes.
- **`method=:mcmle` is an estimator**, not an error; any `method` other than
  `:auto`/`:mple`/`:mcmle` throws an `ArgumentError` naming them.
- **`SumTerm` has a `pow` field** (`SumTerm()` is unchanged in value and
  label).
- **Reference measures follow `ergm.count`.** `GeometricReference()` is the
  counting measure (no `prob`), `BinomialReference(trials)` is `C(trials, y)`
  (no `prob`), Poisson dyads are `Poisson(λ·e^θ)`; the shape parameters are
  absorbed into the `sum` coefficient.
- **The default count support is error-controlled.** With `max_val` unset the
  fit starts at `max(10, 2·largest count)`, refits at twice the bound and
  stops when the estimates move by at most `support_tol` standard errors, at
  most `support_tol` expected dyads lie past the previous bound and no dyad
  puts more than `BOUNDARY_MASS_TOL` on the new top value. `support_control`,
  `support_stable`, `support_delta`, `omitted_tail` and `boundary_mass`
  record the outcome; a support that never settles is `:unconverged` and
  warned about. Fits of unbounded references change numerically (they are
  closer to the untruncated MLE).
- **Boundary statistics follow R.** A statistic at the extreme of its
  attainable range on every dyad has its coefficient fixed at `∓Inf`
  (standard error 0, p 0, excluded from `dof`/AIC/BIC) with R's warning, and
  the rest are estimated on the restricted supports. `method=:mcmle` (the
  default for a dyad-dependent formula) applies R's `drop=TRUE` too: the
  coefficient is fixed and the rest are the MLE with the statistic held at
  its bound by the sampler, with R's warning; it used to refuse the model.
  `drop=false` (R's `control.ergm(drop=FALSE)`) refuses instead, on both
  estimators. `se=:bootstrap`, simulation and `gof` refuse such a fit.
- **Separation and collinearity are detected.** A design with no finite MPLE
  along a combination of statistics is decided exactly by the ecosystem's
  shared verdict (`NetworkCore.clogit_separation`, R's `mple.existence`
  programme) and follows the shared policy: `converged = false`,
  `separated = true`, the separating terms in `separated_terms`, a warning
  naming them with R's "The MPLE does not exist!", `NaN` z values, p-values
  and `confint`, and `se=:bootstrap`/`method=:mcmle` refused. A numerically
  singular pseudo-Hessian is reported through `hessian_cond` and
  `collinear`.
- **Non-convergence is reported, not rescued**: a fit that exhausts `maxiter`
  is `converged = false`, warned, listed in `approximations`, and never
  `is_exact` (nor is a fit with a `∓Inf` coefficient).
- **Counts must be integers under `:weight`.** A network with edges but no
  `:weight`, a partially weighted one, a non-integer weight, or a negative
  weight outside `DiscUnif2Reference(a < 0, b)` is refused, naming the dyad.
  R's `response="w"` is `fit_ergm_count(net, terms; weight=:w)`.
- **Two-mode networks, masked dyads and binary ERGM.jl terms are refused** by
  every entry point, with the count analogue of a binary term named.
- **Directed-only terms** (`CountMutualTerm`, `CyclicalTiesTerm`,
  `NodeOSumTerm`, `NodeISumTerm`) are refused on an undirected network.
- **Terms on negative counts follow R**: `transitiveweights`,
  `cyclicalweights`, `mutual(form="geometric")`, `CMP` and a non-integer
  `sum(pow=)` are refused where R errors or returns `NaN`/`Inf`.
- **`nonzero` counts dyads with a non-zero value** (a 0-valued edge is a zero
  dyad; a negative count is non-zero), and thresholds at or below zero count
  the zero dyads, as in R.
- **`CountMutualTerm()` is labelled `mutual.min`** (was `mutual.count`) and
  takes R's `form`/`threshold`.
- **Types.** `CountERGMModel{T,D,TT,R}` stores its terms as a `Tuple` and has
  no `directed` field (`is_directed(model)`); `CountERGMResult{M}` gained
  fields (positional layout changed) — read results through the StatsAPI
  verbs.
- **Entry point `fit_ergm_count`**; `ergm_count` is an alias of the same
  function. Custom terms implement
  `change_stat_count(term, net, weights, i, j, old, new)`.
- **Sampler defaults and seeds.** `burnin`/`interval` are counted in sweeps
  and default through the dyad-scaled rule shared with ERGM.jl (20 and 1
  beyond ten nodes; 0.1.0 used 1000 and 100); seeded draws differ from 0.1.0.
- **Minimum Julia 1.12**; the package UUID was regenerated.

### Added

- `count_mcmle` / `fit_ergm_count(...; method=:mcmle)`: ERGM.jl's MCMLE
  iteration (Hummel-stepped Monte-Carlo Newton steps from the MPLE, R `ergm`
  4's confidence stopping rule with sample boosting, or `:hotelling`) on the
  Gibbs sampler with ESS-adaptive thinning (`effective_size`); standard
  errors from the Fisher information at the estimate plus the Monte-Carlo
  error; log-likelihood by path sampling (reference measure included). `fit.method`, `fit.mcmc` record the run. A
  dyad-independent model returns its exact fit.
- Terms `AtmostTerm(k)` (`atmost.k`), `CMPTerm()` (`CMP`, `ergm.count`'s
  Conway–Maxwell–Poisson term) and `SumTerm(pow=p)` (`sum<p>`);
  `TransitiveWeightsTerm`, `CyclicalWeightsTerm`, `SmallerthanTerm`,
  `EqualToTerm`, `InIntervalTerm`; the `mutual` forms `:nabsdiff`,
  `:geometric`, `:product`, `:threshold`.
- `se=:bootstrap`: parametric-bootstrap standard errors on the shared
  `NetworkCore.bootstrap_cov` loop (`n_boot`, `boot_burnin`, `boot_interval`,
  `rng`); replicates without a finite MPLE are excluded with one warning and
  kept as `NaN` rows of `fit.boot_replicates`, and the warning, `show` and
  `approximations` say that the standard errors are then biased downward.
  Results are independent of the thread count.
- The full StatsAPI surface (`coef`, `stderror`, `vcov`, `confint`,
  `loglikelihood`, `nobs`, `dof`, `aic`, `bic`, `coeftable`) and the shared
  result metadata (`estimand`, `objective`, `is_exact`, `se_method`,
  `approximations`).
- `gof(fit)` on the shared `NetworkCore.gof` generic: model statistics and the
  dyad count-value distribution.
- Valued covariate terms `CountNodeMatchTerm`, `CountNodeFactorTerm`,
  `CountAbsDiffTerm`, `CountNodeCovTerm`, `CountNodeOCovTerm`,
  `CountNodeICovTerm` and `CountEdgeCovTerm` (R's `nodematch`, `nodefactor`,
  `absdiff`, `nodecov`, `nodeocov`, `nodeicov`, `edgecov` with
  `form="sum"`/`"nonzero"`), with R's labels, pinned by
  `count_covariates.toml`.
- `coefnames(fit)` (StatsAPI) returns the coefficient labels.
- `drop=` keyword on `count_mple`, `count_mcmle` and `fit_ergm_count`.
- `simulate_count_ergm` keyword `max_doublings`; `count_mple` keywords
  `support_tol`, `max_doublings`, `warn`; `rng` keywords throughout.
- `change_stats_support!` (public): a term's change statistics over a dyad's
  whole support in one call; custom terms need only `change_stat_count`.
- Golden fixtures with provenance: `count_mcmle.toml` (the MCMLE of two
  dyad-dependent models, nine `ergm.count` seeds each),
  `count_geometric_drop.toml` (a proper geometric model and a statistic
  dropped at its bound, nine seeds each), `count_covariates.toml` (the
  covariate terms' values and labels, and a fit of them), `zach_poisson.toml`
  (the exact MLE of a dyad-independent model, held at 1e-6) and
  `count_terms.toml` (statistics and R's labels at 1e-9 on four networks).
- A benchmark suite and allocation-regression gates; a precompile workload
  that replays the README's Quick Start (its first call compiled for 4.4 s;
  now 0.1 s).

### Changed

- The MCMLE's building blocks (`mcmle_solve`, `confidence_test`,
  `ess_sample`, `mcmle_covariance`, `bridge_integrate`, `mcmc_defaults`) and
  the dyad count `nobs` reports (`n_observed_dyads`) come from ERGM.jl's
  extension API, `ERGM.Extension`, not from `public` underscore names of
  ERGM, which no longer exist. No behaviour change.

### Fixed

- `Network{Int32}` (any integer vertex type) works with every term, both
  estimators, simulation and `gof`; it threw a `MethodError`.
- The `confint` refusal names the fit it concerns ("an MPLE fit
  (`method=:mple`) of a dyad-dependent model with the default `se`"), not
  "the default count MPLE".
- The Gibbs sampler uses each term's change statistic (0.1.0 ignored them, so
  structural terms had no effect on simulated networks), and keeps negative
  draws under `DiscUnif2Reference(a < 0, b)`.
- Wide supports no longer defeat Newton: every fit climbs a ladder of
  supports, each warm-started from the one below.
- `compute` is type-stable for every term; diagnostic numbers print with
  three significant digits; p-values are floored, never exactly 0.
- Error messages for an observed count outside the support distinguish a
  truncation the caller chose from a bound that is part of the model.

### Performance

- The Gibbs sweep allocates nothing per dyad once warmed up and costs about
  0.6 µs per dyad (3.1 KB and 3.7 µs before); the term profiles are
  O(degree + |support|).
- The MPLE design compresses dyads with identical conditionals into one row:
  a dyad-independent model costs the same at any network size (n = 1000,
  `max_val=30`: 2.9 GiB and 10 s → 0.8 MiB and 0.32 s), and a derivative
  evaluation allocates 192 bytes.

### Known limitations

- **No missing-data MLE, contrastive divergence or stochastic
  approximation.** `method` is `:auto`, `:mple` or `:mcmle`; a network with
  masked dyads is refused by both estimators.
- **`ergm.count`'s MCMC proposals and `control.ergm` tuning are not
  reproduced.** The sampler is a Gibbs sweep over each dyad's enumerated
  conditional; R seeds and tuning constants do not carry over.
- **The reference defaults to `PoissonReference()`**, where R's `ergm`
  requires `reference=` for a valued model.
- **The count MPLE is not an `ergm.count` estimator.** The default uses it
  only where it is the exact MLE (a dyad-independent formula); on a
  dyad-dependent formula `method=:mple` differs from the MLE (1.7–1.8
  standard errors on `zach`), and the default MCMLE reproduces R.
- **Joint normalisability is decided analytically only along five
  configurations**: super-linear growth (`mutual.product`, the squared
  strengths, `sum(pow>1)`, `CMP`) and, under the geometric reference, linear
  growth. Other growth — a user-defined term, a covariate term's
  `form=:sum` — is caught only by the probe at `2^max_doublings` times the
  fitted bound and the sampler's adaptive support; a model improper only
  beyond that bound, and metastable below it, is not detected.
- **An MCMLE with a statistic at its bound has no log-likelihood**
  (`loglikelihood`, AIC and BIC are `NaN`).
- **No `StdNormal` or continuous `Unif` reference**; the estimator enumerates
  integer supports. (`CMP` is a term, as in `ergm.count`.)
- **No `nodecovar`/`nodeocovar`/`nodeicovar`/`nodesqrtcovar`**, no valued
  `nodemix`, `absdiff(pow≠1)` or `nodematch(keep=, levels=)`; the covariate
  terms are one statistic per term (R's `nodefactor` is one term per level).
- **`transitiveweights`/`cyclicalweights` implement the default
  `(min, max, min)` triple only.**
- **`TransitiveTiesTerm`/`CyclicalTiesTerm` and the squared-strength terms
  `NodeOSumTerm`/`NodeISumTerm`/`NodeSumTerm` have no R counterpart** and are
  not validated against R.
- **On negative counts, `sum(pow=)` with a non-integer power and `CMP` are
  refused** (R returns `NaN` and `Inf`).
- **No curved terms, `constraints=`, offsets or two-mode networks**; a
  two-mode network is refused.

## [0.1.0] - 2026-02-09

Development version, never released.
