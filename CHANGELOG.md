# Changelog

All notable changes to ERGMCount.jl are documented in this file. The format
is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
package adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - Unreleased

Release driven by the 2026-07 expert-panel review: reference measures and
estimation are re-specified to match R `ergm.count` (Krivitsky 2012), the
Gibbs sampler now actually uses change statistics, and the package adopts the
ecosystem-wide `fit_*`/StatsAPI/GOF conventions.

### Breaking

- **Canonical entry point renamed to `fit_ergm_count`.** `ergm_count`
  (R-faithful) and `fit_count_ergm` (legacy) remain as `const` aliases, so
  existing calls keep working. *Migration:* none required; prefer
  `fit_ergm_count`.
- **`GeometricReference` lost its `prob` field** — it is now the zero-field
  counting measure `GeometricReference()`, as in ergm.count; the geometric
  shape comes from a negative `SumTerm` coefficient, not the reference.
  *Migration:* drop the probability argument.
- **`BinomialReference(n, prob)` is now `BinomialReference(trials)`** — the
  success probability is absorbed into the `SumTerm` coefficient.
  *Migration:* drop the second argument.
- **Per-term `change_stat` replaced by the uniform
  `change_stat_count(term, net, weights, i, j, old, new)`** generic (the old
  methods had ad-hoc signatures per term). *Migration:* external callers and
  custom count terms must implement/call `change_stat_count`.
- **`CountERGMResult` gained `vcov` and `max_val` fields** (positional
  layout changed); `loglik` is explicitly the pseudo-log-likelihood.
  *Migration:* update positional constructor calls.
- **`NodeOSumTerm`/`NodeISumTerm` return 0 on undirected networks**
  (previously computed as if directed). *Migration:* use `NodeSumTerm` for
  undirected networks.
- **Minimum Julia raised to 1.12**; package UUID regenerated. *Migration:*
  upgrade Julia and re-resolve environments pinning the old UUID.

### Added

- `gof(::CountERGMResult; n_sim, burnin, interval, max_val, rng)` extending
  the ecosystem-wide `Network.gof`, comparing model statistics and the
  dyad count-value distribution in a `Network.GOFResult`.
- StatsAPI accessors: `coef`, `stderror`, `vcov`, `loglikelihood`, `nobs`,
  `dof`; exported `CountERGMModel`/`CountERGMResult`.
- `rng::AbstractRNG` keywords on fitting, simulation, and
  `sample_reference` for reproducible runs; `max_val` truncation keyword;
  `log_reference`/`sample_reference`/`dyad_value` exported.

### Changed

- Reference measures re-specified to Krivitsky (2012)/ergm.count semantics:
  Poisson dyads `Poisson(λ·e^θ)`, Geometric as counting measure, Binomial
  `C(trials, y)`.
- `count_mple` rewritten as a proper pseudo-likelihood: each dyad's full
  conditional is enumerated over the count support (score
  `Σ[Δg(y_obs) − E_θ(Δg)]`, Hessian `−Σ Var_θ(Δg)`), replacing the crude
  logistic `y > 0` approximation.
- `TransitiveTiesTerm`/`CyclicalTiesTerm` normalized to ordered distinct
  triples with correct directed/undirected change statistics.
- Results print through the shared `Network.print_coeftable`; p-values
  computed with `ccdf(Normal(), |z|)` (no underflow to exactly `0.0`).
- Simulation defaults changed: `burnin` 1000 → 100 sweeps, `interval`
  100 → 10 (a sweep already visits every dyad).

### Fixed

- The Gibbs sampler previously ignored change statistics (a placeholder
  `coef·y` conditional), so structural terms (mutuality, transitivity, node
  strength) had no effect on simulated draws; it now uses each term's
  `change_stat_count`.
- `set_edge_attribute!` argument order corrected in the sampler, and edge
  removal now also deletes the stale weight snapshot entry.

### Performance

- A typed `:weight` snapshot (`Dict{Tuple{T,T},Int}` via
  `get_edge_attribute(net, :weight, Int)`) is maintained incrementally,
  removing untyped-`Any` lookups from the innermost Gibbs conditional loop.
- Fitting precomputes the per-dyad change-statistic tensor once and reuses
  it across Newton iterations, delegating optimization to the shared
  `ERGM.newton_fit`.

## [0.1.0] - 2026-02-09

Initial release: count-valued ERGM terms, reference measures, and prototype
estimation/simulation.
