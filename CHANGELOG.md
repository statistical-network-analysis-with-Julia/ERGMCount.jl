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

- **Provenanced golden fixture against `ergm.count` on Zachary's karate club**
  (issue #8). `test/fixtures/zach_poisson.toml`, regenerable with
  `Rscript test/fixtures/r/zach_poisson.R > test/fixtures/zach_poisson.toml`,
  freezes `zach ~ sum + nonzero` under a Poisson reference (ergm.count 4.1.3).

  The model is **dyad-independent on purpose**: each dyad is then an independent
  draw from a two-parameter law with an **exact MLE**, and ERGMCount.jl's
  dyad-conditional enumeration *is* the likelihood — so it computes that exact
  MLE and can be held to 1e-6 rather than to "within Monte-Carlo error".
  ERGMCount.jl reproduces it to **~1e-12** (coefficients and standard errors).

  Two things worth recording:

  - The golden value is solved **analytically** (the score equations collapse to
    one monotone scalar root-find), not with `optim`. `optim(BFGS, reltol=1e-14)`
    — the obvious way to write it — landed **7e-6** from the true optimum, *above*
    the tolerance it was supposed to police. A golden number that is itself only
    good to 7e-6 cannot hold anyone to 1e-6. The frozen residual score is ~1e-14.
  - `ergm.count`'s *own* fit is MCMLE (statnet has no MPLE for valued ERGMs) and
    sits **0.0103** from the exact MLE — further than ERGMCount.jl does. It is
    frozen as a cross-check, not as the reference standard.

  **Truncation is checked, not assumed.** The Poisson reference is unbounded and
  the estimator enumerates `0:max_val`; the exact MLE truncates nothing, so they
  estimate the same thing only if the boundary mass is negligible. The fixture
  freezes `P(y > 30) = 6.9e-23` under the fitted law, against a
  smallest-used conditional probability of 2.3e-3 — nineteen orders of magnitude.
  The testset re-asserts ERGMCount.jl's own reported `boundary_mass` against that
  bound, so a future `max_val` that started to bite goes red rather than quietly
  widening the gap.

- **Robust standard errors: `count_mple(model; se=:bootstrap)`** (also via
  `fit_ergm_count`/`ergm_count`), with the same keywords and semantics as
  `ERGM.mple`'s: `n_boot=100`, `boot_burnin`, `boot_interval`, `rng`. Gibbs-
  simulate `n_boot` count networks at θ̂ with `simulate_count_ergm`, refit the
  count MPLE on each, and report the empirical covariance — on the ONE shared
  `Networks.bootstrap_cov` loop. **The point estimates are unchanged; only the
  covariance is replaced.** Until now the only standard errors available were
  the inverse pseudo-Hessian, which multiplies dyad conditionals as if
  independent and is therefore anticonservative for any dyad-dependent model —
  and they were printed with significance stars (issue #9, ERGMCount#2). On the
  test fixture (`SumTerm` + `CountMutualTerm`) the bootstrap SE of `sum` is ~20%
  **larger** than the Hessian one; that gap is the anticonservatism.
- `se_method(fit)` now reports what was actually used (`:hessian`/`:bootstrap`),
  read off the new `CountERGMResult.se_type` field, and `approximations(fit)`
  and `show` drop the anticonservatism caveat when a bootstrap was used (they
  keep the *point-estimate* pseudo-likelihood caveat, which holds either way).
  `show` now names the standard-error estimator on its own line.

- `gof(::CountERGMResult; n_sim, burnin, interval, max_val, rng)` extending
  the ecosystem-wide `Networks.gof`, comparing model statistics and the
  dyad count-value distribution in a `Networks.GOFResult`.
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
- Results print through the shared `Networks.print_coeftable`; p-values
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

- **The MPLE derivative loop no longer allocates (review finding 15).** The
  conditional moments `E[Δg]` and `E[Δg Δg']` were rebuilt per dyad —
  `zeros(n_terms)`, `zeros(n_terms, n_terms)`, and a fresh `p .* (x * x')` outer
  product for *every support value* — i.e. `n_dyads × (2 + |support|)`
  allocations on every Newton evaluation: **2.1 MB per evaluation** on `zach` at
  `max_val = 30`, and it grew with the support, which is exactly the knob a user
  turns up to make the truncation harmless. They are now filled in place on
  workspaces allocated once (`_count_derivatives`): **192 bytes** per
  evaluation, independent of dyads and support, and **4.8x faster**
  (1.734 ms -> 0.360 ms). The scalar arithmetic is unchanged — each entry is
  still `p * (x[k] * x[l])`, accumulated in the same order — so the fit is
  **bit-for-bit** what it was, which is what the exact-MLE golden fixture
  requires. Pinned by an `@allocated` regression test.
- A typed `:weight` snapshot (`Dict{Tuple{T,T},Int}` via
  `get_edge_attribute(net, :weight, Int)`) is maintained incrementally,
  removing untyped-`Any` lookups from the innermost Gibbs conditional loop.
- Fitting precomputes the per-dyad change-statistic tensor once and reuses
  it across Newton iterations, delegating optimization to the shared
  `ERGM.newton_fit`.

## [0.1.0] - 2026-02-09

Initial release: count-valued ERGM terms, reference measures, and prototype
estimation/simulation.
