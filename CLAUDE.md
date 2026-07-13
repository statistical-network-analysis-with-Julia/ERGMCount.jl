# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ERGMCount.jl is a Julia package for fitting Exponential Random Graph Models (ERGMs) to count-valued (integer-weighted) networks. It is a port of the R `ergm.count` package from the StatNet collection.

## Development Commands

- **Run tests:** `julia --project=. -e 'using Pkg; Pkg.test()'`
- **Load package locally:** `julia --project=.` then `using ERGMCount`
- **Build docs:** `julia --project=docs docs/make.jl`
- **Dev install:** `using Pkg; Pkg.develop(path=".")`

## Architecture

The entire package lives in a single file: `src/ERGMCount.jl`. It is organized into four sections:

1. **Reference Measures** -- Abstract type `AbstractReferenceMeasure` with concrete subtypes (`PoissonReference`, `GeometricReference`, `BinomialReference`, `DiscUnifReference`, `DiscUnif2Reference`). Each implements `log_reference(ref, y)`, `sample_reference(ref)`, and `_support(ref, max_val)`. Following `ergm.count`/Krivitsky (2012), `GeometricReference` is the counting measure h(y)=1 (no parameter) and `BinomialReference(trials)` is h(y)=C(trials,y) (no probability parameter) — the shape/probability is absorbed into the estimated `SumTerm` coefficient.

2. **Count-Specific Terms** -- Subtypes of `AbstractERGMTerm` (`import ERGM: name, compute`). Every term implements `name(t)`, `compute(t, net)`, and `change_stat_count(t, net, weights, i, j, old, new)` — the change in the statistic when dyad (i,j) moves from count `old` to `new`, holding other dyads fixed (implementations must not read the dyad's own current value). Terms: `SumTerm`, `NonzeroTerm`, `GreaterthannTerm`, `CountAtleastnTerm`, `CountMutualTerm`, `TransitiveTiesTerm`, `CyclicalTiesTerm`, `NodeOSumTerm`, `NodeISumTerm`, `NodeSumTerm` (the NodeO/ISum terms are directed-only).

3. **Estimation** -- `ergm_count()` dispatches to `count_mple()`: a proper count MPLE that enumerates each dyad's full conditional P(y_ij=y|rest) ∝ h(y)·exp(θ'Δg(y)) over the (possibly truncated, `max_val`) support — the reference measure enters the estimator directly. Maximized with the shared `ERGM.newton_fit` optimizer (Newton-Raphson with step-halving); `loglik` is the maximized pseudo-log-likelihood. `_count_mple_fit(model, max_val)` is the core (design array + Newton fit), shared by the fit and by the bootstrap's refits.

   **Standard errors** (`se_type` on `CountERGMResult`, reported by `Networks.se_method`): `se=:hessian` (default) is the inverse negative pseudo-Hessian — anticonservative for any dyad-dependent model, because the pseudo-likelihood multiplies the dyad conditionals as if independent. `se=:bootstrap` is a **parametric bootstrap** (Gibbs-simulate `n_boot` networks at θ̂ with `simulate_count_ergm` at the *same* `max_val`, refit the count MPLE on each, empirical covariance), with the same API as `ERGM.mple`'s (`n_boot`/`boot_burnin`/`boot_interval`/`rng`) and running on the ONE shared `Networks.bootstrap_cov` loop — never paste the loop in. The point estimates are identical under both; only the covariance changes. `show` and `approximations(fit)` read `se_type` and drop the anticonservatism caveat when a bootstrap was actually used (that caveat is a claim about the inverse Hessian), while keeping the *point-estimate* pseudo-likelihood caveat, which holds either way.

4. **Simulation** -- `simulate_count_ergm()` Gibbs-samples each dyad from its full conditional using the per-term `change_stat_count`, so structural terms influence draws. Also callable with explicit `(net, terms, coefficients)` without fitting. All sampling entry points (`simulate_count_ergm`, `sample_reference`) take an `rng::AbstractRNG` keyword (default `Random.default_rng()`); the same rng state gives identical draws. The Gibbs loop snapshots the `:weight` attribute once into a typed `Dict{Tuple{T,T},Int}` (Networks.jl's typed `get_edge_attribute(net, :weight, Int)` accessor) and maintains it incrementally alongside the network — the untyped attribute Dict must not be read inside the hot loop.

## The derivative loop allocates O(p²), and the fit is bit-for-bit (review finding 15)

`_count_derivatives(X, log_h, support, y_obs_idx, n_dyads, n_terms)` is the named builder of the `(ll, grad, hess)` closure `newton_fit` maximizes — named, not buried in `_count_mple_fit`, so the `@allocated` regression test measures the code that runs.

The conditional moments `E[Δg]` and `E[Δg Δg']` used to be rebuilt **per dyad**: `zeros(n_terms)`, `zeros(n_terms, n_terms)`, and a fresh `p .* (x * x')` outer product for *every support value* — `n_dyads × (2 + |support|)` allocations on every Newton evaluation, **2.1 MB and 1.73 ms per evaluation** on `zach` at `max_val = 30`. Worse, it grew with the support, which is exactly the knob a user turns up to make the truncation harmless. They are now filled in place on workspaces allocated once: **192 bytes and 0.36 ms** — 4.8x faster, independent of dyads and support.

The scalar arithmetic is **unchanged**: each entry is still `p * (x[k] * x[l])`, accumulated in the same order, and `hess[k,l] -= eXX[k,l] - eX[k]*eX[l]` is the same subtraction. The fit is therefore **bit-for-bit** what it was — which is the standard this package is held to, because `sum + nonzero` under a Poisson reference is dyad-independent and the golden fixture pins the **exact MLE** at 1e-6. A hand loop, not `BLAS.ger!`: at `n_terms` of 2-5 the BLAS call overhead dominates, and `ger!`'s `(α·x[i])·y[j]` rounds differently from `α·(x[i]·y[j])`, which would have cost the bit-identity for nothing.

## Golden fixtures (statnet ergm.count)

`test/fixtures/zach_poisson.toml`, generated by `test/fixtures/r/zach_poisson.R` and loaded with Networks.jl's `load_golden` (which throws without `[provenance]`). Model: `zach ~ sum + nonzero`, `reference = ~Poisson`, on ergm.count's bundled Zachary karate-club network (counts 0-7).

`sum + nonzero` is **dyad-independent on purpose**: each dyad is then an independent draw from a two-parameter law with an **exact MLE**, and `count_mple`'s dyad-conditional enumeration *is* the likelihood — so ERGMCount.jl computes that exact MLE and is held to **1e-6**, not to "within Monte-Carlo error". It reproduces it to ~1e-12.

Two things the fixture had to get right, and both are worth knowing:

- **The golden value is solved analytically**, not with `optim`. The score equations collapse to one monotone scalar root-find. `optim(BFGS, reltol=1e-14)` — the obvious way — landed **7e-6** from the true optimum, *above* the tolerance it was meant to police. A reference number only good to 7e-6 cannot hold anyone to 1e-6. The frozen residual score is ~1e-14.
- **Truncation is measured, not assumed.** The Poisson reference is unbounded; the estimator enumerates `0:max_val` and the exact MLE truncates nothing, so they estimate the same thing only if the boundary mass is negligible. The fixture freezes `P(y > 30) = 6.9e-23` under the fitted law against a smallest-used conditional probability of 2.3e-3 — nineteen orders of magnitude — and the testset re-asserts ERGMCount.jl's own reported `boundary_mass` against that bound, so a future `max_val` that started to bite goes **red** rather than quietly widening the gap. Choosing data with negligible boundary mass is what makes a count-ERGM comparison meaningful at all.

`ergm.count`'s own fit is MCMLE (statnet has no valued MPLE) and sits **0.0103** from the exact MLE — *further* than ERGMCount.jl does. It is frozen as a cross-check on R, not as the reference standard.

## Key Dependencies

- **Networks.jl** -- Provides the `Network` type; edge weights stored via `get_edge_attribute(net, :weight)`
- **ERGM.jl** -- Supplies the `AbstractERGMTerm` base type
- **Distributions.jl** -- Probability distributions for reference measures and sampling
- **StatsBase.jl** -- Weighted sampling utilities (Gibbs draws)

## Conventions

- All term structs subtype `AbstractERGMTerm` and must implement `name()` and `compute()`.
- Edge weights are stored as the `:weight` edge attribute on `Network` objects; absent weight means binary (0/1).
- `fit_ergm_count` is the standardized entry point (ecosystem `fit_<model>` naming); `ergm_count` is the R-faithful alias and `fit_count_ergm` a legacy alias — all the same function.
- Julia >= 1.12 is required (see `[compat]` in Project.toml).
- Tests in `test/runtests.jl` include brute-force verification of every `change_stat_count`, hand-computed term values, Poisson-rate recovery of the MPLE, and estimation↔simulation round trips.
