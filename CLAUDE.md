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

3. **Estimation** -- `ergm_count()` dispatches to `count_mple()`: a proper count MPLE that enumerates each dyad's full conditional P(y_ij=y|rest) ∝ h(y)·exp(θ'Δg(y)) over the (possibly truncated, `max_val`) support — the reference measure enters the estimator directly. Newton-Raphson with step-halving; `loglik` is the maximized pseudo-log-likelihood.

4. **Simulation** -- `simulate_count_ergm()` Gibbs-samples each dyad from its full conditional using the per-term `change_stat_count`, so structural terms influence draws. Also callable with explicit `(net, terms, coefficients)` without fitting.

## Key Dependencies

- **Network.jl** -- Provides the `Network` type; edge weights stored via `get_edge_attribute(net, :weight)`
- **ERGM.jl** -- Supplies the `AbstractERGMTerm` base type
- **Distributions.jl** -- Probability distributions for reference measures and sampling
- **StatsBase.jl** -- Weighted sampling utilities (Gibbs draws)

## Conventions

- All term structs subtype `AbstractERGMTerm` and must implement `name()` and `compute()`.
- Edge weights are stored as the `:weight` edge attribute on `Network` objects; absent weight means binary (0/1).
- The `ergm_count` function is the public API; `fit_count_ergm` is an alias.
- Julia >= 1.12 is required (see `[compat]` in Project.toml).
- Tests in `test/runtests.jl` include brute-force verification of every `change_stat_count`, hand-computed term values, Poisson-rate recovery of the MPLE, and estimation↔simulation round trips.
