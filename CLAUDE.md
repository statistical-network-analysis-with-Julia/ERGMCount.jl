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

1. **Reference Measures** -- Abstract type `AbstractReferenceMeasure` with concrete subtypes (`PoissonReference`, `GeometricReference`, `BinomialReference`, `DiscUnifReference`, `DiscUnif2Reference`). Each implements `log_reference(ref, y)` and `sample_reference(ref)`.

2. **Count-Specific Terms** -- Subtypes of `AbstractERGMTerm` (imported from ERGM.jl). Each term implements `name(t)`, `compute(t, net)`, and optionally `change_stat(...)`. Terms include: `SumTerm`, `NonzeroTerm`, `GreaterthannTerm`, `CountAtleastnTerm`, `CountMutualTerm`, `TransitiveTiesTerm`, `CyclicalTiesTerm`, `NodeOSumTerm`, `NodeISumTerm`, `NodeSumTerm`.

3. **Estimation** -- `ergm_count()` is the main entry point, dispatching to `count_mple()` which performs Newton-Raphson MPLE. Results are stored in `CountERGMResult`.

4. **Simulation** -- `simulate_count_ergm()` generates networks from fitted models via Gibbs sampling over edge values.

## Key Dependencies

- **Network.jl** -- Provides the `Network` type; edge weights stored via `get_edge_attribute(net, :weight)`
- **ERGM.jl** -- Supplies the `AbstractERGMTerm` base type
- **Distributions.jl** -- Probability distributions for reference measures and sampling
- **Optim.jl**, **StatsBase.jl** -- Optimization and weighted sampling utilities

## Conventions

- All term structs subtype `AbstractERGMTerm` and must implement `name()` and `compute()`.
- Edge weights are stored as the `:weight` edge attribute on `Network` objects; absent weight means binary (0/1).
- The `ergm_count` function is the public API; `fit_count_ergm` is an alias.
- Julia >= 1.9 is required (see `[compat]` in Project.toml).
- No test directory exists yet; tests are declared via `[extras]` and `[targets]` in Project.toml.
