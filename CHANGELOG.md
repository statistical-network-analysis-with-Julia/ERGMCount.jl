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

- **Two-mode (bipartite) networks are refused.** `fit_ergm_count`,
  `CountERGMModel`, `simulate_count_ergm` and (through the result's network)
  `gof` throw an `ArgumentError` on `network(n; bipartite=k)` and on a
  `BipartiteNetwork`, with ERGM.jl's reasoning: the estimator enumerates and
  the Gibbs sweep resamples every off-diagonal dyad, so the impossible
  within-mode dyads were counted as observed zeros (0.2.0-dev fitted
  `network(6; bipartite=3)` with `nobs = 15`; the README claimed there was
  no entry point). *Migration:* none — bipartite count ERGMs are not
  implemented; fit the one-mode projection you mean.
- **The MPLE's non-existence by separation is detected** (the count analogue
  of `ERGM._separated`, R's "The MPLE does not exist!"): a design on which
  the pseudo-likelihood is flat along a combination of the statistics —
  `sum + nonzero` when every count is 0 or 1, `sum + atleast(2)` when every
  count is 0 or 2 — used to "converge" in one iteration to θ ≈ (−21, 21)
  with a 2.4e4 standard error and no warning. It is now returned with
  `converged = false` and the new field `separated = true`, warned about,
  listed in `approximations`, printed under `Converged:` by `show`, never
  `is_exact`, and excluded from the bootstrap. The two-signature test (an
  unobserved support value pushed > 18.42 nats down by the tilt alone, and a
  next Newton step > 1e-3·‖θ‖ or an uninvertible Hessian) is documented on
  `count_mple`. *Migration:* a fit that reported large opposite-signed
  coefficients with astronomical standard errors was never a maximum; drop
  or coarsen a term.
- **A numerically singular pseudo-Hessian is loud.** The new field
  `hessian_cond` records the condition number of the negative pseudo-Hessian
  at the estimates; above `ERGMCount._HESSIAN_COND_TOL` (1e8) the fit warns
  naming the statistics that load on the flat direction (`collinear`, e.g.
  `["greaterthan.2", "atleast.3"]` — the same statistic on integer counts —
  or `["sum", "nonzero"]` on a 0/1 network), `approximations` carries the
  caveat and `show` prints it under `Converged:`. The former "not negative
  definite … NaN" warning is subsumed (an exactly singular Hessian has
  `hessian_cond == Inf` and NaN standard errors, said in the same sentence).
- **The counts must be integers under `:weight`.** `CountERGMModel` (hence
  every fit) refuses a network with edges but no `:weight` attribute at all —
  it used to fit as a 0/1 network on which `sum ≡ nonzero`, printing a
  non-identified fit with `Converged: true` — and refuses a non-integer
  (`2.5`, `"3"`) or, under any reference but `DiscUnif2Reference(a < 0, b)`,
  a negative weight, naming the dyad, the value and the rule (a bare
  `InexactError`/`MethodError` before). An integer-valued Float (`2.0`) is a
  count. *Migration:* R's `response="w"` is the new keyword
  `fit_ergm_count(net, terms; weight=:w)`, which fits a `copy` of the network
  with the counts under `:weight` (`fit.model.network` is then the copy); a
  genuinely binary network is fit by setting `:weight` to 1 on every edge.
- **`compute(NonzeroTerm(), net)` counts the dyads with `dyad_value ≠ 0`**,
  not `ne(net)`: an edge whose `:weight` is 0 is a zero dyad (R stores no such
  edge and its `nonzero` says 1 on a two-tie network with one 0-valued tie;
  the estimator's design already read it as 0, so `gof`'s observed statistic
  disagreed with its simulated ones). `GreaterthannTerm(k)` with `k < 0` and
  `CountAtleastnTerm(k)` with `k ≤ 0` count the zero-valued dyads too, as
  R's `greaterthan(-1)`/`atleast(0)` (561 on zach) and as `SmallerthanTerm`
  already did; both rows, and a zero-valued-tie network, are pinned by the
  regenerated `count_terms` fixture. *Migration:* none for networks whose
  edges all carry positive counts.
- **Negative counts under `DiscUnif2Reference(a < 0, b)` survive the Gibbs
  sweep.** `_set_dyad!` treated `y > 0` as "edge present" and removed the
  edge otherwise, so a draw of −1 or −2 became 0: under
  `DiscUnif2Reference(-2, 2)` with θ = 0 the chain produced {0, 1, 2} with
  P(0) = 0.59 instead of the uniform law, and `gof`/`se=:bootstrap` of an
  (exact, correct) fit on negative data were silently wrong. An edge now
  exists for every non-zero count, negative ones stored as negative
  `:weight`s; `gof`'s count-value panel runs over every value the networks
  take (labels `"-2"`, …). *Migration:* none; the estimator was already
  correct, only the simulated side changes.
- **`show(::CountERGMResult)` prints the actual support range** —
  `Support:   1:5  (bounded reference)` for `DiscUnif2Reference(1, 5)`,
  `-2:2` for `DiscUnif2Reference(-2, 2)` — instead of `0:max_val` for every
  reference.
- **`CountERGMModel{T,D,TT,R}` is parameterised on the network's vertex type
  and directedness** (mirroring `ERGM.ERGMModel{T,D}`), on the term tuple type
  and on the reference type; it has **no `directed` field** (use
  `is_directed(model)`), and **`terms` is stored as a `Tuple`**, not a
  `Vector{AbstractERGMTerm}`, so change statistics fold statically. The
  constructor is `CountERGMModel(terms, net, reference=PoissonReference())`
  with `terms` a Tuple, a Vector or a single term. *Migration:*
  `CountERGMModel(AbstractERGMTerm[...], net, ref, is_directed(net))` →
  `CountERGMModel([...], net, ref)`; `model.directed` → `is_directed(model)`;
  code that pushed onto `model.terms` must build a new model.
- **`CountERGMResult{M<:CountERGMModel}`** is parameterised on the model type
  (mirroring `ERGMResult`) and gained the fields `z_values`, `p_values`,
  `iterations`, `gradient_norm`, `support_control`, `support_stable`,
  `support_tol`, `support_delta`, `omitted_tail`, `boot_replicates`,
  `separated`, `hessian_cond` and `collinear`
  (positional layout changed); the 9-argument legacy positional constructor
  (no `se_type`) is removed. *Migration:* read results through the StatsAPI verbs
  (`coef`, `stderror`, `coeftable`, ...) rather than positional construction.
- **`CountMutualTerm` is labelled `mutual.min`, not `mutual.count`, and
  carries R's `form`/`threshold`.** `CountMutualTerm(form=:min;
  threshold=0)` implements `ergm`'s valued `mutual(form=, threshold=)` with
  R's coefficient labels (`mutual.min`, `mutual.nabsdiff`,
  `mutual.geom.mean`, `mutual.product`, `mutual.<threshold>`); the
  zero-argument `CountMutualTerm()` is unchanged in value (`Σ min(y_ij,
  y_ji)`) but its name changed. *Migration:* `coeftable(fit)["mutual.count"]`
  → `coeftable(fit)["mutual.min"]`; code pattern-matching the label must
  follow. The struct now has fields (`form`, `threshold`), so
  `CountMutualTerm` is no longer a singleton type.
- **Directed-only terms are refused on an undirected network.**
  `CountMutualTerm` (every form), `CyclicalTiesTerm`, `NodeOSumTerm` and
  `NodeISumTerm` declare `ERGM.requires_directed`; `CountERGMModel`, `fit_ergm_count` and
  `simulate_count_ergm` throw an `ArgumentError` naming the term and the
  alternative (`NodeSumTerm`, or dropping the term) instead of fitting a
  column that is identically zero. `compute` still returns 0 there.
- **A statistic at the boundary of its attainable range no longer "converges"
  to a large finite number.** When a statistic sits at its smallest (largest)
  attainable value on every dyad's conditional support — `mutual.min` on a
  network without a reciprocated pair, `nonzero` on a complete network, every
  statistic on an empty network — no finite MPLE exists. As in R `ergm`, the
  coefficient is now fixed at `-Inf` (`+Inf`) with standard error 0 and
  p-value 0, the remaining coefficients are estimated with each dyad's support
  restricted to the values the fixed statistic allows (the exact limit of the
  pseudo-likelihood), `dof`/`aic`/`bic` count only the finite coefficients,
  `count_mple` warns with R's sentence, and `show`/`approximations` carry a
  note. Before, the test fixture's separated `mutual.min` was reported as
  θ̂ ≈ −21 with a 1.7e4 standard error and three stars. `se=:bootstrap`,
  `simulate_count_ergm` and `gof` refuse such a fit (there is no finite model
  to simulate from). *Migration:* a `-Inf`/`+Inf` coefficient is a statement
  about the data; drop the term or change the model.
- **The default count support is error-controlled (panel 2026-09 item 17,
  N6).** With `max_val` unset, `count_mple` no longer enumerates a fixed
  `max(10, 2·max count)` — on `zach` that bound moved the estimates by
  1.4e-6/3.9e-6, *above* the golden fixture's own 1e-6 tolerance, which
  passed only because the testset pinned `max_val=30`. It now fits at that
  bound, refits at twice it, and stops at the first doubling at which **all**
  of these hold: every estimate moved by at most `support_tol` (default 1e-3)
  standard errors (absolute scale where the SE is NaN/0), the wider fit leaves
  at most `support_tol` expected dyads past the previous bound, and no dyad
  puts more than `BOUNDARY_MASS_TOL` on the new top value; the fit at the
  larger bound is reported. `fit.support_control` records `:bounded`,
  `:fixed`, `:converged` or `:unconverged` (the doubling hit `max_doublings`,
  default 8, with the estimates still moving — warned with the achieved δ,
  the final `max_val` and the improper-family explanation, listed in
  `approximations`, printed by `show`), `fit.support_stable` is `false` only
  then, and `support_delta`/`omitted_tail` are the achieved bounds (`0.0` for
  a bounded reference, which never loops; `NaN` for a caller-fixed
  `max_val`). `max_val=k` still fixes the bound with the diagnostics computed.
  The default path is pinned against the `zach` fixture at 1e-6 (reproduced
  to ~1e-12 at the doubled bound 28; the old bound 14 is shown to miss it).
  *Migration:* fits of unbounded references with the default support change
  numerically versus 0.1.x (they are closer to the untruncated MLE) and may
  use a larger `max_val` than before (28 rather than 14 on `zach`); pass
  `max_val` to reproduce a fixed-bound fit.
- **`is_exact(fit)` is `false` for an unconverged fit and for a fit with a
  coefficient fixed at `±Inf`** (ERGM.jl's semantics): a limit is not a
  maximizer, and an unconverged iterate is not the MLE of anything. Before,
  a dyad-independent model under a bounded reference was reported exact
  whatever happened to the Newton iteration.
- **Running out of `maxiter` is reported, not rescued.** The coordinate-ascent
  rescue start now runs only when Newton stopped *before* its budget (a
  singular Hessian, a step no halving improves); a fit that exhausts
  `maxiter` is returned with `converged = false`, a warning naming `maxiter`,
  `tol` and the pseudo-score norm, an `approximations` entry ("Newton did not
  converge in N iterations"), and `fit.iterations`/`fit.gradient_norm`
  recording what happened. *Migration:* raise `maxiter` where a tiny budget
  was silently topped up before.
- **Gibbs `burnin`/`interval` default through the dyad-scaled rule shared with
  ERGM.jl.** `simulate_count_ergm`, `gof` and the bootstrap's
  `boot_burnin`/`boot_interval` default to `nothing`, resolved by
  `ERGM._mcmc_defaults(n_dyads)` converted from toggles to sweeps: 20 sweeps of
  burn-in and `cld(max(100, n_dyads ÷ 10), n_dyads)` sweeps between retained
  draws, i.e. 20 / 1 on any network with more than ten nodes (0.1.0: fixed
  1000 and 100 sweeps). The explicit-specification
  form `simulate_count_ergm(net, terms, coefs)` keeps its literal
  `max_val=20`. *Migration:* pass explicit integers to reproduce the old
  spacing.
- **Seeded Gibbs draws differ from 0.1.x.** The sampler draws each dyad by
  one inverse-CDF uniform on the un-normalised conditional weights instead of
  `StatsBase.sample(rng, support, Weights(probs))`, so the same `rng` seed
  yields a different (equally valid) chain. *Migration:* none; results that
  depended on a particular seed's draws must be regenerated.
- **`method=:mcmle` (or any `method` but `:mple`) throws an `ArgumentError`
  explaining what is available** (see "Not implemented" below); before it was
  `Unknown method: mcmle`.
- Removed the private helpers `_z_pvalues`, `_has_dyad_dependent` (the shared
  `Networks.z_pvalues` and a method of `ERGM.has_dyad_dependent` replace them),
  `_dyads` and the dense per-dyad `_count_derivatives(X, log_h, support,
  y_obs_idx, n_dyads, n_terms)` signature (now
  `_count_derivatives(design, cols, mask)` over the compressed design).

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

- **Round 3 (panel 2026-09-08 follow-up) — refusals R makes, said by name:**
  - **Negative counts and the terms R refuses.** On a network with a negative
    dyad weight (admissible only under `DiscUnif2Reference(a < 0, b)`)
    `ergm` 4.12.0 refuses `transitiveweights`/`cyclicalweights` ("Term may
    not be used with networks with negative dyad weights") and its
    `mutual(form="geometric")` returns `NaN`; round 2 returned 3.0/7.0 for
    the weights terms (the `best = 0` floor of the two-path search defining a
    statistic R never produces) and a bare `DomainError` from `sqrt` for the
    geometric mutual. `compute`, `CountERGMModel` (hence every fit) and
    `simulate_count_ergm` (over a support reaching below 0) now throw an
    `ArgumentError` carrying R's sentence for `TransitiveWeightsTerm`,
    `CyclicalWeightsTerm` and `CountMutualTerm(:geometric)`; the other terms
    are unchanged and pinned against R on a seeded 6-actor network with
    values in -2:2 added to `count_terms.toml` (R's two errors and the `NaN`
    are frozen alongside, and the script stops if `ergm` ever starts
    accepting them).
  - **A partially weighted network is refused.** `_validate_count_weights`
    refused only the all-edges-unweighted case; an edge without `:weight`
    among weighted ones (a weight column with NAs, a merge that dropped rows)
    silently counted as 1 — the same `response="w"` mistake on a subset of
    the edges, and R's `response=` never reads a missing value as 1. The
    model constructor (hence `fit_ergm_count`, `count_mple`) and
    `simulate_count_ergm` now refuse it, naming the first bare edge and how
    many there are, with the same `set_edge_attribute!`/`weight=` guidance.
    *Migration:* set the count of every edge (1 if a bare edge really means
    one event).
  - **A binary ERGM.jl term in a count model is refused at construction.**
    `fit_ergm_count(net, [SumTerm(), ERGM.Edges()])` — the natural first
    attempt after `ergm(net ~ edges + mutual, response="w")` — passed
    validation and died in the design build with a `MethodError` on
    `change_stat_count(::Edges, …)`; `CountERGMModel` now throws an
    `ArgumentError` naming the term as a binary ERGM.jl term and its count
    analogue (`NonzeroTerm()` for `edges`, `CountMutualTerm()` for `mutual`,
    `TransitiveWeightsTerm()` for `triangle`, …). `simulate_count_ergm(terms,
    net, coefs)` with the arguments swapped is the same named
    `ArgumentError` `fit_ergm_count(terms, net)` already was.
  - **`PoissonReference(λ)` validates `λ`:** `λ ≤ 0`, `NaN` or `Inf` is an
    `ArgumentError` naming the reference (it used to surface as a
    `DomainError` from `log` or newton_fit's "objective is not finite").
  - **`count_mple(model)` honours the missing-data contract itself:** it
    calls `require_observed(model.network; face_ok=false)` as `fit_ergm_count`
    does — a masked network reached the pseudo-likelihood at its face value
    through this documented entry point (`nobs` counted the masked dyad as an
    observed zero). The refusal names `clear_missing_dyads!`, never `:face`.
  - **`count_mple` and `dyad_value` are exported** (`count_mple` as ERGM.jl
    exports `mple`; both are documented API and the truncation warning names
    `count_mple`), and `change_stats_support!` stays `public`; the test
    suite pins `Base.isexported`/`Base.ispublic` for all three.
  - **Every export's docstring carries a runnable example** (16 were
    prose-only: the reference measures, `SumTerm`, `NonzeroTerm`,
    `GreaterthannTerm`, `CountAtleastnTerm`, the strength terms,
    `log_reference`, `sample_reference`, `is_truncating`, `change_stat_count`,
    `BOUNDARY_MASS_TOL`); a testset executes every ```julia fence of the
    source docstrings and asserts the example is there.

- **R-parity terms (panel 2026-09, item 31; ERGMCount WP4):**
  `TransitiveWeightsTerm()` and `CyclicalWeightsTerm()` — `ergm`'s
  `transitiveweights("min","max","min")` / `cyclicalweights("min","max","min")`,
  each dyad's value capped by its strongest closing two-path, summed over
  ordered pairs on a directed network and unordered pairs on an undirected
  one (R's convention; on an undirected network the two coincide, as in R),
  labels `transitiveweights.min.max.min` / `cyclicalweights.min.max.min`,
  dyad-dependent, with support profiles that register one clamp ramp per
  affected pair (one two-path search each, not one per support value; 0 B,
  pinned) — and the dyad-independent `SmallerthanTerm(k)` (`smallerthan.k`,
  zero dyads counted), `EqualToTerm(v; tolerance=t)` (`equalto.v.pm.t`) and
  `InIntervalTerm(a, b; open=(true, true))` (`ininterval(a,b)`,
  `ininterval[a,b]`, `ininterval(a,b]`, `ininterval[a,b)`; `±Inf` bounds).
  `CountMutualTerm(form; threshold)` gains R's `:nabsdiff`, `:geometric`,
  `:product` and `:threshold` forms (see Breaking). Every new change
  statistic is brute-force tested, every profile held equal to the per-value
  definition, and every term sits in the test suite's `ALL_TERMS`, so the
  Gibbs 0-B pin, the design-build pin and the static-fold pin cover it.
  The non-default weights triples (`geomean`/`sum`) remain unimplemented.
  `TransitiveTiesTerm`/`CyclicalTiesTerm` are kept, with docstrings stating
  that they have **no R counterpart** and are not validated against R.
- **Provenanced term-parity fixture against `ergm` 4.12.0:**
  `test/fixtures/count_terms.toml`, regenerable with `Rscript
  test/fixtures/r/count_terms.R > test/fixtures/count_terms.toml`, freezes
  the summary statistics *and the coefficient names* R prints for `sum`,
  `nonzero`, `greaterthan(2)`, `greaterthan(4)`, `atleast(3)`,
  `transitiveweights`, `cyclicalweights`, `smallerthan(2)`, `equalto(3)`,
  `ininterval(1,3)` (all four bracket forms), `mutual(form =
  min/nabsdiff/geometric/product)` and `greaterthan(-1)`/`atleast(0)` on
  `zach` (undirected, 13 rows) and on a seeded 8-actor directed count matrix
  frozen in the fixture (17 rows), at 1e-9, plus a 3-actor network with a
  0-valued tie (R's `nonzero` is 1 there) and — since round 3 — a seeded
  6-actor directed network with values in -2:2 (11 rows; R's refusal of
  `transitiveweights`/`cyclicalweights` and its `NaN` for
  `mutual(form="geometric")` frozen alongside); the testset rebuilds every
  network from the TOML, compares the names exactly and every statistic row
  by row. `mutual(form="threshold")` is
  **not** pinned: ergm 4.12.0's own `summary()` fails at C model
  initialisation, the fixture records the error, and
  `CountMutualTerm(:threshold)` is tested by hand value (the definition
  R's `emptynwstats` implies: both counts `>= threshold`). The
  `zach_poisson` fixture is untouched.
- **README, docs, CLAUDE.md rewritten to the current behaviour** (panel
  2026-09 items 2 and 8): runnable Quick Start on a real count network, the
  reference-measure comments (`p^y (1-p)^(n-y)` is gone), a Fitting section
  (error-controlled support, `se=:bootstrap`, boundary statistics,
  non-convergence, the StatsAPI surface), Simulation/GOF with the
  dyad-scaled sweep defaults, Missing dyads, Validation against R (both
  fixtures and the exact-MLE argument), and a Not-implemented section
  modelled on ERGM.jl's; the docs gain the real `show` output, a
  Boundary-statistics / Hessian-vs-bootstrap / pseudo-AIC-BIC /
  Non-convergence quartet in the estimation guide, an R-correspondence
  table in the terms guide saying per term whether it is an exact
  counterpart pinned by `count_terms.toml`, and `@docs` entries for every
  new term. Every ```julia block in README and `docs/src` executes under
  the site's `tools/check_snippets.jl`; the Documenter build is strict.
- **The full StatsAPI surface (panel 2026-09, item 15):** `confint(fit;
  level)`, `aic(fit)`, `bic(fit)` and `coeftable(fit)` join `coef`,
  `stderror`, `vcov`, `loglikelihood`, `nobs`, `dof`, all re-exported.
  `aic`/`bic` are **pseudo-likelihood** criteria (`-2·loglik + 2·dof` and
  `-2·loglik + dof·log(nobs)`), comparable only across models on the same
  network, reference and support; `coeftable` returns the shared
  `Networks.CoefficientTable` labelled with the two-argument `name(term, net)`,
  and `show(fit)` prints exactly that table (same vectors), so the two cannot
  disagree. Pinned by `Networks.check_statsapi(fit; strict=true)` on a Poisson,
  a bounded-reference and a bootstrap fit.
- `fit_ergm_count(net, term)` (a single term) and `fit_ergm_count(net,
  terms::Tuple)`; `simulate_count_ergm` accepts the same three forms. Calling
  `fit_ergm_count(terms, net)` with the arguments swapped is an `ArgumentError`
  naming the right order rather than a `MethodError`.
- `has_dyad_dependent(::CountERGMModel)`, a method of the ONE
  `ERGM.has_dyad_dependent` predicate (exported), replacing the same-named
  private.
- **`simulate_count_ergm` refuses a network with masked dyads** through
  `Networks.require_observed` (it used to Gibbs-resample a masked network from
  its face values; `fit_ergm_count` already refused). `gof` inherits the guard.
  `Networks.missing_policies` is declared `(:error,)` for `fit_ergm_count`,
  `count_mple` and `simulate_count_ergm`; the refusal message names
  `clear_missing_dyads!` and never a `missing=:face` keyword that does not
  exist.
- The observed-count-outside-support error branches on
  `is_truncating(reference)`: for a bounded reference (`BinomialReference(3)`
  with an observed 5) it says the bound is part of the model, names the
  reference and the dyad and suggests a wider one instead of advising
  `max_val`, which those references ignore; for an unbounded one it advises
  the `max_val` the caller chose.
- `count_mple` keywords `support_tol` and `max_doublings` (above), and
  `warn::Bool=true`: `warn=false` silences the boundary, truncation,
  non-convergence and undefined-SE warnings without changing what the result
  records (the bootstrap refits use it, as `ERGM.mple`'s do). A non-converged
  Newton iteration warns (naming `maxiter`, `tol`, the pseudo-score norm and
  the two ways out), is listed in `approximations`, printed by `show`, and
  makes `is_exact` false.
- **`change_stats_support!(dest, term, net, weights, i, j, old, support)`**
  (declared `public`, documented): the change-statistic profile of one dyad
  over its whole conditional support, which both the MPLE design build and
  the Gibbs conditional consume. The fallback calls `change_stat_count` per
  value, so a custom count term needs nothing new; the strength terms and
  the triadic terms specialise it to O(degree + |support|) (see
  Performance). A testset holds every specialisation equal to the per-value
  `change_stat_count` for every term, both directednesses, every support
  value and supports that do not start at 0.
- **Benchmark harness and CI regression gates (panel 2026-09, items 7, 24e,
  25).** `benchmark/Project.toml` (`[sources]` at `..`, `../../ERGM.jl`,
  `../../Networks.jl`; instantiates from the committed metadata),
  `benchmark/benchmarks.jl` (BenchmarkTools suite printing the `BENCHJL`
  lines the site's `tools/run_benchmarks.jl` consumes: the zach fixture fit on
  the default path and at `max_val=60`, one derivative evaluation, one Gibbs
  sweep at n = 34 and n = 100 with the same mean degree, and a `SCALING`
  assertion that the per-dyad sweep cost grows by at most 2× between them)
  and `benchmark/regression_tests.jl` (0 B per Gibbs dyad update on a grown
  chain, O(unique slabs) design build, ≤ 512 B per derivative evaluation,
  the scaling ratio). CI runs the gates on the ubuntu/`1` cell, which now
  also runs the test suite with `JULIA_NUM_THREADS=4` so the fresh-process
  thread-count-independence test compares against a different count; the
  workflow's clone list is derived from `[sources]` and says so.
- **PrecompileTools workload (panel 2026-09, item 18).** Two 6-node count
  networks (directed and undirected), the Poisson default path with
  `sum + nonzero`, a `BinomialReference(5)` fit with a dyad-dependent term,
  `show`/`coeftable`/`confint`/`aic`/`bic`, `simulate_count_ergm(n_sim=1)`,
  `gof(n_sim=2)` and a guarded two-replicate `se=:bootstrap`, under a devnull
  logger. Measured in a fresh process (Julia 1.12.6): `using ERGMCount`
  0.8 s → 0.7 s; first `fit_ergm_count` 2.2 s → 0.01 s; first `gof`
  0.5 s → 0.01 s (package precompile 2 s → 8 s).
- A fresh-process **thread-count independence test** for `se=:bootstrap`: the
  replicates are drawn serially from the caller's `rng` and only the refits
  are threaded (in `Networks.bootstrap_cov`), so the standard errors are
  bit-identical under a different `--threads`; and an `@allocated == 0` pin
  on the per-dyad slab fill `_fill_slab!` for every term.
- **Continuation in the support.** A cold Newton start (θ = 0) on a wide
  support overshoots — `GeometricReference()` with `sum + nonzero` on `0:40`
  gave up after two iterations, which the docs' "Comparing Reference
  Measures" example silently exhibited — so every fit now climbs a ladder of
  supports (`0:base`, `0:2·base`, ..., the requested top; `base = max(10,
  2·max count)`), each rung warm-started from the one below; bounded
  references with a wide support (`BinomialReference(100)`) climb the same
  ladder, and the bootstrap refits start at θ̂. When Newton cannot move from
  a start at all (a reference so far from flat on the support that the
  initial conditional is a point mass and the Hessian singular —
  `BinomialReference(100)` on `0:10`), a coordinate-ascent start (bisection
  on each coordinate's score, monotone by concavity) is tried before giving
  up; the common path is untouched. A converged fit whose pseudo-Hessian is
  not negative definite (NaN standard errors) now warns once from
  `count_mple`; `newton_fit`'s own warning is silenced because the rescue
  would otherwise trigger it spuriously. Pinned by regression tests.
- Bootstrap replicates on which the count MPLE does not exist (a boundary
  statistic in the simulated network) or does not converge are **excluded**
  from the covariance with a warning — the NaN-SE hazard — and kept as `NaN`
  rows of the new `fit.boot_replicates`; `show`/`approximations` report the
  count. Fewer than two usable replicates is an error.

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

- **Shared contracts imported, never re-implemented (panel 2026-09, items
  13/14/28):** `newton_fit`, `z_pvalues`, `check_se` and `CoefficientTable`
  come from Networks.jl (`ERGM.newton_fit === Networks.newton_fit`); the
  floored, NaN-aware `Networks.z_pvalues` replaces the unfloored private copy
  (a |z| = 14 coefficient prints `<1e-16`, never `0.0`); the `se=` keyword is
  validated by `Networks.check_se` with the shared message shape
  (`count_mple: se must be one of (:hessian, :bootstrap) (got :sandwich)`).
- **StatsBase dependency dropped.** The Gibbs draw is an inline inverse-CDF
  loop (`_draw_index`) instead of `StatsBase.sample(..., Weights(...))`; the
  bootstrap's exclusion covariance uses stdlib `Statistics.cov`.
- The Gibbs sweep does **not** adopt `ERGM.mh_toggle!`: it redraws each dyad
  from its full conditional (a Gibbs sweep), not a Metropolis toggle, so the
  kernel does not apply. Recorded so nobody tries.

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
- Simulation defaults: the fixed 1000/100-sweep `burnin`/`interval` of 0.1.0
  are replaced by the dyad-scaled `ERGM._mcmc_defaults` rule in sweeps (20
  sweeps of burn-in, 1 between draws beyond ten nodes; see Breaking) — pass
  integers to reproduce the old spacing.
- `CountERGMModel` prints as `ERGM.ERGMModel` does — `CountERGMModel{T,D}:
  n vertices, m edges (directed); terms: sum + nonzero + mutual.min;
  reference: PoissonReference(1.0)` — instead of the raw struct with the
  whole network inside.
- `fit_ergm_count`'s docstring and the `method=:mcmle` error name every
  dyad-independent term (`sum`, `nonzero`, `greaterthan`, `atleast`,
  `smallerthan`, `equalto`, `ininterval`), and the `CountERGMResult`
  docstring states the WP2 stop rule as the conjunction it is (delta AND
  omitted tail AND boundary mass), not the superseded "or".

### Fixed

- **Round 3:**
  - **The adaptive support stops on a rung whose Newton failed.** On an
    exactly collinear design (`sum + greaterthan(2) + atleast(3)` on zach)
    round 2 ran all 8 doublings (`max_val` 14 → 3584, the design growing
    2^8×), each rung's Newton breaking at iteration 1 on the singular
    Hessian, and then reported `support_control = :unconverged` with the
    "not normalisable" sentence — a false diagnosis of a collinear design.
    `_adaptive_support_fit` now returns at the first rung whose fit is not
    usable (`converged == false`, or `hessian_cond > 1e8`) with
    `support_control = :unconverged` and `support_delta = omitted_tail =
    NaN`; `count_mple` says the doubling stopped because that fit did not
    converge / is singular and lets the non-convergence and conditioning
    warnings carry the diagnosis, the "not normalisable" sentence being
    reserved for a doubling that genuinely ran out of `max_doublings`. `show`
    prints "adaptive doubling stopped here because this fit did not
    converge; NOT error-controlled" on the support line. A separated fit
    stops at its first rung the same way.
  - **`compute` is type-stable for every term:** `_get_weights` read the
    untyped `Dict{Tuple{Int,Int},Any}`, so `compute` of `SumTerm`,
    `TransitiveTiesTerm`, `CyclicalTiesTerm`, `TransitiveWeightsTerm` and
    `CyclicalWeightsTerm` inferred `Any` and boxed every `+` (8.4 KB on a
    259-edge network). It now reads the typed snapshot
    `get_edge_attribute(net, :weight, Int)` (an integer-valued Float weight
    converts, as the validator already allowed); `Base.return_types` is
    `[Float64]` for every term in `ALL_TERMS` on both directednesses, pinned.
  - **Diagnostic numbers print with three significant digits:**
    `round(x, sigdigits=3)` printed `Boundary mass: 6.969999999999999e-32`
    (the docs quoted `6.97e-32` under "these are the numbers you get") and
    `condition number 1.6699999999999998e33`; every number in `show`,
    `approximations` and the warnings goes through `_fmt3` (`@sprintf
    "%.3g"`, `Printf` is a new stdlib dependency). The docs' quoted outputs
    are re-checked against a fresh `println(result)`.
  - **Prose that described 0.1.x/round-2 behaviour:** `nonzero` is `Σ I(y ≠
    0)` in the README and the terms guide (negative counts are non-zero); the
    term-parity fixture description names the 13 zach rows, the zero-valued
    tie network and the negative-count network; the CLAUDE.md sentence on
    the bootstrap caveat says what `show` prints (nothing) and what
    `approximations` keeps; `dyad_value`/`compute`/`fit_ergm_count` no
    longer say an edge without `:weight` counts as 1; the landing page
    renders the module docstring once.
  - `[compat] Statistics = "1.11.1"` (the version bundled with Julia 1.12),
    matching Networks/ERGM/ERGMMulti instead of forcing an upgradable-stdlib
    fetch.
- **N6 (panel 2026-09 item 17):** the default count support is
  error-controlled and pinned against the fixture on the default path (see
  Breaking).
- **Silent non-convergence:** `count_mple` returned `converged = false`
  without a word, `is_exact` could still be `true`, and a tiny `maxiter` was
  quietly topped up by the rescue start. Now warned, recorded and never exact
  (see Breaking).
- **Boundary statistics** no longer "converge" to a large finite number with
  a huge standard error (see Breaking); `GreaterthannTerm(k)` with every
  observed count above `k` is `+Inf`, `mutual.min` without a reciprocated
  pair `-Inf`, an empty network all `-Inf`.
- **Bootstrap NaN-SE hazard:** a replicate whose refit is unconverged, hits a
  boundary or returns non-finite coefficients no longer enters the
  covariance; it is a `NaN` row of `boot_replicates`, excluded with one
  aggregate warning (the refits themselves run with `warn=false`).
- The Gibbs sampler previously ignored change statistics (a placeholder
  `coef·y` conditional), so structural terms (mutuality, transitivity, node
  strength) had no effect on simulated draws; it now uses each term's
  `change_stat_count`.
- `set_edge_attribute!` argument order corrected in the sampler; on edge
  removal the typed weight snapshot's entry is zeroed, not deleted (see
  Performance), so `dyad_value` — which reads it behind `has_edge` — never
  sees a stale count.
- A negative observed count under a truncating reference was reported as
  "pass `max_val` ≥ -1"; it is now "counts must be non-negative integers
  under PoissonReference", with the `DiscUnif2Reference` hint (and
  `DiscUnif2Reference(a, b)` with `y < a` keeps the wider-reference hint).
- The test suite pins that the one cross-package `Pkg._name` reach-in in the
  source (`ERGM._mcmc_defaults`) targets a `public` binding
  (`Base.ispublic`), and that no other appears outside comments.

### Performance

- **De-abstracted Gibbs sweep and compressed MPLE design (panel 2026-09, item
  25).** `CountERGMModel.terms` is a `Tuple`, and the per-dyad work is folded
  over it by `@generated` functions (`_accumulate_conditional!` for the
  conditional's linear predictor, `_fill_slab!` for a design slab; 0 B at 36
  terms, above Base's 32-element `map` limit), so nothing dispatches
  dynamically per dyad. Both folds read each term through its **support
  profile** `change_stats_support!` rather than one `change_stat_count` per
  value: the strength terms compute the strength excluding the dyad once and
  fill `(s+y)² − (s+old)²`, the triadic terms collect the third-vertex minima
  once over the neighbour intersections and assemble `Σ_k [min(y, c_k) −
  min(old, c_k)]` for every `y` from a difference array in one pass —
  O(degree + |support|) per dyad instead of O(degree × |support|). The Gibbs
  kernel `_gibbs_update_dyad!` draws by one inverse-CDF uniform on the
  un-normalised weights `exp(η − max η)` (no `probs` vector, no normalising
  pass), touches the network only when the value changes, and zeroes rather
  than deletes the typed weight snapshot's entry on removal (a remove/re-add
  cycle used to rehash the Dict, 32 B amortised); the chain state and every
  retained draw are `copy`s (`Base.copy(::Network)`, the one copier: graph
  and attribute dicts duplicated, vertex attributes preserved, no
  `deepcopy`). Measured per dyad update for `sum + nonzero + mutual.min +
  nodeOSum` on a directed count network at `max_val=10` (Julia 1.12.6, one
  thread, after warm-up):

  | | bytes / dyad | time / dyad, n = 34 | time / dyad, n = 100 |
  |---|---|---|---|
  | committed 0.2.0-dev baseline (`Vector` terms, `sample(Weights)`, `deepcopy`) | 3,108 B | 3.7 µs | 5.4 µs |
  | tuple fold, inverse-CDF draw, `copy` | 0 B (kernel) | 2.1 µs | 4.1 µs |
  | + support profiles (this entry) | **0 B** | **0.6 µs** | **0.8 µs** |

  (The panel measured 1,655 B / 2.7 µs on its own term set.) The 0 B is
  pinned as the **worst over every dyad of a warmed chain, mutation
  included, for every update that does not insert an edge**; an insertion
  is the one mutation that can allocate (Base reallocating Graphs' sorted
  adjacency vector on `insert!` once its front slack is used up — amortised,
  independent of the model, measured 0–35 B per insertion) and is pinned at
  ≤ 128 B per insertion over ten sweeps at 12 and 40 nodes (test suite) and
  five at 34 and 100 nodes (`benchmark/regression_tests.jl`); a sweep in
  which no value changes is exactly 0 B. The benchmark suite's `sweep/n34`
  is 409 µs and 0 allocations for 1,122 dyads and `sweep/n100` 3.0 ms and
  0 allocations for 9,900, a per-dyad ratio of 0.8 with the mean degree held
  fixed, asserted ≤ 2. `_count_mple_fit` compresses dyads with identical (terms × support)
  change-statistic slabs into one row with a per-support-value observation
  count — every dyad of a dyad-independent model shares ONE row whatever the
  network size, so the design is O(unique slabs × support × terms) rather
  than the dense O(dyads × support × terms) tensor (0.5 GB at n=1000,
  max_val=30, p=2 before), and the per-dyad build sweep allocates nothing
  (only a new slab is copied). The derivative closure evaluates in ≤ 512 B
  regardless of rows and support. All three pinned by `@allocated` tests.
  The compressed accumulation reorders floating-point sums relative to the
  0.2.0-dev per-dyad loop, so the fit is no longer claimed bit-for-bit
  identical to it (the "bit-for-bit" claim of the derivative-loop entry
  below is superseded); the exact-MLE golden fixture holds at 1e-6
  (observed ~2e-12).

  Measured, `fit_ergm_count(zach, [SumTerm(), NonzeroTerm()])` (33 actors,
  561 dyads), committed 0.2.0-dev baseline → this release, bytes allocated
  and wall time per fit after warm-up (Julia 1.12.6, one thread):

  | | before | after |
  |---|---|---|
  | `max_val=14` | 1383 KiB, 3.3 ms | 16 KiB, 0.20 ms |
  | `max_val=30` | 2786 KiB, 6.2 ms | 37 KiB, 1.1 ms |
  | `max_val=60` | 5418 KiB, 11.6 ms | 55 KiB, 2.0 ms |
  | default path | 1383 KiB (fixed 14) | 26 KiB (doubled 14 → 28) |
  | n = 1000, `max_val=30` (the panel's case) | 2870 MiB, 10.0 s | 0.8 MiB, 0.32 s |

  The error-controlled default costs one extra fit on the compressed design
  and is still ~50× cheaper in bytes than the old single fixed fit.

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
- Time to first fit: see the PrecompileTools workload under Added.
- Fitting precomputes the per-dyad change-statistic tensor once and reuses
  it across Newton iterations, delegating optimization to the shared
  `ERGM.newton_fit`.

### Not implemented (known limitations, disclosed)

- **Monte-Carlo maximum likelihood (`method=:mcmle`)** — `ergm.count`'s own
  estimator — and with it a missing-data maximum likelihood for masked count
  networks. `fit_ergm_count(...; method=:mcmle)` throws an `ArgumentError`
  saying that ERGMCount.jl fits by maximum pseudo-likelihood only, that the
  MPLE is the exact MLE for a dyad-independent model, and that `se=:bootstrap`
  gives the honest covariance for a dyad-dependent one; a masked network is
  refused. Nothing is silently mis-fit.
- The remaining valued-ERGM surface of `ergm.count`/`ergm` — the `CMP`,
  `StdNormal` and continuous `Unif` references, the
  `nodecovar`/`nodeocovar`/`nodeicovar`/`nodesqrtcovar` family,
  `transitiveweights`/`cyclicalweights` with non-default `(twopath, combine,
  affect)` triples, valued `nodematch`/`nodefactor`/`absdiff`/`edgecov`
  (`form=`), curved and constrained models. There is no term or reference to
  call, so nothing is fit wrongly; `TransitiveWeightsTerm`/
  `CyclicalWeightsTerm` implement the default `(min, max, min)` triple only
  and take no arguments. Two-mode (bipartite) models: a bipartite network is
  refused by every entry point (see Breaking).
- `TransitiveTiesTerm`/`CyclicalTiesTerm` are **not** R terms (a triple-wise
  minimum over ordered triples, kept from 0.1) and are not validated against
  R; their labels say `.count`, not `transitiveweights`.

## [0.1.0] - 2026-02-09

Initial release: count-valued ERGM terms, reference measures, and prototype
estimation/simulation.
