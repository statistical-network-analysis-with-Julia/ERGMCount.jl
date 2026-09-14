# Estimation

ERGMCount.jl estimates count ERGM parameters using Maximum Pseudo-Likelihood Estimation (MPLE). This page covers the estimation procedure, configuration, diagnostics, and best practices.

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
`Networks.newton_fit` optimizer (Newton-Raphson with step-halving). Dyads with
identical conditionals are compressed into one row, so a dyad-independent
model costs the same whatever the network size.

For a **dyad-independent** model (`SumTerm`, `NonzeroTerm`,
`GreaterthannTerm`, `CountAtleastnTerm`, `SmallerthanTerm`, `EqualToTerm`,
`InIntervalTerm` only) each dyad's conditional is its marginal, so the pseudo-likelihood *is* the likelihood and the MPLE is the
exact MLE. `ergm.count`'s Monte-Carlo MLE (`method=:mcmle`) is **not
implemented**; asking for it throws an `ArgumentError` that says so.

## Fitting a Model

The examples below assume the packages are loaded and a small count-valued
network exists:

```julia
using Networks, ERGMCount
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
result = ergm_count(net, terms; reference=PoissonReference(1.0))
```

### Full Options

```julia
result = ergm_count(net, terms;
    reference = PoissonReference(1.0),  # Reference measure
    method = :mple,                     # Estimation method (the only one)
    max_val = nothing,                  # nothing: error-controlled support
    support_tol = 1e-3,                 # stopping rule of the support doubling
    se = :hessian,                      # or :bootstrap (n_boot, rng, ...)
    maxiter = 100,                      # Maximum Newton iterations
    warn = true                         # false silences the fit diagnostics
)
```

### Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `reference` | `AbstractReferenceMeasure` | Baseline distribution for edge values | `PoissonReference()` |
| `method` | `Symbol` | Estimation method; only `:mple` exists | `:mple` |
| `weight` | `Symbol` | The edge attribute holding the counts (R's `response="w"` is `weight=:w`); any name but `:weight` fits a copy of the network | `:weight` |
| `max_val` | `Int` or `nothing` | Fixed top of the enumerated support for an unbounded reference; `nothing` chooses it adaptively | `nothing` |
| `support_tol` | `Float64` | Stop doubling the support when the estimates move by at most this many standard errors *and* at most this many expected dyads sit past the previous bound (and no dyad has appreciable mass at the new bound) | `1e-3` |
| `max_doublings` | `Int` | Cap on the support doublings; hitting it reports `support_stable = false` with a warning | `8` |
| `se` | `Symbol` | `:hessian` (inverse pseudo-Hessian) or `:bootstrap` (parametric bootstrap; `n_boot`, `boot_burnin`, `boot_interval`, `rng`) | `:hessian` |
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
result.support_control   # :converged (or :fixed, :bounded, :unconverged)
result.support_stable    # false only for :unconverged
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
drop the term.

```julia
acyclic = network(5; directed=true)            # 1→2→3→4→5: nothing reciprocated
for (i, j) in ((1, 2), (2, 3), (3, 4), (4, 5))
    add_edge!(acyclic, i, j); set_edge_attribute!(acyclic, :weight, i, j, 2)
end
b = fit_ergm_count(acyclic, [SumTerm(), CountMutualTerm()]; reference=BinomialReference(3))
coef(b)[2]         # -Inf
dof(b)             # 1
is_exact(b)        # false
```

### Standard errors: Hessian vs parametric bootstrap

`se=:hessian` (the default) reports the inverse negative pseudo-Hessian. For
a dyad-independent model the pseudo-likelihood is the likelihood and these
are the usual MLE standard errors. For a dyad-dependent model the
pseudo-likelihood multiplies overlapping conditionals as if independent, so
they are **anticonservative** — `show` says so, and
`Networks.approximations(result)` lists it.

`se=:bootstrap` Gibbs-simulates `n_boot` networks at the estimates with
`simulate_count_ergm`, refits the count MPLE on each, and reports the
empirical covariance of the refits on the shared `Networks.bootstrap_cov`
loop (`n_boot`, `boot_burnin`, `boot_interval`, `rng`, exactly as
`ERGM.mple`). The point estimates are unchanged; only the covariance is
replaced. Replicates without a finite, converged MPLE are excluded with one
aggregate warning and kept as `NaN` rows of `result.boot_replicates`.

```julia
robust = ergm_count(net, terms; reference=PoissonReference(1.0),
                    se=:bootstrap, n_boot=50, rng=Xoshiro(2))
coef(robust) == coef(result)        # true
stderror(robust)                    # typically larger than the Hessian ones for mutual.min
Networks.se_method(robust)          # :bootstrap
```

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
`result.gradient_norm`, listed in `Networks.approximations(result)`, printed
by `show` ("Converged: false" plus the caveat), and never `is_exact`.
Running out of `maxiter` is reported, not silently topped up. `warn=false`
silences the fit diagnostics without changing what the result records.

### Separation: the MPLE does not exist

A boundary statistic is one column at its extreme. The pseudo-likelihood can
also be flat along a *combination* of the statistics — quasi-complete
separation, R's "The MPLE does not exist!" — which no single-column test
sees. The textbook case is `sum + nonzero` on a network whose every count is
0 or 1: `sum − nonzero` is at its minimum on every dyad, so the objective
keeps rising as θ_sum → −∞, θ_nonzero → +∞ and Newton stops somewhere on
that asymptote with coefficients around ±20 and standard errors in the tens
of thousands. `count_mple` tests the two signatures such an asymptote leaves
(as `ERGM.mple` does): an unobserved support value whose exponential tilt is
more than 18.42 nats below the row's largest, *and* a next Newton step still
larger than `1e-3·‖θ‖` (or an uninvertible Hessian). Both together mean the
fit is returned with `converged = false` and `separated = true`, warned
about with R's sentence, listed in `approximations`, printed by `show`, and
never `is_exact`; the bootstrap excludes such replicates.

```julia
binary = network(6; directed=false)          # seven ties, every count 1
for (i, j) in ((1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (1, 3), (2, 5))
    add_edge!(binary, i, j); set_edge_attribute!(binary, :weight, i, j, 1)
end
sep = fit_ergm_count(binary, [SumTerm(), NonzeroTerm()]; warn=false)
sep.separated        # true
sep.converged        # false
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
function, and `fit_count_ergm` is kept for backward compatibility:

```julia
# These are all the same function
result = fit_ergm_count(net, terms; reference=PoissonReference())
result = ergm_count(net, terms; reference=PoissonReference())
result = fit_count_ergm(net, terms; reference=PoissonReference())
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
| `separated` | `Bool` | The MPLE does not exist along a combination of the statistics (then `converged` is `false`) |
| `hessian_cond`, `collinear` | `Float64`, `Vector{String}` | Condition number of the negative pseudo-Hessian, and the statistics loading on its flat direction when it exceeds `1e8` |
| `max_val`, `truncated`, `support_control`, `support_stable`, `support_tol`, `support_delta`, `omitted_tail`, `boundary_mass` | | The support used and how it was chosen |
| `se_type`, `boot_replicates` | | `:hessian`/`:bootstrap`, and the bootstrap refits |

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
Support:   0:40  (TRUNCATED — reference is unbounded; chosen by doubling until max|Δθ| ≤ 0.001·SE and the omitted tail ≤ 0.001; achieved 3.57e-10·SE, tail 1.61e-10)
Boundary mass: 6.97e-32 (max over dyads, at the fitted coefficients)
Pseudo-log-likelihood: -115.8861
AIC: 237.77, BIC: 246.42  (pseudo-likelihood; compare only across models on the same network and support)
Converged: true
Std. errors: inverse pseudo-Hessian

Coefficients:
            Estimate  Std.Error  z value  Pr(>|z|)
sum           1.0529     0.1221   8.6249    <1e-16 ***
nonzero      -4.2415     0.4265  -9.9452    <1e-16 ***
mutual.min    0.2183     0.1729   1.2628    0.2066
---
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Warning: this model contains dyad-dependent terms and was fit by
maximum pseudolikelihood. The standard errors are the inverse
pseudo-Hessian and are expected to be anticonservative; refit with
`se=:bootstrap` for a parametric-bootstrap covariance.
```

(The coefficient table is the shared presentation layer from Networks.jl,
used identically across the ecosystem's model packages; it is exactly
`coeftable(result)`.)

### Accessing Results

```julia
coef(result)             # coefficient vector
stderror(result)         # standard errors
vcov(result)             # covariance matrix
confint(result)          # normal-theory 95% limits, one row per coefficient
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
| Boundary statistic | coefficient `±Inf`, SE 0, a warning naming the statistic | The statistic is at its extreme attainable value on every dyad (e.g. `mutual.min` with no reciprocated pair); no finite MPLE exists — drop the term |
| Separation | `separated = true`, `converged = false`; a warning "the MPLE does not exist (perfect separation)"; coefficients around ±20 with astronomical SEs | A combination of the statistics is at its extreme on every dyad (`sum + nonzero` when every count is 0 or 1); remove or coarsen a term, or collect more varied counts |
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
result1 = ergm_count(net, terms1; reference=PoissonReference())

# Model 2: Add reciprocity
terms2 = [SumTerm(), NonzeroTerm(), CountMutualTerm()]
result2 = ergm_count(net, terms2; reference=PoissonReference())

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
ERGM.jl (`ERGM._mcmc_defaults`, a budget in single-dyad moves, divided by
the number of dyads): **20 sweeps of burn-in and `cld(max(100, n_dyads ÷
10), n_dyads)` sweeps between retained draws** — 1 for any network with more
than ten nodes — and `gof` and the parametric bootstrap (`boot_burnin`,
`boot_interval`) resolve their defaults the same way. The fitted-result
form simulates at the fit's own support; the explicit form
`simulate_count_ergm(net, terms, coefs)` keeps a literal `max_val=20`. Each
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
    max_val=20,
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
3. **Validate with simulation**: Compare observed and simulated network statistics
4. **Start simple**: Begin with basic terms, add complexity gradually
5. **Match reference to data**: Choose a reference that reflects your data's properties
6. **Keep the reference fixed**: the only reference parameter is the Poisson rate $\lambda$ (confounded with the `sum` coefficient — leave it at 1 unless a baseline rate is known) and the binomial `trials`; the geometric shape and the binomial success probability come from the estimated `sum` coefficient
7. **Read a `±Inf` coefficient as a statement about the data**: the statistic is at its boundary and no finite estimate exists; drop the term
8. **Inspect standard errors**: NaN or very large SEs indicate estimation problems
