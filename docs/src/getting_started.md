# Getting Started

Create a count network with an explicit edge-value attribute, choose its admissible support and reference measure, and fit a baseline before adding dependence terms. Inspect both the estimator caveats and the finite-support diagnostics when interpreting a fit.

!!! note "Before you begin"

    Counts must be nonnegative integers on a loop-free, one-mode network; continuous weights are outside this model. Missing dyads are refused. MCMC-MLE is not implemented. Unbounded reference measures use adaptive finite support, so inspect support diagnostics even for a dyad-independent model.

## Installation

```@raw html
<p>Use Julia <strong>1.12 or newer</strong> and the <a href="/getting-started/">shared workspace installation guide</a>. These development packages are not yet registered; the guide prepares the required sibling checkouts and a Julia environment for the examples.</p>
```

Run the blocks below in order in that environment. They build on variables from earlier steps; stochastic examples use seeded random number generators where shown.

## Basic Workflow

The typical ERGMCount.jl workflow consists of four steps:

1. **Create or load a count-valued network** -- Prepare network data with integer edge weights
2. **Choose a reference measure** -- Select the baseline distribution for edge values
3. **Define count terms** -- Specify which network statistics to include
4. **Fit the model** -- Estimate coefficients via MPLE

## Step 1: Create a Count-Valued Network

Count-valued networks use the `Network` type from Networks.jl with the
counts stored as the **`:weight`** edge attribute. That name is fixed: it is
what R's `ergm.count` selects per call with `response="w"`. If your counts
live under another attribute, pass `fit_ergm_count(net, terms; weight=:w)`
(the fit runs on a copy carrying them under `:weight`); a network that has
edges but no `:weight` at all is refused rather than fitted as a 0/1
network, and a non-integer weight (`2.5`, `"3"`) is refused naming the dyad.

```julia
using Networks
using ERGMCount

# Create a directed network with 10 nodes
net = network(10; directed=true)

# Add edges with integer weights
add_edge!(net, 1, 2)
set_edge_attribute!(net, :weight, 1, 2, 3)   # 3 interactions from 1 to 2

add_edge!(net, 2, 1)
set_edge_attribute!(net, :weight, 2, 1, 1)   # 1 interaction from 2 to 1

add_edge!(net, 1, 3)
set_edge_attribute!(net, :weight, 1, 3, 5)   # 5 interactions from 1 to 3

add_edge!(net, 3, 2)
set_edge_attribute!(net, :weight, 3, 2, 2)   # 2 interactions from 3 to 2
```

### Creating from a Matrix

For larger networks, build from a weight matrix:

```julia
# Weight matrix: w[i,j] = number of interactions from i to j
W = [0 3 5 0;
     1 0 0 2;
     0 2 0 4;
     3 0 1 0]

n = size(W, 1)
net = network(n; directed=true)

for i in 1:n, j in 1:n
    if W[i, j] > 0
        add_edge!(net, i, j)
        set_edge_attribute!(net, :weight, i, j, W[i, j])
    end
end
```

### Inspecting Network Properties

```julia
println("Nodes: ", nv(net))          # Number of vertices
println("Edges: ", ne(net))          # Number of non-zero edges
println("Directed: ", is_directed(net))

# Access edge weights
weights = get_edge_attribute(net, :weight)
println("Total weight: ", sum(values(weights)))
```

## Step 2: Choose a Reference Measure

The reference measure determines the baseline distribution for edge values. This is the key distinction from binary ERGMs.

```julia
# Poisson reference: factorial downweighting on nonnegative integers
ref = PoissonReference(1.0)

# Geometric reference: constant reference weight on nonnegative integers
ref = GeometricReference()

# Binomial: for bounded counts (0 to n)
ref = BinomialReference(10)

# Discrete Uniform: equal probability on {0, 1, ..., max}
ref = DiscUnifReference(10)
```

### Choosing the Right Reference

| Reference | Support and baseline weight |
|-----------|-----------------------------|
| `PoissonReference(1.0)` | Nonnegative integers, weight `1/y!` |
| `GeometricReference()` | Nonnegative integers, constant weight; the fitted model must be normalizable |
| `BinomialReference(n)` | Counts from `0` to `n`, with binomial multiplicity weights |
| `DiscUnifReference(b)` | Counts from `0` to `b`, with constant weight |
| `DiscUnif2Reference(a, b)` | Counts from `a` to `b`, with constant weight |

Choose support from what could have been observed, then assess the reference
measure together with the model terms. A uniform reference is still a
modeling assumption; observed variance alone does not select a reference.

See [Reference Measures](guide/references.md) for a detailed comparison.

## Step 3: Define Count Terms

Count terms are sufficient statistics computed from the valued adjacency matrix:

```julia
# Basic terms
terms = [
    SumTerm(),           # Total edge weight: Sigma y_ij
    NonzeroTerm(),       # Number of non-zero dyads: Sigma I(y_ij != 0)
    CountMutualTerm(),   # Weighted mutuality: Sigma min(y_ij, y_ji)
]
```

### Exploring Available Terms

ERGMCount.jl provides terms organized by type:

| Category | Terms | Description |
|----------|-------|-------------|
| **Basic** (dyad-independent) | `SumTerm`, `NonzeroTerm`, `GreaterthannTerm`, `CountAtleastnTerm`, `SmallerthanTerm`, `EqualToTerm`, `InIntervalTerm` | Aggregate edge value statistics, R's `sum`/`nonzero`/`greaterthan`/`atleast`/`smallerthan`/`equalto`/`ininterval` |
| **Structural** | `CountMutualTerm` (forms `:min`, `:nabsdiff`, `:geometric`, `:product`, `:threshold`), `TransitiveWeightsTerm`, `CyclicalWeightsTerm` | R's valued `mutual`, `transitiveweights`, `cyclicalweights` |
| **Structural, no R counterpart** | `TransitiveTiesTerm`, `CyclicalTiesTerm` | Triple-wise minimum sums (not validated against R) |
| **Degree** | `NodeOSumTerm`, `NodeISumTerm`, `NodeSumTerm` | Node strength heterogeneity |

### Example: Comprehensive Model

```julia
terms = [
    # Edge value effects
    SumTerm(),                  # Overall activity level
    NonzeroTerm(),              # Network density
    GreaterthannTerm(3),        # High-frequency ties

    # Structural effects
    CountMutualTerm(),          # Reciprocity in intensity

    # Degree heterogeneity
    NodeOSumTerm(),             # Sender activity variance
    NodeISumTerm(),             # Receiver popularity variance
]
```

## Step 4: Fit the Model

The examples below fit a 12-node seeded network (the 4-node matrix above has
too few dyads to say much):

```julia
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

result = ergm_count(net, terms; reference=PoissonReference(1.0))
```

`fit_ergm_count` is the standardized name of the same function (`ergm_count`
is the R-faithful alias).

### Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `reference` | Reference measure | `PoissonReference()` |
| `method` | Estimation method; only `:mple` exists (`:mcmle` throws) | `:mple` |
| `max_val` | Fixed top of the enumerated support for an unbounded reference; `nothing` chooses it by error control | `nothing` |
| `support_tol` | Stop doubling the support once the estimates move by at most this many SEs and the omitted tail is this small | `1e-3` |
| `se` | `:hessian` (inverse pseudo-Hessian) or `:bootstrap` (parametric bootstrap; `n_boot`, `rng`) | `:hessian` |
| `maxiter` | Maximum Newton iterations | `100` |
| `rng` | Randomness source for the bootstrap | `Random.default_rng()` |
| `warn` | `false` silences the fit diagnostics (still recorded) | `true` |

### Viewing Results

```julia
println(result)
```

Output:

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

The support line says the enumerated count support was chosen by doubling
(`0:40` here) and how far the estimates moved at the last doubling; the
warning identifies the pseudo-likelihood fit with a dyad-dependent term
(`mutual.min`). Inverse pseudo-Hessian standard errors can be too small.
`se=:bootstrap` estimates covariance by simulation and refitting under the
fitted model; it does not remove the approximation in the point estimate.

### Accessing Results Programmatically

Prefer the StatsAPI verbs to the fields:

```julia
coef(result)             # [1.0529, -4.2415, 0.2183]
stderror(result)         # standard errors
vcov(result)             # covariance matrix
confint(result)          # normal-theory 95% limits, one row per coefficient
coeftable(result)        # the table `show` prints
loglikelihood(result)    # pseudo-log-likelihood
nobs(result), dof(result)
aic(result), bic(result) # pseudo-likelihood criteria

result.converged         # true
result.support_control   # :converged — the doubling stopped
result.max_val           # 40
```

## Interpreting Coefficients

Coefficients in count ERGMs relate to the conditional distribution of each edge value given the rest of the network.

| Coefficient | Interpretation |
|-------------|----------------|
| `SumTerm` > 0 | Higher edge values are more likely overall |
| `NonzeroTerm` > 0 | More edges (ties) than expected under reference |
| `CountMutualTerm` > 0 | Tendency for reciprocal intensity matching (`mutual.min`) |
| `TransitiveWeightsTerm` > 0 | Strong two-paths go with strong closing ties |
| `NodeOSumTerm` > 0 | Activity heterogeneity (some nodes send more) |
| `NodeISumTerm` > 0 | Popularity heterogeneity (some nodes receive more) |

**Example interpretations:**

- `SumTerm = 0.2` with Poisson reference: Edge values are shifted upward; the conditional mean for each dyad increases by a factor of `exp(0.2) = 1.22`
- `NonzeroTerm = -2.0`: A strong penalty on having edges, making sparse networks more likely (analogous to the edges term in binary ERGMs)
- `CountMutualTerm = 0.5`: Reciprocal dyads tend to have matched intensities

## Complete Example

```julia
using Networks
using ERGMCount
using Random
using ERGM: compute, name   # generic term interface shared with ERGM.jl
using Statistics: mean

Random.seed!(42)

# Create a count-valued communication network
n = 15
net = network(n; directed=true)

# Simulate edge weights with some structure
for i in 1:n, j in 1:n
    i == j && continue
    # Construct observed zero counts and positive counts
    rate = rand() < 0.3 ? rand(1:8) : 0
    if rate > 0
        add_edge!(net, i, j)
        set_edge_attribute!(net, :weight, i, j, rate)
    end
end

println("Network: $(nv(net)) nodes, $(ne(net)) edges")

# Define model terms
terms = [
    SumTerm(),              # Total communication volume
    NonzeroTerm(),          # Network density (ties vs non-ties)
    CountMutualTerm(),      # Reciprocity in communication frequency
    NodeOSumTerm(),         # Sender activity heterogeneity
]

# Fit with Poisson reference
result = ergm_count(net, terms; reference=PoissonReference(1.0))

# Display results
println(result)

# Check convergence
if result.converged
    println("\nModel converged successfully")
else
    println("\nWarning: Model did not converge")
end
```

## Simulating Count Networks

After fitting a model, simulate new networks from the estimated parameters:

```julia
# Simulate 10 networks from the fitted model. burnin/interval are Gibbs
# sweeps (every dyad redrawn once per sweep); left unset they follow the
# dyad-scaled rule shared with ERGM.jl — 20 sweeps of burn-in, 1 between
# draws on any network with more than ten nodes — and the draws use the
# fit's own support
sim_nets = simulate_count_ergm(result; n_sim=10, rng=Xoshiro(3))

# ... or spell the controls out
sim_nets = simulate_count_ergm(result; n_sim=10, burnin=50, interval=5, rng=Xoshiro(3))

# Compare observed vs simulated statistics
for term in terms
    obs = compute(term, net)
    sim_vals = [compute(term, sn) for sn in sim_nets]
    println("$(name(term)): observed=$(round(obs, digits=2)), " *
            "simulated mean=$(round(mean(sim_vals), digits=2))")
end
```

## Comparing Models

```julia
# Model 1: Basic effects
terms1 = [SumTerm(), NonzeroTerm()]

# Model 2: Add structural effects
terms2 = [SumTerm(), NonzeroTerm(), CountMutualTerm()]

result1 = ergm_count(net, terms1; reference=PoissonReference())
result2 = ergm_count(net, terms2; reference=PoissonReference())

println("Model 1: AIC=", round(aic(result1), digits=2), " BIC=", round(bic(result1), digits=2))
println("Model 2: AIC=", round(aic(result2), digits=2), " BIC=", round(bic(result2), digits=2))
# pseudo-likelihood criteria: same network, reference and support only
```

## Comparing Reference Measures

```julia
# Same terms, different references
terms = [SumTerm(), NonzeroTerm(), CountMutualTerm()]

result_poisson = ergm_count(net, terms; reference=PoissonReference(1.0))
result_geom = ergm_count(net, terms; reference=GeometricReference())

println("Poisson coefficients: ", round.(coef(result_poisson), digits=3))
println("Geometric coefficients: ", round.(coef(result_geom), digits=3))
```

## Best Practices

1. **Start with SumTerm and NonzeroTerm**: These are analogous to intercept terms and control for overall scale and density
2. **Choose reference carefully**: Match the reference measure to your data characteristics (see [Reference Measures](guide/references.md))
3. **Check convergence**: an unconverged fit warns and prints "Converged: false"; verify `result.converged` and `result.support_stable`
4. **Start simple**: Begin with basic terms before adding structural effects
5. **Validate with simulation**: Compare simulated network statistics to observed ones
6. **Scale matters**: the Poisson rate `λ` in `PoissonReference(λ)` is fixed, not estimated; it is confounded with the `sum` coefficient, so keep it at 1 unless a baseline rate is known

## Next Steps

- Learn about [Reference Measures](guide/references.md) in detail
- Explore all [Count Terms](guide/terms.md) available
- Understand the [Estimation](guide/estimation.md) procedure
