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

MPLE maximizes the log-pseudo-likelihood with ERGM.jl's shared `newton_fit`
optimizer (Newton-Raphson with step-halving).

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
    method = :mple,                     # Estimation method
    maxiter = 100                       # Maximum iterations
)
```

### Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `reference` | `AbstractReferenceMeasure` | Baseline distribution for edge values | `PoissonReference()` |
| `method` | `Symbol` | Estimation method (currently `:mple`) | `:mple` |
| `maxiter` | `Int` | Maximum Newton-Raphson iterations | `100` |

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
| `model` | `CountERGMModel` | The fitted model specification |
| `coefficients` | `Vector{Float64}` | Estimated coefficients |
| `std_errors` | `Vector{Float64}` | Standard errors |
| `loglik` | `Float64` | Log-pseudo-likelihood at convergence |
| `converged` | `Bool` | Whether optimization converged |

### Displaying Results

```julia
println(result)
```

Output:

```text
Count ERGM Results
==================
Reference: PoissonReference
Pseudo-log-likelihood: -115.8798
Converged: true

              Estimate  Std.Error  z value  Pr(>|z|)
sum             1.0536     0.1224   8.6097    <1e-16 ***
nonzero        -4.2435     0.4272  -9.9321    <1e-16 ***
mutual.count    0.2182     0.1729   1.2618    0.2070
---
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

(The coefficient table is the shared presentation layer from Networks.jl,
used identically across the ecosystem's model packages.)

### Accessing Results

```julia
# Coefficient vector
result.coefficients

# Standard errors
result.std_errors

# Model specification
result.model.terms      # Vector of terms
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
| Non-convergence | `converged = false` | Increase `maxiter`, simplify model |
| Large coefficients | $|\theta| > 10$ | Possible separation; remove problematic term |
| Large standard errors | SE >> coef | Multicollinearity; remove correlated terms |
| NaN in standard errors | `NaN` values | Singular Hessian; simplify model |

### Handling Non-Convergence

```julia
# Increase iterations
result = ergm_count(net, terms;
    reference=PoissonReference(),
    maxiter=500
)

# Check for problematic terms
for (i, term) in enumerate(terms)
    println("$(name(term)): coef=$(round(result.coefficients[i], digits=4)), " *
            "SE=$(round(result.std_errors[i], digits=4))")
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

println("Model 1: ", length(terms1), " terms, converged=", result1.converged)
println("Model 2: ", length(terms2), " terms, converged=", result2.converged)
```

### Comparing Reference Measures

```julia
terms = [SumTerm(), NonzeroTerm()]

result_pois = ergm_count(net, terms; reference=PoissonReference(1.0))
result_geom = ergm_count(net, terms; reference=GeometricReference())

println("Poisson: ", round.(result_pois.coefficients, digits=3))
println("Geometric: ", round.(result_geom.coefficients, digits=3))
```

## Simulation-Based Validation

After fitting, validate by simulating networks from the estimated model and comparing summary statistics:

```julia
# Simulate from fitted model (pass rng for reproducible draws)
using Random
sim_nets = simulate_count_ergm(result;
    n_sim=100,
    burnin=1000,
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
6. **Fix reference parameters**: Set $\lambda$, $p$, etc. based on domain knowledge, not estimation
7. **Watch for separation**: Very large coefficients indicate model misspecification
8. **Inspect standard errors**: NaN or very large SEs indicate estimation problems
