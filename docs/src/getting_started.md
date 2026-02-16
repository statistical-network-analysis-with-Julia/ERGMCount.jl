# Getting Started

This tutorial walks through common use cases for ERGMCount.jl, from creating count-valued networks to fitting and interpreting models.

## Installation

Install ERGMCount.jl from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/Statistical-network-analysis-with-Julia/ERGMCount.jl")
```

ERGMCount.jl depends on Network.jl and ERGM.jl, which will be installed automatically.

## Basic Workflow

The typical ERGMCount.jl workflow consists of four steps:

1. **Create or load a count-valued network** -- Prepare network data with integer edge weights
2. **Choose a reference measure** -- Select the baseline distribution for edge values
3. **Define count terms** -- Specify which network statistics to include
4. **Fit the model** -- Estimate coefficients via MPLE

## Step 1: Create a Count-Valued Network

Count-valued networks use the `Network` type from Network.jl with edge weights stored as attributes:

```julia
using Network
using ERGMCount

# Create a directed network with 10 nodes
net = Network{Int}(; n=10, directed=true)

# Add edges with integer weights
add_edge!(net, 1, 2)
set_edge_attribute!(net, 1, 2, :weight, 3)   # 3 interactions from 1 to 2

add_edge!(net, 2, 1)
set_edge_attribute!(net, 2, 1, :weight, 1)   # 1 interaction from 2 to 1

add_edge!(net, 1, 3)
set_edge_attribute!(net, 1, 3, :weight, 5)   # 5 interactions from 1 to 3

add_edge!(net, 3, 2)
set_edge_attribute!(net, 3, 2, :weight, 2)   # 2 interactions from 3 to 2
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
net = Network{Int}(; n=n, directed=true)

for i in 1:n, j in 1:n
    if W[i, j] > 0
        add_edge!(net, i, j)
        set_edge_attribute!(net, i, j, :weight, W[i, j])
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
# Poisson: good for unbounded counts (emails, calls)
ref = PoissonReference(1.0)

# Geometric: for counts with high variance
ref = GeometricReference(0.5)

# Binomial: for bounded counts (0 to n)
ref = BinomialReference(10, 0.5)

# Discrete Uniform: equal probability on {0, 1, ..., max}
ref = DiscUnifReference(10)
```

### Choosing the Right Reference

| Data Characteristic | Recommended Reference |
|--------------------|-----------------------|
| Unbounded counts, moderate variance | `PoissonReference` |
| Unbounded counts, high variance | `GeometricReference` |
| Counts bounded by known maximum | `BinomialReference` |
| No prior belief about distribution | `DiscUnifReference` |
| Counts in a known range [a, b] | `DiscUnif2Reference` |

See [Reference Measures](guide/references.md) for a detailed comparison.

## Step 3: Define Count Terms

Count terms are sufficient statistics computed from the valued adjacency matrix:

```julia
# Basic terms
terms = [
    SumTerm(),           # Total edge weight: Sigma y_ij
    NonzeroTerm(),       # Number of non-zero edges: Sigma I(y_ij > 0)
    CountMutualTerm(),   # Weighted mutuality: Sigma min(y_ij, y_ji)
]
```

### Exploring Available Terms

ERGMCount.jl provides terms organized by type:

| Category | Terms | Description |
|----------|-------|-------------|
| **Basic** | `SumTerm`, `NonzeroTerm`, `GreaterthannTerm`, `CountAtleastnTerm` | Aggregate edge value statistics |
| **Structural** | `CountMutualTerm`, `TransitiveTiesTerm`, `CyclicalTiesTerm` | Weighted structural patterns |
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

Use `ergm_count` to estimate model parameters:

```julia
result = ergm_count(net, terms; reference=PoissonReference(1.0))
```

### Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `reference` | Reference measure | `PoissonReference()` |
| `method` | Estimation method | `:mple` |
| `maxiter` | Maximum iterations | `100` |

### Viewing Results

```julia
# Print formatted summary
println(result)

# Output:
# Count ERGM Results
# ==================
# Reference: PoissonReference
# Log-likelihood: -45.6789
# Converged: true
#
# Coefficients:
#   sum                     0.1234 (SE: 0.0456)
#   nonzero                -1.5678 (SE: 0.3210)
#   mutual.count            0.4567 (SE: 0.1234)
```

### Accessing Results Programmatically

```julia
# Coefficient vector
result.coefficients

# Standard errors
result.std_errors

# Check convergence
result.converged
```

## Interpreting Coefficients

Coefficients in count ERGMs relate to the conditional distribution of each edge value given the rest of the network.

| Coefficient | Interpretation |
|-------------|----------------|
| `SumTerm` > 0 | Higher edge values are more likely overall |
| `NonzeroTerm` > 0 | More edges (ties) than expected under reference |
| `CountMutualTerm` > 0 | Tendency for reciprocal intensity matching |
| `NodeOSumTerm` > 0 | Activity heterogeneity (some nodes send more) |
| `NodeISumTerm` > 0 | Popularity heterogeneity (some nodes receive more) |

**Example interpretations:**

- `SumTerm = 0.2` with Poisson reference: Edge values are shifted upward; the conditional mean for each dyad increases by a factor of `exp(0.2) = 1.22`
- `NonzeroTerm = -2.0`: A strong penalty on having edges, making sparse networks more likely (analogous to the edges term in binary ERGMs)
- `CountMutualTerm = 0.5`: Reciprocal dyads tend to have matched intensities

## Complete Example

```julia
using Network
using ERGMCount
using Random

Random.seed!(42)

# Create a count-valued communication network
n = 15
net = Network{Int}(; n=n, directed=true)

# Simulate edge weights with some structure
for i in 1:n, j in 1:n
    i == j && continue
    # Base rate depends on distance
    rate = rand() < 0.3 ? rand(1:8) : 0
    if rate > 0
        add_edge!(net, i, j)
        set_edge_attribute!(net, i, j, :weight, rate)
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
# Simulate 10 networks from the fitted model
sim_nets = simulate_count_ergm(result;
    n_sim=10,
    burnin=1000,
    interval=100,
    max_val=20
)

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

println("Model 1 converged: ", result1.converged)
println("Model 2 converged: ", result2.converged)
```

## Comparing Reference Measures

```julia
# Same terms, different references
terms = [SumTerm(), NonzeroTerm(), CountMutualTerm()]

result_poisson = ergm_count(net, terms; reference=PoissonReference(1.0))
result_geom = ergm_count(net, terms; reference=GeometricReference(0.5))

println("Poisson coefficients: ", round.(result_poisson.coefficients, digits=3))
println("Geometric coefficients: ", round.(result_geom.coefficients, digits=3))
```

## Best Practices

1. **Start with SumTerm and NonzeroTerm**: These are analogous to intercept terms and control for overall scale and density
2. **Choose reference carefully**: Match the reference measure to your data characteristics (see [Reference Measures](guide/references.md))
3. **Check convergence**: Always verify `result.converged == true`
4. **Start simple**: Begin with basic terms before adding structural effects
5. **Validate with simulation**: Compare simulated network statistics to observed ones
6. **Scale matters**: The Poisson parameter lambda in `PoissonReference` affects the baseline rate

## Next Steps

- Learn about [Reference Measures](guide/references.md) in detail
- Explore all [Count Terms](guide/terms.md) available
- Understand the [Estimation](guide/estimation.md) procedure
