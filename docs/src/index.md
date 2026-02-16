# ERGMCount.jl

*Exponential Random Graph Models for Count-Valued Networks in Julia*

A Julia package for statistical modeling of networks with integer-valued edge weights using count-valued ERGMs.

## Overview

Exponential Random Graph Models (ERGMs) for count-valued networks extend classical binary ERGMs to handle networks where edges carry integer-valued weights -- counts, frequencies, or intensities. Instead of modeling the presence or absence of ties, count ERGMs model the magnitude of interactions.

ERGMCount.jl is a port of the R [ergm.count](https://cran.r-project.org/package=ergm.count) package from the [StatNet](https://statnet.org/) collection, providing efficient tools for modeling valued network data.

### What is a Count-Valued Network?

In a count-valued network, each dyad $(i, j)$ carries an integer value $y_{ij} \geq 0$ rather than a binary indicator:

```text
Actor A --[3]--> Actor B   (3 emails sent)
Actor B --[1]--> Actor A   (1 email sent)
Actor A --[0]--> Actor C   (no emails)
```

Examples include:

- Number of emails exchanged between colleagues
- Frequency of trade shipments between countries
- Count of co-authored papers between researchers
- Number of phone calls between friends
- Frequency of animal grooming interactions

### Key Concepts

| Concept | Description |
|---------|-------------|
| **Count-Valued Network** | A network where edges carry non-negative integer values |
| **Reference Measure** | Baseline distribution for edge values (e.g., Poisson, Geometric) |
| **Count Term** | A sufficient statistic computed from the valued adjacency matrix |
| **Edge Weight / Value** | The integer count on a dyad |
| **Strength** | Sum of edge weights incident to a node (weighted analogue of degree) |

### Applications

Count ERGMs are widely used in:

- **Communication networks**: Modeling email, phone call, or message frequencies
- **Trade networks**: Analyzing bilateral trade volumes between countries
- **Collaboration networks**: Studying co-authorship or co-invention counts
- **Animal behavior**: Modeling grooming, mating, or aggression event counts
- **Transportation**: Analyzing passenger or cargo flow frequencies

## Features

- **Multiple reference measures**: Poisson, Geometric, Binomial, Discrete Uniform, and range-bounded Discrete Uniform
- **Count-specific terms**: Sum, Nonzero, Greaterthan-n, Atleast-n, Mutual, Transitive ties, Cyclical ties, and node strength terms
- **MPLE estimation**: Maximum Pseudo-Likelihood Estimation via Newton-Raphson optimization
- **Gibbs sampling simulation**: Generate synthetic count-valued networks from fitted models
- **Integration with Network.jl**: Uses the Network.jl edge attribute system for edge weights

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Statistical-network-analysis-with-Julia/ERGMCount.jl")
```

Or for development:

```julia
using Pkg
Pkg.develop(path="/path/to/ERGMCount.jl")
```

## Quick Start

```julia
using Network
using ERGMCount

# Create a network with edge weights
net = Network{Int}(; n=10, directed=true)
for i in 1:10, j in 1:10
    i == j && continue
    w = rand(0:5)
    if w > 0
        add_edge!(net, i, j)
        set_edge_attribute!(net, i, j, :weight, w)
    end
end

# Define count-specific terms
terms = [
    SumTerm(),           # Total edge weight
    NonzeroTerm(),       # Number of non-zero edges (density)
    CountMutualTerm(),   # Weighted mutuality
]

# Fit with Poisson reference measure
result = ergm_count(net, terms; reference=PoissonReference(1.0))

# View results
println(result)
```

## Choosing a Reference Measure

| Reference | Best For | Edge Support |
|-----------|----------|--------------|
| [`PoissonReference`](@ref) | Unbounded counts (emails, calls) | $\{0, 1, 2, \ldots\}$ |
| [`GeometricReference`](@ref) | High-variance counts | $\{0, 1, 2, \ldots\}$ |
| [`BinomialReference`](@ref) | Bounded counts (meetings out of $n$ days) | $\{0, 1, \ldots, n\}$ |
| [`DiscUnifReference`](@ref) | Equal baseline probability | $\{0, 1, \ldots, \text{max}\}$ |
| [`DiscUnif2Reference`](@ref) | Equal baseline on custom range | $\{a, a+1, \ldots, b\}$ |

## Documentation

```@contents
Pages = [
    "getting_started.md",
    "guide/references.md",
    "guide/terms.md",
    "guide/estimation.md",
    "api/types.md",
    "api/terms.md",
    "api/estimation.md",
]
Depth = 2
```

## Theoretical Background

### The Count ERGM

Count ERGMs model the probability of a valued network as:

$$P(Y = y) \propto h(y) \exp\left(\theta^\top g(y)\right)$$

Where:

- $y$ is the observed valued adjacency matrix with $y_{ij} \in \{0, 1, 2, \ldots\}$
- $h(y) = \prod_{(i,j)} h(y_{ij})$ is the reference measure (baseline distribution)
- $g(y)$ is a vector of sufficient statistics (count terms)
- $\theta$ is the parameter vector to be estimated

The reference measure $h(y)$ generalizes the binary ERGM by specifying a non-uniform baseline over count values. When $h(y_{ij}) = 1$ for $y_{ij} \in \{0, 1\}$, the model reduces to a standard binary ERGM.

### Relationship to Poisson Regression

With a Poisson reference and the `SumTerm` statistic, the count ERGM is closely related to Poisson regression. The coefficient on `SumTerm` shifts the mean of the conditional Poisson distribution for each dyad.

### MPLE for Count ERGMs

Maximum Pseudo-Likelihood Estimation conditions on the rest of the network and maximizes the product of conditional likelihoods for each dyad. This provides a computationally efficient approximation to the full MLE.

## References

1. Krivitsky, P.N. (2012). Exponential-family random graph models for valued networks. *Electronic Journal of Statistics*, 6, 1100-1128.

2. Krivitsky, P.N., Hunter, D.R., Morris, M., Klumb, C. (2023). ergm.count: Fit, Simulate and Diagnose Exponential-Family Models for Networks with Count Edges. R package.

3. Hunter, D.R., Handcock, M.S., Butts, C.T., Goodreau, S.M., Morris, M. (2008). ergm: A package to fit, simulate and diagnose exponential-family models for networks. *Journal of Statistical Software*, 24(3).

4. Desmarais, B.A., Cranmer, S.J. (2012). Statistical mechanics of networks: Estimation and uncertainty. *Physica A: Statistical Mechanics and its Applications*, 391(4), 1865-1876.
