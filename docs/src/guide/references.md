# Reference Measures

Reference measures are the foundation of count ERGMs. They specify the baseline distribution for edge values, playing the same role that the Bernoulli distribution plays in binary ERGMs.

## Overview

In the count ERGM:

$$P(Y = y) \propto h(y) \exp\left(\theta^\top g(y)\right)$$

The reference measure $h(y) = \prod_{(i,j)} h(y_{ij})$ determines the "null model" distribution for each edge value. Without any ERGM terms ($\theta = 0$), the edge values are independent draws from this reference distribution.

## Why Reference Measures Matter

The choice of reference measure affects:

1. **Support**: What values are possible (unbounded vs. bounded)
2. **Baseline shape**: The prior distribution of edge values
3. **Interpretation**: How coefficients relate to changes in edge value probabilities
4. **Convergence**: Some references work better with certain data patterns

## Available Reference Measures

### PoissonReference

The Poisson reference measure is the most common choice for count-valued networks.

```julia
PoissonReference(lambda::Float64=1.0)
```

**Reference function**:

$$h(y_{ij}) = \frac{\lambda^{y_{ij}}}{y_{ij}!}$$

**Support**: $y_{ij} \in \{0, 1, 2, \ldots\}$ (unbounded non-negative integers)

**Properties**:
- Unbounded support means any positive integer is possible
- The parameter $\lambda$ controls the baseline mean
- With $\lambda = 1$ and no ERGM terms, each edge has a Poisson(1) distribution
- Naturally models event counts (emails, calls, meetings)

**When to use**:
- Edge values represent counts of independent events
- No hard upper bound on edge values
- Moderate variance (mean approximately equals variance)

```julia
# Default: lambda = 1.0
ref = PoissonReference()

# Higher baseline rate
ref = PoissonReference(2.0)

# Lower baseline rate (sparser networks)
ref = PoissonReference(0.5)
```

**Log-reference computation**:

```julia
log_reference(ref::PoissonReference, y) = y * log(ref.lambda) - logfactorial(y)
```

### GeometricReference

The Geometric reference is useful for over-dispersed count data.

```julia
GeometricReference(prob::Float64=0.5)
```

**Reference function**:

$$h(y_{ij}) = (1 - p)^{y_{ij}}$$

**Support**: $y_{ij} \in \{0, 1, 2, \ldots\}$ (unbounded)

**Properties**:
- Unbounded support like Poisson
- Has heavier tails than Poisson -- allows more extreme values
- The parameter $p$ controls the decay rate of probabilities
- Smaller $p$ means heavier tails

**When to use**:
- Edge values show more variance than a Poisson would predict
- Some dyads have very high counts while most are zero or low
- Heavy-tailed distributions are expected

```julia
# Default: prob = 0.5
ref = GeometricReference()

# Heavier tails (more high-count edges)
ref = GeometricReference(0.2)

# Lighter tails
ref = GeometricReference(0.8)
```

### BinomialReference

The Binomial reference is appropriate when edge values have a known upper bound.

```julia
BinomialReference(n::Int, prob::Float64=0.5)
```

**Reference function**:

$$h(y_{ij}) = \binom{n}{y_{ij}} p^{y_{ij}} (1-p)^{n-y_{ij}}$$

**Support**: $y_{ij} \in \{0, 1, \ldots, n\}$ (bounded)

**Properties**:
- Bounded support -- edge values cannot exceed $n$
- Natural when each edge represents "successes out of trials"
- The parameter $n$ is the number of trials (maximum edge value)
- The parameter $p$ is the baseline success probability

**When to use**:
- Edge values have a natural upper bound
- Counts represent "k out of n" scenarios
- Examples: meetings out of possible workdays, co-occurrences in fixed-size groups

```julia
# 10 possible meetings, baseline 50% attendance
ref = BinomialReference(10, 0.5)

# 5 weekdays, low baseline
ref = BinomialReference(5, 0.2)
```

### DiscUnifReference

The Discrete Uniform reference assigns equal probability to all values in the support.

```julia
DiscUnifReference(max::Int)
```

**Reference function**:

$$h(y_{ij}) = \frac{1}{\text{max} + 1}$$

**Support**: $y_{ij} \in \{0, 1, \ldots, \text{max}\}$

**Properties**:
- All values are equally likely under the null model
- No prior preference for any particular edge value
- The constant reference cancels in the normalizing constant, simplifying computation

**When to use**:
- No strong prior belief about the baseline distribution
- Data is bounded and you want a "non-informative" reference
- Likert-scale or ordinal data treated as counts

```julia
# Values 0 through 10
ref = DiscUnifReference(10)

# Values 0 through 5
ref = DiscUnifReference(5)
```

### DiscUnif2Reference

The range-bounded Discrete Uniform reference allows specifying both lower and upper bounds.

```julia
DiscUnif2Reference(a::Int, b::Int)
```

**Reference function**:

$$h(y_{ij}) = \frac{1}{b - a + 1}$$

**Support**: $y_{ij} \in \{a, a+1, \ldots, b\}$

**Properties**:
- Uniform distribution over an arbitrary integer range
- Allows non-zero lower bounds
- Useful when zero is not a valid edge value

**When to use**:
- Edge values must be at least some minimum (e.g., rating scales starting at 1)
- Data has known lower and upper bounds

```julia
# Rating scale 1 to 5
ref = DiscUnif2Reference(1, 5)

# Bounded range 0 to 100
ref = DiscUnif2Reference(0, 100)
```

## Comparing Reference Measures

### Decision Tree

```text
Is the count bounded?
|
+-- Yes --> Is there a natural "trials" interpretation?
|           |
|           +-- Yes --> BinomialReference
|           |
|           +-- No  --> Is there a non-zero lower bound?
|                       |
|                       +-- Yes --> DiscUnif2Reference
|                       +-- No  --> DiscUnifReference
|
+-- No  --> Is the variance close to the mean?
            |
            +-- Yes --> PoissonReference
            +-- No  --> GeometricReference
```

### Comparison Table

| Property | Poisson | Geometric | Binomial | DiscUnif |
|----------|---------|-----------|----------|----------|
| Support | $[0, \infty)$ | $[0, \infty)$ | $[0, n]$ | $[0, \text{max}]$ |
| Tail weight | Light | Heavy | Light | None |
| Parameters | $\lambda$ | $p$ | $n, p$ | max |
| Variance | $= \lambda$ | $(1-p)/p^2$ | $np(1-p)$ | $(m+1)^2/12$ |
| Best for | Event counts | Over-dispersed | Bounded trials | No prior |

### Sensitivity Analysis

When unsure which reference to use, compare results across multiple references:

```julia
terms = [SumTerm(), NonzeroTerm(), CountMutualTerm()]

refs = [
    ("Poisson(1.0)", PoissonReference(1.0)),
    ("Geometric(0.5)", GeometricReference(0.5)),
    ("DiscUnif(10)", DiscUnifReference(10)),
]

for (name, ref) in refs
    result = ergm_count(net, terms; reference=ref)
    println("$name: ", round.(result.coefficients, digits=3),
            " converged=", result.converged)
end
```

If coefficients are qualitatively similar across references (same signs, similar magnitudes), the results are robust. Large differences suggest that the reference choice matters for your data, and domain knowledge should guide the selection.

## Sampling from References

Each reference measure supports random sampling, which is used in Gibbs simulation:

```julia
# Draw a random value from the reference
sample_reference(PoissonReference(2.0))    # Random Poisson(2) draw
sample_reference(GeometricReference(0.3))  # Random Geometric(0.3) draw
sample_reference(BinomialReference(10, 0.5)) # Random Binomial(10, 0.5) draw
```

## Log-Reference Values

For numerical stability, ERGMCount.jl works with log-reference values internally:

```julia
# Log-probability under reference
log_reference(PoissonReference(1.0), 3)     # log(1^3 / 3!)
log_reference(GeometricReference(0.5), 3)   # 3 * log(0.5)
log_reference(BinomialReference(10, 0.5), 3) # log(C(10,3) * 0.5^10)
```

## Mathematical Details

### Relationship Between Reference and Conditional Distribution

Given parameters $\theta$ and the rest of the network $y_{-ij}$, the conditional distribution of $y_{ij}$ is:

$$P(Y_{ij} = y \mid Y_{-ij} = y_{-ij}) \propto h(y) \exp\left(\theta^\top \Delta g(y_{ij})\right)$$

Where $\Delta g(y_{ij})$ is the change in sufficient statistics when $y_{ij}$ changes.

### Poisson-SumTerm Connection

With Poisson reference and SumTerm, the conditional distribution of $y_{ij}$ is:

$$P(Y_{ij} = y) \propto \frac{\lambda^y}{y!} \exp(\theta \cdot y) = \frac{(\lambda e^\theta)^y}{y!}$$

This is a Poisson distribution with rate $\lambda \exp(\theta)$. The coefficient $\theta$ on SumTerm multiplicatively adjusts the baseline Poisson rate.

### Identifiability

The reference measure and the ERGM terms must be identifiable together. For example, with a Poisson reference, the SumTerm coefficient and the Poisson $\lambda$ are confounded -- changing one can be compensated by changing the other. In practice, $\lambda$ is fixed (often at 1.0) and only the ERGM coefficients are estimated.
