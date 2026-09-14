# ERGMCount.jl

Model nonnegative integer counts attached to network dyads: for example, the number of messages sent from one actor to another during a fixed observation window. ERGMCount.jl estimates count-valued ERGMs by maximum pseudo-likelihood and simulates from their reference measures and sufficient statistics.

| First analysis | Learn the model or data | Reference and detail |
|:--|:--|:--|
| [Fit a count network](getting_started.md) | [Choose the count support](guide/references.md) | [Inspect estimation and uncertainty](guide/estimation.md) |

!!! note "Supported scope"

    Counts must be nonnegative integers on a loop-free, one-mode network; continuous weights are outside this model. Missing dyads are refused. MCMC-MLE is not implemented. Unbounded reference measures use adaptive finite support, so inspect support diagnostics even for a dyad-independent model.

## Installation

```@raw html
<p>Use Julia <strong>1.12 or newer</strong> and the <a href="/getting-started/">shared workspace installation guide</a>. These development packages are not yet registered; the guide prepares the required sibling checkouts and a Julia environment for the examples.</p>
```

## Quick Start

Construct counts for six actors and fit a common Poisson baseline. An absent edge represents an observed zero count:

```julia
using Networks, ERGMCount

net = network(6; directed=true)
for i in 1:6, j in 1:6
    i == j && continue
    count = mod(i + 2j, 4)
    count == 0 && continue
    add_edge!(net, i, j)
    set_edge_attribute!(net, :weight, i, j, count)
end
fit = fit_ergm_count(net, [SumTerm()]; reference=PoissonReference(1.0))
display(fit)
```

```@raw html
<p>The sum coefficient shifts the Poisson mean for each dyad. This constructed example is a baseline, with no reciprocity or transitivity term. Counts aggregate an observation window; use <a href="/REM.jl/dev/">REM.jl</a> when event timing itself is the outcome.</p>
```

## Choosing a Reference Measure

| Reference | Baseline weight | Edge Support |
|-----------|----------|--------------|
| [`PoissonReference`](@ref ERGMCount.PoissonReference) | Factorial downweighting | $\{0, 1, 2, \ldots\}$ |
| [`GeometricReference`](@ref ERGMCount.GeometricReference) | Constant reference weight | $\{0, 1, 2, \ldots\}$ |
| [`BinomialReference`](@ref ERGMCount.BinomialReference) | Binomial multiplicity | $\{0, 1, \ldots, n\}$ |
| [`DiscUnifReference`](@ref ERGMCount.DiscUnifReference) | Equal baseline probability | $\{0, 1, \ldots, \text{max}\}$ |
| [`DiscUnif2Reference`](@ref ERGMCount.DiscUnif2Reference) | Equal baseline on custom range | $\{a, a+1, \ldots, b\}$ |

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

Maximum Pseudo-Likelihood Estimation conditions on the rest of the network and maximizes the product of conditional likelihoods for each dyad. For a dyad-independent model this objective equals the likelihood, but unbounded reference measures are still evaluated on adaptive finite support. For dependent models it remains a pseudo-likelihood approximation. `se=:bootstrap` estimates uncertainty by simulation and refitting under the fitted model; its usefulness depends on that model, adequate simulation, and successful refits.

## Not implemented

- **Monte-Carlo maximum likelihood** (`method=:mcmle`, `ergm.count`'s own estimator) and a missing-data MLE for masked count networks: `fit_ergm_count(...; method=:mcmle)` throws an `ArgumentError` explaining that ERGMCount fits by maximum pseudo-likelihood only; masked networks are refused everywhere.
- The remaining valued-ERGM surface of `ergm.count`/`ergm` (the `CMP`, `StdNormal`, continuous `Unif` references; `nodecovar`-family terms; non-default `transitiveweights`/`cyclicalweights` triples; valued `nodematch`/`nodefactor`/`absdiff`/`edgecov`; curved and constrained models). There is no term or reference to call, so nothing is silently mis-fit. A two-mode (bipartite) network is **refused** with an `ArgumentError` by `fit_ergm_count`, `CountERGMModel`, `simulate_count_ergm` and `gof` — enumerating the one-mode dyads would count the impossible within-mode dyads as observed zeros.
- `TransitiveTiesTerm`/`CyclicalTiesTerm` are **not** R terms (a triple-wise minimum over ordered triples) and are not validated against R; `TransitiveWeightsTerm`/`CyclicalWeightsTerm` are the `ergm` counterparts.

## Module Reference

```@docs
ERGMCount.ERGMCount
```

## References

1. Krivitsky, P.N. (2012). Exponential-family random graph models for valued networks. *Electronic Journal of Statistics*, 6, 1100-1128.

2. Krivitsky, P.N., Hunter, D.R., Morris, M., Klumb, C. (2023). ergm.count: Fit, Simulate and Diagnose Exponential-Family Models for Networks with Count Edges. R package.

3. Hunter, D.R., Handcock, M.S., Butts, C.T., Goodreau, S.M., Morris, M. (2008). ergm: A package to fit, simulate and diagnose exponential-family models for networks. *Journal of Statistical Software*, 24(3).

4. Desmarais, B.A., Cranmer, S.J. (2012). Statistical mechanics of networks: Estimation and uncertainty. *Physica A: Statistical Mechanics and its Applications*, 391(4), 1865-1876.


## Citation

If you use ERGMCount.jl in your work, please cite it using the entry in
[`CITATION.bib`](https://github.com/statistical-network-analysis-with-Julia/ERGMCount.jl/blob/main/CITATION.bib):

```biblatex
@misc{SNWJERGMCountJL,
  author = {{Statistical Network Analysis with Julia}},
  title = {ERGMCount.jl: Exponential Random Graph Models for Count-Valued Networks in Julia},
  year = {2026},
  url = {https://github.com/statistical-network-analysis-with-Julia/ERGMCount.jl},
  note = {Homepage: https://statistical-network-analysis-with-Julia.github.io/ERGMCount.jl; GitHub: https://github.com/statistical-network-analysis-with-Julia}
}
```
