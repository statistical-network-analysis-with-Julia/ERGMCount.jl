# Count Terms

Count terms are sufficient statistics computed from the valued adjacency matrix. They generalize binary ERGM terms to handle integer-valued edges.

## Terms Interface

All count terms are subtypes of `AbstractERGMTerm` and implement:

<!-- skip-check -->
```julia
compute(term, net) -> Float64       # Full network statistic
change_stat_count(term, net, weights, i, j, old, new) -> Float64  # Change when dyad (i,j) moves from count `old` to `new`
name(term) -> String                 # Human-readable name
```

Both the estimator (one row of the compressed MPLE design per dyad) and the
Gibbs sampler (one full conditional per dyad) need a dyad's change statistic
at *every* value of its support at once. They read it through the support
profile `ERGMCount.change_stats_support!(dest, term, net, weights, i, j,
old, support)`, whose fallback calls `change_stat_count` once per value —
so a custom term that implements the three functions above works
everywhere. The strength terms (`NodeOSumTerm`, `NodeISumTerm`,
`NodeSumTerm`) and the triadic terms (`TransitiveTiesTerm`,
`CyclicalTiesTerm`, `TransitiveWeightsTerm`, `CyclicalWeightsTerm`)
specialise the profile to compute the part that does not depend on the
value once per dyad; the test suite holds every
specialisation equal to the per-value definition.

## Basic Terms

### SumTerm

The most fundamental count term -- the sum of all edge values.

<!-- skip-check -->
```julia
SumTerm()
```

**Statistic**:

$$g(y) = \sum_{(i,j)} y_{ij}$$

**Change statistic**: When edge $(i,j)$ changes by $\delta$:

$$\Delta g = \delta$$

**Interpretation**:
- Analogous to "intercept" for the count part of the model
- Controls the overall magnitude of edge values
- With Poisson reference, shifts the conditional mean
- Positive coefficient: higher edge values are more likely
- Negative coefficient: lower edge values are more likely

**Example**:

```julia
using Networks, ERGMCount
using ERGM: compute, name   # generic term interface shared with ERGM.jl

net = network(5; directed=true)
add_edge!(net, 1, 2)
set_edge_attribute!(net, :weight, 1, 2, 3)
add_edge!(net, 2, 3)
set_edge_attribute!(net, :weight, 2, 3, 2)

term = SumTerm()
println(compute(term, net))  # 5.0 (3 + 2)
println(name(term))          # "sum"
```

### NonzeroTerm

Counts the number of non-zero dyads (an edge whose `:weight` is 0 is a zero dyad, as in R; a negative count, admissible under `DiscUnif2Reference(a < 0, b)`, is non-zero).

```julia
NonzeroTerm()
```

**Statistic**:

$$g(y) = \sum_{(i,j)} \mathbb{I}(y_{ij} \neq 0)$$

**Change statistic**: When edge $(i,j)$ changes from old value to new value:

$$\Delta g = \begin{cases} 1 & \text{if old} = 0 \text{ and new} \neq 0 \\ -1 & \text{if old} \neq 0 \text{ and new} = 0 \\ 0 & \text{otherwise} \end{cases}$$

**Interpretation**:
- Analogous to the edges term in binary ERGMs
- Controls network density (proportion of non-zero dyads)
- Negative coefficient: sparse networks preferred (fewer ties)
- Positive coefficient: dense networks preferred (more ties)

**Example**:

```julia
term = NonzeroTerm()
println(compute(term, net))  # 2.0 (two non-zero edges)
println(name(term))          # "nonzero"
```

### GreaterthannTerm

Counts edges with values exceeding a threshold.

<!-- skip-check -->
```julia
GreaterthannTerm(threshold::Int)
```

**Statistic**:

$$g(y) = \sum_{(i,j)} \mathbb{I}(y_{ij} > n)$$

**Interpretation**:
- Captures the prevalence of high-intensity ties
- Useful for testing whether strong ties are over-represented
- Multiple thresholds can be included for a more flexible model

**Example**:

```julia
# Count edges with weight > 2
term = GreaterthannTerm(2)
println(name(term))          # "greaterthan.2"

# Use multiple thresholds
terms = [
    GreaterthannTerm(1),   # More than 1 interaction
    GreaterthannTerm(3),   # More than 3 interactions
    GreaterthannTerm(5),   # More than 5 interactions
]
```

### CountAtleastnTerm

Counts edges with values at or above a threshold.

<!-- skip-check -->
```julia
CountAtleastnTerm(threshold::Int)
```

**Statistic**:

$$g(y) = \sum_{(i,j)} \mathbb{I}(y_{ij} \geq n)$$

**Interpretation**:
- Similar to `GreaterthannTerm` but uses $\geq$ instead of $>$
- `CountAtleastnTerm(1)` is equivalent to `NonzeroTerm()`
- Useful for modeling different intensity thresholds

**Example**:

```julia
# Count edges with weight >= 3
term = CountAtleastnTerm(3)
println(name(term))          # "atleast.3"
```

### SmallerthanTerm

`ergm`'s `smallerthan(threshold)`: the number of dyads whose value is
**below** the threshold — zero-valued dyads included.

<!-- skip-check -->
```julia
SmallerthanTerm(threshold::Int)
```

**Statistic**:

$$g(y) = \sum_{(i,j)} \mathbb{I}(y_{ij} < n)$$

**Example**:

```julia
term = SmallerthanTerm(2)
println(name(term))          # "smallerthan.2"
```

### EqualToTerm

`ergm`'s `equalto(value, tolerance)`: the number of dyads whose value lies
within `tolerance` of `value`, inclusive — zero-valued dyads included.

<!-- skip-check -->
```julia
EqualToTerm(value::Int; tolerance::Int=0)
```

**Statistic**:

$$g(y) = \sum_{(i,j)} \mathbb{I}(\lvert y_{ij} - v \rvert \le t)$$

**Example**:

```julia
term = EqualToTerm(3)
println(name(term))          # "equalto.3.pm.0"
```

### InIntervalTerm

`ergm`'s `ininterval(lower, upper, open)`: the number of dyads whose value
lies between the bounds, each end open (exclusive, R's default) or closed.

<!-- skip-check -->
```julia
InIntervalTerm(lower, upper; open=(true, true))
```

**Statistic** (for the default open interval):

$$g(y) = \sum_{(i,j)} \mathbb{I}(a < y_{ij} < b)$$

**Example**:

```julia
println(name(InIntervalTerm(1, 3)))                        # "ininterval(1,3)"
println(name(InIntervalTerm(1, 3; open=(false, false))))   # "ininterval[1,3]"
println(name(InIntervalTerm(0, Inf; open=(false, true))))  # "ininterval[0,Inf)"
```

Together with `SumTerm`, `NonzeroTerm`, `GreaterthannTerm` and
`CountAtleastnTerm` these are the **dyad-independent** terms: a model made
of them only has an exact MLE (see [Estimation](estimation.md)).

## Structural Terms

### CountMutualTerm

Valued reciprocity — `ergm`'s `mutual(form=, threshold=)` for valued
networks.

<!-- skip-check -->
```julia
CountMutualTerm(form=:min; threshold=0)
```

**Statistic**: over the unordered pairs of a directed network,

$$g(y) = \sum_{i < j} m(y_{ij}, y_{ji})$$

| `form` | $m(a, b)$ | R label |
|---|---|---|
| `:min` (default) | $\min(a, b)$ | `mutual.min` |
| `:nabsdiff` | $-\lvert a - b \rvert$ | `mutual.nabsdiff` |
| `:geometric` | $\sqrt{ab}$ | `mutual.geom.mean` |
| `:product` | $ab$ | `mutual.product` |
| `:threshold` | $\mathbb{I}(a \ge t)\,\mathbb{I}(b \ge t)$ — binary mutuality after thresholding | `mutual.<t>` |

**Interpretation**:
- Generalizes binary mutuality to counts; `:min` captures the matched
  component of reciprocal ties, `:nabsdiff` penalises unmatched intensity
- Positive coefficient: dyads tend to reciprocate at similar intensities
- Directed networks only — an undirected `CountERGMModel` refuses it

**Example**:

```julia
net = network(3; directed=true)

# Dyad (1,2): 3 from 1->2, 2 from 2->1
add_edge!(net, 1, 2)
set_edge_attribute!(net, :weight, 1, 2, 3)
add_edge!(net, 2, 1)
set_edge_attribute!(net, :weight, 2, 1, 2)

# Dyad (1,3): 1 from 1->3, 4 from 3->1
add_edge!(net, 1, 3)
set_edge_attribute!(net, :weight, 1, 3, 1)
add_edge!(net, 3, 1)
set_edge_attribute!(net, :weight, 3, 1, 4)

term = CountMutualTerm()
# min(3,2) + min(1,4) = 2 + 1 = 3
println(compute(term, net))  # 3.0
println(name(term))          # "mutual.min"

println(compute(CountMutualTerm(:nabsdiff), net))              # -4.0 = -(1 + 3)
println(compute(CountMutualTerm(:product), net))               # 10.0 = 6 + 4
println(compute(CountMutualTerm(:threshold; threshold=2), net)) # 1.0 — only {1,2}
```

`:min`, `:nabsdiff`, `:geometric` and `:product` are pinned against `ergm`
4.12 by the `count_terms` fixture; `:threshold` follows the documented
definition but cannot be pinned, because `ergm` 4.12's own
`mutual(form="threshold")` fails at C model initialisation (the fixture
records the error).

### TransitiveWeightsTerm

`ergm`'s `transitiveweights("min", "max", "min")` (Krivitsky 2012, eq. 13):
each dyad's value capped by the strongest two-path that closes it.

```julia
TransitiveWeightsTerm()
```

**Statistic**:

$$g(y) = \sum_{(i,j)} \min\Big( y_{ij},\ \max_k \min(y_{ik}, y_{kj}) \Big)$$

over every ordered pair of a directed network and every unordered pair of
an undirected one (R's convention).

**Interpretation**:
- A valued transitivity: a strong `i → k → j` path makes a strong `i → j` tie
  more likely
- The default `(min, max, min)` triple is the "stable" one in R; the
  `geomean`/`sum` alternatives are not implemented
- **Negative counts are refused**, exactly as `ergm` refuses the term: on a
  network with a negative dyad weight (admissible only under
  `DiscUnif2Reference(a < 0, b)`) `compute`, `CountERGMModel` and
  `simulate_count_ergm` throw an `ArgumentError` — "`TransitiveWeightsTerm`
  (`transitiveweights.min.max.min`) may not be used with networks with
  negative dyad weights" — rather than define a statistic R never produces.
  The same holds for `CyclicalWeightsTerm` and for
  `CountMutualTerm(:geometric)` (R returns `NaN` there).

**Example**:

```julia
net = network(3; directed=true)
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 3)
add_edge!(net, 2, 3); set_edge_attribute!(net, :weight, 2, 3, 2)
add_edge!(net, 1, 3); set_edge_attribute!(net, :weight, 1, 3, 1)

term = TransitiveWeightsTerm()
println(compute(term, net))  # 1.0 — pair (1,3): min(y_13 = 1, min(3, 2) = 2)
println(name(term))          # "transitiveweights.min.max.min"
```

### CyclicalWeightsTerm

`ergm`'s `cyclicalweights("min", "max", "min")`: each dyad's value capped by
the strongest two-path that closes a directed 3-cycle through it.

```julia
CyclicalWeightsTerm()
```

**Statistic**:

$$g(y) = \sum_{(i,j)} \min\Big( y_{ij},\ \max_k \min(y_{jk}, y_{ki}) \Big)$$

On an undirected network the cycle and the transitive two-path coincide,
so the statistic equals `TransitiveWeightsTerm` there (R allows both).

**Example**:

```julia
net = network(3; directed=true)
add_edge!(net, 1, 2); set_edge_attribute!(net, :weight, 1, 2, 3)
add_edge!(net, 2, 3); set_edge_attribute!(net, :weight, 2, 3, 2)
add_edge!(net, 3, 1); set_edge_attribute!(net, :weight, 3, 1, 2)

term = CyclicalWeightsTerm()
println(compute(term, net))  # 6.0 — each dyad of the cycle capped at 2
println(name(term))          # "cyclicalweights.min.max.min"
```

### TransitiveTiesTerm

A simple valued transitivity over ordered triples. **Not an R term**: it is
neither `ergm`'s `transitiveweights` (use `TransitiveWeightsTerm`) nor the
binary `transitiveties`, and it is not validated against R.

```julia
TransitiveTiesTerm()
```

**Statistic**:

$$g(y) = \sum_{i,j,k} \min(y_{ij}, y_{jk}, y_{ik})$$

Where the sum is over distinct triples with $i \neq j \neq k$.

**Interpretation**:
- Generalizes binary transitivity to counts
- Captures the tendency for transitive closure to carry weight
- A positive coefficient indicates that when $i \to j$ and $j \to k$ are strong, $i \to k$ also tends to be strong

**Example**:

```julia
term = TransitiveTiesTerm()
println(name(term))  # "transitiveties.count"
```

### CyclicalTiesTerm

A simple valued cyclicality over directed 3-cycles. **Not an R term**: it is
neither `ergm`'s `cyclicalweights` (use `CyclicalWeightsTerm`) nor the
binary `cyclicalties`, and it is not validated against R.

```julia
CyclicalTiesTerm()
```

**Statistic**:

$$g(y) = \frac{1}{3}\sum_{i,j,k} \min(y_{ij}, y_{jk}, y_{ki})$$

**Interpretation**:
- Generalizes binary cyclical closure to counts
- Captures the tendency for directed cycles to carry weight
- The $1/3$ factor corrects for each cycle being counted three times
- Only meaningful for directed networks

**Example**:

```julia
term = CyclicalTiesTerm()
println(name(term))  # "cyclicalties.count"
```

## Degree / Strength Terms

These terms capture heterogeneity in node-level activity and popularity, measured by edge weight sums (strength) rather than edge counts (degree).

### NodeOSumTerm

Out-strength heterogeneity: measures whether some actors send more total weight than others.

```julia
NodeOSumTerm()
```

**Statistic**:

$$g(y) = \sum_i \left(\sum_j y_{ij}\right)^2$$

The sum of squared out-strengths.

**Interpretation**:
- Positive coefficient: heterogeneous sending activity (some actors send much more)
- Analogous to out-degree heterogeneity in binary ERGMs
- Captures "high-volume senders"

**Example**:

```julia
term = NodeOSumTerm()
println(name(term))  # "nodeOSum"
```

### NodeISumTerm

In-strength heterogeneity: measures whether some actors receive more total weight.

```julia
NodeISumTerm()
```

**Statistic**:

$$g(y) = \sum_j \left(\sum_i y_{ij}\right)^2$$

The sum of squared in-strengths.

**Interpretation**:
- Positive coefficient: heterogeneous receiving (some actors receive much more)
- Captures "popular" or "attractive" targets
- Analogous to in-degree heterogeneity in binary ERGMs

**Example**:

```julia
term = NodeISumTerm()
println(name(term))  # "nodeISum"
```

### NodeSumTerm

Total strength heterogeneity: combines in-strength and out-strength.

```julia
NodeSumTerm()
```

**Statistic**:

$$g(y) = \sum_i \left(\sum_j y_{ij} + \sum_j y_{ji}\right)^2$$

The sum of squared total strengths.

**Interpretation**:
- Captures overall centrality heterogeneity
- Useful for undirected networks or when in/out distinction is not needed
- Positive coefficient: some nodes are much more active overall

**Example**:

```julia
term = NodeSumTerm()
println(name(term))  # "nodeSum"
```

## Choosing Terms

### By Research Question

| Question | Recommended Terms |
|----------|-------------------|
| What is the overall interaction level? | `SumTerm` |
| How dense is the network? | `NonzeroTerm` |
| Are strong ties over-represented? | `GreaterthannTerm(k)` |
| Is communication reciprocal? | `CountMutualTerm` |
| Are there transitive groups? | `TransitiveWeightsTerm` |
| How many dyads sit in an intensity band? | `InIntervalTerm(a, b)`, `EqualToTerm(v)`, `SmallerthanTerm(k)` |
| Do some people send much more? | `NodeOSumTerm` |
| Do some people receive much more? | `NodeISumTerm` |

### Model Building Strategy

1. **Start with SumTerm + NonzeroTerm**: These control the baseline scale and density
2. **Add structural terms**: CountMutualTerm for reciprocity, TransitiveWeightsTerm for clustering
3. **Add heterogeneity terms**: NodeOSumTerm and/or NodeISumTerm if needed
4. **Add threshold terms**: GreaterthannTerm for specific intensity effects

### Example: Full Model

A model this rich needs a network with enough dyads to identify it (the
3-node example above has six); build one with counts first:

```julia
using Random
rng = Xoshiro(7)
net = network(15; directed=true)
for i in 1:15, j in 1:15
    if i != j && rand(rng) < 0.3
        add_edge!(net, i, j)
        set_edge_attribute!(net, :weight, i, j, rand(rng, 1:6))
    end
end

terms = [
    # Baseline
    SumTerm(),
    NonzeroTerm(),

    # Structural
    CountMutualTerm(),
    TransitiveWeightsTerm(),

    # Heterogeneity
    NodeOSumTerm(),
    NodeISumTerm(),

    # Intensity thresholds
    GreaterthannTerm(3),
    GreaterthannTerm(5),
]

result = ergm_count(net, terms; reference=PoissonReference())
```

## Computing Terms Manually

You can compute any term on a network without fitting a model:

```julia
net = network(10; directed=true)
# ... add edges with weights ...

for term in [SumTerm(), NonzeroTerm(), CountMutualTerm()]
    val = compute(term, net)
    println("$(name(term)) = $val")
end
```

## Correspondence with R

Every term below is labelled with the string R's `summary()`/`coef()`
prints, so a by-name comparison with a statnet fit works. "Pinned" means
the statistic *and* the name are frozen from `ergm` 4.12 in
`test/fixtures/count_terms.toml` (generated by
`test/fixtures/r/count_terms.R`) and compared at 1e-9 on `zach` and on a
directed count network.

| ERGMCount term | R term (`ergm`/`ergm.count`) | Label | Exact counterpart, pinned? | Binary analogue |
|---|---|---|---|---|
| `SumTerm()` | `sum` | `sum` | yes | (none) |
| `NonzeroTerm()` | `nonzero` | `nonzero` | yes | `edges` |
| `GreaterthannTerm(k)` | `greaterthan(k)` | `greaterthan.k` | yes | — |
| `CountAtleastnTerm(k)` | `atleast(k)` | `atleast.k` | yes | — |
| `SmallerthanTerm(k)` | `smallerthan(k)` | `smallerthan.k` | yes | — |
| `EqualToTerm(v; tolerance=t)` | `equalto(v, t)` | `equalto.v.pm.t` | yes | — |
| `InIntervalTerm(a, b; open)` | `ininterval(a, b, open)` | `ininterval(a,b)`, `[a,b]`, ... | yes (all four bracket forms) | — |
| `CountMutualTerm()` | `mutual(form="min")` | `mutual.min` | yes | `mutual` |
| `CountMutualTerm(:nabsdiff)` | `mutual(form="nabsdiff")` | `mutual.nabsdiff` | yes | — |
| `CountMutualTerm(:geometric)` | `mutual(form="geometric")` | `mutual.geom.mean` | yes | — |
| `CountMutualTerm(:product)` | `mutual(form="product")` | `mutual.product` | yes | — |
| `CountMutualTerm(:threshold; threshold=t)` | `mutual(form="threshold", threshold=t)` | `mutual.t` | **no** — documented definition; `ergm` 4.12 errors on its own term (recorded in the fixture) | `mutual` |
| `TransitiveWeightsTerm()` | `transitiveweights("min","max","min")` | `transitiveweights.min.max.min` | yes (default triple only) | `ttriple`-like |
| `CyclicalWeightsTerm()` | `cyclicalweights("min","max","min")` | `cyclicalweights.min.max.min` | yes (default triple only) | `ctriple`-like |
| `TransitiveTiesTerm()` | **none** | `transitiveties.count` | no R counterpart; not validated against R | — |
| `CyclicalTiesTerm()` | **none** | `cyclicalties.count` | no R counterpart; not validated against R | — |
| `NodeOSumTerm()` / `NodeISumTerm()` / `NodeSumTerm()` | **none** (`ergm`'s `nodeocovar`/`nodeicovar`/`nodecovar` are covariance statistics, not squared strengths) | `nodeOSum` / `nodeISum` / `nodeSum` | no R counterpart; tested by hand value and brute force | `ostar(2)`/`istar(2)`-like |

Not implemented (there is no term to call): the non-default
`transitiveweights`/`cyclicalweights` triples (`twopath="geomean"`,
`combine="sum"`, `affect="geomean"`), the
`nodecovar`/`nodeocovar`/`nodeicovar`/`nodesqrtcovar` family, and the valued
forms of `nodematch`/`nodefactor`/`absdiff`/`edgecov`.
