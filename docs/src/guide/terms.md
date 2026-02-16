# Count Terms

Count terms are sufficient statistics computed from the valued adjacency matrix. They generalize binary ERGM terms to handle integer-valued edges.

## Terms Interface

All count terms are subtypes of `AbstractERGMTerm` and implement:

```julia
compute(term, net) -> Float64       # Full network statistic
change_stat(term, net, i, j, ...) -> Float64  # Change when edge value changes
name(term) -> String                 # Human-readable name
```

## Basic Terms

### SumTerm

The most fundamental count term -- the sum of all edge values.

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
using Network, ERGMCount

net = Network{Int}(; n=5, directed=true)
add_edge!(net, 1, 2)
set_edge_attribute!(net, 1, 2, :weight, 3)
add_edge!(net, 2, 3)
set_edge_attribute!(net, 2, 3, :weight, 2)

term = SumTerm()
println(compute(term, net))  # 5.0 (3 + 2)
println(name(term))          # "sum"
```

### NonzeroTerm

Counts the number of non-zero edges (dyads with at least one event).

```julia
NonzeroTerm()
```

**Statistic**:

$$g(y) = \sum_{(i,j)} \mathbb{I}(y_{ij} > 0)$$

**Change statistic**: When edge $(i,j)$ changes from old value to new value:

$$\Delta g = \begin{cases} 1 & \text{if old} = 0 \text{ and new} > 0 \\ -1 & \text{if old} > 0 \text{ and new} = 0 \\ 0 & \text{otherwise} \end{cases}$$

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

## Structural Terms

### CountMutualTerm

Measures reciprocity in edge weights -- the tendency for dyads to have matched intensities.

```julia
CountMutualTerm()
```

**Statistic**:

$$g(y) = \sum_{i < j} \min(y_{ij}, y_{ji})$$

**Interpretation**:
- Generalizes binary mutuality to counts
- Captures the matched component of reciprocal ties
- Positive coefficient: dyads tend to reciprocate at similar intensities
- Only meaningful for directed networks

**Example**:

```julia
net = Network{Int}(; n=3, directed=true)

# Dyad (1,2): 3 from 1->2, 2 from 2->1
add_edge!(net, 1, 2)
set_edge_attribute!(net, 1, 2, :weight, 3)
add_edge!(net, 2, 1)
set_edge_attribute!(net, 2, 1, :weight, 2)

# Dyad (1,3): 1 from 1->3, 4 from 3->1
add_edge!(net, 1, 3)
set_edge_attribute!(net, 1, 3, :weight, 1)
add_edge!(net, 3, 1)
set_edge_attribute!(net, 3, 1, :weight, 4)

term = CountMutualTerm()
# min(3,2) + min(1,4) = 2 + 1 = 3
println(compute(term, net))  # 3.0
```

### TransitiveTiesTerm

Measures weighted transitivity in the network.

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

Measures weighted cyclicality in the network.

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
| Are there transitive groups? | `TransitiveTiesTerm` |
| Do some people send much more? | `NodeOSumTerm` |
| Do some people receive much more? | `NodeISumTerm` |

### Model Building Strategy

1. **Start with SumTerm + NonzeroTerm**: These control the baseline scale and density
2. **Add structural terms**: CountMutualTerm for reciprocity, TransitiveTiesTerm for clustering
3. **Add heterogeneity terms**: NodeOSumTerm and/or NodeISumTerm if needed
4. **Add threshold terms**: GreaterthannTerm for specific intensity effects

### Example: Full Model

```julia
terms = [
    # Baseline
    SumTerm(),
    NonzeroTerm(),

    # Structural
    CountMutualTerm(),
    TransitiveTiesTerm(),

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
net = Network{Int}(; n=10, directed=true)
# ... add edges with weights ...

for term in [SumTerm(), NonzeroTerm(), CountMutualTerm()]
    val = compute(term, net)
    println("$(name(term)) = $val")
end
```

## Comparison with Binary ERGM Terms

| Binary ERGM Term | Count ERGM Analogue | Key Difference |
|------------------|-----------------------|----------------|
| Edges | NonzeroTerm | Counts non-zero vs. all edges |
| (none) | SumTerm | New: total weight |
| Mutual | CountMutualTerm | min(y_ij, y_ji) vs. I(both present) |
| Triangle | TransitiveTiesTerm | min of weights vs. I(all present) |
| Out-degree | NodeOSumTerm | Squared out-strength vs. squared out-degree |
| In-degree | NodeISumTerm | Squared in-strength vs. squared in-degree |
