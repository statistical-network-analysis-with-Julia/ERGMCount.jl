# Terms API Reference

```@meta
CurrentModule = ERGMCount
```

This page documents all count-specific ERGM terms available in ERGMCount.jl.

## Basic Terms

### SumTerm

```@docs
SumTerm
```

### NonzeroTerm

```@docs
NonzeroTerm
```

### GreaterthannTerm

```@docs
GreaterthannTerm
```

### CountAtleastnTerm

```@docs
CountAtleastnTerm
```

### SmallerthanTerm

```@docs
SmallerthanTerm
```

### EqualToTerm

```@docs
EqualToTerm
```

### InIntervalTerm

```@docs
InIntervalTerm
```

## Structural Terms

### CountMutualTerm

```@docs
CountMutualTerm
```

### TransitiveWeightsTerm

```@docs
TransitiveWeightsTerm
```

### CyclicalWeightsTerm

```@docs
CyclicalWeightsTerm
```

### TransitiveTiesTerm

```@docs
TransitiveTiesTerm
```

### CyclicalTiesTerm

```@docs
CyclicalTiesTerm
```

## Degree / Strength Terms

### NodeOSumTerm

```@docs
NodeOSumTerm
```

### NodeISumTerm

```@docs
NodeISumTerm
```

### NodeSumTerm

```@docs
NodeSumTerm
```

## Interface Functions

Count terms implement the shared ERGM.jl term interface (`compute`,
`name`) plus the count-specific change statistic `change_stat_count`,
which replaces the binary `change_stat` for dyads that move between
arbitrary count values. The estimator and the Gibbs sampler consume the
change statistics one dyad at a time over the whole count support, through
the *support profile* `change_stats_support!`; its fallback calls
`change_stat_count` per value, so a custom term needs only the latter.

### compute

<!-- skip-check -->
```julia
compute(term::AbstractERGMTerm, net) -> Float64
```

Compute the full-network value of the term statistic. This is the shared
ERGM.jl term interface (its docstring lives in the ERGM.jl manual); every
count term in this package implements a method that reads the `:weight`
edge attribute through the same dyad-value view the estimator uses: an edge
whose `:weight` is 0 is a zero
dyad (so `nonzero` does not count it, exactly as R's valued terms), and the
threshold terms (`greaterthan`, `atleast`, `smallerthan`, `equalto`,
`ininterval`) count the zero-valued dyads whenever 0 meets the threshold —
`GreaterthannTerm(-1)` and `CountAtleastnTerm(0)` are the number of dyads.
A `CountERGMModel` refuses a network with edges but no `:weight` at all, and
one where only some edges carry a `:weight` (a bare edge is a data gap, not a
count of 1); `TransitiveWeightsTerm`, `CyclicalWeightsTerm` and
`CountMutualTerm(:geometric)` refuse a network with a negative count, as
`ergm` does.

### change_stat_count

```@docs
change_stat_count
```

### change_stats_support!

```@docs
change_stats_support!
```

### name

<!-- skip-check -->
```julia
name(term::AbstractERGMTerm) -> String
```

Return the term's coefficient label, the string R prints for the same term (e.g. `"sum"`, `"mutual.min"`, `"transitiveweights.min.max.min"`),
used to label coefficients in fitted results. This is the shared ERGM.jl
term interface (its docstring lives in the ERGM.jl manual); every count
term implements a method.

### dyad_value

```@docs
dyad_value
```
