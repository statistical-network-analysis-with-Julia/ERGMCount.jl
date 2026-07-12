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

## Structural Terms

### CountMutualTerm

```@docs
CountMutualTerm
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
arbitrary count values.

### compute

```julia
compute(term::AbstractERGMTerm, net) -> Float64
```

Compute the full-network value of the term statistic. This is the shared
ERGM.jl term interface (its docstring lives in the ERGM.jl manual); every
count term in this package implements a method that reads the `:weight`
edge attribute (edges without the attribute count as 1).

### change_stat_count

```@docs
change_stat_count
```

### name

```julia
name(term::AbstractERGMTerm) -> String
```

Return the descriptive name of the term (e.g. `"sum"`, `"mutual.count"`),
used to label coefficients in fitted results. This is the shared ERGM.jl
term interface (its docstring lives in the ERGM.jl manual); every count
term implements a method.

### dyad_value

```@docs
dyad_value
```
