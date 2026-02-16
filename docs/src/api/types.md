# Types API Reference

This page documents the core data types in ERGMCount.jl.

## Reference Measures

### AbstractReferenceMeasure

```@docs
AbstractReferenceMeasure
```

### PoissonReference

```@docs
PoissonReference
```

### GeometricReference

```@docs
GeometricReference
```

### BinomialReference

```@docs
BinomialReference
```

### DiscUnifReference

```@docs
DiscUnifReference
```

### DiscUnif2Reference

```@docs
DiscUnif2Reference
```

## Model Types

### CountERGMModel

```@docs
CountERGMModel
```

### CountERGMResult

```@docs
CountERGMResult
```

## Reference Measure Functions

These functions operate on reference measures for computing log-probabilities and sampling.

### log_reference

```@docs
log_reference
```

### sample_reference

```@docs
sample_reference
```
