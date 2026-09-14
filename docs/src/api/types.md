# Types API Reference

```@meta
CurrentModule = ERGMCount
```

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

### Reference Measure Traits

A reference measure with unbounded support (Poisson, geometric) must be
truncated to be normalizable, so the fitted model is exact only if the mass
beyond the truncation point is negligible. `is_truncating` declares which
measures are affected and `BOUNDARY_MASS_TOL` is the tolerance above which the
boundary mass is reported as an approximation.

```@docs
is_truncating
BOUNDARY_MASS_TOL
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
