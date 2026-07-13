# Estimation API Reference

```@meta
CurrentModule = ERGMCount
```

This page documents the functions for fitting, simulating, and assessing count ERGMs.

## Model Fitting

### fit_ergm_count

```@docs
fit_ergm_count
```

### ergm_count

```@docs
ergm_count
```

### fit_count_ergm

```@docs
fit_count_ergm
```

### count_mple

```@docs
count_mple
```

## Simulation

### simulate_count_ergm

```@docs
simulate_count_ergm
```

## Goodness of Fit

### gof

Goodness-of-fit assessment is provided as a method of the shared
`Networks.gof` generic, so the same `gof(result)` call works across the
model packages of the ecosystem.

```@docs
gof(::CountERGMResult)
```
