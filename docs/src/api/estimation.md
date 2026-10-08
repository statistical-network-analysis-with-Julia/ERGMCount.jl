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

### count_mple

```@docs
count_mple
```

### count_mcmle

```@docs
count_mcmle
```

## StatsAPI surface

Every fit answers the ecosystem's StatsAPI surface (`coef`, `stderror`,
`vcov`, `confint`, `loglikelihood`, `nobs`, `dof`, `aic`, `bic`, `coeftable`,
`coefnames`);
the methods below carry ERGMCount-specific semantics.

```@docs
aic(::CountERGMResult)
bic(::CountERGMResult)
confint(::CountERGMResult)
coeftable(::CountERGMResult)
coefnames(::CountERGMResult)
```

## Model classification

```@docs
has_dyad_dependent(::CountERGMModel)
is_exact(::CountERGMResult)
se_method(::CountERGMResult)
```

## Simulation

### simulate_count_ergm

```@docs
simulate_count_ergm
```

## Goodness of Fit

### gof

Goodness-of-fit assessment is provided as a method of the shared
`NetworkCore.gof` generic, so the same `gof(result)` call works across the
model packages of the ecosystem.

```@docs
gof(::CountERGMResult)
```

## Renamed and removed names

Names that existed during development but were never released. They are not
kept as deprecated aliases; calling one is an `UndefVarError` (a removed
keyword is an `ArgumentError` naming the estimator that takes it).

| Development name | Use instead |
|:--|:--|
| `fit_count_ergm(net, terms; ...)` | `fit_ergm_count(net, terms; ...)` or its R-named alias `ergm_count` |
| `method=:mple` as the implicit default for every formula | `method=:auto` (the default) follows R: the MPLE for a dyad-independent formula, the MCMLE otherwise; pass `method=:mple` for the count MPLE of a dyad-dependent one |
