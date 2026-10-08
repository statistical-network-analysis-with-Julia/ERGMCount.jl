# Golden fixture: statnet `ergm`'s VALUED covariate terms (form = "sum" and
# form = "nonzero"), their labels, and one ergm.count fit that uses them.
#
# Regenerate from the package root (about a minute):
#
#   Rscript test/fixtures/r/count_covariates.R > test/fixtures/count_covariates.toml
#
# WHAT THIS FIXTURE PINS
#
# (a) Summary statistics and the labels R prints for nodematch (uniform, per
#     level, nonzero), nodefactor (each level but the first, sum and
#     nonzero), absdiff, nodecov, nodeocov, nodeicov and edgecov (a network
#     attribute "dist"), on a seeded 12-actor DIRECTED count network and a
#     seeded 12-actor UNDIRECTED one, both with a categorical vertex attribute
#     "g" (levels a, b, c) and a numeric one "x" (all frozen below). These
#     are deterministic functions of the network: 1e-9.
# (b) The dyad-independent model  ~ sum + nonzero + nodematch("g") +
#     absdiff("x") + edgecov("dist")  on the directed network, reference =
#     ~Poisson, fit by ergm.count (MCMLE, its only estimator) under five
#     seeds. ERGMCount.jl's count MPLE of this model is its exact MLE, so a
#     Julia fit must lie within R's Monte-Carlo resolution of R's mean (the
#     band of [tolerance]).

suppressMessages({
  .libPaths(c(path.expand("~/R/library"), .libPaths()))
  library(ergm.count)
  library(parallel)
})

seed <- 20261008
set.seed(seed)

num <- function(x) paste(sprintf("%.17g", x), collapse = ", ")
ints <- function(x) paste(sprintf("%d", as.integer(x)), collapse = ", ")
strs <- function(x) paste(sprintf('"%s"', x), collapse = ", ")

n <- 12L
g <- sample(c("a", "b", "c"), n, replace = TRUE)
x <- round(runif(n, 0, 4), 2)

# Directed
Md <- matrix(rpois(n * n, 1), n, n); diag(Md) <- 0L
Wd <- matrix(round(runif(n * n), 2), n, n); diag(Wd) <- 0
dnw <- network(Md, directed = TRUE, matrix.type = "adjacency",
               ignore.eval = FALSE, names.eval = "w")
dnw %v% "g" <- g
dnw %v% "x" <- x
dnw %n% "dist" <- Wd
di <- which(Md > 0, arr.ind = TRUE)

# Undirected (symmetric counts and covariate)
Mu <- matrix(rpois(n * n, 1), n, n); Mu[lower.tri(Mu, diag = TRUE)] <- 0L
Mu <- Mu + t(Mu)
Wu <- matrix(round(runif(n * n), 2), n, n); Wu[lower.tri(Wu, diag = TRUE)] <- 0
Wu <- Wu + t(Wu)
unw <- network(Mu, directed = FALSE, matrix.type = "adjacency",
               ignore.eval = FALSE, names.eval = "w")
unw %v% "g" <- g
unw %v% "x" <- x
unw %n% "dist" <- Wu
ui <- which(Mu > 0 & upper.tri(Mu), arr.ind = TRUE)

f_d <- dnw ~ nodematch("g") + nodematch("g", diff = TRUE) +
  nodematch("g", form = "nonzero") + nodefactor("g") +
  nodefactor("g", form = "nonzero") + absdiff("x") +
  absdiff("x", form = "nonzero") + nodecov("x") + nodecov("x", form = "nonzero") +
  nodeocov("x") + nodeicov("x") + edgecov("dist") + edgecov("dist", form = "nonzero")
f_u <- unw ~ nodematch("g") + nodematch("g", diff = TRUE) + nodefactor("g") +
  absdiff("x") + nodecov("x") + edgecov("dist") + edgecov("dist", form = "nonzero")
sd_ <- summary(f_d, response = "w")
su <- summary(f_u, response = "w")

# (b) the fit
f_fit <- dnw ~ sum + nonzero + nodematch("g") + absdiff("x") + edgecov("dist")
seeds <- 1:5
fit_one <- function(s) {
  set.seed(s)
  fit <- NULL
  invisible(capture.output(
    fit <- suppressWarnings(suppressMessages(
      ergm(f_fit, response = "w", reference = ~Poisson,
           control = control.ergm(seed = s)))), type = "output"))
  list(coef = as.numeric(coef(fit)), se = as.numeric(sqrt(diag(vcov(fit)))),
       names = names(coef(fit)), failed = isTRUE(fit$failure))
}
fits <- mclapply(seeds, fit_one, mc.cores = 3)
stopifnot(!any(sapply(fits, `[[`, "failed")))
coefs <- t(sapply(fits, `[[`, "coef"))
ses <- t(sapply(fits, `[[`, "se"))

cat('name = "count_covariates"\n\n')
cat("[provenance]\n")
cat(sprintf('r_version = "%s"\n', as.character(getRversion())))
cat(sprintf('ergm_count_version = "%s"\n', as.character(packageVersion("ergm.count"))))
cat(sprintf('ergm_version = "%s"\n', as.character(packageVersion("ergm"))))
cat(sprintf('network_version = "%s"\n', as.character(packageVersion("network"))))
cat(sprintf("seed = %d\n", seed))
cat('script = "test/fixtures/r/count_covariates.R"\n')
cat(sprintf('date = "%s"\n', format(Sys.Date())))
cat('models = "(a) summary() of the valued nodematch/nodefactor/absdiff/nodecov/nodeocov/nodeicov/edgecov terms (form sum and nonzero) on a seeded 12-actor directed and a seeded 12-actor undirected count network, response=w; (b) ~ sum + nonzero + nodematch(g) + absdiff(x) + edgecov(dist) on the directed one, reference=~Poisson, ergm.count MCMLE at control.ergm defaults, seeds 1..5"\n\n')

cat("[tolerance]\n")
cat("# Summary statistics are deterministic functions of the network.\n")
cat("summary_statistics = 1e-9\n")
cat("# (b): the exact MLE against R's five-seed Monte-Carlo mean, within\n")
cat("# band_sd * sd_k * sqrt(1 + 1/n_seeds), floored at floor_se_fraction of\n")
cat("# R's mean standard error (the band of count_mcmle.toml); standard errors\n")
cat("# likewise against R's spread of standard errors.\n")
cat("band_sd = 4.0\n")
cat("floor_se_fraction = 0.1\n\n")

cat("[values]\n")
cat(sprintf("n = %d\n", n))
cat(sprintf("g = [%s]\n", strs(g)))
cat(sprintf("x = [%s]\n", num(x)))
cat(sprintf("directed_edge_src = [%s]\n", ints(di[, 1])))
cat(sprintf("directed_edge_dst = [%s]\n", ints(di[, 2])))
cat(sprintf("directed_edge_weight = [%s]\n", ints(Md[di])))
cat(sprintf("directed_dist = [%s]\n", num(as.numeric(t(Wd)))))
cat(sprintf("directed_names = [%s]\n", strs(names(sd_))))
cat(sprintf("directed_summary = [%s]\n", num(as.numeric(sd_))))
cat(sprintf("undirected_edge_src = [%s]\n", ints(ui[, 1])))
cat(sprintf("undirected_edge_dst = [%s]\n", ints(ui[, 2])))
cat(sprintf("undirected_edge_weight = [%s]\n", ints(Mu[ui])))
cat(sprintf("undirected_dist = [%s]\n", num(as.numeric(t(Wu)))))
cat(sprintf("undirected_names = [%s]\n", strs(names(su))))
cat(sprintf("undirected_summary = [%s]\n", num(as.numeric(su))))
cat("# (b); `*_dist` above are row-major (row i is W[i, ])\n")
cat(sprintf("fit_names = [%s]\n", strs(fits[[1]]$names)))
cat(sprintf("fit_seeds = [%s]\n", ints(seeds)))
cat(sprintf("fit_coefficients_mean = [%s]\n", num(colMeans(coefs))))
cat(sprintf("fit_coefficients_sd = [%s]\n", num(apply(coefs, 2, sd))))
cat(sprintf("fit_std_errors_mean = [%s]\n", num(colMeans(ses))))
cat(sprintf("fit_std_errors_sd = [%s]\n", num(apply(ses, 2, sd))))
