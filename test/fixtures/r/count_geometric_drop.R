# Golden fixture: statnet `ergm.count` MCMLE on two models at the edges of the
# count sample space -- a proper geometric-reference model with strong
# reciprocity, and a statistic at the boundary of its attainable range.
#
# Regenerate from the package root (a few minutes: 2 models x 9 default fits,
# on 3 cores):
#
#   Rscript test/fixtures/r/count_geometric_drop.R > test/fixtures/count_geometric_drop.toml
#
# WHAT THIS FIXTURE PINS
#
#   (c) a seeded 20-actor DIRECTED count network with strong reciprocity
#       (frozen as an edge list below) ~ sum + mutual(form = "min"),
#       reference = ~Geometric. At its MLE theta_sum + theta_mutual > 0 >
#       2 theta_sum + theta_mutual: a dyad whose reciprocal sits at the top of
#       a truncated support is pulled to the top too, so the all-top network
#       is a coordinate-wise mode of the truncated joint, yet the weight of
#       a reciprocated pair at counts (Y, Y) falls like
#       exp((2 theta_sum + theta_mutual) Y): the model is proper, R fits it,
#       and R's simulations at its fit reproduce the observed mean count.
#       The script checks both inequalities at R's mean and freezes the
#       simulated mean. ERGMCount.jl must fit it (it used to refuse it as
#       having "a mode on the truncation bound").
#   (d) a seeded 16-actor UNDIRECTED count network that is a perfect matching
#       (frozen as an edge list) ~ sum + nonzero + transitiveweights("min",
#       "max", "min"), reference = ~Poisson. No two ties share an actor, so
#       transitiveweights is 0 -- its smallest attainable value (R's minval)
#       -- and its change statistic is 0 on every dyad. R fixes its
#       coefficient at -Inf (`drop = TRUE`, "Observed statistic(s)
#       transitiveweights.min.max.min are at their smallest attainable
#       values") and estimates sum and nonzero with the statistic held at 0:
#       the MLE of the model restricted to the networks without a weighted
#       two-path closing a tie. ERGMCount.jl's `method=:mcmle` must do the
#       same (it used to refuse the model).
#
# Both sides are Monte-Carlo MLEs; the frozen nine-seed mean and spread give
# the band of [tolerance], the same rule as count_mcmle.toml.

suppressMessages({
  .libPaths(c(path.expand("~/R/library"), .libPaths()))
  library(ergm.count)
  library(parallel)
})

seed <- 20261007
set.seed(seed)

num <- function(x) paste(sprintf("%.17g", x), collapse = ", ")
ints <- function(x) paste(sprintf("%d", as.integer(x)), collapse = ", ")
strs <- function(x) paste(sprintf('"%s"', x), collapse = ", ")

# --- (c) 20 actors, directed, strong reciprocity -------------------------------
# Each unordered pair is reciprocated (one count on both arcs) with
# probability 0.6, and each arc gets two further Bernoulli(0.05) counts
nc <- 20L
Mc <- matrix(0L, nc, nc)
for (i in 1:(nc - 1)) for (j in (i + 1):nc) {
  s <- rbinom(1, 1, 0.6)
  Mc[i, j] <- s + rbinom(1, 2, 0.05)
  Mc[j, i] <- s + rbinom(1, 2, 0.05)
}
cnw <- network(Mc, directed = TRUE, matrix.type = "adjacency",
               ignore.eval = FALSE, names.eval = "w")
ci <- which(Mc > 0, arr.ind = TRUE)

# --- (d) 16 actors, undirected perfect matching with Poisson+1 counts ----------
nd <- 16L
Md <- matrix(0L, nd, nd)
for (k in seq(1, nd, by = 2)) {
  v <- rpois(1, 1.5) + 1L
  Md[k, k + 1] <- v
  Md[k + 1, k] <- v
}
dnw <- network(Md, directed = FALSE, matrix.type = "adjacency",
               ignore.eval = FALSE, names.eval = "w")
di <- which(Md > 0 & upper.tri(Md), arr.ind = TRUE)

f_c <- cnw ~ sum + mutual(form = "min")
f_d <- dnw ~ sum + nonzero + transitiveweights("min", "max", "min")

fit_one <- function(f, ref, s) {
  set.seed(s)
  fit <- NULL
  warned <- character(0)
  msgs <- character(0)
  invisible(capture.output(
    fit <- withCallingHandlers(
      ergm(f, response = "w", reference = ref, control = control.ergm(seed = s)),
      warning = function(w) { warned <<- c(warned, conditionMessage(w))
                              invokeRestart("muffleWarning") },
      message = function(m) { msgs <<- c(msgs, conditionMessage(m))
                              invokeRestart("muffleMessage") }),
    type = "output"))
  list(fit = fit, coef = as.numeric(coef(fit)),
       se = as.numeric(sqrt(diag(vcov(fit)))), names = names(coef(fit)),
       warned = warned, msgs = msgs, failed = isTRUE(fit$failure))
}

seeds <- 1:9
run_model <- function(f, ref) {
  fits <- mclapply(seeds, function(s) fit_one(f, ref, s), mc.cores = 3)
  stopifnot(!any(sapply(fits, `[[`, "failed")))
  coefs <- t(sapply(fits, `[[`, "coef"))
  ses <- t(sapply(fits, `[[`, "se"))
  # A coefficient fixed at -Inf has the same value in every seed: its sd is 0
  fin <- apply(coefs, 2, function(x) all(is.finite(x)))
  sdv <- rep(0, ncol(coefs)); sdv[fin] <- apply(coefs[, fin, drop = FALSE], 2, sd)
  sev <- rep(0, ncol(ses)); sev[fin] <- apply(ses[, fin, drop = FALSE], 2, sd)
  list(fits = fits, names = fits[[1]]$names, coefs = coefs, ses = ses,
       coef_mean = colMeans(coefs), coef_sd = sdv,
       se_mean = colMeans(ses), se_sd = sev)
}

rc <- run_model(f_c, ~Geometric)
rd <- run_model(f_d, ~Poisson)

# (c): the regime that matters, at R's mean
th <- rc$coef_mean
stopifnot(th[1] + th[2] > 0, 2 * th[1] + th[2] < 0)
# ... and R's simulations at its own seed-1 fit reproduce the observed mean
sim <- simulate(rc$fits[[1]]$fit, nsim = 100, output = "stats", seed = 1)
n_arcs <- nc * (nc - 1)
sim_mean <- mean(sim[, 1]) / n_arcs
obs_mean <- sum(Mc) / n_arcs
stopifnot(abs(sim_mean - obs_mean) < 0.1)

# (d): R drops transitiveweights in every seed, and says so
stopifnot(all(rd$coefs[, 3] == -Inf), all(is.finite(rd$coefs[, 1:2])))
dropped_msg <- any(grepl("smallest attainable", unlist(lapply(rd$fits, `[[`, "msgs"))))
stopifnot(dropped_msg)

sc <- summary(f_c, response = "w")
sd_ <- summary(f_d, response = "w")

emit <- function(tag, r, s) {
  cat(sprintf("%s_term_names = [%s]\n", tag, strs(r$names)))
  cat(sprintf("%s_summary_statistics = [%s]\n", tag, num(as.numeric(s))))
  cat(sprintf("%s_seeds = [%s]\n", tag, ints(seeds)))
  for (k in seq_along(seeds)) {
    cat(sprintf("%s_coefficients_seed%d = [%s]\n", tag, seeds[k],
                gsub("-Inf", "-inf", num(r$coefs[k, ]))))
    cat(sprintf("%s_std_errors_seed%d = [%s]\n", tag, seeds[k], num(r$ses[k, ])))
  }
  cat(sprintf("%s_coefficients_mean = [%s]\n", tag, gsub("-Inf", "-inf", num(r$coef_mean))))
  cat(sprintf("%s_coefficients_sd = [%s]\n", tag, num(r$coef_sd)))
  cat(sprintf("%s_std_errors_mean = [%s]\n", tag, num(r$se_mean)))
  cat(sprintf("%s_std_errors_sd = [%s]\n", tag, num(r$se_sd)))
}

cat('name = "count_geometric_drop"\n\n')

cat("[provenance]\n")
cat(sprintf('r_version = "%s"\n', as.character(getRversion())))
cat(sprintf('ergm_count_version = "%s"\n', as.character(packageVersion("ergm.count"))))
cat(sprintf('ergm_version = "%s"\n', as.character(packageVersion("ergm"))))
cat(sprintf('network_version = "%s"\n', as.character(packageVersion("network"))))
cat(sprintf("seed = %d\n", seed))
cat('script = "test/fixtures/r/count_geometric_drop.R"\n')
cat(sprintf('date = "%s"\n', format(Sys.Date())))
cat('models = "(c) a seeded 20-actor directed count network with strong reciprocity ~ sum + mutual(form=min), response=w, reference=~Geometric; (d) a seeded 16-actor undirected perfect matching ~ sum + nonzero + transitiveweights(min,max,min), response=w, reference=~Poisson"\n')
cat('estimator = "ergm.count MCMLE at control.ergm defaults, nine seeds (1..9) per model"\n')
cat("\n")

cat("[tolerance]\n")
cat("# Observed statistics are deterministic functions of the graph.\n")
cat("summary_statistics = 1e-9\n")
cat("# The band of count_mcmle.toml: a Julia fit within\n")
cat("# band_sd * sqrt(sd_k^2 + sd_k^2 / n_seeds) of R's nine-seed mean, floored\n")
cat("# at floor_se_fraction of R's mean standard error; standard errors likewise.\n")
cat("band_sd = 4.0\n")
cat("floor_se_fraction = 0.1\n")
cat("# (c) R's simulated mean count per arc at its seed-1 fit (100 networks)\n")
cat("# against the observed mean, |difference| < this:\n")
cat("sim_mean_tolerance = 0.1\n")
cat("\n")

cat("[values]\n")
cat("# --- (c) the seeded directed network, frozen as a valued edge list ------\n")
cat(sprintf("geometric_n = %d\n", nc))
cat(sprintf("geometric_edge_src = [%s]\n", ints(ci[, 1])))
cat(sprintf("geometric_edge_dst = [%s]\n", ints(ci[, 2])))
cat(sprintf("geometric_edge_weight = [%s]\n", ints(Mc[ci])))
emit("geometric", rc, sc)
cat(sprintf("geometric_r_sim_mean = %.17g\n", sim_mean))
cat(sprintf("geometric_observed_mean = %.17g\n", obs_mean))
cat("\n# --- (d) the seeded matching, frozen as a valued edge list --------------\n")
cat(sprintf("drop_n = %d\n", nd))
cat(sprintf("drop_edge_src = [%s]\n", ints(di[, 1])))
cat(sprintf("drop_edge_dst = [%s]\n", ints(di[, 2])))
cat(sprintf("drop_edge_weight = [%s]\n", ints(Md[di])))
emit("drop", rd, sd_)
