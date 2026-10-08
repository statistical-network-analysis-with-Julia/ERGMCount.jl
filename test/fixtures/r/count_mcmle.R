# Golden fixture: statnet `ergm.count` MCMLE of DYAD-DEPENDENT count ERGMs.
#
# Regenerate from the package root (about ten minutes: 2 models x (9 default
# fits + 1 long-chain fit); the fits run on 3 cores):
#
#   Rscript test/fixtures/r/count_mcmle.R > test/fixtures/count_mcmle.toml
#
# WHAT THIS FIXTURE PINS
#
# ERGMCount.jl's `method=:mcmle` -- Monte-Carlo maximum likelihood on the exact
# Gibbs sampler -- against ergm.count's own estimator, which is MCMLE (statnet
# has no valued MPLE). Both sides carry Monte-Carlo error, so the comparison is
# made at the resolution R itself has: every model is refit under nine seeds
# and the seed-to-seed mean and standard deviation of every coefficient and
# standard error are frozen. The Julia testset compares its own fits with R's
# mean inside a band derived from R's spread (see [tolerance]).
#
#   (a) zach ~ sum + nonzero + transitiveweights("min","max","min"),
#       reference = ~Poisson  (ergm.count's bundled karate-club counts, the
#       model of the ergm.count vignette);
#   (b) a seeded 14-actor DIRECTED count matrix with built-in reciprocity
#       (frozen as an edge list below):  ~ sum + nonzero + mutual(form="min"),
#       reference = ~Poisson.
#
# A tenth, long-chain fit per model (`hp_*`: MCMC.samplesize 16384, interval
# 2048) is frozen as the best available value of the MLE itself; the default
# fits' mean must agree with it, which the script checks.
#
# The log-likelihood R reports for a valued ERGM is relative to the reference
# measure (`logLik` at theta = 0 is defined as 0), estimated by bridge
# sampling; its seed-to-seed spread is frozen too.

suppressMessages({
  .libPaths(c(path.expand("~/R/library"), .libPaths()))
  library(ergm.count)
  library(parallel)
})

seed <- 20261002
set.seed(seed)

num <- function(x) paste(sprintf("%.17g", x), collapse = ", ")
ints <- function(x) paste(sprintf("%d", as.integer(x)), collapse = ", ")
strs <- function(x) paste(sprintf('"%s"', x), collapse = ", ")

# --- (a) zach ----------------------------------------------------------------
data(zach)
W <- as.matrix(zach, attrname = "contexts")
zi <- which(W > 0 & upper.tri(W), arr.ind = TRUE)

# --- (b) a seeded 14-actor directed count matrix with reciprocity -------------
nd <- 14L
A <- matrix(rpois(nd * nd, 0.45), nd, nd)
B <- matrix(rpois(nd * nd, 0.55), nd, nd)
B[lower.tri(B)] <- t(B)[lower.tri(B)]          # a symmetric component
M <- A + B
diag(M) <- 0
dnw <- network(M, directed = TRUE, matrix.type = "adjacency",
               ignore.eval = FALSE, names.eval = "w")
di <- which(M > 0, arr.ind = TRUE)

f_a <- zach ~ sum + nonzero + transitiveweights("min", "max", "min")
f_b <- dnw ~ sum + nonzero + mutual(form = "min")

fit_one <- function(f, response, s, ctrl = control.ergm(seed = s)) {
  set.seed(s)
  fit <- NULL
  warned <- character(0)
  invisible(capture.output(
    fit <- withCallingHandlers(
      suppressMessages(ergm(f, response = response, reference = ~Poisson,
                            control = ctrl)),
      warning = function(w) { warned <<- c(warned, conditionMessage(w))
                              invokeRestart("muffleWarning") }),
    type = "output"))
  ll <- NULL
  invisible(capture.output(ll <- suppressMessages(as.numeric(logLik(fit))),
                           type = "output"))
  list(coef = as.numeric(coef(fit)), se = as.numeric(sqrt(diag(vcov(fit)))),
       loglik = ll, names = names(coef(fit)), warned = warned,
       failed = isTRUE(fit$failure))
}

seeds <- 1:9
run_model <- function(f, response) {
  fits <- mclapply(seeds, function(s) fit_one(f, response, s), mc.cores = 3)
  hp <- fit_one(f, response, 4242,
                control.ergm(seed = 4242, MCMC.samplesize = 16384,
                             MCMC.interval = 2048, MCMLE.maxit = 60))
  coefs <- t(sapply(fits, `[[`, "coef"))
  ses <- t(sapply(fits, `[[`, "se"))
  lls <- sapply(fits, `[[`, "loglik")
  stopifnot(!any(sapply(fits, `[[`, "failed")), !hp$failed)
  out <- list(names = fits[[1]]$names, coefs = coefs, ses = ses, lls = lls,
              coef_mean = colMeans(coefs), coef_sd = apply(coefs, 2, sd),
              se_mean = colMeans(ses), se_sd = apply(ses, 2, sd),
              ll_mean = mean(lls), ll_sd = sd(lls), hp = hp)
  # The default fits' mean must agree with the long-chain fit within the
  # default fits' own spread (4 sd of a mean of 9 plus the long fit's own
  # error, bounded by one default sd)
  stopifnot(all(abs(out$coef_mean - hp$coef) <=
                4 * out$coef_sd * sqrt(1 / length(seeds) + 1)))
  out
}

ra <- run_model(f_a, "contexts")
rb <- run_model(f_b, "w")

# Observed statistics: deterministic
sa <- summary(f_a, response = "contexts")
sb <- summary(f_b, response = "w")

emit <- function(tag, r, s) {
  cat(sprintf("%s_term_names = [%s]\n", tag, strs(r$names)))
  cat(sprintf("%s_summary_statistics = [%s]\n", tag, num(as.numeric(s))))
  cat(sprintf("%s_seeds = [%s]\n", tag, ints(seeds)))
  for (k in seq_along(seeds)) {
    cat(sprintf("%s_coefficients_seed%d = [%s]\n", tag, seeds[k], num(r$coefs[k, ])))
    cat(sprintf("%s_std_errors_seed%d = [%s]\n", tag, seeds[k], num(r$ses[k, ])))
  }
  cat(sprintf("%s_coefficients_mean = [%s]\n", tag, num(r$coef_mean)))
  cat(sprintf("%s_coefficients_sd = [%s]\n", tag, num(r$coef_sd)))
  cat(sprintf("%s_std_errors_mean = [%s]\n", tag, num(r$se_mean)))
  cat(sprintf("%s_std_errors_sd = [%s]\n", tag, num(r$se_sd)))
  cat(sprintf("%s_logliks = [%s]\n", tag, num(r$lls)))
  cat(sprintf("%s_loglik_mean = %.17g\n", tag, r$ll_mean))
  cat(sprintf("%s_loglik_sd = %.17g\n", tag, r$ll_sd))
  cat(sprintf("%s_hp_coefficients = [%s]\n", tag, num(r$hp$coef)))
  cat(sprintf("%s_hp_std_errors = [%s]\n", tag, num(r$hp$se)))
  cat(sprintf("%s_hp_loglik = %.17g\n", tag, r$hp$loglik))
}

cat('name = "count_mcmle"\n\n')

cat("[provenance]\n")
cat(sprintf('r_version = "%s"\n', as.character(getRversion())))
cat(sprintf('ergm_count_version = "%s"\n', as.character(packageVersion("ergm.count"))))
cat(sprintf('ergm_version = "%s"\n', as.character(packageVersion("ergm"))))
cat(sprintf('network_version = "%s"\n', as.character(packageVersion("network"))))
cat(sprintf("seed = %d\n", seed))
cat('script = "test/fixtures/r/count_mcmle.R"\n')
cat(sprintf('date = "%s"\n', format(Sys.Date())))
cat('models = "(a) zach ~ sum + nonzero + transitiveweights(min,max,min), response=contexts, reference=~Poisson; (b) a seeded 14-actor directed count matrix ~ sum + nonzero + mutual(form=min), response=w, reference=~Poisson"\n')
cat('estimator = "ergm.count MCMLE at control.ergm defaults, nine seeds (1..9) per model; hp_* = one fit at MCMC.samplesize=16384, MCMC.interval=2048 (seed 4242)"\n')
cat("\n")

cat("[tolerance]\n")
cat("# Observed statistics are deterministic functions of the graph.\n")
cat("summary_statistics = 1e-9\n")
cat("#\n")
cat("# Both estimators are Monte-Carlo MLEs of the same likelihood, so they are\n")
cat("# compared at R's own resolution. For coefficient k a Julia fit must lie\n")
cat("# within  band_sd * sqrt(sd_k^2 + sd_k^2 / n_seeds)  of R's nine-seed mean,\n")
cat("# where sd_k is R's seed-to-seed standard deviation (`*_coefficients_sd`):\n")
cat("# the first term is one fit's own Monte-Carlo error (taken equal to R's),\n")
cat("# the second the error of R's mean. A band floored at\n")
cat("# `floor_se_fraction` of R's mean standard error keeps a coefficient whose\n")
cat("# nine R fits happened to agree unusually closely from demanding more than\n")
cat("# either sampler's O(1/ESS) bias allows. Standard errors likewise, against\n")
cat("# `*_std_errors_sd`, floored at `floor_se_fraction` of the standard error.\n")
cat("band_sd = 4.0\n")
cat("floor_se_fraction = 0.1\n")
cat("# The relative log-likelihood (bridge sampling on both sides): within\n")
cat("# band_sd * sqrt(2) * R's seed-to-seed sd, floored at 1.0 log unit.\n")
cat("loglik_floor = 1.0\n")
cat("\n")

cat("[values]\n")
cat("# --- (a) zach, frozen as a valued edge list -----------------------------\n")
cat(sprintf("zach_n = %d\n", nrow(W)))
cat(sprintf("zach_edge_src = [%s]\n", ints(zi[, 1])))
cat(sprintf("zach_edge_dst = [%s]\n", ints(zi[, 2])))
cat(sprintf("zach_edge_weight = [%s]\n", ints(W[zi])))
emit("zach", ra, sa)
cat("\n# --- (b) the seeded directed matrix, frozen as a valued edge list -------\n")
cat(sprintf("directed_n = %d\n", nd))
cat(sprintf("directed_edge_src = [%s]\n", ints(di[, 1])))
cat(sprintf("directed_edge_dst = [%s]\n", ints(di[, 2])))
cat(sprintf("directed_edge_weight = [%s]\n", ints(M[di])))
emit("directed", rb, sb)
