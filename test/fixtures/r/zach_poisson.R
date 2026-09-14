# Golden fixture: statnet `ergm.count` on Zachary's karate club (valued).
#
# Regenerate from the package root (~90s: the replication study refits the MCMLE
# under five further seeds):
#
#   Rscript test/fixtures/r/zach_poisson.R > test/fixtures/zach_poisson.toml
#
# THE MODEL, AND WHY IT IS DYAD-INDEPENDENT ON PURPOSE
#
# `zach` ships with ergm.count: Zachary's 34 monks^W karate-club members, with an
# edge attribute `contexts` counting in how many of 8 social contexts each pair
# interacted (0-7). The model is
#
#     zach ~ sum + nonzero,  response = "contexts",  reference = ~Poisson
#
# Both terms are functions of a single dyad's own count, so the model is
# DYAD-INDEPENDENT: under the Poisson reference h(y) = 1/y!, every dyad is an
# independent draw from the two-parameter law
#
#     P(y) = exp(theta_sum * y + theta_nonzero * 1{y>0}) / (y! * Z),
#     Z(theta) = 1 + exp(theta_nonzero) * (exp(exp(theta_sum)) - 1).
#
# That has an EXACT maximum likelihood estimator, and this script solves for it
# ANALYTICALLY (the score equations collapse to one monotone scalar root-find;
# see below) and freezes it as `exact_*`. It is the point of the fixture.
#
# It is worth saying why it is not just `optim`-ed, because `optim` was the first
# thing tried: BFGS at reltol = 1e-14 landed ~7e-6 from the true optimum, which is
# ABOVE the 1e-6 tolerance this fixture wants to hold ERGMCount.jl to. A golden
# value that is itself only good to 7e-6 cannot police anyone at 1e-6. The frozen
# `exact_score` is now ~1e-14, and that is what licenses the tolerance.
#
# ERGMCount.jl's `count_mple` enumerates each dyad's full conditional over the
# count support. For a dyad-independent model the conditional IS the marginal, so
# the pseudo-likelihood IS the likelihood: ERGMCount.jl is computing the exact
# MLE, and it must agree with R's exact MLE to optimizer precision (1e-6). No
# Monte Carlo, no excuses. A dyad-dependent model (transitiveties, mutual) would
# have been a prettier demonstration and would have made every number a
# Monte-Carlo blur on both sides, testing nothing sharply.
#
# `ergm.count`'s own fit is MCMLE (there is no MPLE for valued ERGMs in statnet),
# so it is frozen SEPARATELY, at a tolerance derived from its measured
# seed-to-seed spread. It is a check on the exact answer, not the reference for
# it: R's MCMLE and R's exact MLE differ by more than ERGMCount.jl and R's exact
# MLE do. See [values].mcmle_vs_exact_max_abs_diff.
#
# TRUNCATION -- THE THING THAT COULD SILENTLY INVALIDATE THIS COMPARISON
#
# The Poisson reference has UNBOUNDED support, and ERGMCount.jl cannot enumerate
# an infinite one: it truncates each dyad's conditional at 0:max_val. The exact
# MLE above does not truncate. So the two are estimating the same thing ONLY IF
# the probability mass past max_val is negligible; otherwise ERGMCount.jl is
# fitting a different (truncated) model and any agreement is luck.
#
# It is negligible here, and the fixture proves it rather than asserting it.
# `boundary_tail_mass` in [values] is the exact P(y > max_val) under the fitted
# law at max_val = 30, computed in log space from the same normalizing constant
# (`1 - sum(head)` underflows to exactly 0 in double precision, which would look
# like a suspiciously round number rather than the real bound it is). It comes
# out around 1e-22, against a smallest-used conditional probability of ~1e-3
# (`min_used_conditional_prob`) -- nineteen orders of magnitude. The truncation
# moves the likelihood far below the last bit of a Float64; it cannot move a
# coefficient at the 1e-6 level being asserted. The Julia testset re-checks
# ERGMCount.jl's own reported `boundary_mass` against the same bound, so a future
# change to max_val that started to bite would turn this test red rather than
# quietly widen the gap.
#
# This choice of data is exactly why the model is `sum + nonzero` on zach and not
# something with a fat tail: the counts run 0-7, the fitted Poisson rate is ~2.8,
# and nothing is anywhere near the boundary.

suppressMessages({
  .libPaths(c(path.expand("~/R/library"), .libPaths()))
  library(ergm.count)
})

seed <- 20260713
set.seed(seed)

data(zach)
W <- as.matrix(zach, attrname = "contexts")
n <- nrow(W)
stopifnot(isSymmetric(W))
y <- W[upper.tri(W)]              # the 561 undirected dyads
N <- length(y)
max_val <- 30L

S_sum <- sum(y)
S_nonzero <- sum(y > 0)

# --- the EXACT MLE of the dyad-independent Poisson-reference model ----------
# Untruncated: Z is summed in closed form, Z = 1 + e^{t2}(e^{e^{t1}} - 1).
#
# This is solved ANALYTICALLY, not by `optim`. It was written with optim(BFGS,
# reltol=1e-14) first, and that turned out to be the loose side of the
# comparison: it left the coefficients ~7e-6 from the true optimum, which is
# ABOVE the 1e-6 this fixture wants to assert. A golden value that is itself only
# good to 7e-6 cannot be used to hold anyone to 1e-6, so it is not used.
#
# The score equations for a dyad-independent exponential family are just moment
# matching, E[T] = T_obs, and with lambda = e^{t1}, b = e^{t2} they collapse:
#
#   dlogZ/dt2 = b(e^lam - 1)/Z = p            (p   = nonzero / N)
#   dlogZ/dt1 = b lam e^lam / Z = m           (m   = sum / N)
#
# The first gives Z = 1/(1-p) directly, hence b(e^lam - 1) = p/(1-p); dividing
# the second by the first eliminates b entirely and leaves ONE scalar equation,
#
#   lam * e^lam / (e^lam - 1) = m / p,
#
# monotone in lam, solved by uniroot to 1e-14. b, and therefore both
# coefficients, then follow in closed form. This is the exact MLE to ~1e-15.
logZ <- function(th) log1p(exp(th[2]) * expm1(exp(th[1])))

m <- S_sum / N
p <- S_nonzero / N
g_lam <- function(lam) lam * exp(lam) / expm1(lam) - m / p
lam_hat <- uniroot(g_lam, c(1e-8, 50), tol = 1e-14)$root
b_hat <- (p / (1 - p)) / expm1(lam_hat)
exact_coef <- c(log(lam_hat), log(b_hat))

# Fisher information = N * Cov(T) under the fitted law, T = (y, 1{y>0}). Summed
# exactly over the support in log space (the tail is negligible -- see below --
# but it costs nothing to carry 400 terms).
supp <- 0:400
lp <- exact_coef[1] * supp + exact_coef[2] * (supp > 0) -
      lfactorial(supp) - logZ(exact_coef)
pr <- exp(lp)
E1 <- sum(pr * supp); E2 <- sum(pr * (supp > 0))
V11 <- sum(pr * supp^2) - E1^2
V22 <- sum(pr * (supp > 0)) - E2^2
V12 <- sum(pr * supp * (supp > 0)) - E1 * E2
FI <- N * matrix(c(V11, V12, V12, V22), 2, 2)
exact_se <- sqrt(diag(solve(FI)))
exact_loglik <- exact_coef[1] * S_sum + exact_coef[2] * S_nonzero -
                N * logZ(exact_coef) - sum(lfactorial(y))

# The score at the analytic solution: this is HOW EXACT the golden value is, and
# it is the number that licenses asserting 1e-6 against it.
exact_score <- c(S_sum - N * E1, S_nonzero - N * E2)

# --- the truncation check ---------------------------------------------------
# P(y > max_val) under the fitted law. If this were not negligible, ERGMCount.jl
# (which enumerates 0:max_val) and the exact MLE above would be fitting different
# models and the whole comparison would be void.
# Computed in LOG space: the tail is ~1e-22 and 1 - sum(head) underflows to
# exactly 0 in double precision, which would look like a suspiciously round
# number rather than the real bound it is.
log_unnorm <- function(yv, th) th[1] * yv + th[2] * (yv > 0) - lfactorial(yv)
lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }
tail_ks <- (max_val + 1):(max_val + 400)
log_tail_mass <- lse(log_unnorm(tail_ks, exact_coef)) - logZ(exact_coef)
tail_mass <- exp(log_tail_mass)

# For scale: the smallest conditional probability the estimator actually uses on
# the observed support (0..7). The truncation is negligible only RELATIVE to this.
log_used <- log_unnorm(0:max(y), exact_coef) - logZ(exact_coef)
min_used_prob <- exp(min(log_used))

# --- ergm.count's own MCMLE -------------------------------------------------
# ergm() prints iteration chatter to stdout, which IS this fixture. Suppress it.
quiet_ergm <- function(s) {
  set.seed(s)
  fit <- NULL
  invisible(capture.output(
    fit <- ergm(zach ~ sum + nonzero, response = "contexts",
                reference = ~Poisson,
                control = control.ergm(seed = s, MCMC.samplesize = 4096,
                                       MCMC.burnin = 16384,
                                       MCMC.interval = 1024)),
    type = "output"))
  fit
}
fit_mc <- quiet_ergm(seed)
mcmle_coef <- coef(fit_mc)
mcmle_se <- sqrt(diag(vcov(fit_mc)))

rep_seeds <- c(101, 202, 303, 404, 505)
reps <- t(sapply(rep_seeds, function(s) coef(quiet_ergm(s))))
mcmle_sd <- apply(reps, 2, sd)
mcmle_vs_exact <- max(abs(as.numeric(mcmle_coef) - exact_coef))

num <- function(x) paste(sprintf("%.17g", x), collapse = ", ")
strs <- function(x) paste(sprintf('"%s"', x), collapse = ", ")
ei <- which(W > 0 & upper.tri(W), arr.ind = TRUE)

cat('name = "zach_poisson"\n\n')

cat("[provenance]\n")
cat(sprintf('r_version = "%s"\n', as.character(getRversion())))
cat(sprintf('ergm_count_version = "%s"\n', as.character(packageVersion("ergm.count"))))
cat(sprintf('ergm_version = "%s"\n', as.character(packageVersion("ergm"))))
cat(sprintf('network_version = "%s"\n', as.character(packageVersion("network"))))
cat(sprintf("seed = %d\n", seed))
cat('script = "test/fixtures/r/zach_poisson.R"\n')
cat(sprintf('date = "%s"\n', format(Sys.Date())))
cat('dataset = "ergm.count::zach (Zachary 1977): 34 karate-club members, undirected, `contexts` edge counts 0-7"\n')
cat('model = "zach ~ sum + nonzero, response=\\"contexts\\", reference=~Poisson -- DYAD-INDEPENDENT, so an exact MLE exists and is computed below"\n')
cat('exact_estimator = "analytic: the score equations reduce to lam*e^lam/(e^lam-1) = mean(y)/P(y>0), solved by uniroot to 1e-14; NOT optim, which was only good to 7e-6 and could not license the 1e-6 tolerance below"\n')
cat('mcmle_control = "control.ergm(MCMC.samplesize=4096, MCMC.burnin=16384, MCMC.interval=1024)"\n')
cat(sprintf('replication_seeds = "%s"\n', paste(rep_seeds, collapse = ",")))
cat(sprintf("max_val = %d\n", max_val))
cat("\n")

cat("[tolerance]\n")
cat("# Observed sufficient statistics: a deterministic function of the data.\n")
cat("# Machine precision; a disagreement is a bug in a term formula.\n")
cat("summary_statistics = 1e-9\n")
cat("#\n")
cat("# EXACT MLE. sum + nonzero under a Poisson reference is dyad-independent, so\n")
cat("# ERGMCount.jl's dyad-conditional enumeration IS the likelihood and it is\n")
cat("# computing the same exact MLE that `optim` computes here. Both are solving a\n")
cat("# two-parameter concave problem to convergence, with no Monte Carlo on either\n")
cat("# side. 1e-6 is optimizer precision and ~1e-5 of the smallest standard error\n")
cat("# in the model. If this fails it is a BUG in the estimator, the reference\n")
cat("# measure, or a change statistic. DO NOT LOOSEN IT.\n")
cat("#\n")
cat("# The one thing that could make this comparison meaningless is truncation:\n")
cat("# ERGMCount.jl enumerates 0:max_val and the exact MLE does not truncate at\n")
cat("# all. `boundary_tail_mass` below is the exact P(y > 30) under the fitted\n")
cat("# law -- see `boundary_tail_mass_log10` -- and it is ~19 orders of magnitude\n")
cat("# below the smallest conditional probability the estimator actually uses\n")
cat("# (`min_used_conditional_prob`). The truncation cannot\n")
cat("# reach the 1e-6 being asserted. That is a property of THIS dataset (counts\n")
cat("# 0-7, fitted rate ~2.8), not of the estimator, and the Julia testset\n")
cat("# re-checks ERGMCount.jl's own reported `boundary_mass` so that a future\n")
cat("# max_val that started to bite goes red instead of quietly drifting.\n")
cat("exact_coefficients = 1e-6\n")
cat("exact_std_errors = 1e-6\n")
cat("#\n")
cat("# ergm.count's OWN fit is MCMLE -- statnet has no MPLE for valued ERGMs -- so\n")
cat("# it carries Monte-Carlo error that the exact MLE does not. Read `mcmle_sd`\n")
cat("# in [values]: that is ergm.count refitting this model under five further\n")
cat("# seeds and disagreeing with itself. It is the floor here. The frozen MCMLE\n")
cat("# fit sits `mcmle_vs_exact_max_abs_diff` away from the exact answer -- which\n")
cat("# is FURTHER than ERGMCount.jl sits from it. So this comparison is a\n")
cat("# consistency check on ergm.count, not the reference standard; the reference\n")
cat("# standard is `exact_coefficients` above. Tolerance = 6x the largest observed\n")
cat("# MCMLE seed-to-seed sd, which comfortably covers the single frozen draw.\n")
cat(sprintf("mcmle_coefficients = %.4g\n", max(0.05, 6 * max(mcmle_sd))))
cat("\n")

cat("[values]\n")
cat("# --- the network, frozen. Julia rebuilds it exactly. --------------------\n")
cat(sprintf("n_actors = %d\n", n))
cat("directed = false\n")
cat(sprintf("edge_src = [%s]\n", paste(ei[, 1], collapse = ", ")))
cat(sprintf("edge_dst = [%s]\n", paste(ei[, 2], collapse = ", ")))
cat(sprintf("edge_weight = [%s]\n", paste(W[ei], collapse = ", ")))
cat(sprintf("n_dyads = %d\n", N))
cat(sprintf("max_val = %d\n", max_val))
cat("\n# --- observed sufficient statistics, deterministic ----------------------\n")
cat('summary_statistic_names = ["sum", "nonzero"]\n')
cat(sprintf("summary_statistics = [%s]\n", num(c(S_sum, S_nonzero))))
cat("\n# --- the EXACT MLE. This is the assertion with teeth. -------------------\n")
cat('term_names = ["sum", "nonzero"]\n')
cat(sprintf("exact_coefficients = [%s]\n", num(exact_coef)))
cat(sprintf("exact_std_errors = [%s]\n", num(exact_se)))
cat(sprintf("exact_loglik = %.17g\n", exact_loglik))
cat("# The score (= observed minus expected sufficient statistics) at the frozen\n")
cat("# solution. This is HOW EXACT the golden value is, and it is what licenses\n")
cat("# asserting 1e-6 against it. Compare: `optim(BFGS, reltol=1e-14)`, the\n")
cat("# obvious way to write this, left a score of order 1e-3 and coefficients\n")
cat("# ~7e-6 off -- above the tolerance it was supposed to police.\n")
cat(sprintf("exact_score = [%s]\n", num(exact_score)))
cat("\n# P(y > max_val) under the fitted law: the entire error the 0:max_val\n")
cat("# truncation introduces into ERGMCount.jl's likelihood. If this were not\n")
cat("# negligible the comparison above would be void.\n")
cat(sprintf("boundary_tail_mass = %.17g\n", tail_mass))
cat(sprintf("boundary_tail_mass_log10 = %.17g\n", log_tail_mass / log(10)))
cat("# For scale: the SMALLEST conditional probability the estimator actually\n")
cat("# uses on the observed support (counts 0-7). The truncation error is\n")
cat("# negligible only relative to this, and it is ~19 orders of magnitude below.\n")
cat(sprintf("min_used_conditional_prob = %.17g\n", min_used_prob))
cat("\n# --- ergm.count's own MCMLE, for cross-checking the exact answer --------\n")
cat(sprintf("mcmle_coefficients = [%s]\n", num(as.numeric(mcmle_coef))))
cat(sprintf("mcmle_std_errors = [%s]\n", num(as.numeric(mcmle_se))))
cat("# ergm.count disagreeing with ITSELF over five further seeds.\n")
cat(sprintf("mcmle_sd = [%s]\n", num(mcmle_sd)))
cat("# ...and how far its frozen fit sits from the exact MLE. Compare this with\n")
cat("# the gap the Julia testset reports between ERGMCount.jl and the exact MLE.\n")
cat(sprintf("mcmle_vs_exact_max_abs_diff = %.17g\n", mcmle_vs_exact))
