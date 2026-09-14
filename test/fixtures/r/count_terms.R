# Golden fixture: statnet `ergm`/`ergm.count` VALUED TERM PARITY — summary
# statistics and the coefficient NAMES R emits — on two networks.
#
# Regenerate from the package root (a second or two; nothing here is Monte Carlo):
#
#   Rscript test/fixtures/r/count_terms.R > test/fixtures/count_terms.toml
#
# WHAT THIS FIXTURE PINS, AND WHY
#
# The 2026-09 panel (item 31, and the ERGMCount reviews) found that the
# package's valued term layer was validated against hand literals only, and
# that two of its terms carried names suggesting R parity they did not have:
# `TransitiveTiesTerm`/`CyclicalTiesTerm` (a triple-wise minimum over ordered
# triples) are NOT `ergm`'s `transitiveweights`/`cyclicalweights` (a dyad's
# value capped by its strongest two-path, Krivitsky 2012 eq. 13). This fixture
# freezes, FROM R, every statistic ERGMCount.jl claims an exact counterpart for:
#
#   (a) on ergm.count's bundled `zach` (undirected, `contexts` counts 0-7):
#       sum, nonzero, greaterthan(2), greaterthan(4), atleast(3),
#       transitiveweights("min","max","min"), cyclicalweights("min","max","min"),
#       smallerthan(2), equalto(3), ininterval(1,3) (R's default open ends) and
#       ininterval(1,3,open=c(FALSE,FALSE)), greaterthan(-1) and atleast(0) —
#       every one a deterministic function of the observed graph, so agreement
#       is at machine precision and any disagreement is a bug in a term formula,
#       full stop. The undirected rows pin R's convention that an undirected
#       network's transitive/cyclical weights are summed over UNORDERED pairs
#       (and coincide); the last two pin that a threshold at or below zero
#       counts the ZERO dyads (561 on zach: every dyad), as smallerthan does.
#   (b) on a seeded 8-actor DIRECTED count matrix generated below (rpois(64, 1.2)
#       with the diagonal zeroed; frozen as an edge list in [values] so the Julia
#       test rebuilds it from the TOML alone): the same terms plus the valued
#       `mutual` forms `min`, `nabsdiff`, `geometric`, `product`, and
#       `ininterval(1,3)` with the four bracket combinations, so the directed
#       (ORDERED-pair) sums, the reciprocity forms and the interval endpoints are
#       all R numbers.
#   (c) the NAMES `summary()` prints for the same formulas, compared exactly:
#       `mutual.min`, `mutual.nabsdiff`, `mutual.geom.mean`, `mutual.product`,
#       `transitiveweights.min.max.min`, `equalto.3.pm.0`, `ininterval(1,3)`,
#       `ininterval[1,3]`, ... ERGMCount.jl labels its coefficients with exactly
#       these strings, so a by-name comparison with a statnet fit works.
#   (e) a 3-actor directed network with an edge whose value is 0 (set with
#       set.edge.value on an existing tie; frozen WITH the zero edge in the
#       edge list): R's valued terms read a 0-valued edge as a zero dyad —
#       nonzero = 1, not the 2 ties the network object holds — and the
#       threshold terms count it with the empty dyads. ERGMCount.jl's
#       `compute` must agree (the 2026-09 round-1 review found `nonzero`
#       counting `ne(net)`), and this is what makes `gof`'s observed
#       statistics comparable with the estimator's `dyad_value` view.
#   (f) a seeded 6-actor DIRECTED network with values in -2:2 (sample(-2:2, 36)
#       with the diagonal zeroed, drawn AFTER the 8-actor matrix so its frozen
#       values are unchanged): `ergm` REFUSES `transitiveweights` and
#       `cyclicalweights` on a network with a negative dyad weight ("Term may
#       not be used with networks with negative dyad weights"; both error
#       strings frozen) and its `mutual(form="geometric")` returns NaN (frozen
#       as the string "NaN"), while sum, nonzero, greaterthan(0), atleast(0),
#       smallerthan(0), equalto(-1), ininterval(-2,1) (both bracket forms) and
#       mutual(min/nabsdiff/product) are ordinary numbers. ERGMCount.jl refuses
#       the three (at `compute`, at model construction and at simulation over a
#       negative support) and must reproduce the rest exactly — the round-2
#       review found `compute(TransitiveWeightsTerm(), neg)` returning 3.0, a
#       statistic R never produces, and `mutual(:geometric)` throwing a bare
#       DomainError from `sqrt`.
#   (d) `mutual(form="threshold", threshold=2)`: ergm 4.12.0's OWN summary fails
#       at C model initialisation ("term with functions ergm::mutual_wt_threshold
#       is declared to have statistics but does not appear to have a change, a
#       difference, or a summary function"), so the statistic CANNOT be pinned
#       against R today. The error string is frozen so the fixture documents why;
#       ERGMCount.jl's `CountMutualTerm(:threshold)` follows the documented
#       definition (binary mutuality after thresholding at >= threshold, which is
#       what R's `emptynwstats` rule — every pair mutual on an empty network when
#       threshold <= 0 — implies) and is tested by hand value only.
#
# There is no estimation here. `ergm.count` has no MPLE for valued models, and
# an MCMLE comparison on a dyad-dependent model would be a Monte-Carlo blur on
# both sides; the exact-MLE fixture (zach_poisson.toml) is where estimation is
# pinned. This file is the term layer's fixture.

suppressMessages({
  .libPaths(c(path.expand("~/R/library"), .libPaths()))
  library(ergm.count)
})

seed <- 20260912   # only the 8-actor matrix is random; recorded for provenance
set.seed(seed)

num <- function(x) paste(ifelse(is.finite(x), sprintf("%.17g", x),
                                ifelse(x > 0, "inf", "-inf")), collapse = ", ")
ints <- function(x) paste(sprintf("%d", as.integer(x)), collapse = ", ")
strs <- function(x) paste(sprintf('"%s"', x), collapse = ", ")
# A TOML basic string: escape backslashes, double quotes and newlines
tstr <- function(x) sprintf('"%s"', gsub("\n", "\\\\n", gsub('"', '\\\\"', gsub("\\\\", "\\\\\\\\", x))))
err_of <- function(expr) tryCatch({ expr; "" }, error = function(e) conditionMessage(e))

# --- (a) zach, undirected --------------------------------------------------
data(zach)
W <- as.matrix(zach, attrname = "contexts")
stopifnot(isSymmetric(W))
ei <- which(W > 0 & upper.tri(W), arr.ind = TRUE)
f_zach <- zach ~ sum + nonzero + greaterthan(2) + greaterthan(4) + atleast(3) +
  transitiveweights("min", "max", "min") + cyclicalweights("min", "max", "min") +
  smallerthan(2) + equalto(3) + ininterval(1, 3) + ininterval(1, 3, open = c(FALSE, FALSE)) +
  greaterthan(-1) + atleast(0)
sum_zach <- summary(f_zach, response = "contexts")

# --- (e) an edge whose value is 0 -------------------------------------------
zw <- network.initialize(3, directed = TRUE)
zw[1, 2] <- 1
zw[2, 3] <- 1
Wz <- matrix(0, 3, 3)
Wz[1, 2] <- 2
Wz[2, 3] <- 0
set.edge.value(zw, "w", Wz)
stopifnot(network.edgecount(zw) == 2)          # the tie exists, its value is 0
zwi <- which(as.matrix(zw) > 0, arr.ind = TRUE)  # the TIES, zero-valued one included
f_zw <- zw ~ nonzero + sum + greaterthan(-1) + atleast(0) + greaterthan(0) + atleast(1) +
  smallerthan(1) + equalto(0) + ininterval(-1, 1)
sum_zw <- summary(f_zw, response = "w")

# --- (b) a seeded 8-actor directed count matrix ----------------------------
M <- matrix(rpois(64, 1.2), 8, 8)
diag(M) <- 0
nw <- network(M, directed = TRUE, matrix.type = "adjacency",
              ignore.eval = FALSE, names.eval = "w")
di <- which(M > 0, arr.ind = TRUE)
f_dir <- nw ~ sum + nonzero + greaterthan(2) + greaterthan(4) + atleast(3) +
  transitiveweights("min", "max", "min") + cyclicalweights("min", "max", "min") +
  mutual(form = "min") + mutual(form = "nabsdiff") + mutual(form = "geometric") +
  mutual(form = "product") + smallerthan(2) + equalto(3) +
  ininterval(1, 3) + ininterval(1, 3, open = c(FALSE, FALSE)) +
  ininterval(1, 3, open = c(TRUE, FALSE)) + ininterval(1, 3, open = c(FALSE, TRUE))
sum_dir <- summary(f_dir, response = "w")

# --- (d) what R cannot compute today ---------------------------------------
err_threshold <- err_of(summary(nw ~ mutual(form = "threshold", threshold = 2), response = "w"))
if (err_threshold == "") stop("ergm's mutual(form=\"threshold\") now works: pin it (add a row above) and drop this branch")

# --- (f) a seeded 6-actor directed network with NEGATIVE counts --------------
N <- matrix(sample(-2:2, 36, replace = TRUE), 6, 6)
diag(N) <- 0
nn <- network(N, directed = TRUE, matrix.type = "adjacency",
              ignore.eval = FALSE, names.eval = "w")
ni <- which(N != 0, arr.ind = TRUE)               # every non-zero dyad is a tie
stopifnot(any(N < 0))
f_neg <- nn ~ sum + nonzero + greaterthan(0) + atleast(0) + smallerthan(0) + equalto(-1) +
  ininterval(-2, 1) + ininterval(-2, 1, open = c(FALSE, FALSE)) +
  mutual(form = "min") + mutual(form = "nabsdiff") + mutual(form = "product")
sum_neg <- summary(f_neg, response = "w")
err_tw_neg <- err_of(summary(nn ~ transitiveweights("min", "max", "min"), response = "w"))
err_cw_neg <- err_of(summary(nn ~ cyclicalweights("min", "max", "min"), response = "w"))
if (err_tw_neg == "" || err_cw_neg == "") stop("ergm now accepts transitiveweights/cyclicalweights on negative weights: pin the values and drop ERGMCount's refusal")
geom_neg <- as.numeric(summary(nn ~ mutual(form = "geometric"), response = "w"))
if (!is.nan(geom_neg)) stop("ergm's mutual(form=\"geometric\") no longer returns NaN on negative weights: pin the value and drop ERGMCount's refusal")

cat('name = "count_terms"\n\n')

cat("[provenance]\n")
cat(sprintf('r_version = "%s"\n', as.character(getRversion())))
cat(sprintf('ergm_count_version = "%s"\n', as.character(packageVersion("ergm.count"))))
cat(sprintf('ergm_version = "%s"\n', as.character(packageVersion("ergm"))))
cat(sprintf('network_version = "%s"\n', as.character(packageVersion("network"))))
cat(sprintf("seed = %d\n", seed))
cat('script = "test/fixtures/r/count_terms.R"\n')
cat(sprintf('date = "%s"\n', format(Sys.Date())))
cat('datasets = "ergm.count::zach (Zachary 1977): 34 karate-club members, undirected, `contexts` edge counts 0-7; an 8-actor directed count matrix rpois(64, 1.2) with zero diagonal under the seed above, frozen below"\n')
cat('model_zach = "zach ~ sum + nonzero + greaterthan(2) + greaterthan(4) + atleast(3) + transitiveweights(min,max,min) + cyclicalweights(min,max,min) + smallerthan(2) + equalto(3) + ininterval(1,3) + ininterval(1,3,open=c(FALSE,FALSE)) + greaterthan(-1) + atleast(0), response=\\"contexts\\" -- summary() only"\n')
cat('model_zero_weight = "zw ~ nonzero + sum + greaterthan(-1) + atleast(0) + greaterthan(0) + atleast(1) + smallerthan(1) + equalto(0) + ininterval(-1,1), response=\\"w\\" -- summary() only; zw is a 3-actor directed network with ties 1->2 (value 2) and 2->3 (value 0)"\n')
cat('model_directed = "nw ~ sum + nonzero + greaterthan(2) + greaterthan(4) + atleast(3) + transitiveweights(min,max,min) + cyclicalweights(min,max,min) + mutual(form=min/nabsdiff/geometric/product) + smallerthan(2) + equalto(3) + ininterval(1,3) with open = (TRUE,TRUE), (FALSE,FALSE), (TRUE,FALSE), (FALSE,TRUE), response=\\"w\\" -- summary() only"\n')
cat('model_negative = "nn ~ sum + nonzero + greaterthan(0) + atleast(0) + smallerthan(0) + equalto(-1) + ininterval(-2,1) + ininterval(-2,1,open=c(FALSE,FALSE)) + mutual(form=min/nabsdiff/product), response=\\"w\\" -- summary() only; nn is a 6-actor directed network with values in -2:2 frozen below; transitiveweights, cyclicalweights are REFUSED by ergm there (errors frozen) and mutual(form=geometric) is NaN"\n')
cat('not_pinned = "mutual(form=threshold): ergm 4.12.0 errors at C model initialisation (r_error_mutual_threshold below); ERGMCount.jl implements the documented definition and tests it by hand value"\n')
cat("\n")

cat("[tolerance]\n")
cat("# Every value is a DETERMINISTIC function of the observed graph -- no\n")
cat("# estimator, no simulation. Machine precision; a disagreement is a bug in a\n")
cat("# term formula (or in the ordered/unordered pair convention). 1e-9 leaves\n")
cat("# room only for the sqrt in mutual.geom.mean. DO NOT LOOSEN.\n")
cat("zach_summary = 1e-9\n")
cat("directed_summary = 1e-9\n")
cat("zero_weight_summary = 1e-9\n")
cat("negative_summary = 1e-9\n\n")

cat("[values]\n")
cat("# --- zach (undirected), frozen as a valued edge list; the Julia test\n")
cat("# rebuilds it exactly -------------------------------------------------\n")
cat(sprintf("zach_n = %d\n", nrow(W)))
cat(sprintf("zach_edge_src = [%s]\n", ints(ei[, 1])))
cat(sprintf("zach_edge_dst = [%s]\n", ints(ei[, 2])))
cat(sprintf("zach_edge_weight = [%s]\n", ints(W[ei])))
cat(sprintf("zach_summary_names = [%s]\n", strs(names(sum_zach))))
cat(sprintf("zach_summary = [%s]\n", num(as.numeric(sum_zach))))
cat("\n# --- the 8-actor DIRECTED count network (tail, head, count) and its\n")
cat("# statistics; every sum is over ORDERED pairs in R's C code --------------\n")
cat(sprintf("directed_n = %d\n", nrow(M)))
cat(sprintf("directed_edge_src = [%s]\n", ints(di[, 1])))
cat(sprintf("directed_edge_dst = [%s]\n", ints(di[, 2])))
cat(sprintf("directed_edge_weight = [%s]\n", ints(M[di])))
cat(sprintf("directed_summary_names = [%s]\n", strs(names(sum_dir))))
cat(sprintf("directed_summary = [%s]\n", num(as.numeric(sum_dir))))
cat("\n# --- a 3-actor DIRECTED network with a ZERO-valued tie (tail, head, value;\n")
cat("# the zero-valued tie is in the list): R reads it as a zero dyad ------------\n")
cat(sprintf("zero_weight_n = %d\n", network.size(zw)))
cat(sprintf("zero_weight_edge_src = [%s]\n", ints(zwi[, 1])))
cat(sprintf("zero_weight_edge_dst = [%s]\n", ints(zwi[, 2])))
cat(sprintf("zero_weight_edge_weight = [%s]\n", ints(Wz[zwi])))
cat(sprintf("zero_weight_summary_names = [%s]\n", strs(names(sum_zw))))
cat(sprintf("zero_weight_summary = [%s]\n", num(as.numeric(sum_zw))))
cat("\n# --- a 6-actor DIRECTED network with NEGATIVE counts (tail, head, value):\n")
cat("# the terms R computes there, and the ones it refuses -------------------\n")
cat(sprintf("negative_n = %d\n", nrow(N)))
cat(sprintf("negative_edge_src = [%s]\n", ints(ni[, 1])))
cat(sprintf("negative_edge_dst = [%s]\n", ints(ni[, 2])))
cat(sprintf("negative_edge_weight = [%s]\n", ints(N[ni])))
cat(sprintf("negative_summary_names = [%s]\n", strs(names(sum_neg))))
cat(sprintf("negative_summary = [%s]\n", num(as.numeric(sum_neg))))
cat(sprintf("r_error_transitiveweights_negative = %s\n", tstr(err_tw_neg)))
cat(sprintf("r_error_cyclicalweights_negative = %s\n", tstr(err_cw_neg)))
cat('r_mutual_geometric_negative = "NaN"\n')
cat("\n# --- what R refuses (conditionMessage of the error): mutual(form=threshold)\n")
cat("# cannot be pinned against ergm 4.12.0 -----------------------------------\n")
cat(sprintf("r_error_mutual_threshold = %s\n", tstr(err_threshold)))
