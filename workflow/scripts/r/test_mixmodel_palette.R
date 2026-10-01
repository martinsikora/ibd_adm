## Unit tests for the palette scale helpers. Run: Rscript workflow/scripts/r/test_mixmodel_palette.R
suppressMessages(library(lsei))
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "mixmodel_palette.R"))
ok <- function(cond, msg) { if (!isTRUE(cond)) stop("FAILED: ", msg); cat("ok:", msg, "\n") }

## model_pred
S <- cbind(a = c(2, 1, 1), b = c(1, 3, 0))
ok(all.equal(model_pred(S, c(.5, .5), raw = FALSE), c(1.5, 2, .5)), "normalised mode returns S p unchanged")
ok(all.equal(sum(model_pred(S, c(.5, .5), raw = TRUE)), 1), "raw mode returns a proportion vector")

## raw_source_matrix: per-individual mean of the summed source palettes
sum_mat <- cbind(A = c(40, 10, 0), B = c(5, 60, 20))
rownames(sum_mat) <- c("A", "B", "C")
m <- raw_source_matrix(sum_mat, n_src = c(A = 4, B = 5))
ok(all.equal(m["A", "A"], 40 / 4), "entries are the sum divided by the number of source individuals")
ok(all.equal(m["C", "B"], 20 / 5), "every entry is a plain per-individual mean (no donor-count rescaling)")

## the point of the option: sources with unequal total IBD, exact mixture of the raw palettes
S1 <- c(60, 20, 10, 10); S2 <- c(10, 10, 20, 10) # totals 100 and 50
Sraw <- cbind(S1, S2)
for (alpha in c(.2, .5, .8)) {
  y <- alpha * S1 + (1 - alpha) * S2
  w <- pnnls(Sraw, y)$x; w <- w / sum(w)
  Sn <- Sraw / rep(colSums(Sraw), each = nrow(Sraw))
  pn <- pnnls(Sn, y / sum(y), sum = 1)$x
  ok(abs(w[1] - alpha) < 1e-8, sprintf("raw fit recovers alpha = %.1f exactly", alpha))
  ok(pn[1] > alpha + 0.01, sprintf("normalised fit over-credits the larger source at alpha = %.1f (%.3f)", alpha, pn[1]))
}
## with equal totals the two scales agree
S3 <- c(40, 30, 20, 10); S4 <- c(10, 20, 30, 40)
Se <- cbind(S3, S4); y <- .3 * S3 + .7 * S4
w <- pnnls(Se, y)$x; w <- w / sum(w)
pn <- pnnls(Se / 100, y / 100, sum = 1)$x
ok(max(abs(w - pn)) < 1e-8, "raw and normalised agree when the sources have equal totals")
cat("all palette tests passed\n")
