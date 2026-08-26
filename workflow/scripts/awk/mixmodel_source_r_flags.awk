# Copyright 2025 Martin Sikora <martin.sikora@sund.ku.dk>. GPL v2+ (see repo).
#
# QC flag for mixture sources based on genome-wide per-individual IBD emission R.
#
# WHY. The fitted p is an IBD *share*, not a genome share: p_i ~ f_i * R_i, where
# R_i is the total IBD individual i emits into the donor palette. normalize_matrix_
# cols() in mixmodel_ibd.R scales every palette column to sum 1 and discards those
# totals, so the fit never sees R. A source whose R sits far below the panel's
# typical R therefore has its proportion systematically DEFLATED (and one far above,
# inflated). Because R is a property of the source and the panel, not of any target,
# it can be flagged before looking at a single fit.
#
# WHAT THIS IS NOT. It flags risk and DIRECTION. It is not a correction and does not
# give a magnitude. Validated against published per-individual qpAdm (Harney 2023
# Catoctin, Brielle 2023 Swahili, Simoes 2023 / Serrano 2023 North Africa, Patterson
# 2021 Britain): the multiplicative offset for one and the same flagged source ranges
# from 1.3x to 6.8x across cohorts, and no per-source transform, emission-budget
# residual, or segment-length stratification predicted which applies. See
# docs/ for the full negative result; do not reintroduce a correction here.
#
# A SEVERE/WARN TIER DOES NOT MEAN "BAD SOURCE". It means the source's proportions
# sit on a different scale from the rest of the panel. Low R and low informativeness
# are unrelated. Worked example: X01:Ethiopia_Neolithic has the lowest R in
# world_base_2 (589, n=1) yet is the dominant source for the Hadza at p = 0.46-0.83
# with p/se = 66-241 -- among the most confident assignments in the whole run -- and
# correctly orders Sandawe/Somali/Masai/Mbuti/Dinka above Yemeni/Egyptian. Its
# RANKING is excellent while its SCALE is off ~400x. It is also the panel's only
# proxy for deep East African ancestry. Never drop a source on the strength of R.
#
# TIERS ARE PANEL-RELATIVE. R is a donor-panel-size-weighted average of coalescent
# rates, so the median is "a typical source in THIS panel". Flags are not comparable
# between panels -- adding donors of one ancestry moves every other source's tier
# without anything about those sources changing. The panel and median are stamped
# into the header for exactly this reason.
#
# INPUTS
#   src_prof.tsv  pop_id1, pop_id2, sum(ibd) over source samples  (= n_i * R_i)
#   SRCNF         pop_id1, n_source_individuals
# Both produced by ibd_residual_profiles.awk in the same streaming pass.
#
# Usage:
#   gawk -v NF_FILE=src_n.tsv -v PANEL=<name> [-v WARN=2.5] [-v SEVERE=10] \
#        -f mixmodel_source_r_flags.awk src_prof.tsv
#
# R_excl_self is the REPORTED statistic (fold_vs_median is computed from it):
# within-component sharing is 8% of the total for Morocco but 21-33% for
# YRI/JuHoan/ShumLaka, i.e. it inflates precisely the sources most at risk of being
# flagged. R_genomewide is emitted alongside for reference. For n=1 components the
# two are identical by construction.

BEGIN {
  FS = OFS = "\t"
  if (WARN   == "") WARN   = 2.5
  if (SEVERE == "") SEVERE = 10
  if (PANEL  == "") PANEL  = "NA"
  while ((getline l < NF_FILE) > 0) { split(l, nf, "\t"); nsrc[nf[1]] = nf[2] + 0 }
}
{
  p1 = $1; p2 = $2; v = $3 + 0
  tot[p1] += v
  if (p1 != p2) ext[p1] += v            # exclude within-component sharing
}
END {
  n = 0
  for (c in tot) {
    k = (c in nsrc && nsrc[c] > 0) ? nsrc[c] : 1
    n++
    comp[n] = c; cnt[n] = k
    Rg[n] = tot[c] / k
    Rx[n] = (ext[c] + 0) / k
    rv[n] = Rx[n]
  }
  if (n == 0) { print "mixmodel_source_r_flags: no source rows" > "/dev/stderr"; exit 1 }
  asort(rv)
  med = (n % 2) ? rv[(n + 1) / 2] : (rv[n / 2] + rv[n / 2 + 1]) / 2

  for (i = 1; i <= n; i++) idx[i] = i    # ascending R: most-deflated first
  for (i = 1; i <= n; i++)
    for (j = i + 1; j <= n; j++)
      if (Rx[idx[j]] < Rx[idx[i]]) { t = idx[i]; idx[i] = idx[j]; idx[j] = t }

  printf "# panel=%s  n_sources=%d  median_R_excl_self=%.1f  WARN=%gx  SEVERE=%gx\n",
         PANEL, n, med, WARN, SEVERE
  printf "# tier = scale offset vs this panel's median; NOT a statement about source quality\n"
  print "component", "n", "R_excl_self", "R_genomewide", "fold_vs_median", "direction", "flag"
  for (m = 1; m <= n; m++) {
    i = idx[m]
    if (med <= 0 || Rx[i] <= 0) { f = "NA"; dir = "NA"; fold = 0 }
    else {
      r = Rx[i] / med
      fold = (r < 1) ? 1 / r : r
      dir  = (r < 1) ? "underestimated" : "overestimated"
      if      (fold >= SEVERE) f = "SEVERE"
      else if (fold >= WARN)   f = "WARN"
      else                   { f = "ok"; dir = "-" }
    }
    printf "%s\t%d\t%.1f\t%.1f\t%.2f\t%s\t%s\n", comp[i], cnt[i], Rx[i], Rg[i], fold, dir, f
  }
}
