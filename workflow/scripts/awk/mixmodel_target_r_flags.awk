# Copyright 2025 Martin Sikora <martin.sikora@sund.ku.dk>. GPL v2+ (see repo).
#
# Per-target companion to mixmodel_source_r_flags.awk: how much of each target's estimate
# rests on sources whose R marks them as offset-prone.
#
# The fitted p is an IBD share, roughly f_i * R_i. A source with a low R has a deflated p,
# so tiering on the raw p under-calls exactly the worst cases. The tier is therefore
# computed on the R-corrected share q = (p / R) renormalised over the sources. 1/R
# over-corrects as an estimator, which is the wanted behaviour for a screen: q_flagged is
# an upper bound on what the flagged sources could contribute. It is not a corrected
# proportion, and no per-source correction was found that generalises across cohorts.
# Raw p is reported alongside for reference.
#
# Usage:
#   gawk -v fQC=<source_R_flags.tsv> [-v PMIN=0.002] [-v HIGH=0.02] [-v MODERATE=0.05] \
#        [-v LOW=0.01] -f mixmodel_target_r_flags.awk <source_R_flags.tsv> <mixmodel_*.tsv>
#
# A source with raw p below PMIN cannot escalate the tier: dividing a p of 1e-4 by the
# panel's smallest R manufactures a sizeable share out of nothing (47 Patterson targets
# on an axis with no real offset in the first version).
#
# Tiers: HIGH if the R-corrected share on SEVERE sources is >= HIGH; MODERATE if the share
# on WARN sources is >= MODERATE; LOW if the total flagged share is >= LOW; else none.

BEGIN {
  FS = OFS = "\t"
  if (PMIN == "") PMIN = 0.002
  if (HIGH == "") HIGH = 0.02
  if (MODERATE == "") MODERATE = 0.05
  if (LOW == "") LOW = 0.01
  PROCINFO["sorted_in"] = "@ind_str_asc"
}
FILENAME == fQC {
  if ($1 ~ /^#/) next
  if ($1 == "component") {
    for (i = 1; i <= NF; i++) h[$i] = i
    if (!("R_excl_self" in h) || !("flag" in h)) {
      print "mixmodel_target_r_flags: unexpected source flag header" > "/dev/stderr"; exit 1
    }
    next
  }
  R[$1] = $(h["R_excl_self"]) + 0
  flag[$1] = $(h["flag"])
  next
}
FNR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
$(c["group"]) != "target" { next }
{
  t = $(c["sample_id"]); s = $(c["source_pop"]); p = $(c["p"]) + 0
  if (p <= 0 || !(s in R) || R[s] <= 0) next
  q = p / R[s]
  esc = (p >= PMIN)
  tot[t] += p; qtot[t] += q; lab[t] = $(c["label"])
  if (flag[s] == "SEVERE") {
    sev[t] += p
    if (esc) { qsev[t] += q; if (q > wsq[t] + 0) { wsq[t] = q; wss[t] = s } }
  } else if (flag[s] == "WARN") {
    warn[t] += p
    if (esc) { qwarn[t] += q; if (q > wwq[t] + 0) { wwq[t] = q; wws[t] = s } }
  }
}
END {
  print "sample_id", "label", "p_severe", "p_warn", "q_severe", "q_warn", "q_flagged", "top_flagged_source", "risk"
  for (t in tot) {
    if (tot[t] <= 0 || qtot[t] <= 0) continue
    a = (sev[t] + 0) / tot[t]; b = (warn[t] + 0) / tot[t]
    qa = (qsev[t] + 0) / qtot[t]; qb = (qwarn[t] + 0) / qtot[t]; qf = qa + qb
    top = (wsq[t] + 0 >= wwq[t] + 0) ? wss[t] : wws[t]
    if (top == "") top = "-"
    if      (qa >= HIGH)     r = "HIGH"
    else if (qb >= MODERATE) r = "MODERATE"
    else if (qf >= LOW)      r = "LOW"
    else                     r = "none"
    printf "%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%s\t%s\n", t, lab[t], a, b, qa, qb, qf, top, r
  }
}
