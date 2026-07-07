# Copyright 2025 Martin Sikora <martin.sikora@sund.ku.dk>. GPL v2+ (see repo).
#
# One streaming pass over the concatenated per-chromosome ibd_pop tables
# (columns: chrom, sample1, pop_id1, pop_id2, ibd, n_inds) producing the
# genome-total profiles needed by mixmodel_residual_diagnostic.R:
#   POPF  pop_id1, pop_id2, sum(ibd)   over ALL samples   -> palette P
#   SRCPF pop_id1, pop_id2, sum(ibd)   over SOURCE samples -> source_mat S
#   VALPF sample1, pop_id2, sum(ibd)   over VAL samples    -> validation y
# SRCF / VALF list the source / validation sample ids (one per line).
BEGIN {
  FS = OFS = "\t"
  while ((getline l < SRCF) > 0) src[l] = 1
  while ((getline l < VALF) > 0) val[l] = 1
}
$1 == "chrom" { next }                       # skip repeated headers
{
  s = $2; p1 = $3; p2 = $4; v = $5 + 0
  pop[p1 SUBSEP p2] += v
  if (s in src) sp[p1 SUBSEP p2] += v
  if (s in val) vp[s SUBSEP p2] += v
}
END {
  for (k in pop) { split(k, a, SUBSEP); print a[1], a[2], pop[k] > POPF }
  for (k in sp)  { split(k, a, SUBSEP); print a[1], a[2], sp[k]  > SRCPF }
  for (k in vp)  { split(k, a, SUBSEP); print a[1], a[2], vp[k]  > VALPF }
}
