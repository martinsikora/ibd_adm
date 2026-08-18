#!/usr/bin/env python3
# Copyright 2025 Martin Sikora <martin.sikora@sund.ku.dk>
#
#  This file is free software: you may copy, redistribute and/or modify it
#  under the terms of the GNU General Public License as published by the
#  Free Software Foundation, either version 2 of the License, or (at your
#  option) any later version.
#
#  This file is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Which panel pop_ids are treated as full clusters rather than recipients.

A cluster whose samples are all `recipient` still has a well-defined profile
over the donor palette -- aggregate_ibd.py restricts only the donor side -- so
it can legitimately appear as a population in its own right. By default it does
not: tvd_matrix.py drops it from the TVD matrix, and the mixmodel sample_map
gives it an "_r" suffix that pca_ibd.R and plot_mixmodel.R read as
"cluster_min_dist", which keeps it out of the plots under its own name.

Manually curated single-individual clusters (ancestry anchors split out of a
larger cluster) are recipient-only by construction and are exactly the case
where that default is wrong. These helpers resolve one list of pop_ids to
exempt, so the TVD matrix and the sample_map agree on it.
"""

import argparse


def add_args(p):
    """Attach the shared pop-selection options to an ArgumentParser."""
    p.add_argument('--include_pops', default=None,
                   help='Comma-separated pop_ids to treat as full clusters even '
                        'when their samples are group == recipient')
    p.add_argument('--include_recipient_only_pops', action='store_true',
                   help='Treat every pop_id that has no donor_recipient sample at '
                        'all as a full cluster (i.e. recipient-only by construction)')
    p.add_argument('--exclude_pops', default=None,
                   help='Comma-separated pop_ids to drop, applied after the include '
                        'options; use for catch-all bins such as "unassigned"')
    return p


def csv_arg(v):
    return [x.strip() for x in v.split(',') if x.strip()] if v else []


def resolve(pop_ids, groups, args):
    """Return (forced, dropped) sets of pop_ids.

    pop_ids/groups are parallel sequences over the panel rows. `forced` are
    pop_ids to treat as full clusters despite being recipients; `dropped` are
    pop_ids to remove outright. Exclusions win over inclusions.
    """
    dropped = set(csv_arg(args.exclude_pops))
    forced = set()
    if getattr(args, 'include_recipient_only_pops', False):
        rows = list(zip(pop_ids, groups))
        donor_pops = {p for p, g in rows if g == 'donor_recipient'}
        forced |= {p for p, _ in rows} - donor_pops
    forced |= set(csv_arg(args.include_pops))
    forced -= dropped
    return forced, dropped


def report(forced, dropped, what):
    if forced:
        print(f'__ treating {len(forced)} recipient-only pop(s) as full clusters '
              f'for {what}: {", ".join(sorted(forced))} __')
    if dropped:
        print(f'__ excluding pop(s) from {what}: {", ".join(sorted(dropped))} __')
