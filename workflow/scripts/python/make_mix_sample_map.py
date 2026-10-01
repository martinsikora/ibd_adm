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

"""Build the mixmodel sample_map (sample_id, pop_id) from an aggregation panel.

Recipient samples get an "_r" suffix on their pop_id. That suffix is a labelling
device -- mixmodel_ibd.R strips it when building the donor palette and the source
populations -- but it is also the marker pca_ibd.R uses to decide which samples
shape the PC axes, and the one plot_mixmodel.R uses to give a cluster a separate
"_r" entry in the plots.

Pops named by the shared --include_* options are exempted, so a manually curated
recipient-only cluster keeps its own name throughout the plots instead of being
rendered as a recipient of a cluster it no longer belongs to. Keep these options
in step with the ones passed to tvd_matrix.py, or a pop that appears in the TVD
tree will still show up "_r"-suffixed downstream.

Replaces the awk one-liner previously inlined in the rules; with no options it
produces byte-identical output.
"""

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import panel_pops


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--panel', required=True,
                   help='Aggregation panel TSV (sample_id, pop_id, group)')
    p.add_argument('-o', '--out', dest='out_file', required=True)
    panel_pops.add_args(p)
    args = p.parse_args()

    with open(args.panel, newline='') as fh:
        rows = [[c.replace('\r', '') for c in r]
                for r in csv.reader(fh, delimiter='\t') if r]
    if not rows:
        raise ValueError(f'empty panel file: {args.panel}')
    header, body = rows[0], rows[1:]

    pop_ids = [r[1] for r in body]
    groups = [r[2] if len(r) >= 3 else '' for r in body]
    # `dropped` only subtracts from `forced` here (panel_pops.resolve does that)
    # -- it must NOT remove rows. The sample_map defines which samples the model
    # fits at all, so excluding a pop from the TVD tree should not drop its
    # samples from the mixmodel; they keep the default "_r".
    forced, _ = panel_pops.resolve(pop_ids, groups, args)
    panel_pops.report(forced, set(), 'the mixmodel sample_map')

    n_sfx = 0
    out = [header[:2]]
    for r, pop, grp in zip(body, pop_ids, groups):
        if grp == 'recipient' and pop not in forced:
            pop += '_r'
            n_sfx += 1
        out.append([r[0], pop])

    with open(args.out_file, 'w', newline='') as fh:
        csv.writer(fh, delimiter='\t', lineterminator='\n').writerows(out)
    print(f'__ wrote {len(out) - 1} samples, {n_sfx} suffixed "_r" __')


if __name__ == '__main__':
    main()
