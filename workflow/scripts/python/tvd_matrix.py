#!/usr/bin/env python3
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import panel_pops


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('files', nargs='+', help='Files with ibd sharing data for each chromosome')
    p.add_argument('-s', '--sample_file', default=None,
                   help='Panel file (sample_id, pop_id, group); restricts TVD recipients '
                        'to group == donor_recipient (cluster_full) individuals')
    p.add_argument('-o', '--out', dest='out_file', required=True)
    panel_pops.add_args(p)
    return p.parse_args()


def main():
    args = parse_args()

    print('__ reading data __')
    dfs = [pd.read_table(f) for f in args.files]
    ibd_pop = pd.concat(dfs, ignore_index=True)

    req = {'pop_id1', 'pop_id2', 'ibd'}
    if not req.issubset(ibd_pop.columns):
        raise ValueError('Input files must include columns: pop_id1, pop_id2, ibd')

    # Restrict recipient side to donor_recipient (cluster_full) individuals so that
    # only cluster_full samples contribute to the per-cluster TVD profiles. The donor
    # side (pop_id2) is already cluster_full-only from aggregate_ibd.py.
    if args.sample_file is not None:
        panel = pd.read_table(args.sample_file)
        if {'sample_id', 'group'}.issubset(panel.columns):
            if 'sample1' not in ibd_pop.columns:
                raise ValueError(
                    'Input files must include a sample1 column to restrict TVD '
                    'recipients by group'
                )
            keep = panel['group'] == 'donor_recipient'
            forced, dropped = panel_pops.resolve(
                panel['pop_id'], panel['group'], args
            )
            panel_pops.report(forced, dropped, 'the TVD matrix')
            if forced:
                keep = keep | panel['pop_id'].isin(forced)
            if dropped:
                keep = keep & ~panel['pop_id'].isin(dropped)

            full_ids = set(panel.loc[keep, 'sample_id'])
            n_before = ibd_pop['sample1'].nunique()
            ibd_pop = ibd_pop[ibd_pop['sample1'].isin(full_ids)]
            print(
                f'__ restricted recipients to {len(full_ids)} samples '
                f'({ibd_pop["sample1"].nunique()}/{n_before} present in IBD data) __'
            )

    print('__ calculating TVD __')
    d = ibd_pop.groupby(['pop_id1', 'pop_id2'], as_index=False)['ibd'].sum()
    group_totals = d.groupby('pop_id1')['ibd'].transform('sum')
    d['p_ibd'] = np.where(group_totals > 0, d['ibd'] / group_totals, 0.0)

    wide = d.pivot(index='pop_id1', columns='pop_id2', values='p_ibd').fillna(0.0)
    m = wide.to_numpy(dtype=float).T
    col_ids = list(wide.index)

    out_parts = []
    for i, pid in enumerate(col_ids):
        r = np.sum(np.abs(m[:, i][:, None] - m) / 2.0, axis=0)
        out_parts.append(pd.DataFrame({'pop_id1': pid, 'pop_id2': col_ids, 'tvd': r}))

    tvd = pd.concat(out_parts, ignore_index=True)

    print('__ writing output __')
    tvd.to_csv(args.out_file, sep='\t', index=False)
    print('__ done! __')


if __name__ == '__main__':
    main()
