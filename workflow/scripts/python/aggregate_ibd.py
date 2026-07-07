#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('-f', '--first', action='store_true', default=False,
                   help='Remove first sample from between-pop donor set instead of random subset.')
    p.add_argument('--seed', type=int, default=None,
                   help='Optional RNG seed for random between-pop donor subset.')
    p.add_argument('-i', '--in', dest='ibd_file', required=True)
    p.add_argument('-o', '--out', dest='out_file', required=True)
    p.add_argument('-s', '--sample_file', required=True)
    return p.parse_args()


def main():
    args = parse_args()

    print('__ reading IBD data __')
    ibd = pd.read_table(
        args.ibd_file,
        header=None,
        names=['sample1', 'sample2', 'chrom', 'ibd'],
        dtype={'sample1': str, 'sample2': str, 'chrom': str, 'ibd': float}
    )

    print('__ reading metadata __')
    sample_map = pd.read_table(args.sample_file)
    req = {'sample_id', 'pop_id', 'group'}
    if not req.issubset(sample_map.columns):
        raise ValueError('sample_file must include columns: sample_id, pop_id, group')

    print('__ aggregating data __')
    ibd_rev = ibd.rename(columns={'sample1': 'sample2', 'sample2': 'sample1'})
    ibd = pd.concat([ibd, ibd_rev], ignore_index=True)

    excl = set(sample_map.loc[sample_map['pop_id'] == 'exclude', 'sample_id'])
    ibd = ibd[~ibd['sample1'].isin(excl) & ~ibd['sample2'].isin(excl)].copy()

    sm1 = sample_map.rename(columns={'sample_id': 'sample1', 'pop_id': 'pop_id1', 'group': 'group1'})
    sm2 = sample_map.rename(columns={'sample_id': 'sample2', 'pop_id': 'pop_id2', 'group': 'group2'})
    ibd = ibd.merge(sm1[['sample1', 'pop_id1', 'group1']], on='sample1', how='left')
    ibd = ibd.merge(sm2[['sample2', 'pop_id2', 'group2']], on='sample2', how='left')

    donors = sample_map[sample_map['group'] == 'donor_recipient'].copy()
    pop_size = donors.groupby('pop_id', as_index=False).size().rename(columns={'size': 'n'})

    donors_multi = donors[donors['pop_id'].isin(pop_size.loc[pop_size['n'] > 1, 'pop_id'])]

    if args.first:
        donors_between = donors_multi[donors_multi.groupby('pop_id', sort=False).cumcount() >= 1]
    else:
        rng = np.random.default_rng(args.seed)
        idx = []
        for _, g in donors_multi.groupby('pop_id'):
            n = len(g)
            keep = n - 1
            take = rng.choice(g.index.to_numpy(), size=keep, replace=False)
            idx.extend(take.tolist())
        donors_between = donors_multi.loc[idx]

    between_ids = set(donors_between['sample_id'])

    grp_cols = ['chrom', 'sample1', 'pop_id1', 'pop_id2']

    ibd_within = ibd[
        (ibd['group2'] == 'donor_recipient') &
        (ibd['pop_id1'] == ibd['pop_id2'])
    ].groupby(grp_cols, as_index=False)['ibd'].sum()
    ibd_within = ibd_within.merge(pop_size.rename(columns={'pop_id': 'pop_id2'}), on='pop_id2', how='left')
    ibd_within['n_inds'] = ibd_within['n'] - 1
    ibd_within = ibd_within.drop(columns=['n'])

    ibd_between = ibd[
        (ibd['group2'] == 'donor_recipient') &
        (ibd['pop_id1'] != ibd['pop_id2']) &
        (ibd['sample2'].isin(between_ids))
    ].groupby(grp_cols, as_index=False)['ibd'].sum()
    ibd_between = ibd_between.merge(pop_size.rename(columns={'pop_id': 'pop_id2'}), on='pop_id2', how='left')
    ibd_between['n_inds'] = ibd_between['n']
    ibd_between = ibd_between.drop(columns=['n'])

    out = pd.concat([ibd_within, ibd_between], ignore_index=True)
    out = out.sort_values(['sample1', 'pop_id2'])

    print('__ writing output __')
    out.to_csv(args.out_file, sep='\t', index=False, compression='gzip' if args.out_file.endswith('.gz') else None)
    print('__ done! __')


if __name__ == '__main__':
    main()
