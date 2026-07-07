#!/usr/bin/env python3
import argparse
import csv


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Create an aggregation panel from clustering output at a given height."
        )
    )
    parser.add_argument(
        "--clusters",
        required=True,
        help="Clustering TSV from the cluster pipeline (cluster_cut.R)",
    )
    parser.add_argument(
        "--height",
        required=True,
        help="Cluster cut height to select",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output panel TSV with columns sample_id, pop_id, group",
    )
    return parser.parse_args()


def height_match(val, target):
    try:
        return abs(float(val) - float(target)) < 1e-9
    except ValueError:
        return str(val) == str(target)


def map_group(group):
    if group == "cluster_full":
        return "donor_recipient"
    if group == "cluster_min_dist":
        return "recipient"
    return group


def main():
    args = parse_args()
    rows = []
    with open(args.clusters, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if not height_match(row.get("cut_height", ""), args.height):
                continue
            rows.append(row)

    with open(args.out, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "pop_id", "group"])
        for row in rows:
            sample_id = row.get("sample_id", "")
            cluster_label = row.get("cluster_label", "")
            raw_group = row.get("group", "")
            group = map_group(row.get("group", ""))
            pop_id = cluster_label
            if not sample_id or not pop_id:
                continue
            writer.writerow([sample_id, pop_id, group])


if __name__ == "__main__":
    main()
