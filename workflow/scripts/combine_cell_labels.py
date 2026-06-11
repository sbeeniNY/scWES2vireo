#!/usr/bin/env python3
"""
Combine every pool's Vireo per-cell assignment into one flat table for Seurat.

One row per cell across all pools:
    pool, cell, donor, sample

- `cell`   : the 10x barcode as written by Vireo (e.g. AAACCTGAGACTACAA-1).
             Barcodes repeat across pools, so join on (pool, cell) in R.
- `donor`  : Vireo donor_id — a real donor name, or the literal 'doublet' /
             'unassigned'.
- `sample` : human-readable label from config `sample_labels[pool][donor]`.
             For 'doublet'/'unassigned' the literal token is carried through so
             those cells stay explicitly labelled in Seurat. Falls back to the
             donor id if no label is configured.

Runs as a stdlib-only CLI under a `shell:` rule (rule combine_cell_labels) with
`unset PYTHONPATH` — no pandas / pyyaml (the rule's conda env ships neither).
The `--summary` input (demux_summary.tsv) is read only to verify the pool set
matches before writing, so this step runs after demux_summary.
"""

import argparse
import csv
import os
import re
import sys

SPECIAL = ("doublet", "unassigned")


def load_sample_labels(config_path: str):
    """
    Minimal stdlib parser for the fixed 2-space `sample_labels:` block of
    config.yaml (the cellsnp conda env has no pyyaml). Returns
    {pool: {donor_id: label}}; {} if absent.
    """
    labels: dict = {}
    in_block = False
    cur_pool = None
    with open(config_path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            indent = len(line) - len(line.lstrip(" "))
            stripped = line.strip()
            if indent == 0:  # a top-level key ends the sample_labels block
                in_block = stripped.rstrip(":").strip() == "sample_labels"
                cur_pool = None
                continue
            if not in_block:
                continue
            if indent == 2 and stripped.endswith(":"):
                cur_pool = stripped[:-1].strip()
                labels[cur_pool] = {}
            elif indent >= 4 and cur_pool is not None and ":" in stripped:
                k, v = stripped.split(":", 1)
                labels[cur_pool][k.strip()] = v.strip().strip('"').strip("'")
    return labels


def pool_from_path(path: str) -> str:
    """Extract pool from .../demux/{pool}/vireo/donor_ids.tsv."""
    parts = os.path.normpath(path).split(os.sep)
    try:
        idx = parts.index("vireo")
        return parts[idx - 1]
    except ValueError:
        m = re.search(r"/demux/([^/]+)/vireo/donor_ids\.tsv$", path)
        if not m:
            raise ValueError(f"Cannot infer pool name from path: {path}")
        return m.group(1)


def pools_in_summary(summary_path: str):
    """Distinct pool names present in demux_summary.tsv."""
    pools = set()
    with open(summary_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            p = (row.get("pool") or "").strip()
            if p:
                pools.add(p)
    return pools


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donor-ids", nargs="+", required=True,
                    help="per-pool vireo/donor_ids.tsv paths")
    ap.add_argument("--config", required=True, help="config/config.yaml")
    ap.add_argument("--summary", required=True,
                    help="demux_summary.tsv (for pool-count verification)")
    ap.add_argument("--out", required=True, help="output per-cell TSV")
    args = ap.parse_args()

    sample_labels = load_sample_labels(args.config)

    rows = []
    pools_seen = []
    per_pool_counts = {}
    for tsv_path in args.donor_ids:
        pool = pool_from_path(tsv_path)
        pools_seen.append(pool)
        label_map = sample_labels.get(pool, {}) or {}
        n = 0
        with open(tsv_path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if reader.fieldnames is None or "donor_id" not in reader.fieldnames:
                print(f"WARNING: {tsv_path} missing 'donor_id'; "
                      f"columns={reader.fieldnames}", file=sys.stderr)
                per_pool_counts[pool] = 0
                continue
            for row in reader:
                cell = (row.get("cell") or "").strip()
                if not cell:
                    continue
                donor = (row.get("donor_id") or "").strip()
                if donor in SPECIAL:
                    sample = donor
                else:
                    sample = label_map.get(donor, donor)
                rows.append({"pool": pool, "cell": cell,
                             "donor": donor, "sample": sample})
                n += 1
        per_pool_counts[pool] = n

    # verify the pool set lines up with demux_summary
    seen = set(pools_seen)
    sum_pools = pools_in_summary(args.summary)
    if sum_pools and sum_pools != seen:
        print(f"WARNING: pool-set mismatch vs demux_summary — "
              f"donor_ids={sorted(seen)} summary={sorted(sum_pools)}",
              file=sys.stderr)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    cols = ["pool", "cell", "donor", "sample"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # report
    print("=" * 56, file=sys.stderr)
    print(f"Combined {len(pools_seen)} pools, {len(rows)} cells -> {args.out}",
          file=sys.stderr)
    for p in pools_seen:
        print(f"  {p:12s} {per_pool_counts[p]:>9d} cells", file=sys.stderr)
    print("=" * 56, file=sys.stderr)


if __name__ == "__main__":
    main()
