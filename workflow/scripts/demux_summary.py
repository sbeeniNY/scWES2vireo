#!/usr/bin/env python3
"""
Aggregate per-pool Vireo `donor_ids.tsv` outputs into one cohort-wide summary.

Runs as a standalone CLI (invoked from `rules/demux.smk :: rule demux_summary`
via a `shell:` directive — NOT a Snakemake `script:` directive). The `script:`
directive is avoided because Snakemake can inject the host snakemake env's
package path into the rule's conda env, so `import pandas` may load a mismatched
build and crash with a numpy ABI error. Using a plain CLI under `shell:` (with
`unset PYTHONPATH`) keeps imports inside the rule's own env, and this script
depends only on the standard library (the cellsnp env has no pyyaml, so the
small fixed `sample_labels:` config block is parsed by hand).

Output columns:
    pool, donor_id, mapped_label, n_cells, doublet_rate, unassigned_rate, total_cells
"""

import argparse
import csv
import os
import re
import sys
from collections import Counter


def load_sample_labels(config_path: str):
    """
    Minimal stdlib parser for the `sample_labels:` block of config.yaml.

    Avoids a pyyaml dependency: the rule's cellsnp conda env does not ship
    pyyaml, and pulling it from the host snakemake env (PYTHONPATH leak) is what
    we deliberately stopped doing. The block is a controlled, fixed 2-space
    structure:
        sample_labels:
          <pool>:
            <donor_id>: <label>
    Returns {pool: {donor_id: label}}; {} if the block is absent.
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
            if indent == 0:  # a top-level key ends any sample_labels block
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


def count_donor_ids(tsv_path: str):
    """Return (Counter of donor_id values, total_rows). Empty/headerless -> (Counter(), 0)."""
    counts: Counter = Counter()
    total = 0
    with open(tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or "donor_id" not in reader.fieldnames:
            print(f"WARNING: {tsv_path} missing 'donor_id' column; "
                  f"columns={reader.fieldnames}", file=sys.stderr)
            return counts, total
        for row in reader:
            total += 1
            counts[(row.get("donor_id") or "").strip()] += 1
    return counts, total


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donor-ids", nargs="+", required=True,
                    help="per-pool donor_ids.tsv paths")
    ap.add_argument("--config", required=True, help="config/config.yaml")
    ap.add_argument("--out", required=True, help="output summary TSV")
    args = ap.parse_args()

    sample_labels = load_sample_labels(args.config)

    rows = []
    for tsv_path in args.donor_ids:
        pool = pool_from_path(tsv_path)
        print(f"Reading {pool}: {tsv_path}", file=sys.stderr)
        counts, total = count_donor_ids(tsv_path)

        n_doublet = counts.get("doublet", 0)
        n_unassigned = counts.get("unassigned", 0)
        doublet_rate = (n_doublet / total) if total else 0.0
        unassigned_rate = (n_unassigned / total) if total else 0.0
        pool_label_map = sample_labels.get(pool, {}) or {}

        # real donors first (exclude special tokens and blanks)
        for donor_id, n_cells in counts.items():
            if donor_id in ("doublet", "unassigned", ""):
                continue
            rows.append({
                "pool": pool,
                "donor_id": donor_id,
                "mapped_label": pool_label_map.get(donor_id, ""),
                "n_cells": int(n_cells),
                "doublet_rate": round(doublet_rate, 6),
                "unassigned_rate": round(unassigned_rate, 6),
                "total_cells": total,
            })
        # then doublet / unassigned rows for completeness
        for special in ("doublet", "unassigned"):
            if counts.get(special, 0) > 0:
                rows.append({
                    "pool": pool,
                    "donor_id": special,
                    "mapped_label": "",
                    "n_cells": int(counts[special]),
                    "doublet_rate": round(doublet_rate, 6),
                    "unassigned_rate": round(unassigned_rate, 6),
                    "total_cells": total,
                })

    rows.sort(key=lambda r: (r["pool"], r["donor_id"]))

    cols = ["pool", "donor_id", "mapped_label", "n_cells",
            "doublet_rate", "unassigned_rate", "total_cells"]
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"[demux_summary] {len(rows)} rows -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
