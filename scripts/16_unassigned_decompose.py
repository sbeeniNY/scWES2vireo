#!/usr/bin/env python3
"""
Decompose Vireo "unassigned" cells into recoverable vs genuinely ambiguous.

Motivation
----------
Vireo assigns a cell to a donor only when the best-singlet posterior >= 0.9
(internal cutoff); otherwise the cell is labelled "unassigned". A high
unassigned rate can therefore mean two very different things:

  1. data_limited   - the cell has too few informative SNPs covered, so no
                      method could confidently assign it (honest "unassigned").
  2. recoverable    - the cell HAS enough SNPs and a fairly confident best
                      singlet, but sits just under the 0.9 cutoff. Relaxing the
                      threshold (or relaxing cellSNP coverage so it gets more
                      SNPs) would recover it.
  3. true_ambiguous - enough SNPs but the posterior is genuinely split between
                      donors -> likely a doublet / high ambient RNA, correctly
                      left unassigned.

This script reads each pool's `vireo/donor_ids.tsv`, classifies every
unassigned cell, and writes per-pool + cohort summaries plus a per-cell table.
It is NON-DESTRUCTIVE: it never edits Vireo output. It additionally reports how
many cells WOULD be recovered (relabelled to their best singlet) if the
threshold were relaxed, and writes that proposed relabelling as a separate
column so you can decide whether to adopt it.

Vireo donor_ids.tsv columns used:
  cell, donor_id, prob_max, prob_doublet, n_vars, best_singlet
(missing columns are handled gracefully).

Usage
-----
  python scripts/16_unassigned_decompose.py \
      --demux-dir <OUTDIR>/demux \
      --out-dir   <OUTDIR>/demux/qc_unassigned \
      --min-vars 10 --recover-prob 0.7 --assign-thr 0.9

  # or restrict to specific pools
  python scripts/16_unassigned_decompose.py --demux-dir <OUTDIR>/demux \
      --pools Pool1 Pool8 Pool9
"""

import argparse
import csv
import glob
import os
import sys


def find_pools(demux_dir, pools):
    """Return [(pool, donor_ids_path), ...] for pools that have Vireo output."""
    found = []
    if pools:
        candidates = pools
    else:
        candidates = sorted(
            os.path.basename(os.path.dirname(p))
            for p in glob.glob(os.path.join(demux_dir, "*", "vireo"))
        )
    for pool in candidates:
        path = os.path.join(demux_dir, pool, "vireo", "donor_ids.tsv")
        if os.path.isfile(path):
            found.append((pool, path))
        else:
            print(f"WARNING: no donor_ids.tsv for {pool} ({path})", file=sys.stderr)
    return found


def read_donor_ids(path):
    """Yield dict rows from a Vireo donor_ids.tsv."""
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            yield row


def to_float(x, default=float("nan")):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def to_int(x, default=0):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return default


def classify(row, min_vars, recover_prob):
    """Classify a single unassigned cell. Returns (class, proposed_donor)."""
    n_vars = to_int(row.get("n_vars"))
    prob_max = to_float(row.get("prob_max"))
    prob_doublet = to_float(row.get("prob_doublet"), 0.0)
    best = row.get("best_singlet") or row.get("best_doublet") or ""

    if n_vars < min_vars:
        return "data_limited", ""
    # enough SNPs from here on
    if prob_doublet >= recover_prob:
        return "true_ambiguous", ""
    if prob_max >= recover_prob:
        return "recoverable", best
    return "true_ambiguous", ""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--demux-dir", required=True,
                    help="<OUTDIR>/demux (contains <pool>/vireo/donor_ids.tsv)")
    ap.add_argument("--out-dir", default=None,
                    help="output dir (default: <demux-dir>/qc_unassigned)")
    ap.add_argument("--pools", nargs="*", default=None,
                    help="subset of pools (default: auto-discover)")
    ap.add_argument("--min-vars", type=int, default=10,
                    help="min informative SNPs for a cell to be assignable [10]")
    ap.add_argument("--recover-prob", type=float, default=0.7,
                    help="best-singlet posterior to call an unassigned cell recoverable [0.7]")
    ap.add_argument("--assign-thr", type=float, default=0.9,
                    help="Vireo's assignment cutoff, for reporting only [0.9]")
    args = ap.parse_args()

    out_dir = args.out_dir or os.path.join(args.demux_dir, "qc_unassigned")
    os.makedirs(out_dir, exist_ok=True)

    pools = find_pools(args.demux_dir, args.pools)
    if not pools:
        print("ERROR: no pools with donor_ids.tsv found.", file=sys.stderr)
        sys.exit(1)

    classes = ["data_limited", "recoverable", "true_ambiguous"]
    per_pool = []          # summary rows
    cell_rows = []         # per-unassigned-cell detail

    for pool, path in pools:
        n_total = 0
        n_unassigned = 0
        counts = {c: 0 for c in classes}

        for row in read_donor_ids(path):
            n_total += 1
            if (row.get("donor_id") or "").strip() != "unassigned":
                continue
            n_unassigned += 1
            cls, proposed = classify(row, args.min_vars, args.recover_prob)
            counts[cls] += 1
            cell_rows.append({
                "pool": pool,
                "cell": row.get("cell", ""),
                "n_vars": to_int(row.get("n_vars")),
                "prob_max": row.get("prob_max", ""),
                "prob_doublet": row.get("prob_doublet", ""),
                "best_singlet": row.get("best_singlet", ""),
                "class": cls,
                "proposed_donor": proposed,
            })

        rec = counts["recoverable"]
        row_out = {
            "pool": pool,
            "n_cells": n_total,
            "n_unassigned": n_unassigned,
            "unassigned_pct": round(100.0 * n_unassigned / n_total, 2) if n_total else 0.0,
            "data_limited": counts["data_limited"],
            "recoverable": rec,
            "true_ambiguous": counts["true_ambiguous"],
            "recoverable_pct_of_unassigned":
                round(100.0 * rec / n_unassigned, 2) if n_unassigned else 0.0,
            "recoverable_pct_of_cells":
                round(100.0 * rec / n_total, 2) if n_total else 0.0,
        }
        per_pool.append(row_out)

    # cohort totals
    tot = {"pool": "TOTAL"}
    for k in ["n_cells", "n_unassigned", "data_limited", "recoverable", "true_ambiguous"]:
        tot[k] = sum(r[k] for r in per_pool)
    tot["unassigned_pct"] = round(100.0 * tot["n_unassigned"] / tot["n_cells"], 2) if tot["n_cells"] else 0.0
    tot["recoverable_pct_of_unassigned"] = round(100.0 * tot["recoverable"] / tot["n_unassigned"], 2) if tot["n_unassigned"] else 0.0
    tot["recoverable_pct_of_cells"] = round(100.0 * tot["recoverable"] / tot["n_cells"], 2) if tot["n_cells"] else 0.0
    per_pool.append(tot)

    summary_path = os.path.join(out_dir, "unassigned_decomposition_per_pool.tsv")
    cols = ["pool", "n_cells", "n_unassigned", "unassigned_pct",
            "data_limited", "recoverable", "true_ambiguous",
            "recoverable_pct_of_unassigned", "recoverable_pct_of_cells"]
    with open(summary_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(per_pool)

    cells_path = os.path.join(out_dir, "unassigned_decomposition_cells.tsv")
    cell_cols = ["pool", "cell", "n_vars", "prob_max", "prob_doublet",
                 "best_singlet", "class", "proposed_donor"]
    with open(cells_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cell_cols, delimiter="\t")
        w.writeheader()
        w.writerows(cell_rows)

    # console report
    print("=" * 64)
    print("Unassigned decomposition")
    print(f"  thresholds: min_vars={args.min_vars}  recover_prob={args.recover_prob}"
          f"  (Vireo assign_thr={args.assign_thr})")
    print("=" * 64)
    hdr = f"{'pool':10s} {'cells':>8s} {'unassgn':>8s} {'data_lim':>9s} {'recover':>8s} {'true_amb':>9s}"
    print(hdr)
    for r in per_pool:
        print(f"{r['pool']:10s} {r['n_cells']:8d} {r['n_unassigned']:8d} "
              f"{r['data_limited']:9d} {r['recoverable']:8d} {r['true_ambiguous']:9d}")
    print("-" * 64)
    print(f"Of {tot['n_unassigned']} unassigned cells: "
          f"{tot['recoverable']} recoverable "
          f"({tot['recoverable_pct_of_unassigned']}%), "
          f"{tot['data_limited']} data-limited, "
          f"{tot['true_ambiguous']} truly ambiguous.")
    print(f"\nSaved: {summary_path}")
    print(f"Saved: {cells_path}")

    # optional scatter plot (prob_max vs n_vars, coloured by class) if matplotlib present
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        color = {"data_limited": "#999999", "recoverable": "#4DAF4A",
                 "true_ambiguous": "#E41A1C"}
        fig, ax = plt.subplots(figsize=(7, 5))
        for cls in classes:
            xs = [to_int(r["n_vars"]) for r in cell_rows if r["class"] == cls]
            ys = [to_float(r["prob_max"]) for r in cell_rows if r["class"] == cls]
            ax.scatter(xs, ys, s=4, alpha=0.3, c=color[cls], label=cls)
        ax.axvline(args.min_vars, ls="--", c="grey", lw=0.8)
        ax.axhline(args.recover_prob, ls="--", c="grey", lw=0.8)
        ax.set_xscale("symlog")
        ax.set_xlabel("n_vars (informative SNPs per cell, symlog)")
        ax.set_ylabel("prob_max (best-singlet posterior)")
        ax.set_title("Unassigned cells: recoverable vs ambiguous")
        ax.legend(markerscale=3, fontsize=8)
        fig.tight_layout()
        plot_path = os.path.join(out_dir, "unassigned_decomposition_scatter.png")
        fig.savefig(plot_path, dpi=150)
        print(f"Saved: {plot_path}")
    except Exception as e:  # noqa: BLE001 - plotting is best-effort
        print(f"(plot skipped: {e})", file=sys.stderr)


if __name__ == "__main__":
    main()
