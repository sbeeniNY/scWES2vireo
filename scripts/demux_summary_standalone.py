#!/usr/bin/env python3
"""
Standalone cohort-wide demux summary (no Snakemake).

Reads per-pool Vireo donor_ids.tsv files and config.yaml sample_labels,
writes a single TSV matching workflow/scripts/demux_summary.py output.

Example:
  python3 scripts/demux_summary_standalone.py \\
    --config config/config.yaml \\
    --pools Pool1 Pool3 \\
    --input-dir /path/to/demux \\
    --output /path/to/demux_summary.tsv
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import sys

import pandas as pd
import yaml

log = logging.getLogger("demux_summary_standalone")


def _pool_from_path(path: str) -> str:
    parts = os.path.normpath(path).split(os.sep)
    try:
        idx = parts.index("vireo")
        return parts[idx - 1]
    except ValueError:
        m = re.search(r"/demux/([^/]+)/vireo/donor_ids\.tsv$", path)
        if not m:
            raise ValueError(f"Cannot infer pool name from path: {path}")
        return m.group(1)


def aggregate(config_path: str, donor_tsv_paths: list[str], out_tsv: str) -> None:
    with open(config_path, "r") as fh:
        cfg = yaml.safe_load(fh)

    sample_labels = cfg.get("sample_labels", {}) or {}
    rows: list[dict] = []

    for tsv_path in donor_tsv_paths:
        pool = _pool_from_path(tsv_path)
        log.info("Reading %s: %s", pool, tsv_path)

        if not os.path.isfile(tsv_path):
            log.warning("Missing donor_ids.tsv (skip): %s", tsv_path)
            continue

        df = pd.read_csv(tsv_path, sep="\t")

        if "donor_id" not in df.columns:
            log.error("%s missing 'donor_id'; columns=%s", tsv_path, list(df.columns))
            continue

        total = int(len(df))
        counts = df["donor_id"].value_counts(dropna=False).to_dict()

        n_doublet = int(counts.get("doublet", 0))
        n_unassigned = int(counts.get("unassigned", 0))
        doublet_rate = (n_doublet / total) if total else 0.0
        unassigned_rate = (n_unassigned / total) if total else 0.0

        pool_label_map = sample_labels.get(pool, {}) or {}

        for donor_id, n_cells in counts.items():
            if donor_id in ("doublet", "unassigned"):
                continue
            if pd.isna(donor_id):
                continue
            rows.append(
                {
                    "pool": pool,
                    "donor_id": donor_id,
                    "mapped_label": pool_label_map.get(donor_id, ""),
                    "n_cells": int(n_cells),
                    "doublet_rate": round(doublet_rate, 6),
                    "unassigned_rate": round(unassigned_rate, 6),
                    "total_cells": total,
                }
            )

        for special in ("doublet", "unassigned"):
            if counts.get(special, 0) > 0:
                rows.append(
                    {
                        "pool": pool,
                        "donor_id": special,
                        "mapped_label": "",
                        "n_cells": int(counts[special]),
                        "doublet_rate": round(doublet_rate, 6),
                        "unassigned_rate": round(unassigned_rate, 6),
                        "total_cells": total,
                    }
                )

    summary = pd.DataFrame(
        rows,
        columns=[
            "pool",
            "donor_id",
            "mapped_label",
            "n_cells",
            "doublet_rate",
            "unassigned_rate",
            "total_cells",
        ],
    ).sort_values(["pool", "donor_id"]).reset_index(drop=True)

    os.makedirs(os.path.dirname(out_tsv) or ".", exist_ok=True)
    summary.to_csv(out_tsv, sep="\t", index=False)
    log.info("Wrote %s rows to %s", len(summary), out_tsv)
    print(f"[demux_summary_standalone] {len(summary)} rows -> {out_tsv}", file=sys.stderr)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build demux_summary.tsv from Vireo donor_ids.tsv files.")
    parser.add_argument("--config", required=True, help="Path to config.yaml (for sample_labels).")
    parser.add_argument(
        "--pools",
        nargs="+",
        required=True,
        help="Pool names (subset or full list), e.g. Pool1 Pool3 ...",
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Demux root containing <pool>/vireo/donor_ids.tsv",
    )
    parser.add_argument("--output", required=True, help="Output TSV path.")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    paths = [os.path.join(args.input_dir, p, "vireo", "donor_ids.tsv") for p in args.pools]
    aggregate(args.config, paths, args.output)


if __name__ == "__main__":
    main()
