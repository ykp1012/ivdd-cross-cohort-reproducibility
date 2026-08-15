"""Merge per-archive 10x audits into an analysis-ready discovery raw-data ledger."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("audits", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = [row for audit in args.audits for row in read_csv(audit)]
    gsm_counts: dict[str, int] = {}
    for row in rows:
        gsm_counts[row["gsm"]] = gsm_counts.get(row["gsm"], 0) + 1
    duplicate_gsms = sorted(gsm for gsm, count in gsm_counts.items() if count != 1)
    if duplicate_gsms:
        raise ValueError(f"Expected exactly one audit row per GSM: {duplicate_gsms}")

    rows.sort(key=lambda row: (row["dataset"], row["gsm"]))
    for row in rows:
        row["raw_count_eligible"] = row["matrix_status"].startswith("raw 10x") and row["dimension_check_pass"] == "True"
        row["analysis_status"] = "eligible for pre-QC raw-count ingestion" if row["raw_count_eligible"] else "blocked"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["gsm"])
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
