"""Validate and hash the pre-specified IVDD program definitions."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def normalized_genes(genes: list[str]) -> list[str]:
    values = [gene.strip().upper() for gene in genes if gene.strip()]
    if len(values) != len(set(values)):
        raise ValueError("Duplicate gene symbol in a module")
    return sorted(values)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("definition", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    definition = json.loads(args.definition.read_text(encoding="utf-8"))
    rows: list[dict[str, object]] = []
    locked_at = datetime.now(timezone.utc).isoformat()
    for module in definition["modules"]:
        genes = normalized_genes(module["genes"])
        payload = "\n".join(genes).encode("utf-8")
        digest = hashlib.sha256(payload).hexdigest()
        rows.append(
            {
                "module_id": module["module_id"],
                "label": module["label"],
                "source_class": module["source_class"],
                "source_ids": ";".join(module["source_ids"]),
                "gene_count": len(genes),
                "gene_symbols_sorted": ";".join(genes),
                "gene_list_sha256": digest,
                "score_direction": definition["score_direction"],
                "locked_at_utc": locked_at,
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
