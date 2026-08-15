"""Merge audited GSE230809 child-series metadata into a discovery sample ledger."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()

    rows = [row for path in args.metadata for row in read_csv(path)]
    rows.sort(key=lambda row: (row["dataset"], row["gsm"]))
    donors_by_compartment: Counter[tuple[str, str]] = Counter()
    for row in rows:
        patient_id = row["patient_id"]
        if not patient_id:
            raise ValueError(f"Missing patient id for {row['dataset']} {row['gsm']}")
        row["discovery_project"] = "GSE230809"
        row["donor_id"] = patient_id
        row["compartment"] = {"annulus fibrosus": "AF", "nucleus pulposus": "NP"}.get(row["tissue"].lower(), "")
        row["analysis_status"] = "pending raw archive audit"
        row["confounding_note"] = "age and disease state are fully confounded across the discovery contrast"
        donors_by_compartment[(patient_id, row["compartment"])] += 1

    for row in rows:
        if donors_by_compartment[(row["donor_id"], row["compartment"])] != 1:
            raise ValueError(f"Expected one library per donor/compartment: {row['gsm']}")

    donor_to_compartments: dict[str, set[str]] = {}
    for row in rows:
        donor_to_compartments.setdefault(row["donor_id"], set()).add(row["compartment"])
    for row in rows:
        row["pair_status"] = "AF_NP_paired" if donor_to_compartments[row["donor_id"]] == {"AF", "NP"} else f"{row['compartment']}_only"

    fields = list(rows[0]) if rows else ["dataset", "gsm"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    controls = {row["donor_id"] for row in rows if row["disease_state"].lower() == "healthy"}
    advanced = {row["donor_id"] for row in rows if row["disease_state"].lower() != "healthy"}
    paired = {donor for donor, compartments in donor_to_compartments.items() if compartments == {"AF", "NP"}}
    summary = [
        {"metric": "samples", "value": len(rows)},
        {"metric": "biological_donors", "value": len(donor_to_compartments)},
        {"metric": "healthy_donors", "value": len(controls)},
        {"metric": "advanced_donors", "value": len(advanced)},
        {"metric": "AF_NP_paired_donors", "value": len(paired)},
        {"metric": "healthy_paired_donors", "value": len(controls & paired)},
        {"metric": "advanced_paired_donors", "value": len(advanced & paired)},
        {"metric": "AF_libraries", "value": sum(row["compartment"] == "AF" for row in rows)},
        {"metric": "NP_libraries", "value": sum(row["compartment"] == "NP" for row in rows)},
        {"metric": "sexes", "value": ";".join(sorted({row["sex"] for row in rows if row["sex"]}))},
        {"metric": "age_disease_relation", "value": "fully confounded in discovery contrast"},
    ]
    with args.summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value"])
        writer.writeheader()
        writer.writerows(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
