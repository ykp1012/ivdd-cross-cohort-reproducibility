"""Create a transparent GSE165722 donor ledger from the audited raw archive.

The archive is read in place and never extracted.  This ledger preserves the
GEO-versus-source-publication grade discrepancy and records the analysis limit
caused by GEO's "normalized counts" description.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import tarfile
from pathlib import Path

from audit_gse165722 import SAMPLES


def count_lines(member: object) -> int:
    """Count data lines in a compressed cell-name member without extraction."""
    with gzip.open(member, "rt", encoding="utf-8", errors="replace") as handle:
        next(handle, None)
        return sum(1 for line in handle if line.strip())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    fields = [
        "dataset", "gsm", "sample", "donor_id", "compartment",
        "geo_reported_grade", "source_publication_grade",
        "source_publication_severity_group", "cells_in_supplied_matrix",
        "matrix_status", "eligible_for_raw_count_pseudobulk",
        "eligible_for_score_level_validation", "metadata_note",
    ]
    rows: list[dict[str, object]] = []
    with tarfile.open(args.archive, "r") as tar:
        for sample in SAMPLES:
            name = f"{sample.gsm}_{sample.sample}.cellname.txt.gz"
            member = tar.extractfile(name)
            if member is None:
                raise FileNotFoundError(name)
            rows.append(
                {
                    "dataset": "GSE165722",
                    "gsm": sample.gsm,
                    "sample": sample.sample,
                    "donor_id": sample.sample,
                    "compartment": "NP",
                    "geo_reported_grade": sample.geo_grade,
                    "source_publication_grade": sample.source_grade,
                    "source_publication_severity_group": sample.severity_group,
                    "cells_in_supplied_matrix": count_lines(member),
                    "matrix_status": "integer count-like; GEO calls supplementary values normalized counts",
                    "eligible_for_raw_count_pseudobulk": False,
                    "eligible_for_score_level_validation": True,
                    "metadata_note": "Sample title is used as a presumed donor-level sample ID; GEO does not expose a patient identifier, age, sex, disc level, or batch. GEO grades I-IV conflict with Tu et al. PMID 34825784 Table 1 grades II-V; source-publication group used.",
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
