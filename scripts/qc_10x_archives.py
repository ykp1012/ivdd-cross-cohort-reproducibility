"""Stream cell-level QC metrics from GEO 10x TAR archives.

The archive is never extracted.  Each GSM is read as a barcode/features/
Matrix Market triple and summarized independently.  The output is deliberately
limited to technical QC; it does not assign cell identities or treat cells as
independent biological replicates.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import tarfile
from collections import defaultdict
from pathlib import Path
from typing import BinaryIO, Iterable


SUFFIXES = {
    "barcodes": "_barcodes.tsv.gz",
    "features": "_features.tsv.gz",
    "matrix": "_matrix.mtx.gz",
}
GSM_RE = re.compile(r"^(GSM\d+)_")


def read_tsv_lines(member: tarfile.TarInfo, tar: tarfile.TarFile) -> list[str]:
    stream = tar.extractfile(member)
    if stream is None:
        raise OSError(f"Unable to open TAR member: {member.name}")
    with stream, gzip.open(stream, "rt", encoding="utf-8", errors="replace") as handle:
        return [line.rstrip("\r\n") for line in handle if line.strip()]


def read_barcodes(member: tarfile.TarInfo, tar: tarfile.TarFile) -> list[str]:
    return [line.split("\t", 1)[0] for line in read_tsv_lines(member, tar)]


def read_features(member: tarfile.TarInfo, tar: tarfile.TarFile) -> tuple[list[str], list[bool]]:
    names: list[str] = []
    mitochondrial: list[bool] = []
    for line in read_tsv_lines(member, tar):
        fields = line.split("\t")
        if len(fields) >= 3:
            name = fields[1]
        elif len(fields) == 2:
            name = fields[1]
        elif fields:
            name = fields[0]
        else:
            continue
        names.append(name)
        mitochondrial.append(name.upper().startswith("MT-"))
    return names, mitochondrial


def parse_matrix(
    member: tarfile.TarInfo,
    tar: tarfile.TarFile,
    n_features: int,
    n_cells: int,
    mitochondrial: list[bool],
) -> tuple[list[int], list[int], list[int], int, int]:
    """Return total UMI, detected genes and mitochondrial UMI per cell.

    Cell Ranger Matrix Market files are column-oriented (gene, cell, value).
    Counts are accumulated in fixed-size arrays, so memory scales with cells,
    not with the number of non-zero entries in the matrix.
    """
    total = [0] * n_cells
    detected = [0] * n_cells
    mt_total = [0] * n_cells
    nnz = 0
    duplicate_coordinates = 0
    seen: set[tuple[int, int]] | None = None

    stream = tar.extractfile(member)
    if stream is None:
        raise OSError(f"Unable to open TAR member: {member.name}")
    with stream, gzip.open(stream, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().strip()
        if not header.startswith("%%MatrixMarket matrix coordinate"):
            raise ValueError(f"Unexpected Matrix Market header in {member.name}: {header!r}")
        line = handle.readline().strip()
        while line.startswith("%"):
            line = handle.readline().strip()
        shape = line.split()
        if len(shape) != 3:
            raise ValueError(f"Malformed Matrix Market dimensions in {member.name}: {line!r}")
        rows, columns, declared_nnz = (int(value) for value in shape)
        if (rows, columns) != (n_features, n_cells):
            raise ValueError(
                f"Dimension mismatch for {member.name}: matrix {rows}x{columns}, "
                f"features/barcodes {n_features}x{n_cells}"
            )

        for raw in handle:
            if not raw.strip() or raw.startswith("%"):
                continue
            fields = raw.split()
            if len(fields) < 3:
                raise ValueError(f"Malformed Matrix Market entry in {member.name}: {raw!r}")
            feature = int(fields[0]) - 1
            cell = int(fields[1]) - 1
            raw_value = fields[2]
            try:
                value = int(raw_value)
            except ValueError as exc:
                raise ValueError(f"Non-integer count in {member.name}: {raw_value!r}") from exc
            if not (0 <= feature < n_features and 0 <= cell < n_cells):
                raise ValueError(f"Out-of-range coordinate in {member.name}: {feature + 1}, {cell + 1}")
            if value < 0:
                raise ValueError(f"Negative count in {member.name}: {value}")
            if value == 0:
                continue
            total[cell] += value
            detected[cell] += 1
            if mitochondrial[feature]:
                mt_total[cell] += value
            nnz += 1

        if declared_nnz != nnz:
            # A zero-valued entry is legal in Matrix Market but unusual for
            # 10x.  Keep the discrepancy auditable rather than failing silently.
            duplicate_coordinates = declared_nnz - nnz
    return total, detected, mt_total, declared_nnz, duplicate_coordinates


def member_map(tar: tarfile.TarFile) -> dict[str, tarfile.TarInfo]:
    return {member.name: member for member in tar.getmembers() if member.isfile()}


def gsm_members(members: dict[str, tarfile.TarInfo], gsm: str) -> dict[str, tarfile.TarInfo]:
    found: dict[str, tarfile.TarInfo] = {}
    for label, suffix in SUFFIXES.items():
        matches = [member for name, member in members.items() if name.startswith(f"{gsm}_") and name.endswith(suffix)]
        if len(matches) != 1:
            raise ValueError(f"{gsm}: expected one {label} member, found {[m.name for m in matches]}")
        found[label] = matches[0]
    return found


def load_ledger(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"Empty ledger: {path}")
    required = {"gsm", "dataset"}
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"Ledger missing required columns: {sorted(missing)}")
    return {row["gsm"]: row for row in rows}


def write_csv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("ledger", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--min-umi", type=int, default=500)
    parser.add_argument("--max-mt-pct", type=float, default=20.0)
    parser.add_argument("--cell-output", type=Path, default=None)
    args = parser.parse_args()
    if args.min_genes < 0 or args.min_umi < 0 or args.max_mt_pct < 0:
        raise ValueError("QC thresholds must be non-negative")

    ledger = load_ledger(args.ledger)
    cell_output = args.cell_output or args.output_dir / f"{args.archive.stem}_cell_qc.tsv.gz"
    cell_output.parent.mkdir(parents=True, exist_ok=True)
    summary_rows: list[dict[str, object]] = []
    donor_counts: defaultdict[tuple[str, str, str, str], dict[str, int | str]] = defaultdict(
        lambda: {"total_cells": 0, "qc_cells": 0}
    )

    with tarfile.open(args.archive, "r") as tar, cell_output.open("wb") as raw_out:
        members = member_map(tar)
        archive_gsms = {match.group(1) for name in members if (match := GSM_RE.match(name))}
        gsms = sorted(set(ledger) & archive_gsms)
        if not gsms:
            raise ValueError(f"No ledger GSMs found in {args.archive}")
        unexpected = archive_gsms - set(gsms)
        if unexpected:
            raise ValueError(f"Archive has GSMs absent from ledger: {sorted(unexpected)}")
        with gzip.GzipFile(fileobj=raw_out, mode="wb") as compressed:
            header = [
                "dataset", "gsm", "donor_id", "compartment", "disease_state", "barcode",
                "total_umi", "detected_genes", "mt_umi", "pct_mt", "qc_pass", "qc_reason",
                "min_genes", "min_umi", "max_mt_pct",
            ]
            compressed.write(("\t".join(header) + "\n").encode("utf-8"))
            for gsm in gsms:
                row = ledger[gsm]
                parts = gsm_members(members, gsm)
                barcodes = read_barcodes(parts["barcodes"], tar)
                feature_names, mitochondrial = read_features(parts["features"], tar)
                totals, detected, mt_total, declared_nnz, nnz_difference = parse_matrix(
                    parts["matrix"], tar, len(feature_names), len(barcodes), mitochondrial
                )
                if len(mitochondrial) != len(feature_names):
                    raise AssertionError("Feature name and mitochondrial arrays differ")
                qc_cells = 0
                reason_counts: defaultdict[str, int] = defaultdict(int)
                donor_key = (
                    row.get("dataset", ""), row.get("donor_id", "") or row.get("patient_id", "") or gsm,
                    row.get("compartment", ""), row.get("disease_state", ""),
                )
                donor_counts[donor_key]["total_cells"] = int(donor_counts[donor_key]["total_cells"]) + len(barcodes)
                for index, barcode in enumerate(barcodes):
                    pct_mt = 0.0 if totals[index] == 0 else 100.0 * mt_total[index] / totals[index]
                    reasons: list[str] = []
                    if detected[index] < args.min_genes:
                        reasons.append("detected_genes_below_min")
                    if totals[index] < args.min_umi:
                        reasons.append("umi_below_min")
                    if pct_mt >= args.max_mt_pct:
                        reasons.append("mt_pct_at_or_above_max")
                    passed = not reasons
                    if passed:
                        qc_cells += 1
                    else:
                        reason_counts[";".join(reasons)] += 1
                    values = [
                        row.get("dataset", ""), gsm, donor_key[1], row.get("compartment", ""),
                        row.get("disease_state", ""), barcode, str(totals[index]), str(detected[index]),
                        str(mt_total[index]), f"{pct_mt:.6f}", "True" if passed else "False",
                        ";".join(reasons) if reasons else "pass", str(args.min_genes), str(args.min_umi),
                        f"{args.max_mt_pct:g}",
                    ]
                    compressed.write(("\t".join(values) + "\n").encode("utf-8"))
                donor_counts[donor_key]["qc_cells"] = int(donor_counts[donor_key]["qc_cells"]) + qc_cells
                summary_rows.append(
                    {
                        "dataset": row.get("dataset", ""), "gsm": gsm, "donor_id": donor_key[1],
                        "compartment": row.get("compartment", ""), "disease_state": row.get("disease_state", ""),
                        "barcodes": len(barcodes), "features": len(feature_names), "matrix_declared_nnz": declared_nnz,
                        "matrix_nonzero_entries_read": declared_nnz - nnz_difference,
                        "zero_or_duplicate_entry_difference": nnz_difference,
                        "qc_cells": qc_cells, "excluded_cells": len(barcodes) - qc_cells,
                        "qc_fraction": f"{qc_cells / len(barcodes):.6f}" if barcodes else "0",
                        "min_genes": args.min_genes, "min_umi": args.min_umi, "max_mt_pct": args.max_mt_pct,
                        "exclusion_reason_counts": ";".join(f"{key}={value}" for key, value in sorted(reason_counts.items())),
                        "status": "technical QC only; anatomical source and resident annotation pending",
                    }
                )

    donor_rows: list[dict[str, object]] = []
    for key, counts in sorted(donor_counts.items()):
        dataset, donor_id, compartment, disease_state = key
        qc_cells = int(counts["qc_cells"])
        donor_rows.append(
            {
                "dataset": dataset, "donor_id": donor_id, "compartment": compartment,
                "disease_state": disease_state, "total_cells": counts["total_cells"], "qc_cells": qc_cells,
                "resident_threshold_20_pass": qc_cells >= 20,
                "resident_threshold_30_pass": qc_cells >= 30,
                "resident_threshold_50_pass": qc_cells >= 50,
                "interpretation": "pre-annotation QC cells; not resident-cell count",
            }
        )

    fields_summary = list(summary_rows[0]) if summary_rows else ["gsm"]
    write_csv(args.output_dir / f"{args.archive.stem}_library_qc_summary.csv", summary_rows, fields_summary)
    write_csv(
        args.output_dir / f"{args.archive.stem}_donor_qc_sensitivity.csv",
        donor_rows,
        list(donor_rows[0]) if donor_rows else ["donor_id"],
    )
    write_csv(
        args.output_dir / f"{args.archive.stem}_qc_parameters.csv",
        [{
            "archive": str(args.archive), "ledger": str(args.ledger), "min_genes": args.min_genes,
            "min_umi": args.min_umi, "max_mt_pct": args.max_mt_pct,
            "cell_output": str(cell_output),
            "scope": "streaming technical QC before anatomical/resident annotation",
        }],
        ["archive", "ledger", "min_genes", "min_umi", "max_mt_pct", "cell_output", "scope"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
