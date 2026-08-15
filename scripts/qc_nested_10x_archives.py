"""Stream QC for GSE251686-style nested per-GSM 10x-like archives."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import tarfile
from collections import defaultdict
from pathlib import Path


def read_lines(inner: tarfile.TarFile, name: str) -> list[str]:
    stream = inner.extractfile(name)
    if stream is None:
        raise OSError(name)
    with stream:
        return [line.decode("utf-8", errors="replace").rstrip("\r\n") for line in stream if line.strip()]


def nested_members(outer: tarfile.TarFile, member: tarfile.TarInfo) -> tarfile.TarFile:
    stream = outer.extractfile(member)
    if stream is None:
        raise OSError(member.name)
    with stream:
        payload = gzip.decompress(stream.read())
    return tarfile.open(fileobj=io.BytesIO(payload), mode="r:")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("ledger", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--min-umi", type=int, default=500)
    parser.add_argument("--max-mt-pct", type=float, default=20.0)
    args = parser.parse_args()
    with args.ledger.open(newline="", encoding="utf-8") as handle:
        ledger = {row["gsm"]: row for row in csv.DictReader(handle)}
    cell_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    donor_counts: defaultdict[tuple[str, str, str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    with tarfile.open(args.archive, "r") as outer:
        for member in sorted((m for m in outer.getmembers() if m.isfile()), key=lambda m: m.name):
            gsm = member.name.split("_", 1)[0]
            if gsm not in ledger:
                raise ValueError(f"{gsm} absent from ledger")
            inner = nested_members(outer, member)
            with inner:
                names = {Path(m.name).name: m.name for m in inner.getmembers() if m.isfile()}
                genes = [line.split("\t")[-1] for line in read_lines(inner, names["genes.tsv"])]
                barcodes = [line.split("\t", 1)[0] for line in read_lines(inner, names["barcodes.tsv"])]
                mt = [gene.upper().startswith("MT-") for gene in genes]
                total = [0] * len(barcodes)
                detected = [0] * len(barcodes)
                mt_total = [0] * len(barcodes)
                matrix_stream = inner.extractfile(names["matrix.mtx"])
                if matrix_stream is None:
                    raise OSError(names["matrix.mtx"])
                with matrix_stream, io.TextIOWrapper(matrix_stream, encoding="utf-8", errors="replace") as handle:
                    first = handle.readline().strip()
                    if not first.startswith("%%MatrixMarket matrix coordinate"):
                        raise ValueError(f"Unexpected matrix header: {member.name}")
                    line = handle.readline().strip()
                    while line.startswith("%"):
                        line = handle.readline().strip()
                    rows, cols, _ = (int(value) for value in line.split())
                    if (rows, cols) != (len(genes), len(barcodes)):
                        raise ValueError(f"Dimension mismatch: {member.name}")
                    for raw in handle:
                        if not raw.strip() or raw.startswith("%"):
                            continue
                        feature, cell, value = (int(value) for value in raw.split()[:3])
                        feature -= 1
                        cell -= 1
                        if value <= 0:
                            continue
                        total[cell] += value
                        detected[cell] += 1
                        if mt[feature]:
                            mt_total[cell] += value
                metadata = ledger[gsm]
                donor = metadata.get("presumed_donor_or_library_key", gsm)
                key = (metadata.get("dataset", ""), donor, metadata.get("compartment", "NP"), metadata.get("severity_group", ""))
                passed = 0
                reasons_count: defaultdict[str, int] = defaultdict(int)
                for index, barcode in enumerate(barcodes):
                    pct = 0.0 if total[index] == 0 else 100.0 * mt_total[index] / total[index]
                    reasons: list[str] = []
                    if detected[index] < args.min_genes:
                        reasons.append("detected_genes_below_min")
                    if total[index] < args.min_umi:
                        reasons.append("umi_below_min")
                    if pct >= args.max_mt_pct:
                        reasons.append("mt_pct_at_or_above_max")
                    ok = not reasons
                    passed += int(ok)
                    reasons_count[";".join(reasons) if reasons else "pass"] += 1
                    cell_rows.append({
                        "dataset": metadata.get("dataset", ""), "gsm": gsm, "donor_id": donor,
                        "compartment": metadata.get("compartment", "NP"), "severity_group": metadata.get("severity_group", ""),
                        "barcode": barcode, "total_umi": total[index], "detected_genes": detected[index],
                        "mt_umi": mt_total[index], "pct_mt": f"{pct:.6f}", "qc_pass": ok,
                        "qc_reason": ";".join(reasons) if reasons else "pass",
                    })
                donor_counts[key]["total_cells"] += len(barcodes)
                donor_counts[key]["qc_cells"] += passed
                summary_rows.append({
                    "dataset": metadata.get("dataset", ""), "gsm": gsm, "donor_id": donor,
                    "compartment": metadata.get("compartment", "NP"), "severity_group": metadata.get("severity_group", ""),
                    "barcodes": len(barcodes), "features": len(genes), "qc_cells": passed,
                    "excluded_cells": len(barcodes) - passed, "qc_fraction": f"{passed / len(barcodes):.6f}",
                    "min_genes": args.min_genes, "min_umi": args.min_umi, "max_mt_pct": args.max_mt_pct,
                    "exclusion_reason_counts": ";".join(f"{k}={v}" for k, v in sorted(reasons_count.items())),
                    "status": "technical QC only; source is GEO NP",
                })
    args.output_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.archive.stem
    with gzip.open(args.output_dir / f"{prefix}_cell_qc.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(cell_rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(cell_rows)
    with (args.output_dir / f"{prefix}_library_qc_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0])); writer.writeheader(); writer.writerows(summary_rows)
    donor_rows = []
    for (dataset, donor, compartment, severity), counts in sorted(donor_counts.items()):
        qc = counts["qc_cells"]
        donor_rows.append({"dataset": dataset, "donor_id": donor, "compartment": compartment, "severity_group": severity,
                           "total_cells": counts["total_cells"], "qc_cells": qc,
                           "qc_threshold_20_pass": qc >= 20, "qc_threshold_30_pass": qc >= 30, "qc_threshold_50_pass": qc >= 50,
                           "interpretation": "pre-annotation QC cells; presumed sample/library key"})
    with (args.output_dir / f"{prefix}_donor_qc_sensitivity.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(donor_rows[0])); writer.writeheader(); writer.writerows(donor_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
