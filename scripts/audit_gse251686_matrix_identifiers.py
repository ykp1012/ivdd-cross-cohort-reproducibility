"""Audit identifier uniqueness in GSE251686 nested CeleScope matrices.

The archive is read without extraction.  Each library's `genes.tsv` is
checked as a two-column Ensembl-ID/gene-symbol table and each `barcodes.tsv`
is checked for within-library duplicate barcodes.  These checks establish the
appropriate feature key for later cross-library mapping; they do not perform
expression analysis.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import re
import tarfile
from pathlib import Path


GSM_RE = re.compile(r"^(GSM\d+)_")


def member_by_basename(inner: tarfile.TarFile, basename: str) -> tarfile.TarInfo:
    matches = [item for item in inner.getmembers() if item.isfile() and Path(item.name).name == basename]
    if len(matches) != 1:
        raise ValueError(f"Expected one {basename!r}; found {[item.name for item in matches]}")
    return matches[0]


def read_nonempty_lines(stream: io.BufferedReader) -> list[str]:
    with io.TextIOWrapper(stream, encoding="utf-8", errors="replace") as handle:
        return [line.rstrip("\n") for line in handle if line.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows: list[dict[str, object]] = []
    with tarfile.open(args.archive, "r") as outer:
        for outer_member in sorted((item for item in outer.getmembers() if item.isfile()), key=lambda item: item.name):
            gsm_match = GSM_RE.match(outer_member.name)
            if gsm_match is None:
                raise ValueError(f"Archive member lacks GSM prefix: {outer_member.name}")
            source = outer.extractfile(outer_member)
            if source is None:
                raise OSError(outer_member.name)
            with source:
                compressed = source.read()
            with gzip.GzipFile(fileobj=io.BytesIO(compressed), mode="rb") as unpacked:
                nested = unpacked.read()
            with tarfile.open(fileobj=io.BytesIO(nested), mode="r:") as inner:
                gene_stream = inner.extractfile(member_by_basename(inner, "genes.tsv"))
                barcode_stream = inner.extractfile(member_by_basename(inner, "barcodes.tsv"))
                if gene_stream is None or barcode_stream is None:
                    raise OSError(outer_member.name)
                gene_lines = read_nonempty_lines(gene_stream)
                barcode_lines = read_nonempty_lines(barcode_stream)

            parsed_genes = [line.split("\t") for line in gene_lines]
            malformed_gene_rows = sum(len(fields) != 2 for fields in parsed_genes)
            if malformed_gene_rows:
                raise ValueError(f"{outer_member.name}: {malformed_gene_rows} genes.tsv rows are not two-column")
            feature_ids = [fields[0] for fields in parsed_genes]
            gene_symbols = [fields[1] for fields in parsed_genes]
            if any(not value for value in feature_ids):
                raise ValueError(f"{outer_member.name}: empty Ensembl feature ID")
            if any(not value for value in barcode_lines):
                raise ValueError(f"{outer_member.name}: empty barcode")
            unique_feature_ids = len(set(feature_ids))
            unique_gene_symbols = len(set(gene_symbols))
            unique_barcodes = len(set(barcode_lines))
            rows.append(
                {
                    "gsm": gsm_match.group(1),
                    "outer_member": outer_member.name,
                    "genes_tsv_rows": len(gene_lines),
                    "genes_tsv_two_column_check_pass": True,
                    "unique_ensembl_feature_ids": unique_feature_ids,
                    "duplicate_ensembl_feature_ids": len(feature_ids) - unique_feature_ids,
                    "unique_gene_symbols": unique_gene_symbols,
                    "duplicate_gene_symbols": len(gene_symbols) - unique_gene_symbols,
                    "barcodes_tsv_rows": len(barcode_lines),
                    "unique_barcodes": unique_barcodes,
                    "duplicate_barcodes": len(barcode_lines) - unique_barcodes,
                    "recommended_cross_library_feature_key": "Ensembl feature ID",
                    "identifier_check_pass": unique_feature_ids == len(feature_ids) and unique_barcodes == len(barcode_lines),
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["gsm"])
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
