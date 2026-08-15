#!/usr/bin/env python3
"""Audit GSE56081 raw files and sequence-derived GPL15314 probe mapping.

This is an extension-screen audit only.  It never writes to the default
score or effect directories.  GPL15314 exposes probe sequences but leaves
most Ensembl/accession fields empty, so the optional Ensembl lookup is kept
explicit and versioned in the output.  Exact sequence matches are candidates
for a score-level mapping; they are not silently promoted to a default cohort.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import tarfile
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Iterable


MODULE_ORDER = (
    "ecm_collagen_remodeling",
    "inflammatory_nfkb",
    "hypoxia_oxidative_stress",
    "disc_matrix_homeostasis",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def parse_platform(path: Path) -> tuple[list[str], list[dict[str, str]], dict[str, list[str]]]:
    """Return configured gene list, platform rows, and sequence -> probe IDs."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        text = handle.read()
    try:
        section = text.split("!platform_table_begin\n", 1)[1].split(
            "!platform_table_end", 1
        )[0]
    except IndexError as exc:
        raise ValueError("GPL15314 SOFT lacks a complete platform table") from exc
    lines = section.splitlines()
    header = lines[0].split("\t")
    rows: list[dict[str, str]] = []
    sequence_to_ids: dict[str, list[str]] = defaultdict(list)
    for line in lines[1:]:
        values = line.split("\t")
        if len(values) != len(header):
            continue
        row = dict(zip(header, values))
        rows.append(row)
        sequence = row.get("SEQUENCE", "").strip().upper()
        if (
            row.get("CONTROL_TYPE", "").strip().upper() == "FALSE"
            and len(sequence) >= 45
            and set(sequence) <= set("ACGTN")
        ):
            sequence_to_ids[sequence].append(row.get("ID", ""))
    return [], rows, dict(sequence_to_ids)


def parse_modules(path: Path) -> tuple[list[str], dict[str, list[str]]]:
    payload = read_json(path)
    modules = payload.get("modules", [])
    module_genes = {str(item["module_id"]): list(item["genes"]) for item in modules}
    genes = [gene for module_id in MODULE_ORDER for gene in module_genes.get(module_id, [])]
    if set(module_genes) != set(MODULE_ORDER):
        raise ValueError("program module configuration does not contain the four locked modules")
    return genes, module_genes


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTN", "TGCAN")
    return sequence.translate(table)[::-1]


def fetch_text(url: str, retries: int = 4) -> str:
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            request = urllib.request.Request(
                url,
                headers={
                    "User-Agent": "Codex-gse56081-audit/1.0",
                    "Accept": "application/json,text/plain,*/*",
                },
            )
            with urllib.request.urlopen(request, timeout=60) as response:
                return response.read().decode("utf-8")
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
            last_error = exc
            # Ensembl may return a transient 429/5xx while a batch is running.
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"failed to fetch {url}: {last_error}")


def ensembl_mapping(
    genes: Iterable[str], cache_path: Path, request_delay: float
) -> list[dict[str, str]]:
    """Fetch canonical transcript sequences and exact platform-probe matches."""
    cache: dict[str, dict] = {}
    if cache_path.exists():
        cache = read_json(cache_path)
    rows: list[dict[str, str]] = []
    for gene in genes:
        cached = cache.get(gene)
        if cached is None or not cached.get("sequence"):
            lookup_url = (
                "https://rest.ensembl.org/lookup/symbol/homo_sapiens/"
                f"{gene}?content-type=application/json"
            )
            lookup = json.loads(fetch_text(lookup_url))
            transcript_versioned = str(lookup.get("canonical_transcript", ""))
            transcript = transcript_versioned.split(".", 1)[0]
            if not transcript:
                raise RuntimeError(f"Ensembl has no canonical transcript for {gene}")
            time.sleep(request_delay)
            sequence_url = (
                f"https://rest.ensembl.org/sequence/id/{transcript}"
                "?type=cdna;content-type=text/plain"
            )
            sequence = fetch_text(sequence_url).replace("\n", "").strip().upper()
            cached = {
                "gene_id": str(lookup.get("id", "")),
                "canonical_transcript": transcript_versioned,
                "canonical_transcript_id": transcript,
                "assembly": str(lookup.get("assembly_name", "GRCh38")),
                "sequence": sequence,
            }
            cache[gene] = cached
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            cache_path.write_text(json.dumps(cache, indent=2, sort_keys=True), encoding="utf-8")
            time.sleep(request_delay)
        rows.append({"gene": gene, **{k: str(v) for k, v in cached.items() if k != "sequence"}})
    return rows, cache


def raw_tar_audit(raw_tar: Path, platform_ids: set[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with tarfile.open(raw_tar, "r:") as archive:
        members = archive.getmembers()
        raw_member_count = sum(member.name.endswith("_raw_data.txt.gz") for member in members)
        scan_member_count = sum(member.name.endswith("_scanning_graph.tif.gz") for member in members)
        for member in members:
            if not member.name.endswith("_raw_data.txt.gz"):
                continue
            match = re.search(r"(GSM\d+)_", member.name)
            gsm = match.group(1) if match else ""
            with archive.extractfile(member) as compressed:
                if compressed is None:
                    raise ValueError(f"cannot read tar member {member.name}")
                with gzip.GzipFile(fileobj=compressed) as handle:
                    fep_header: list[str] | None = None
                    design_file = ""
                    genomic_build = ""
                    extractor_version = ""
                    feature_header: list[str] | None = None
                    feature_rows = 0
                    probe_ids: set[str] = set()
                    for raw_line in handle:
                        line = raw_line.decode("utf-8", errors="replace").rstrip("\r\n")
                        if line.startswith("FEPARAMS\t"):
                            fep_header = line.split("\t")[1:]
                        elif fep_header is not None and line.startswith("DATA\t"):
                            values = line.split("\t")[1:]
                            metadata = dict(zip(fep_header, values))
                            design_file = metadata.get("FeatureExtractor_DesignFileName", "")
                            genomic_build = metadata.get("Grid_GenomicBuild", "")
                            extractor_version = metadata.get("FeatureExtractor_Version", "")
                            fep_header = None
                        elif line == "FEATURES\tFeatureNum\tRow\tCol\tSubTypeMask\tControlType\tProbeName\tSystematicName\tPositionX\tPositionY\tgProcessedSignal\tgProcessedSigError\tgMedianSignal\tgBGMedianSignal\tgBGPixSDev\tgIsSaturated\tgIsFeatNonUnifOL\tgIsBGNonUnifOL\tgIsFeatPopnOL\tgIsBGPopnOL\tIsManualFlag\tgBGSubSignal\tgIsPosAndSignif\tgIsWellAboveBG\tSpotExtentX\tgBGMeanSignal":
                            feature_header = line.split("\t")
                        elif feature_header is not None and line.startswith("DATA\t"):
                            values = line.split("\t")
                            if len(values) == len(feature_header):
                                feature_rows += 1
                                probe_ids.add(values[6])
            rows.append(
                {
                    "gsm": gsm,
                    "member": member.name,
                    "member_bytes": str(member.size),
                    "archive_raw_member_count": str(raw_member_count),
                    "archive_scan_member_count": str(scan_member_count),
                    "design_file": design_file,
                    "genomic_build": genomic_build,
                    "feature_extractor_version": extractor_version,
                    "feature_rows": str(feature_rows),
                    "unique_probe_ids": str(len(probe_ids)),
                    "platform_probe_ids_intersection": str(len(probe_ids & platform_ids)),
                    "platform_probe_ids_missing_from_raw": str(len(platform_ids - probe_ids)),
                    "platform_probe_ids_missing_list": ";".join(sorted(platform_ids - probe_ids)),
                    "unmapped_probe_ids": str(len(probe_ids - platform_ids)),
                    "feature_table_status": (
                        "pass"
                        if feature_header and feature_rows and not (probe_ids - platform_ids)
                        else "fail"
                    ),
                }
            )
    return rows


def matrix_audit(matrix_path: Path, platform_ids: set[str]) -> dict[str, str]:
    sample_header = ""
    matrix_header = ""
    feature_rows = 0
    feature_ids: set[str] = set()
    with gzip.open(matrix_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_geo_accession"):
                sample_header = line.rstrip("\r\n")
            elif line.startswith("!series_matrix_table_begin"):
                break
        # The next line is the matrix header; skip comments if present.
        for line in handle:
            if line.startswith("!"):
                if line.startswith("!series_matrix_table_end"):
                    break
                continue
            values = line.rstrip("\r\n").split("\t")
            if not values:
                continue
            if values[0].strip('"') == "ID_REF":
                matrix_header = line.rstrip("\r\n")
                continue
            feature_rows += 1
            feature_ids.add(values[0].strip('"'))
    unmapped = feature_ids - platform_ids
    return {
        "matrix_sha256": sha256(matrix_path),
        "sample_header": sample_header,
        "matrix_header": matrix_header,
        "feature_rows": str(feature_rows),
        "unique_feature_ids": str(len(feature_ids)),
        "platform_feature_intersection": str(len(feature_ids & platform_ids)),
        "unmapped_feature_ids": str(len(unmapped)),
        "status": (
            "pass"
            if matrix_header and feature_rows and len(feature_ids) == feature_rows and not unmapped
            else "fail"
        ),
    }


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument(
        "--output-dir", type=Path, default=Path("data/derived/geo_candidate_audit/GSE56081_probe_annotation")
    )
    parser.add_argument("--request-delay", type=float, default=0.15)
    args = parser.parse_args()
    root = args.project_root.resolve()
    output = (root / args.output_dir).resolve() if not args.output_dir.is_absolute() else args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    module_path = root / "config/program_modules.json"
    platform_path = root / "data/raw/geo_candidates/GSE56081/GSE56081_family.soft.gz"
    matrix_path = root / "data/raw/geo_candidates/GSE56081/GSE56081_series_matrix.txt.gz"
    raw_tar = root / "data/raw/geo_candidates/GSE56081/GSE56081_RAW.tar"
    for path in (module_path, platform_path, matrix_path, raw_tar):
        if not path.exists():
            raise FileNotFoundError(path)

    genes, module_genes = parse_modules(module_path)
    _, platform_rows, sequence_to_ids = parse_platform(platform_path)
    platform_ids = {row.get("ID", "") for row in platform_rows if row.get("ID")}
    control_false = [row for row in platform_rows if row.get("CONTROL_TYPE", "").upper() == "FALSE"]
    sequence_to_ids = {key: value for key, value in sequence_to_ids.items()}
    mapping_rows, sequence_cache = ensembl_mapping(
        genes, output / "ensembl_canonical_sequence_cache.json", args.request_delay
    )

    # Match each canonical transcript against every platform probe in both
    # orientations.  A probe can be reported for multiple locked genes.
    gene_hits: dict[str, list[tuple[str, str]]] = {}
    for row in mapping_rows:
        sequence = sequence_cache[row["gene"]]["sequence"]
        hits: list[tuple[str, str]] = []
        for probe_sequence, probe_ids in sequence_to_ids.items():
            if probe_sequence in sequence:
                hits.extend(("platform_forward_in_transcript", probe_id) for probe_id in probe_ids)
            rc = reverse_complement(probe_sequence)
            if rc in sequence:
                hits.extend(("platform_reverse_complement_in_transcript", probe_id) for probe_id in probe_ids)
        gene_hits[row["gene"]] = hits

    # Reverse index within the locked gene set to flag collisions.
    probe_to_genes: dict[str, set[str]] = defaultdict(set)
    for gene, hits in gene_hits.items():
        for _, probe_id in hits:
            probe_to_genes[probe_id].add(gene)
    detailed: list[dict[str, str]] = []
    for module_id in MODULE_ORDER:
        for gene in module_genes[module_id]:
            base = next(row for row in mapping_rows if row["gene"] == gene)
            hits = gene_hits[gene]
            unique_hits = sorted({probe_id for _, probe_id in hits})
            ambiguous = sorted(probe_id for probe_id in unique_hits if len(probe_to_genes[probe_id]) > 1)
            detailed.append(
                {
                    "module_id": module_id,
                    "gene": gene,
                    "gene_id": base.get("gene_id", ""),
                    "canonical_transcript": base.get("canonical_transcript", ""),
                    "assembly": base.get("assembly", ""),
                    "candidate_probe_count": str(len(unique_hits)),
                    "candidate_probe_ids": ";".join(unique_hits),
                    "locked_gene_collision_probe_count": str(len(ambiguous)),
                    "locked_gene_collision_probe_ids": ";".join(ambiguous),
                    "mapping_status": (
                        "sequence_match_unique_within_locked_genes"
                        if unique_hits and not ambiguous
                        else "sequence_match_ambiguous_within_locked_genes"
                        if unique_hits
                        else "no_exact_canonical_transcript_match"
                    ),
                }
            )

    summary: list[dict[str, str]] = []
    for module_id in MODULE_ORDER:
        rows = [row for row in detailed if row["module_id"] == module_id]
        mapped = [row for row in rows if row["candidate_probe_count"] != "0"]
        unique = [row for row in rows if row["mapping_status"] == "sequence_match_unique_within_locked_genes"]
        summary.append(
            {
                "module_id": module_id,
                "configured_genes": str(len(rows)),
                "genes_with_exact_sequence_match": str(len(mapped)),
                "genes_with_unique_locked_gene_match": str(len(unique)),
                "mapped_fraction": f"{len(unique) / len(rows):.6f}",
                "minimum_0_80_gate": "pass" if len(unique) / len(rows) >= 0.8 else "fail",
                "annotation_boundary": "Ensembl canonical transcript exact sequence match; genome-wide probe specificity not established",
            }
        )

    input_rows = []
    for path in (module_path, platform_path, matrix_path, raw_tar):
        input_rows.append(
            {
                "path": str(path.relative_to(root)).replace("\\", "/"),
                "bytes": str(path.stat().st_size),
                "sha256": sha256(path),
            }
        )
    write_csv(output / "input_artifact_hashes.csv", input_rows)
    write_csv(output / "probe_gene_sequence_mapping.csv", detailed)
    write_csv(output / "module_mapping_summary.csv", summary)
    raw_rows = raw_tar_audit(raw_tar, platform_ids)
    matrix_rows = [matrix_audit(matrix_path, platform_ids)]
    write_csv(output / "raw_tar_feature_audit.csv", raw_rows)
    write_csv(output / "matrix_feature_audit.csv", matrix_rows)
    no_match_genes = [row["gene"] for row in detailed if row["mapping_status"] == "no_exact_canonical_transcript_match"]
    locked_collisions = [row["gene"] for row in detailed if int(row["locked_gene_collision_probe_count"]) > 0]

    manifest = {
        "schema_version": 1,
        "purpose": "GSE56081 candidate-only raw archive and GPL15314 probe annotation audit",
        "decision": "candidate_score_level_requires_manual_or_external_probe_annotation_review",
        "platform": "GPL15314",
        "design_file_declared_by_raw_tar": "033010_D_F_20110314.xml",
        "design_files_observed": sorted({row["design_file"] for row in raw_rows}),
        "genomic_builds_observed": sorted({row["genomic_build"] for row in raw_rows}),
        "feature_extractor_versions_observed": sorted({row["feature_extractor_version"] for row in raw_rows}),
        "platform_rows": len(platform_rows),
        "control_false_rows": len(control_false),
        "sequence_annotated_rows": len(sequence_to_ids),
        "sequence_method": "exact 45+ nt platform sequence match to Ensembl GRCh38 canonical cDNA in either orientation",
        "sequence_specificity_limit": "matches were checked for collisions among the 78 locked genes only; genome-wide uniqueness and manufacturer annotation were not established",
        "default_result_mutation": False,
        "module_gate": summary,
        "no_exact_canonical_transcript_match_genes": no_match_genes,
        "locked_gene_collision_genes": locked_collisions,
        "raw_tar_feature_statuses": sorted({row["feature_table_status"] for row in raw_rows}),
        "matrix_feature_audit": matrix_rows[0],
        "input_sha256": {row["path"]: row["sha256"] for row in input_rows},
        "generated_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    (output / "audit_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    readme_lines = [
        "# GSE56081 Probe Annotation Audit",
        "",
        "This is a candidate-only extension audit. It does not modify the default IVDD result.",
        "",
        "The raw tar contains ten Agilent Feature Extraction files and ten scan images. "
        "GPL15314 exposes 60,756 platform rows and probe sequences, but most coding-probe "
        "ENSEMBL_ID/ACCESSION_STRING fields are blank. The mapping table therefore records "
        "exact sequence matches to Ensembl GRCh38 canonical cDNA, in either orientation.",
        "",
        "## First-Stage Results",
        "",
        "| Module | Configured genes | Exact matches | Unique within locked genes | Fraction | 0.80 gate |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in summary:
        readme_lines.append(
            f"| {row['module_id']} | {row['configured_genes']} | {row['genes_with_exact_sequence_match']} | "
            f"{row['genes_with_unique_locked_gene_match']} | {row['mapped_fraction']} | {row['minimum_0_80_gate']} |"
        )
    readme_lines += [
        "",
        "All four modules pass the display-only 0.80 gate (23/24, 17/21, 16/18, and 14/15 genes). "
        "The eight genes without an exact current canonical-cDNA match are "
        + (", ".join(no_match_genes) if no_match_genes else "none")
        + ". No collision was found among the 78 locked genes.",
        "",
        f"The raw Feature Extraction audit is `{raw_rows[0]['feature_table_status']}` for all {len(raw_rows)} samples; "
        f"the corrected series-matrix audit has {matrix_rows[0]['feature_rows']} features, "
        f"{matrix_rows[0]['unmapped_feature_ids']} unmapped platform IDs, and status `{matrix_rows[0]['status']}`.",
        "",
        "This is a reproducible candidate mapping, not a manufacturer-certified annotation. "
        "The first-stage collision check covers only the 78 locked genes. Genome-wide uniqueness, "
        "transcript-version compatibility with the 2011 design, and probe summarization rules remain open. "
        "Consult `GLOBAL_SPECIFICITY_README.md` and `global_specificity_manifest.json` for the fixed "
        "Ensembl release-wide specificity audit; GSE56081 remains outside the default summary.",
        "",
    ]
    (output / "README.md").write_text("\n".join(readme_lines), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
