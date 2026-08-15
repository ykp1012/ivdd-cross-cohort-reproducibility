#!/usr/bin/env python3
"""Audit candidate GPL15314 probes against a fixed Ensembl transcriptome/genome.

The first-stage GSE56081 audit anchors probes to the four locked modules by
exact matching against current Ensembl canonical cDNA.  This second-stage
audit repeats the search against every transcript in a fixed Ensembl release
and against the GRCh38 primary assembly.  A small dependency-free
Aho-Corasick matcher is used so overlapping exact matches are retained.

This remains an extension-screen artifact.  It never writes the default
20-effect result and it does not imply that the 2011 manufacturer annotation
has been recovered.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import time
from collections import defaultdict, deque
from pathlib import Path
from typing import Iterable, Iterator


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


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


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


def parse_platform(path: Path) -> dict[str, dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        text = handle.read()
    section = text.split("!platform_table_begin\n", 1)[1].split(
        "!platform_table_end", 1
    )[0]
    lines = section.splitlines()
    header = lines[0].split("\t")
    rows: dict[str, dict[str, str]] = {}
    for line in lines[1:]:
        values = line.split("\t")
        if len(values) == len(header):
            row = dict(zip(header, values))
            if row.get("ID"):
                rows[row["ID"]] = row
    return rows


def parse_modules(path: Path) -> dict[str, list[str]]:
    with path.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    modules = {str(item["module_id"]): list(item["genes"]) for item in payload["modules"]}
    if set(modules) != set(MODULE_ORDER):
        raise ValueError("program module configuration does not contain the four locked modules")
    return modules


def parse_probe_queries(
    mapping_path: Path, platform_path: Path
) -> tuple[list[dict[str, str]], dict[str, dict[str, str]]]:
    platform = parse_platform(platform_path)
    queries: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    with mapping_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            probe_ids = [item for item in row["candidate_probe_ids"].split(";") if item]
            for probe_id in probe_ids:
                key = (row["gene"], probe_id)
                if key in seen:
                    continue
                seen.add(key)
                platform_row = platform.get(probe_id)
                if platform_row is None:
                    raise ValueError(f"candidate probe {probe_id} is absent from GPL15314")
                sequence = platform_row.get("SEQUENCE", "").strip().upper()
                if len(sequence) < 45 or set(sequence) - set("ACGT"):
                    raise ValueError(f"candidate probe {probe_id} has no usable A/C/G/T sequence")
                queries.append(
                    {
                        "module_id": row["module_id"],
                        "gene": row["gene"],
                        "expected_gene_id": row.get("gene_id", "").split(".", 1)[0],
                        "probe_id": probe_id,
                        "sequence": sequence,
                        "reverse_complement": reverse_complement(sequence),
                        "length": str(len(sequence)),
                    }
                )
    queries.sort(key=lambda row: (MODULE_ORDER.index(row["module_id"]), row["gene"], row["probe_id"]))
    return queries, platform


class ExactMatcher:
    """A dependency-free Aho-Corasick matcher for A/C/G/T patterns."""

    _code = [-1] * 256
    for _letter, _value in (("A", 0), ("C", 1), ("G", 2), ("T", 3)):
        _code[ord(_letter)] = _value
        _code[ord(_letter.lower())] = _value

    def __init__(self, patterns: Iterable[str]):
        self.patterns = list(patterns)
        self.next: list[list[int]] = [[-1, -1, -1, -1]]
        self.fail: list[int] = [0]
        self.outputs: list[list[int]] = [[]]
        for pattern_id, pattern in enumerate(self.patterns):
            state = 0
            for letter in pattern:
                code = self._code[ord(letter)]
                if code < 0:
                    raise ValueError(f"pattern contains unsupported base: {pattern}")
                child = self.next[state][code]
                if child < 0:
                    child = len(self.next)
                    self.next[state][code] = child
                    self.next.append([-1, -1, -1, -1])
                    self.fail.append(0)
                    self.outputs.append([])
                state = child
            self.outputs[state].append(pattern_id)

        queue: deque[int] = deque()
        for code in range(4):
            child = self.next[0][code]
            if child < 0:
                self.next[0][code] = 0
            else:
                queue.append(child)
        while queue:
            state = queue.popleft()
            failure = self.fail[state]
            if self.outputs[failure]:
                self.outputs[state].extend(self.outputs[failure])
            for code in range(4):
                child = self.next[state][code]
                if child < 0:
                    self.next[state][code] = self.next[failure][code]
                else:
                    self.fail[child] = self.next[failure][code]
                    queue.append(child)

    def scan(self, sequence: bytes) -> Iterator[tuple[int, int]]:
        state = 0
        code_table = self._code
        transitions = self.next
        outputs = self.outputs
        for position, value in enumerate(sequence):
            code = code_table[value]
            if code < 0:
                state = 0
                continue
            state = transitions[state][code]
            for pattern_id in outputs[state]:
                yield pattern_id, position - len(self.patterns[pattern_id]) + 1


def fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    header = ""
    chunks: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield header, "".join(chunks).upper()
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
        if header:
            yield header, "".join(chunks).upper()


def header_gene_id(header: str) -> str:
    match = re.search(r"(?:^|\s)gene:([^\s]+)", header)
    return match.group(1).split(".", 1)[0] if match else ""


def parse_gtf_gene_spans(path: Path) -> dict[str, list[dict[str, str]]]:
    spans: dict[str, list[dict[str, str]]] = defaultdict(list)
    attr_re = re.compile(r'([A-Za-z_]+) "([^"]*)"')
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = dict(attr_re.findall(fields[8]))
            gene_id = attrs.get("gene_id", "").split(".", 1)[0]
            if not gene_id:
                continue
            spans[fields[0]].append(
                {
                    "gene_id": gene_id,
                    "gene_name": attrs.get("gene_name", ""),
                    "start": fields[3],
                    "end": fields[4],
                }
            )
    for values in spans.values():
        values.sort(key=lambda row: int(row["start"]))
    return dict(spans)


def scan_transcriptome(
    paths: Iterable[Path], matcher: ExactMatcher, pattern_to_queries: dict[int, list[dict[str, str]]],
    output_path: Path,
) -> tuple[dict[str, list[dict[str, str]]], dict[str, str]]:
    hits: dict[str, list[dict[str, str]]] = defaultdict(list)
    stats = {"files": "0", "records": "0", "bases": "0", "pattern_hits": "0"}
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        fields = [
            "source_file", "record_id", "gene_id", "gene_symbol", "pattern_id",
            "probe_id", "gene", "orientation", "start_0based", "end_0based_inclusive",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for path in paths:
            stats["files"] = str(int(stats["files"]) + 1)
            for header, sequence in fasta_records(path):
                stats["records"] = str(int(stats["records"]) + 1)
                stats["bases"] = str(int(stats["bases"]) + len(sequence))
                record_id = header.split(None, 1)[0]
                gene_id = header_gene_id(header)
                symbol_match = re.search(r"(?:^|\s)gene_symbol:([^\s]+)", header)
                gene_symbol = symbol_match.group(1) if symbol_match else ""
                for pattern_id, start in matcher.scan(sequence.encode("ascii", "ignore")):
                    for query in pattern_to_queries[pattern_id]:
                        orientation = query["orientation"]
                        row = {
                            "source_file": path.name,
                            "record_id": record_id,
                            "gene_id": gene_id,
                            "gene_symbol": gene_symbol,
                            "pattern_id": str(pattern_id),
                            "probe_id": query["probe_id"],
                            "gene": query["gene"],
                            "orientation": orientation,
                            "start_0based": str(start),
                            "end_0based_inclusive": str(start + int(query["length"]) - 1),
                        }
                        writer.writerow(row)
                        hits[query["probe_id"]].append(row)
                        stats["pattern_hits"] = str(int(stats["pattern_hits"]) + 1)
    return dict(hits), stats


def scan_genome(
    path: Path, matcher: ExactMatcher, pattern_to_queries: dict[int, list[dict[str, str]]],
    output_path: Path,
) -> tuple[dict[str, list[dict[str, str]]], dict[str, str]]:
    hits: dict[str, list[dict[str, str]]] = defaultdict(list)
    stats = {"files": "1", "records": "0", "bases": "0", "pattern_hits": "0"}
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        fields = [
            "source_file", "record_id", "pattern_id", "probe_id", "gene", "orientation",
            "start_1based", "end_1based", "sequence_length",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for header, sequence in fasta_records(path):
            stats["records"] = str(int(stats["records"]) + 1)
            stats["bases"] = str(int(stats["bases"]) + len(sequence))
            record_id = header.split(None, 1)[0]
            for pattern_id, start in matcher.scan(sequence.encode("ascii", "ignore")):
                for query in pattern_to_queries[pattern_id]:
                    row = {
                        "source_file": path.name,
                        "record_id": record_id,
                        "pattern_id": str(pattern_id),
                        "probe_id": query["probe_id"],
                        "gene": query["gene"],
                        "orientation": query["orientation"],
                        "start_1based": str(start + 1),
                        "end_1based": str(start + int(query["length"])),
                        "sequence_length": query["length"],
                    }
                    writer.writerow(row)
                    hits[query["probe_id"]].append(row)
                    stats["pattern_hits"] = str(int(stats["pattern_hits"]) + 1)
    return dict(hits), stats


def reference_manifest(root: Path, reference_dir: Path) -> dict:
    urls = {
        "cdna": "https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
        "ncrna": "https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz",
        "gtf": "https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz",
        "primary_assembly": "https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
    }
    names = {
        "cdna": "Homo_sapiens.GRCh38.cdna.all.fa.gz",
        "ncrna": "Homo_sapiens.GRCh38.ncrna.fa.gz",
        "gtf": "Homo_sapiens.GRCh38.113.gtf.gz",
        "primary_assembly": "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
    }
    artifacts = {}
    for key, name in names.items():
        path = reference_dir / name
        if not path.exists():
            raise FileNotFoundError(path)
        artifacts[key] = {
            "path": str(path.relative_to(root)).replace("\\", "/"),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            "url": urls[key],
        }
    return {
        "ensembl_release": 113,
        "assembly": "GRCh38",
        "primary_assembly_scope": "Ensembl GRCh38 primary assembly FASTA; alternate loci are not searched",
        "downloaded_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "artifacts": artifacts,
    }


def classify_probes(
    queries: list[dict[str, str]],
    transcript_hits: dict[str, list[dict[str, str]]],
    genome_hits: dict[str, list[dict[str, str]]],
    gene_spans: dict[str, list[dict[str, str]]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for query in queries:
        probe_id = query["probe_id"]
        tx = transcript_hits.get(probe_id, [])
        tx_genes = sorted({row["gene_id"] for row in tx if row["gene_id"]})
        expected = query["expected_gene_id"]
        genome = genome_hits.get(probe_id, [])
        loci = sorted(
            {
                (row["record_id"], int(row["start_1based"]), int(row["end_1based"]))
                for row in genome
            }
        )
        overlap_ids: set[str] = set()
        overlap_names: set[str] = set()
        for record_id, start, end in loci:
            for span in overlapping_genes(gene_spans, record_id, start, end):
                overlap_ids.add(span["gene_id"])
                if span["gene_name"]:
                    overlap_names.add(span["gene_name"])
        reasons: list[str] = []
        if expected not in tx_genes:
            reasons.append("no_expected_release113_transcript_match")
        if any(gene_id != expected for gene_id in tx_genes):
            reasons.append("cross_gene_transcript_match")
        if len(loci) > 1:
            reasons.append("multiple_primary_assembly_loci")
        if overlap_ids - {expected}:
            reasons.append("cross_gene_primary_assembly_overlap")
        if loci and not overlap_ids:
            reasons.append("primary_assembly_hit_outside_annotated_gene")
        if not loci and expected in tx_genes:
            boundary = "transcript_only_exact_match_no_contiguous_primary_assembly_hit"
        elif loci:
            boundary = "primary_assembly_locus_checked_against_release113_gtf"
        else:
            boundary = "no_exact_reference_hit"
        status = "pass_unique_expected_gene" if not reasons else "fail_specificity_gate"
        rows.append(
            {
                "module_id": query["module_id"],
                "gene": query["gene"],
                "expected_gene_id": expected,
                "probe_id": probe_id,
                "probe_length": query["length"],
                "platform_sequence": query["sequence"],
                "reverse_complement_sequence": query["reverse_complement"],
                "release113_transcript_hit_count": str(len(tx)),
                "release113_transcript_gene_count": str(len(tx_genes)),
                "release113_transcript_gene_ids": ";".join(tx_genes),
                "primary_assembly_locus_count": str(len(loci)),
                "primary_assembly_loci": ";".join(
                    f"{chrom}:{start}-{end}" for chrom, start, end in loci
                ),
                "primary_assembly_overlap_gene_ids": ";".join(sorted(overlap_ids)),
                "primary_assembly_overlap_gene_names": ";".join(sorted(overlap_names)),
                "specificity_status": status,
                "failure_reasons": ";".join(reasons),
                "specificity_boundary": boundary,
            }
        )
    return rows


def overlapping_genes(
    spans: dict[str, list[dict[str, str]]], chromosome: str, start: int, end: int
) -> list[dict[str, str]]:
    # Candidate hit counts are small; a linear scan is transparent and avoids
    # hiding interval semantics in an external package.
    return [
        row
        for row in spans.get(chromosome, [])
        if int(row["start"]) <= end and int(row["end"]) >= start
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/derived/geo_candidate_audit/GSE56081_probe_annotation"),
    )
    args = parser.parse_args()
    root = args.project_root.resolve()
    output = (root / args.output_dir).resolve() if not args.output_dir.is_absolute() else args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    reference_dir = root / "data/raw/reference/ensembl_release_113"
    mapping_path = output / "probe_gene_sequence_mapping.csv"
    platform_path = root / "data/raw/geo_candidates/GSE56081/GSE56081_family.soft.gz"
    module_path = root / "config/program_modules.json"
    for path in (reference_dir, mapping_path, platform_path, module_path):
        if not path.exists():
            raise FileNotFoundError(path)

    modules = parse_modules(module_path)
    queries, _ = parse_probe_queries(mapping_path, platform_path)
    # A single sequence is searched once; duplicate platform probes retain
    # separate ledger rows through pattern_to_queries.
    pattern_to_queries: dict[int, list[dict[str, str]]] = defaultdict(list)
    pattern_ids: dict[str, int] = {}
    for query in queries:
        for orientation, sequence in (
            ("platform_forward_in_reference", query["sequence"]),
            ("platform_reverse_complement_in_reference", query["reverse_complement"]),
        ):
            pattern_id = pattern_ids.setdefault(sequence, len(pattern_ids))
            item = dict(query)
            item["orientation"] = orientation
            pattern_to_queries[pattern_id].append(item)
    patterns = [""] * len(pattern_ids)
    for sequence, pattern_id in pattern_ids.items():
        patterns[pattern_id] = sequence
    matcher = ExactMatcher(patterns)
    write_csv(
        output / "probe_specificity_queries.csv",
        [
            {"pattern_id": str(pattern_ids[q["sequence"]]), **q, "orientation": "platform_forward_in_reference"}
            for q in queries
        ]
        + [
            {"pattern_id": str(pattern_ids[q["reverse_complement"]]), **q, "orientation": "platform_reverse_complement_in_reference"}
            for q in queries
        ],
    )

    transcript_paths = [
        reference_dir / "Homo_sapiens.GRCh38.cdna.all.fa.gz",
        reference_dir / "Homo_sapiens.GRCh38.ncrna.fa.gz",
    ]
    transcript_hits, transcript_stats = scan_transcriptome(
        transcript_paths, matcher, pattern_to_queries, output / "release113_transcriptome_probe_hits.csv"
    )
    genome_hits, genome_stats = scan_genome(
        reference_dir / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
        matcher,
        pattern_to_queries,
        output / "release113_primary_assembly_probe_hits.csv",
    )
    spans = parse_gtf_gene_spans(reference_dir / "Homo_sapiens.GRCh38.113.gtf.gz")
    detailed = classify_probes(queries, transcript_hits, genome_hits, spans)
    write_csv(output / "probe_global_specificity_ledger.csv", detailed)

    summary: list[dict[str, str]] = []
    for module_id in MODULE_ORDER:
        module_rows = [row for row in detailed if row["module_id"] == module_id]
        genes = modules[module_id]
        specific_genes = {
            row["gene"] for row in module_rows if row["specificity_status"] == "pass_unique_expected_gene"
        }
        failed = [row for row in module_rows if row["specificity_status"] != "pass_unique_expected_gene"]
        summary.append(
            {
                "module_id": module_id,
                "configured_genes": str(len(genes)),
                "candidate_probes": str(len(module_rows)),
                "globally_specific_probes": str(len(module_rows) - len(failed)),
                "genes_with_globally_specific_probe": str(len(specific_genes)),
                "specific_gene_fraction": f"{len(specific_genes) / len(genes):.6f}",
                "minimum_0_80_gate": "pass" if len(specific_genes) / len(genes) >= 0.8 else "fail",
                "failed_probe_count": str(len(failed)),
                "cross_gene_or_multi_locus_exclusions": str(
                    sum(
                        1
                        for row in failed
                        if "cross_gene" in row["failure_reasons"] or "multiple_primary_assembly_loci" in row["failure_reasons"]
                    )
                ),
            }
        )
    write_csv(output / "module_global_specificity_summary.csv", summary)

    ref_manifest = reference_manifest(root, reference_dir)
    (output / "ensembl_release_113_reference_manifest.json").write_text(
        json.dumps(ref_manifest, indent=2), encoding="utf-8"
    )
    all_failed = [row for row in detailed if row["specificity_status"] != "pass_unique_expected_gene"]
    all_cross = [row for row in all_failed if row["failure_reasons"]]
    decision = (
        "candidate_score_level_eligible_after_excluding_failed_probes"
        if not all_cross and all(row["minimum_0_80_gate"] == "pass" for row in summary)
        else "candidate_only_blocked_global_probe_specificity_not_resolved"
    )
    manifest = {
        "schema_version": 1,
        "purpose": "GSE56081 candidate-only whole-transcriptome and primary-assembly exact probe specificity audit",
        "decision": decision,
        "default_result_mutation": False,
        "platform": "GPL15314",
        "matching": {
            "rule": "complete exact A/C/G/T sequence matching in both platform orientation and reverse complement",
            "minimum_probe_length": min(int(row["length"]) for row in queries),
            "maximum_probe_length": max(int(row["length"]) for row in queries),
            "overlapping_matches_retained": True,
            "multi_hit_handling": "all transcript records, gene IDs, and primary-assembly loci are retained; failed probes are explicitly excluded from any exploratory score",
        },
        "reference": ref_manifest,
        "queries": {
            "candidate_probe_count": len(queries),
            "unique_search_pattern_count": len(patterns),
        },
        "transcriptome_scan": transcript_stats,
        "primary_assembly_scan": genome_stats,
        "module_gate": summary,
        "failed_probe_count": len(all_failed),
        "failed_probe_ids": [row["probe_id"] for row in all_failed],
        "generated_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    (output / "global_specificity_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    lines = [
        "# GSE56081 Global Probe Specificity Audit",
        "",
        "This is a candidate-only extension audit. It does not modify the frozen default IVDD result.",
        "",
        "GPL15314 candidate probes were searched exactly (A/C/G/T, both orientations) against every",
        "Ensembl release-113 human cDNA and ncRNA record and the GRCh38 primary assembly. The",
        "dependency-free matcher retains overlapping hits. The ledger records all transcript gene IDs",
        "and genomic loci; probes with cross-gene transcript hits, multiple genomic loci, or an",
        "unexpected genomic overlap are excluded from any exploratory score.",
        "",
        f"Reference release: Ensembl 113 / GRCh38; cDNA records={transcript_stats['records']}; "
        f"primary-assembly records={genome_stats['records']}.",
        f"Candidate probes searched: {len(queries)}; failed specificity probes: {len(all_failed)}.",
        "",
        "## Module Gate",
        "",
        "| Module | Candidate probes | Globally specific probes | Specific genes | Fraction | 0.80 gate |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in summary:
        lines.append(
            f"| {row['module_id']} | {row['candidate_probes']} | {row['globally_specific_probes']} | "
            f"{row['genes_with_globally_specific_probe']} | {row['specific_gene_fraction']} | {row['minimum_0_80_gate']} |"
        )
    lines += [
        "",
        f"Decision: `{decision}`.",
        "A pass here is still sequence-anchored evidence rather than manufacturer-certified",
        "annotation; the 2011 Arraystar design/transcript-version and probe summarization boundaries",
        "remain documented in the first-stage audit.",
    ]
    (output / "GLOBAL_SPECIFICITY_README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
