#!/usr/bin/env python3
"""Build a sanitized, publication-ready snapshot outside the working tree.

The project workspace can contain large third-party inputs and machine-specific
provenance. This utility copies only the public release surface, replaces the
current project-root path in text artifacts with repository-relative paths, and
writes a checksum manifest for the snapshot. It never deletes an existing
output directory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sys
from pathlib import Path


INCLUDED_PATHS = (
    ".github",
    ".gitattributes",
    ".gitignore",
    ".zenodo.json",
    ".zenodo.json.template",
    "CHANGELOG.md",
    "CITATION.cff",
    "CITATION.cff.template",
    "LICENSE",
    "README.md",
    "RELEASE_METADATA.md",
    "REPRODUCIBILITY.md",
    "config",
    "data/DATA_ACCESS.md",
    "data/public_data_manifest.tsv",
    "data/derived",
    "docs",
    "ivdd_results_viewer.html",
    "manuscript",
    "results",
    "scripts",
    "tools/python/requirements-lock.txt",
    "tools/r",
)

EXCLUDED_PREFIXES = (
    Path("data/derived/_tmp_GSE251686_soft.csv"),
    Path("data/derived/donor_module_effect_summary_canonical_audit"),
    Path("data/derived/donor_module_effect_summary/resolved_contrast_spec_before_GSE165722.csv"),
    Path("data/derived/geo_candidate_audit/GSE56081_probe_annotation/tmp_transcript_hits.csv"),
    Path("data/derived/geo_candidate_audit/GSE56081_probe_annotation/ensembl_canonical_sequence_cache.json"),
    Path("data/derived/np_post_hoc_external_expansion_meta_analysis_rerun_20260814"),
    Path("manuscript/formal_submission/pdf_export.log"),
    Path("manuscript/05_revision_notes.md"),
    Path("results/reproducibility_rerun_20260814"),
)

EXCLUDED_DIRECTORY_NAMES = {"__pycache__", ".pytest_cache", ".git", ".venv", "venv"}
EXCLUDED_SUFFIXES = {".partial", ".pyc", ".pyo", ".log"}
TEXT_SUFFIXES = {
    ".bib",
    ".cff",
    ".csv",
    ".html",
    ".json",
    ".md",
    ".ps1",
    ".py",
    ".r",
    ".rst",
    ".toml",
    ".tsv",
    ".txt",
    ".yaml",
    ".yml",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def should_exclude(relative_path: Path) -> bool:
    if any(part in EXCLUDED_DIRECTORY_NAMES for part in relative_path.parts):
        return True
    if any(
        relative_path == prefix or prefix in relative_path.parents
        for prefix in EXCLUDED_PREFIXES
    ):
        return True
    if relative_path.suffix.lower() in EXCLUDED_SUFFIXES:
        return True
    return relative_path.name.endswith(".exporting.pdf")


def iter_included_files(project_root: Path) -> list[Path]:
    files: list[Path] = []
    for relative in INCLUDED_PATHS:
        source = project_root / relative
        if not source.exists():
            continue
        if source.is_file():
            files.append(source)
            continue
        files.extend(path for path in source.rglob("*") if path.is_file())
    return sorted(set(files), key=lambda item: item.as_posix())


def sanitized_text(text: str, project_root: Path) -> tuple[str, bool]:
    native_root = str(project_root)
    slash_root = project_root.as_posix()
    json_root = json.dumps(native_root, ensure_ascii=True)[1:-1]
    replacements = (
        (json_root, "."),
        (native_root.replace("\\", "\\\\"), "."),
        (native_root, "."),
        (slash_root, "."),
    )
    sanitized = text
    for source, replacement in replacements:
        sanitized = sanitized.replace(source, replacement)
    return sanitized, sanitized != text


def copy_file(source: Path, destination: Path, project_root: Path) -> bool:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.suffix.lower() not in TEXT_SUFFIXES:
        shutil.copy2(source, destination)
        return False

    try:
        original = source.read_bytes().decode("utf-8")
    except UnicodeDecodeError:
        shutil.copy2(source, destination)
        return False

    normalized = original.replace("\r\n", "\n").replace("\r", "\n")
    content, sanitized = sanitized_text(normalized, project_root)
    if content != original:
        destination.write_text(content, encoding="utf-8", newline="\n")
        shutil.copystat(source, destination)
    else:
        shutil.copy2(source, destination)
    return sanitized


def write_manifest(output_root: Path, records: list[dict[str, str]]) -> None:
    manifest_path = output_root / "RELEASE_FILE_MANIFEST.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("path", "bytes", "sha256", "path_sanitized"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(records)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a public IVDD release snapshot without raw data or local environments."
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="New output directory outside the project root.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    project_root = Path(__file__).resolve().parents[1]
    output_root = args.output.resolve()

    if output_root == project_root or project_root in output_root.parents:
        raise ValueError("The release snapshot must be created outside the project root.")
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output directory: {output_root}")
    output_root.mkdir(parents=True, exist_ok=True)

    manifest_records: list[dict[str, str]] = []
    excluded_count = 0
    sanitized_count = 0
    for source in iter_included_files(project_root):
        relative = source.relative_to(project_root)
        if should_exclude(relative):
            excluded_count += 1
            continue
        destination = output_root / relative
        sanitized = copy_file(source, destination, project_root)
        sanitized_count += int(sanitized)
        manifest_records.append(
            {
                "path": relative.as_posix(),
                "bytes": str(destination.stat().st_size),
                "sha256": sha256_file(destination),
                "path_sanitized": str(sanitized).lower(),
            }
        )

    write_manifest(output_root, manifest_records)
    print(
        f"Created {len(manifest_records)} release files at {output_root}; "
        f"excluded {excluded_count} files and sanitized {sanitized_count} text artifacts."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
