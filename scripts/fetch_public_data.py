#!/usr/bin/env python3
"""Download and verify public source data declared in data/public_data_manifest.tsv.

The release repository intentionally excludes the large third-party GEO and
Ensembl files. This standard-library-only utility restores them at the exact
paths expected by the analysis scripts and verifies every downloaded file.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sys
from datetime import UTC, datetime
from pathlib import Path, PurePosixPath
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


REQUIRED_COLUMNS = {
    "release_groups",
    "local_path",
    "accession",
    "source",
    "url",
    "bytes",
    "sha256",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_relative_path(value: str) -> Path:
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"Manifest path is not repository-relative: {value!r}")
    return Path(*path.parts)


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Manifest has no header: {path}")
        missing = REQUIRED_COLUMNS.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"Manifest is missing columns: {', '.join(sorted(missing))}")
        records = list(reader)

    if not records:
        raise ValueError(f"Manifest has no records: {path}")
    for record in records:
        safe_relative_path(record["local_path"])
        if len(record["sha256"]) != 64:
            raise ValueError(f"Invalid SHA-256 for {record['local_path']}")
        int(record["bytes"])
    return records


def groups_for(record: dict[str, str]) -> set[str]:
    return {item.strip() for item in record["release_groups"].split(";") if item.strip()}


def select_records(
    records: list[dict[str, str]], selected_groups: set[str], include_all: bool
) -> list[dict[str, str]]:
    if include_all:
        return records
    return [record for record in records if groups_for(record).intersection(selected_groups)]


def append_retrieval_log(project_root: Path, record: dict[str, str], destination: Path) -> None:
    log_path = project_root / "data" / "raw" / "retrieval_manifest.ndjson"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "accession": record["accession"],
        "asset": record["local_path"],
        "url": record["url"],
        "retrieved_at_utc": datetime.now(UTC).isoformat(),
        "bytes": destination.stat().st_size,
        "sha256": record["sha256"],
    }
    with log_path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(payload, separators=(",", ":")) + "\n")


def verify_existing(destination: Path, record: dict[str, str]) -> tuple[bool, str]:
    if not destination.is_file():
        return False, "missing"
    expected_size = int(record["bytes"])
    if destination.stat().st_size != expected_size:
        return False, f"size {destination.stat().st_size} != {expected_size}"
    actual_hash = sha256_file(destination)
    if actual_hash != record["sha256"].lower():
        return False, "SHA-256 mismatch"
    return True, "verified"


def download(destination: Path, record: dict[str, str], timeout: int) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".partial")
    resume_at = partial.stat().st_size if partial.exists() else 0
    expected_size = int(record["bytes"])
    expected_hash = record["sha256"].lower()
    if resume_at > expected_size:
        raise RuntimeError(
            f"Partial file is larger than the expected size; remove {partial} and retry."
        )
    if resume_at == expected_size:
        if sha256_file(partial) != expected_hash:
            raise RuntimeError(
                f"Complete partial file has an invalid SHA-256; remove {partial} and retry."
            )
        partial.replace(destination)
        return

    request = Request(record["url"], headers={"User-Agent": "ivdd-cross-cohort-reproducibility/1.0"})
    if resume_at:
        request.add_header("Range", f"bytes={resume_at}-")

    try:
        response = urlopen(request, timeout=timeout)
    except (HTTPError, URLError) as error:
        raise RuntimeError(f"Download failed for {record['url']}: {error}") from error

    with response:
        status = getattr(response, "status", response.getcode())
        if resume_at and status != 206:
            raise RuntimeError(
                f"Server did not honor resume for {record['url']}; remove {partial} and retry."
            )
        if not resume_at and status not in (200, 206):
            raise RuntimeError(f"Unexpected HTTP status {status} for {record['url']}")
        mode = "ab" if resume_at else "wb"
        with partial.open(mode) as handle:
            shutil.copyfileobj(response, handle, length=1024 * 1024)

    if partial.stat().st_size != expected_size:
        raise RuntimeError(
            f"Download is incomplete ({partial.stat().st_size} of {expected_size} bytes); "
            f"partial file retained at {partial}."
        )
    if sha256_file(partial) != expected_hash:
        raise RuntimeError(f"Downloaded SHA-256 does not match; partial file retained at {partial}.")
    partial.replace(destination)


def format_target(project_root: Path, record: dict[str, str]) -> Path:
    return project_root / safe_relative_path(record["local_path"])


def parse_args() -> argparse.Namespace:
    project_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=project_root / "data" / "public_data_manifest.tsv",
        help="TSV manifest to use (default: data/public_data_manifest.tsv).",
    )
    parser.add_argument(
        "--group",
        action="append",
        default=[],
        help="Release group to fetch; may be repeated. Default: default.",
    )
    parser.add_argument("--all", action="store_true", help="Select every manifest record.")
    parser.add_argument("--dry-run", action="store_true", help="Print selected records without reading or downloading files.")
    parser.add_argument("--verify-only", action="store_true", help="Verify selected local files without downloading.")
    parser.add_argument("--timeout", type=int, default=120, help="HTTP timeout in seconds (default: 120).")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.timeout <= 0:
        raise ValueError("--timeout must be positive")

    project_root = Path(__file__).resolve().parents[1]
    manifest = args.manifest if args.manifest.is_absolute() else project_root / args.manifest
    records = read_manifest(manifest)
    available_groups = {group for record in records for group in groups_for(record)}
    requested_groups = set(args.group) or {"default"}
    unknown_groups = requested_groups.difference(available_groups)
    if unknown_groups and not args.all:
        raise ValueError(
            f"Unknown group(s): {', '.join(sorted(unknown_groups))}. "
            f"Available groups: {', '.join(sorted(available_groups))}"
        )

    selected = select_records(records, requested_groups, args.all)
    if not selected:
        raise ValueError("No manifest records selected")

    print(f"Selected {len(selected)} source asset(s).")
    if args.dry_run:
        for record in selected:
            print(f"{record['local_path']} <- {record['url']}")
        return 0

    failures = 0
    for record in selected:
        destination = format_target(project_root, record)
        valid, status = verify_existing(destination, record)
        if valid:
            print(f"[verified] {record['local_path']}")
            continue
        if args.verify_only:
            failures += 1
            print(f"[failed] {record['local_path']}: {status}", file=sys.stderr)
            continue
        if destination.exists():
            failures += 1
            print(
                f"[failed] {record['local_path']}: existing file is invalid ({status}); "
                "it was not overwritten.",
                file=sys.stderr,
            )
            continue
        print(f"[download] {record['local_path']}")
        try:
            download(destination, record, args.timeout)
            valid, status = verify_existing(destination, record)
            if not valid:
                failures += 1
                print(f"[failed] {record['local_path']}: {status}", file=sys.stderr)
                continue
            append_retrieval_log(project_root, record, destination)
            print(f"[verified] {record['local_path']}")
        except RuntimeError as error:
            failures += 1
            print(f"[failed] {record['local_path']}: {error}", file=sys.stderr)

    return 1 if failures else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)
