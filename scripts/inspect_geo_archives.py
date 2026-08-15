"""Inspect GEO raw archives without modifying them.

The output is an auditable file inventory that informs later parsing. It does
not extract the archive and therefore preserves the raw-data boundary.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import tarfile
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    archive = args.archive.resolve()
    if not archive.is_file():
        raise FileNotFoundError(archive)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r") as handle, args.output.open("w", newline="", encoding="utf-8") as out:
        fields = ["archive", "archive_sha256", "archive_bytes", "inspected_at_utc", "member", "bytes", "is_file"]
        writer = csv.DictWriter(out, fieldnames=fields)
        writer.writeheader()
        archive_hash = sha256(archive)
        inspected_at = datetime.now(timezone.utc).isoformat()
        for member in handle.getmembers():
            writer.writerow(
                {
                    "archive": archive.name,
                    "archive_sha256": archive_hash,
                    "archive_bytes": archive.stat().st_size,
                    "inspected_at_utc": inspected_at,
                    "member": member.name,
                    "bytes": member.size,
                    "is_file": member.isfile(),
                }
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
