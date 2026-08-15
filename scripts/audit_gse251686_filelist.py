"""Cross-check GSE251686 GEO filelist entries against the downloaded TAR.

This is an archive-level provenance check.  It compares the GEO-published
outer archive size and every nested per-GSM archive name/size without
extracting any matrix payload.
"""

from __future__ import annotations

import argparse
import csv
import tarfile
from pathlib import Path


def read_filelist(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines or lines[0] != "#Archive/File\tName\tTime\tSize\tType":
        raise ValueError(f"Unexpected GEO filelist header in {path}")
    rows: list[dict[str, str]] = []
    for line in lines[1:]:
        values = line.split("\t")
        if len(values) != 5:
            raise ValueError(f"Malformed GEO filelist row: {line!r}")
        rows.append(
            {
                "filelist_entry_type": values[0],
                "name": values[1],
                "published_time": values[2],
                "published_bytes": values[3],
                "published_type": values[4],
            }
        )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("filelist", type=Path)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = read_filelist(args.filelist)
    archive_rows = [row for row in rows if row["filelist_entry_type"] == "Archive"]
    member_rows = [row for row in rows if row["filelist_entry_type"] == "File"]
    if len(archive_rows) != 1:
        raise ValueError(f"Expected one Archive row, found {len(archive_rows)}")
    if archive_rows[0]["name"] != args.archive.name:
        raise ValueError(f"GEO archive name {archive_rows[0]['name']!r} != {args.archive.name!r}")
    if int(archive_rows[0]["published_bytes"]) != args.archive.stat().st_size:
        raise ValueError("Downloaded archive byte size does not match GEO filelist")

    with tarfile.open(args.archive, "r") as handle:
        members = {member.name: member.size for member in handle.getmembers() if member.isfile()}
    published_names = {row["name"] for row in member_rows}
    if set(members) != published_names:
        raise ValueError(
            "GEO filelist/archive member mismatch: "
            f"missing_in_archive={sorted(published_names - set(members))}, "
            f"unexpected_in_archive={sorted(set(members) - published_names)}"
        )

    output_rows: list[dict[str, object]] = [
        {
            **archive_rows[0],
            "observed_bytes": args.archive.stat().st_size,
            "exists_in_downloaded_archive": True,
            "byte_size_match": True,
        }
    ]
    for row in member_rows:
        observed = members[row["name"]]
        expected = int(row["published_bytes"])
        if observed != expected:
            raise ValueError(f"Size mismatch for {row['name']}: GEO={expected}, TAR={observed}")
        output_rows.append(
            {
                **row,
                "observed_bytes": observed,
                "exists_in_downloaded_archive": True,
                "byte_size_match": True,
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
