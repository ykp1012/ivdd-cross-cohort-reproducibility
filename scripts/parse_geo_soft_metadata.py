"""Parse selected GEO SOFT sample metadata into a long-form, auditable ledger.

The parser is deliberately conservative: it retains raw characteristic text,
normalizes only known key/value fields, and never invents patient metadata.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path


SAMPLE_START = re.compile(r"^\^SAMPLE\s*=\s*(?P<gsm>\S+)")
FIELD = re.compile(r"^!(?P<field>Sample_[^=]+)\s*=\s*(?P<value>.*)$")


def normalize_key(value: str) -> str:
    return re.sub(r"\s+", " ", value.strip().lower()).rstrip(":")


def split_characteristic(value: str) -> tuple[str, str] | None:
    if ":" not in value:
        return None
    key, content = value.split(":", maxsplit=1)
    return normalize_key(key), content.strip()


def get_first(values: dict[str, list[str]], key: str) -> str:
    return values.get(key, [""])[0]


def grade_from_title(title: str) -> str:
    match = re.search(r"Thompson\s+grade\s+([A-Za-z0-9-]+)", title, flags=re.IGNORECASE)
    return match.group(1) if match else ""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("soft", type=Path)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    current: dict[str, object] | None = None
    samples: list[dict[str, object]] = []
    with gzip.open(args.soft, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            start = SAMPLE_START.match(line)
            if start:
                if current is not None:
                    samples.append(current)
                current = {"gsm": start.group("gsm"), "fields": {}}
                continue
            if current is None:
                continue
            match = FIELD.match(line)
            if not match:
                continue
            fields = current["fields"]
            assert isinstance(fields, dict)
            fields.setdefault(match.group("field").strip(), []).append(match.group("value"))
    if current is not None:
        samples.append(current)

    rows: list[dict[str, str]] = []
    for sample in samples:
        fields = sample["fields"]
        assert isinstance(fields, dict)
        characteristics: dict[str, list[str]] = {}
        for value in fields.get("Sample_characteristics_ch1", []):
            parsed = split_characteristic(value)
            if parsed is not None:
                characteristics.setdefault(parsed[0], []).append(parsed[1])
        title = get_first(fields, "Sample_title")
        rows.append(
            {
                "dataset": args.dataset,
                "gsm": str(sample["gsm"]),
                "sample_title": title,
                "source_name": get_first(fields, "Sample_source_name_ch1"),
                "tissue": get_first(characteristics, "tissue"),
                "disease_state": get_first(characteristics, "disease state") or get_first(characteristics, "status"),
                "patient_id": get_first(characteristics, "patient id"),
                "age": get_first(characteristics, "age"),
                "sex": get_first(characteristics, "sex") or get_first(characteristics, "gender"),
                "disc": get_first(characteristics, "disc"),
                "degeneration_grade": (
                    get_first(characteristics, "thompson grade")
                    or get_first(characteristics, "grade")
                    or grade_from_title(title)
                ),
                "sra": next((value.removeprefix("SRA: ") for value in fields.get("Sample_relation", []) if value.startswith("SRA:")), ""),
                "data_processing": " | ".join(fields.get("Sample_data_processing", [])),
                "supplementary_format": " | ".join(
                    value for value in fields.get("Sample_data_processing", []) if "Supplementary files format" in value
                ),
                "raw_characteristics": " | ".join(fields.get("Sample_characteristics_ch1", [])),
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["dataset", "gsm"])
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
