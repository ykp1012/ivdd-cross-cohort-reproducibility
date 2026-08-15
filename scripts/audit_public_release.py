#!/usr/bin/env python3
"""Audit a staged public release before it is committed or archived."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


REQUIRED_PATHS = (
    ".gitattributes",
    ".gitignore",
    "README.md",
    "CHANGELOG.md",
    "REPRODUCIBILITY.md",
    "data/DATA_ACCESS.md",
    "data/public_data_manifest.tsv",
    "scripts/fetch_public_data.py",
    "tools/python/requirements-lock.txt",
    "tools/r/requirements.tsv",
)
REQUIRED_RELEASE_METADATA = ("LICENSE", "CITATION.cff", ".zenodo.json")
FORBIDDEN_PATH_PREFIXES = (
    Path("data/raw"),
    Path("tools/python/venv"),
    Path("r_library"),
    Path("logs"),
    Path("sources"),
)
FORBIDDEN_DIRECTORY_NAMES = {"__pycache__", ".pytest_cache", ".git", ".venv", "venv"}
FORBIDDEN_SUFFIXES = {".partial", ".pyc", ".pyo", ".log"}
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
PATH_PATTERNS = (
    re.compile(r"[A-Za-z]:\\\\(?:Users|home)\\\\", re.IGNORECASE),
    re.compile(r"[A-Za-z]:\\\\\\\\(?:Users|home)\\\\\\\\", re.IGNORECASE),
    re.compile(r"[A-Za-z]:/(?:Users|home)/", re.IGNORECASE),
    re.compile(r"/(?:Users|home)/"),
    re.compile("file:" + "/" * 3),
)
SECRET_PATTERNS = (
    re.compile(r"\bghp_[A-Za-z0-9]{30,}\b"),
    re.compile(r"\bgithub_pat_[A-Za-z0-9_]{20,}\b"),
    re.compile(r"\bsk-[A-Za-z0-9]{20,}\b"),
    re.compile(r"\bAKIA[0-9A-Z]{16}\b"),
)
MAX_GITHUB_FILE_BYTES = 100 * 1024 * 1024


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit a staged public release directory.")
    parser.add_argument(
        "release_root",
        nargs="?",
        type=Path,
        default=Path.cwd(),
        help="Release root to audit (default: current directory).",
    )
    parser.add_argument(
        "--require-release-metadata",
        action="store_true",
        help="Require completed LICENSE, CITATION.cff, and .zenodo.json files.",
    )
    return parser.parse_args()


def is_forbidden_path(relative: Path) -> bool:
    if any(part in FORBIDDEN_DIRECTORY_NAMES for part in relative.parts):
        return True
    if any(
        relative == prefix or prefix in relative.parents
        for prefix in FORBIDDEN_PATH_PREFIXES
    ):
        return True
    if relative.suffix.lower() in FORBIDDEN_SUFFIXES:
        return True
    return relative.name.endswith(".exporting.pdf")


def read_text_if_applicable(path: Path) -> str | None:
    if path.suffix.lower() not in TEXT_SUFFIXES:
        return None
    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return None


def main() -> int:
    args = parse_args()
    release_root = args.release_root.resolve()
    if not release_root.is_dir():
        print(f"Release root does not exist: {release_root}", file=sys.stderr)
        return 2

    errors: list[str] = []
    required_paths = list(REQUIRED_PATHS)
    if args.require_release_metadata:
        required_paths.extend(REQUIRED_RELEASE_METADATA)
    for relative in required_paths:
        if not (release_root / relative).is_file():
            errors.append(f"Missing required file: {relative}")

    files = sorted(path for path in release_root.rglob("*") if path.is_file())
    public_files: list[Path] = []
    for path in files:
        relative = path.relative_to(release_root)
        if ".git" in relative.parts:
            continue
        public_files.append(path)
        if is_forbidden_path(relative):
            errors.append(f"Forbidden release path: {relative.as_posix()}")
        if path.stat().st_size >= MAX_GITHUB_FILE_BYTES:
            errors.append(
                f"GitHub file-size limit exceeded: {relative.as_posix()} ({path.stat().st_size} bytes)"
            )

        content = read_text_if_applicable(path)
        if content is None:
            continue
        if relative.name not in {"CITATION.cff.template", ".zenodo.json.template"}:
            for pattern in PATH_PATTERNS:
                if pattern.search(content):
                    errors.append(f"Machine-specific path in: {relative.as_posix()}")
                    break
        for pattern in SECRET_PATTERNS:
            if pattern.search(content):
                errors.append(f"Possible credential in: {relative.as_posix()}")
                break

    if args.require_release_metadata:
        for relative in REQUIRED_RELEASE_METADATA:
            path = release_root / relative
            if path.is_file() and "REPLACE_WITH_" in path.read_text(encoding="utf-8"):
                errors.append(f"Unresolved metadata placeholder in: {relative}")

    if errors:
        print("Public release audit failed:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 1

    total_bytes = sum(path.stat().st_size for path in public_files)
    print(f"Public release audit passed: {len(public_files)} files, {total_bytes} bytes.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
