"""Metadata helpers for files passed between workflow stages."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Sequence


def file_artifact(
    path: Path,
    logical_name: str,
    header: Sequence[str],
    row_count: int,
) -> dict[str, object]:
    """Describe the exact bytes and tabular contract of one file."""
    if not logical_name:
        raise ValueError("Artifact logical_name must not be empty")
    if not header:
        raise ValueError("Artifact header must not be empty")
    if row_count < 0:
        raise ValueError("Artifact row_count must not be negative")

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return {
        "logical_name": logical_name,
        "header": list(header),
        "row_count": row_count,
        "size_bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
    }
