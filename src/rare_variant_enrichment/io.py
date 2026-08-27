import gzip
import json
from io import TextIOWrapper
from pathlib import Path
from typing import Mapping, Sequence, TextIO


def open_text(path: Path) -> TextIO:
    raw = path.open("rb")
    if raw.read(2) == b"\x1f\x8b":
        raw.seek(0)
        return TextIOWrapper(gzip.GzipFile(fileobj=raw), encoding="utf-8")
    raw.seek(0)
    return TextIOWrapper(raw, encoding="utf-8")


def read_nonempty_lines(path: Path) -> list[str]:
    with open_text(path) as handle:
        return [line.strip() for line in handle if line.strip()]


def write_json(
    path: Path,
    payload: Mapping[str, object] | Sequence[object],
    *,
    sort_keys: bool = True,
) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=sort_keys) + "\n")
