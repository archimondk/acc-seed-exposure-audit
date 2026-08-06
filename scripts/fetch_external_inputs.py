"""Download and verify the two STRING v12.0 files required by the core run."""

from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, BinaryIO, Callable, Sequence


@dataclass(frozen=True)
class ExternalFile:
    filename: str
    url: str
    sha256: str
    bytes: int


STRING_INPUTS = (
    ExternalFile(
        filename="9606.protein.info.v12.0.txt.gz",
        url=(
            "https://stringdb-downloads.org/download/protein.info.v12.0/"
            "9606.protein.info.v12.0.txt.gz"
        ),
        sha256="144de4b0d98c6a7dfde6ddc2591cf88657f27b989eadff4f501450c3ed1f0f1c",
        bytes=1_970_090,
    ),
    ExternalFile(
        filename="9606.protein.links.v12.0.txt.gz",
        url=(
            "https://stringdb-downloads.org/download/protein.links.v12.0/"
            "9606.protein.links.v12.0.txt.gz"
        ),
        sha256="3e22f32572211aa341d5b4bd08d30c32e693e294603202120936872f87719d4f",
        bytes=83_164_437,
    ),
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ensure_external_file(
    spec: ExternalFile,
    destination: Path,
    opener: Callable[[str], BinaryIO] = urllib.request.urlopen,
) -> dict[str, Any]:
    destination.mkdir(parents=True, exist_ok=True)
    target = destination / spec.filename
    if target.is_file():
        actual_hash = _sha256(target)
        if actual_hash != spec.sha256:
            raise ValueError(
                f"existing file hash mismatch for {target}: "
                f"expected {spec.sha256}, found {actual_hash}"
            )
        if target.stat().st_size != spec.bytes:
            raise ValueError(
                f"existing file size mismatch for {target}: "
                f"expected {spec.bytes}, found {target.stat().st_size}"
            )
        return {"file": str(target), "status": "verified", "sha256": actual_hash}

    temporary = destination / f".{spec.filename}.part"
    if temporary.exists():
        temporary.unlink()
    try:
        with opener(spec.url) as response, temporary.open("wb") as output:
            while True:
                block = response.read(1024 * 1024)
                if not block:
                    break
                output.write(block)
        actual_size = temporary.stat().st_size
        actual_hash = _sha256(temporary)
        if actual_size != spec.bytes or actual_hash != spec.sha256:
            raise ValueError(
                f"download verification failed for {spec.filename}: "
                f"size {actual_size}/{spec.bytes}, hash {actual_hash}/{spec.sha256}"
            )
        temporary.replace(target)
    except Exception:
        if temporary.exists():
            temporary.unlink()
        raise
    return {"file": str(target), "status": "downloaded", "sha256": spec.sha256}


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    results = [
        ensure_external_file(spec, args.project_root.resolve())
        for spec in STRING_INPUTS
    ]
    print(json.dumps({"status": "ok", "files": results}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
