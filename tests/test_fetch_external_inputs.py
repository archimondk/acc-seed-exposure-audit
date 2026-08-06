import hashlib
import io
from pathlib import Path

import pytest

from scripts.fetch_external_inputs import ExternalFile, ensure_external_file


def test_fetcher_writes_and_verifies_download(tmp_path: Path) -> None:
    payload = b"locked external input"
    spec = ExternalFile(
        filename="example.bin",
        url="https://example.invalid/example.bin",
        sha256=hashlib.sha256(payload).hexdigest(),
        bytes=len(payload),
    )
    result = ensure_external_file(spec, tmp_path, opener=lambda _: io.BytesIO(payload))
    assert result["status"] == "downloaded"
    assert (tmp_path / spec.filename).read_bytes() == payload


def test_fetcher_refuses_mismatched_existing_file(tmp_path: Path) -> None:
    target = tmp_path / "example.bin"
    target.write_bytes(b"wrong")
    spec = ExternalFile(
        filename=target.name,
        url="https://example.invalid/example.bin",
        sha256=hashlib.sha256(b"right").hexdigest(),
        bytes=5,
    )
    with pytest.raises(ValueError, match="existing file hash mismatch"):
        ensure_external_file(spec, tmp_path, opener=lambda _: io.BytesIO(b"right"))
