"""Tests for analysis.human_rereview_comparison (run from the repository root)."""
import csv
import hashlib
import json

import pytest

from analysis import human_rereview_comparison as h


def _write(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def _claude_as_human(tmp_path):
    rows = h.read_csv(h.IR / "claude_blinded_evidence_review_v1" / "claude_review_locked.csv")
    b02 = h.read_csv(h.IR / "cabozantinib_b02_rerate" / "claude_b02_rerated.csv")[0]
    rows = [b02 if r["record_id"] == "B02" else r for r in rows]
    path = tmp_path / "human_review_locked.csv"
    _write(path, rows)
    return path, hashlib.sha256(path.read_bytes()).hexdigest()


def test_reproduces_published_model_human_agreement(tmp_path):
    """Feeding the locked Claude output must reproduce the archived Claude-human agreement
    (results/evidence_audit/second_reviewer_agreement.md)."""
    path, digest = _claude_as_human(tmp_path)
    assert h.main(["--human", str(path), "--expected-sha256", digest, "--out", str(tmp_path / "out")]) == 0
    rep = json.loads((tmp_path / "out" / "human_rereview_agreement.json").read_text())
    p = rep["primary"]
    assert p["strict_candidate_eligible"]["agree"] == 19 and p["strict_binary_label"]["agree"] == 19
    assert p["drug_specificity"]["agree"] == 18 and p["direction"]["agree"] == 14
    assert rep["final_adjudicated"]["direction"]["agree"] == 18
    assert all(v["agree"] == 19 for v in rep["claude_model_audit"].values())
    assert rep["human_strict_positive"] == ["Cabozantinib", "Mitotane"]
    assert p["strict_binary_label"]["exact_95ci"][0] == pytest.approx(0.8235, abs=1e-3)


def test_refuses_changed_file(tmp_path):
    path, digest = _claude_as_human(tmp_path)
    with pytest.raises(RuntimeError, match="differs from the hash"):
        h.main(["--human", str(path), "--expected-sha256", "0" * 64, "--out", str(tmp_path / "o1")])


def test_refuses_reordered_rows(tmp_path):
    path, _ = _claude_as_human(tmp_path)
    rows = h.read_csv(path)
    rows[0], rows[1] = rows[1], rows[0]
    _write(path, rows)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    with pytest.raises(RuntimeError, match="identity column"):
        h.main(["--human", str(path), "--expected-sha256", digest, "--out", str(tmp_path / "o2")])


def test_vocabulary_and_consistency_flags():
    rows = [{"record_id": "B01", "drug": "x", "evidence_domain": "clinic", "strict_candidate_eligible": "no",
             "strict_binary_label": "positive"}]
    issues = h.vocabulary_issues(rows)
    kinds = {(i["field"], i["issue"]) for i in issues}
    assert ("evidence_domain", "not_in_vocabulary") in kinds
    assert ("strict_binary_label", "label_without_eligibility") in kinds


def test_exact_interval():
    assert h.exact_ci(19, 19) == [pytest.approx(0.8235, abs=1e-4), 1.0]
    lo, hi = h.exact_ci(0, 19)
    assert lo == 0.0 and hi == pytest.approx(0.1765, abs=1e-4)
