"""Compare the locked independent human classification of the 19 evidence records.

Plan (independent_review/human_blinded_evidence_review_protocol_rev12.md, executed in R1):
  1. check the locked file: SHA-256 equals the hash recorded on receipt, 19 rows in template
     order, identity columns unchanged, controlled vocabulary, required fields present;
  2. primary comparison: human vs the primary human curator (data/evidence/evidence_labels_v2.csv)
     on strict_candidate_eligible and strict_binary_label: raw agreement with exact
     (Clopper-Pearson) 95% intervals; field-level agreement for all seven categorical fields;
  3. secondary: human vs the final adjudicated labels, and human vs the locked Claude output
     (with the B02 access-completion record), reported as a model audit;
  4. no Cohen's kappa (sparse, structurally unbalanced strict endpoint).
The locked file is never modified.

Usage (from the repository root):
  python -m analysis.human_rereview_comparison --human <human_review_locked.csv> \
      --expected-sha256 <hash recorded on receipt> --out <new_dir>
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from scipy.stats import beta

REPO = Path(__file__).resolve().parents[1]
IR = REPO / "independent_review"
FIELDS = ["evidence_domain", "evidence_design", "drug_specificity", "exposure_context", "direction",
          "strict_candidate_eligible", "strict_binary_label"]
PRIMARY_FIELD = {f: f for f in FIELDS} | {"direction": "direction_v2"}
VOCAB = {
    "evidence_domain": {"clinical", "preclinical", "class_extrapolation", "uncertain"},
    "evidence_design": {"guideline_and_clinical_standard", "prospective_single_arm_phase_2",
                        "randomized_phase_3_regimen", "retrospective_multicenter_regimen",
                        "retrospective_case_series", "multicenter_case_series", "single_pediatric_case_report",
                        "phase_1_expansion", "in_vitro_ACC_cell_lines", "no_direct_ACC_study", "other", "uncertain"},
    "drug_specificity": {"direct", "class_only", "none", "uncertain"},
    "exposure_context": {"monotherapy", "monotherapy_or_standard", "combination_regimen", "none", "uncertain"},
    "direction": {"positive", "negative", "mixed_or_limited", "regimen_effect", "positive_preclinical",
                  "contextual_case", "class_extrapolated", "no_direct_evidence", "uncertain"},
    "strict_candidate_eligible": {"yes", "no"},
    "strict_binary_label": {"positive", "negative", "NA"},
    "confidence": {"high", "medium", "low"},
    "source_accessed": {"yes", "no"},
    "access_level": {"full_text", "abstract_only", "title_only", "none"},
}
IDENTITY = ["record_id", "source_id", "drug", "doi_or_identifier"]
# final adjudication of the first comparison (results/evidence_audit/second_reviewer_agreement.md)
ADOPTED_DIRECTION = {"B11": "regimen_effect", "B14": "regimen_effect", "B15": "class_extrapolated",
                     "B19": "regimen_effect"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def norm(v: str | None) -> str:
    v = (v or "").strip()
    return "NA" if v in ("", "NA", "N/A", "na") else v


def exact_ci(k: int, n: int) -> list[float]:
    lo = 0.0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(0.975, k + 1, n - k))
    return [lo, hi]


def check_locked(human: list[dict[str, str]], template: list[dict[str, str]]) -> list[str]:
    problems = []
    if len(human) != len(template):
        problems.append(f"expected {len(template)} rows, found {len(human)}")
    for h, t in zip(human, template):
        for f in IDENTITY:
            if h.get(f, "").strip() != t[f].strip():
                problems.append(f"{t['record_id']}: identity column {f} changed")
    if human and list(human[0]) != list(template[0]):
        problems.append("column names or order differ from the template")
    return problems


def vocabulary_issues(human: list[dict[str, str]]) -> list[dict[str, str]]:
    issues = []
    for row in human:
        for f, allowed in VOCAB.items():
            v = norm(row.get(f))
            if v not in allowed and not (f == "strict_binary_label" and v == "NA"):
                issues.append({"record_id": row["record_id"], "field": f, "value": v, "issue": "not_in_vocabulary"})
        elig, lab = norm(row.get("strict_candidate_eligible")), norm(row.get("strict_binary_label"))
        if elig != "yes" and lab != "NA":
            issues.append({"record_id": row["record_id"], "field": "strict_binary_label", "value": lab,
                           "issue": "label_without_eligibility"})
        if elig == "yes" and lab == "NA":
            issues.append({"record_id": row["record_id"], "field": "strict_binary_label", "value": lab,
                           "issue": "eligible_without_binary_label"})
    return issues


def compare(human, reference, ref_field, label):
    rows, summary = [], {}
    for f in FIELDS:
        agree = 0
        for h, r in zip(human, reference):
            a, b = norm(h.get(f)), norm(r.get(ref_field(f)))
            if a == b:
                agree += 1
            else:
                rows.append({"comparison": label, "record_id": h["record_id"], "drug": h["drug"], "field": f,
                             "human": a, "reference": b})
        summary[f] = {"agree": agree, "n": len(human), "proportion": agree / len(human),
                      "exact_95ci": exact_ci(agree, len(human))}
    return summary, rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--human", type=Path, required=True)
    ap.add_argument("--expected-sha256", required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    if args.out.exists() and any(args.out.iterdir()):
        raise RuntimeError(f"refusing to overwrite non-empty {args.out}")
    digest = hashlib.sha256(args.human.read_bytes()).hexdigest()
    if digest != args.expected_sha256.strip().lower():
        raise RuntimeError(f"locked file hash {digest} differs from the hash recorded on receipt")
    template = read_csv(IR / "claude_blinded_evidence_review_v1" / "output_template.csv")
    human = read_csv(args.human)
    problems = check_locked(human, template)
    if problems:
        raise RuntimeError("locked file does not match the template: " + "; ".join(problems))
    issues = vocabulary_issues(human)

    by_drug = {r["drug"]: r for r in read_csv(REPO / "data" / "evidence" / "evidence_labels_v2.csv")}
    primary = [by_drug[h["drug"]] for h in human]
    final = [dict(p, direction_v2=ADOPTED_DIRECTION.get(h["record_id"], p["direction_v2"]))
             for h, p in zip(human, primary)]
    claude = read_csv(IR / "claude_blinded_evidence_review_v1" / "claude_review_locked.csv")
    b02 = read_csv(IR / "cabozantinib_b02_rerate" / "claude_b02_rerated.csv")[0]
    claude = [b02 if r["record_id"] == "B02" else r for r in claude]
    assert [r["record_id"] for r in claude] == [h["record_id"] for h in human]

    s_primary, d_primary = compare(human, primary, lambda f: PRIMARY_FIELD[f], "human_vs_primary_curator")
    s_final, d_final = compare(human, final, lambda f: PRIMARY_FIELD[f], "human_vs_final_adjudicated")
    s_claude, d_claude = compare(human, claude, lambda f: f, "human_vs_claude_model_audit")
    strict_positive = sorted(h["drug"] for h in human if norm(h["strict_binary_label"]) == "positive")
    strict_negative = sorted(h["drug"] for h in human if norm(h["strict_binary_label"]) == "negative")
    report = {
        "locked_file_sha256": digest, "rows": len(human), "vocabulary_or_consistency_issues": issues,
        "access": {lv: sum(norm(h.get("access_level")) == lv for h in human) for lv in ("full_text", "abstract_only", "title_only", "none")},
        "human_strict_positive": strict_positive, "human_strict_negative": strict_negative,
        "primary": s_primary, "final_adjudicated": s_final, "claude_model_audit": s_claude,
        "kappa": "not reported (protocol): sparse and structurally unbalanced strict endpoint",
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "human_rereview_agreement.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    with (args.out / "human_rereview_disagreements.csv").open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=["comparison", "record_id", "drug", "field", "human", "reference"])
        w.writeheader()
        w.writerows(d_primary + d_final + d_claude)
    lines = ["# Independent human classification of the 19 evidence records", "",
             f"Locked file SHA-256: `{digest}`. Vocabulary/consistency issues: {len(issues)}.", "",
             "| Field | vs primary curator | vs final adjudicated | vs Claude (model audit) |", "|---|---|---|---|"]
    for f in FIELDS:
        cell = lambda s: f"{s[f]['agree']}/{s[f]['n']} ({s[f]['exact_95ci'][0]:.3f}–{s[f]['exact_95ci'][1]:.3f})"
        lines.append(f"| {f} | {cell(s_primary)} | {cell(s_final)} | {cell(s_claude)} |")
    lines += ["", f"Human strict positives: {', '.join(strict_positive) or 'none'}; strict negatives: {', '.join(strict_negative) or 'none'}.",
              "", "Intervals: exact binomial 95%. Field-level disagreements: human_rereview_disagreements.csv."]
    (args.out / "human_rereview_agreement.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
