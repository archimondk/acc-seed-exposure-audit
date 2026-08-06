"""Validate formal leakage-audit outputs and apply frozen decision rules."""

import json
from pathlib import Path

from analysis.leakage_audit import evaluate_formal_outputs


if __name__ == "__main__":
    result = evaluate_formal_outputs(Path(__file__).resolve().parent)
    print(json.dumps(result, indent=2, ensure_ascii=False))
