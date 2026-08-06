"""Create the frozen seed files for the leakage-audit arms."""

from pathlib import Path

from analysis.leakage_audit import write_frozen_arm_inputs


if __name__ == "__main__":
    write_frozen_arm_inputs(Path(__file__).resolve().parent)
