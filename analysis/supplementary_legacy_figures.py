"""Regenerate exploratory Supplementary Figures S1 and S2.

Figure S1 reproduces the prespecified four-component weight grid in the locked
108-drug complete-case universe. Figure S2 reproduces the leave-one-mechanism-
class-out (LOMCO) diagnostic for the exploratory drug-neighbour component.
Neither analysis is used as a primary validation result.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 9,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
    }
)


COMPONENT_COLUMNS = (
    "C_ACC_pct",
    "ACCrelative_resid_pct",
    "S_Bneighbor_pct",
    "S_external",
)
SELECTED_DRUGS = (
    "Mitotane",
    "Abemaciclib",
    "Palbociclib",
    "Ribociclib",
    "Ixazomib",
    "Olaparib",
    "Cobimetinib",
    "Actinomycin D",
    "Doxorubicin",
    "Afatinib",
)
PALETTE = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "purple": "#CC79A7",
    "red": "#D55E00",
    "gray": "#6B7280",
    "light": "#E2E8F0",
    "dark": "#1A202C",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_weight_grid() -> list[tuple[float, float, float, float]]:
    """Return the locked 0.1-step simplex used in the legacy analysis."""
    settings: list[tuple[float, float, float, float]] = []
    values = [index / 10 for index in range(7)]
    for weight_c in values:
        for weight_r in values:
            for weight_b in values:
                weight_e = round(1.0 - weight_c - weight_r - weight_b, 10)
                if weight_e < 0 or weight_e > 0.6:
                    continue
                if weight_c == 0 and weight_r == 0:
                    continue
                settings.append((weight_c, weight_r, weight_b, weight_e))
    return settings


def _load_locked_components(project_root: Path) -> pd.DataFrame:
    path = project_root / "data/bindex_network/ADRS_final_fullSTRING_ranked.csv"
    frame = pd.read_csv(path)
    frame = frame.loc[frame["ACCrelative_resid_pct"].notna()].copy()
    if len(frame) != 108:
        raise ValueError(f"Expected locked 108-drug universe, found {len(frame)}")
    if frame["drug"].duplicated().any():
        raise ValueError("Drug names must be unique")
    return frame


def compute_weight_grid(project_root: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    frame = _load_locked_components(project_root)
    settings = build_weight_grid()
    records: list[dict[str, Any]] = []

    for setting_id, weights in enumerate(settings, start=1):
        values = frame.loc[:, COMPONENT_COLUMNS].to_numpy(dtype=float)
        weight_vector = np.asarray(weights, dtype=float)
        present = np.isfinite(values)
        numerator = np.nansum(values * weight_vector, axis=1)
        denominator = np.sum(present * weight_vector, axis=1)
        if np.any(denominator <= 0):
            raise ValueError("A weight setting left one or more drugs unscored")
        scores = numerator / denominator
        order = sorted(
            range(len(frame)),
            key=lambda index: (-float(scores[index]), str(frame.iloc[index]["drug"])),
        )
        ranks = np.empty(len(frame), dtype=int)
        for rank, index in enumerate(order, start=1):
            ranks[index] = rank
        for index, row in frame.reset_index(drop=True).iterrows():
            records.append(
                {
                    "setting_id": setting_id,
                    "w_C_ACC": weights[0],
                    "w_residual": weights[1],
                    "w_B_neighbor": weights[2],
                    "w_external": weights[3],
                    "drug": row["drug"],
                    "score": float(scores[index]),
                    "rank": int(ranks[index]),
                }
            )

    results = pd.DataFrame.from_records(records)
    summaries = (
        results.groupby("drug", sort=False)["rank"]
        .agg(
            median_rank="median",
            q1_rank=lambda values: float(np.percentile(values, 25)),
            q3_rank=lambda values: float(np.percentile(values, 75)),
            minimum_rank="min",
            maximum_rank="max",
            top10_rate=lambda values: float(np.mean(values <= 10)),
            top20_rate=lambda values: float(np.mean(values <= 20)),
        )
        .reset_index()
    )
    metrics = {
        "analysis_status": "legacy exploratory; excluded from primary robustness claims",
        "universe_n": 108,
        "n_weight_settings": len(settings),
        "grid_step": 0.1,
        "maximum_component_weight": 0.6,
        "requires_C_ACC_or_residual_positive": True,
        "selected_drug_summary": summaries.loc[
            summaries["drug"].isin(SELECTED_DRUGS)
        ].to_dict(orient="records"),
    }
    return results.merge(summaries, on="drug", how="left"), metrics


def _mechanism_class(drug: str, moa: str) -> str:
    name = drug.lower()
    moa_lower = (moa or "").lower()
    if name in {"ribociclib", "palbociclib", "abemaciclib", "trilaciclib"} or "cdk4" in moa_lower:
        return "CDK4/6"
    if "mek" in moa_lower:
        return "MEK"
    if "proteasome" in moa_lower:
        return "Proteasome"
    if "hdac" in moa_lower:
        return "HDAC"
    if "parp" in moa_lower:
        return "PARP"
    if "egfr" in moa_lower or name in {
        "afatinib",
        "erlotinib",
        "osimertinib",
        "neratinib",
        "dacomitinib",
        "tucatinib",
    }:
        return "EGFR"
    if "btk" in moa_lower:
        return "BTK"
    if name in {"crizotinib", "ceritinib", "brigatinib", "tepotinib"} or (
        "alk" in moa_lower.replace("alkyl", "")
    ):
        return "ALK/MET"
    if "topoisomerase" in moa_lower or name in {
        "doxorubicin",
        "daunorubicin",
        "idarubicin",
        "epirubicin",
        "etoposide",
        "irinotecan",
        "topotecan",
    }:
        return "Topoisomerase"
    if "tubulin" in moa_lower or "microtubule" in moa_lower or name in {
        "docetaxel",
        "paclitaxel",
        "vinblastine",
        "vinorelbine",
        "ixabepilone",
    }:
        return "Tubulin"
    if "alkylat" in moa_lower or name in {
        "cisplatin",
        "carboplatin",
        "oxaliplatin",
        "mitomycin",
        "melphalan",
    }:
        return "Alkylator"
    if "antimetabol" in moa_lower or name in {
        "gemcitabine",
        "cytarabine",
        "fluorouracil",
        "pemetrexed",
        "pralatrexate",
    }:
        return "Antimetabolite"
    return "Other"


def _load_associations(path: Path) -> dict[str, set[str]]:
    associations: dict[str, set[str]] = defaultdict(set)
    with path.open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            associations[row["drug"].strip()].add(row["gene"].strip())
    return dict(associations)


def _load_racc(path: Path) -> dict[str, float]:
    with path.open(encoding="utf-8", newline="") as stream:
        return {
            row["gene"]: float(row["rACC_full"])
            for row in csv.DictReader(stream)
        }


def _percentile_by_order(values: dict[str, float], keys: list[str]) -> dict[str, float]:
    ordered = sorted(keys, key=lambda key: values[key])
    denominator = len(ordered) - 1
    return {
        key: index / denominator
        for index, key in enumerate(ordered)
    }


def _full_precision_residuals(
    project_root: Path,
    ordered_drugs: list[str],
) -> dict[str, float]:
    activity: dict[str, float | None] = {}
    with (project_root / "data/bindex_network/Sactivity_124_v1.csv").open(
        encoding="utf-8", newline=""
    ) as stream:
        for row in csv.DictReader(stream):
            activity[row["drug"]] = (
                float(row["mean_ZAUC"]) if row["mean_ZAUC"] else None
            )
    with (project_root / "data/bindex_network/NCI60_potency_124.csv").open(
        encoding="utf-8", newline=""
    ) as stream:
        nci = {
            row["drug"]: float(row["NCI60_mean_neglogGI50"])
            for row in csv.DictReader(stream)
        }
    common = [drug for drug in ordered_drugs if activity.get(drug) is not None]
    acc_percentile = _percentile_by_order(
        {drug: -float(activity[drug]) for drug in common}, common
    )
    nci_percentile = _percentile_by_order(nci, common)
    x = np.asarray([nci_percentile[drug] for drug in common], dtype=float)
    y = np.asarray([acc_percentile[drug] for drug in common], dtype=float)
    slope = float(np.cov(x, y, bias=True)[0, 1] / np.var(x))
    intercept = float(y.mean() - slope * x.mean())
    residual = {
        drug: acc_percentile[drug] - (intercept + slope * nci_percentile[drug])
        for drug in common
    }
    return _percentile_by_order(residual, common)


def compute_lomco(project_root: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    components = _load_locked_components(project_root).set_index("drug")
    associations = _load_associations(
        project_root / "data/bindex_network/bindex_edges_1304.csv"
    )
    racc = _load_racc(
        project_root / "data/bindex_network/rACC_399_fullSTRING.csv"
    )
    drugs = sorted(associations)
    observed = _full_precision_residuals(project_root, drugs)
    drugs = [drug for drug in drugs if drug in observed]
    if set(drugs) != set(components.index):
        raise ValueError("Full-precision residual universe differs from locked 108 drugs")
    missing = sorted(set(drugs) - associations.keys())
    if missing:
        raise ValueError(f"Missing association sets: {missing}")

    gene_weights = {gene: 1.0 + value for gene, value in racc.items()}
    sum_weights = {
        drug: sum(gene_weights[gene] for gene in associations[drug])
        for drug in drugs
    }
    with (project_root / "data/bindex_network/ADRS_v2_ranked.csv").open(
        encoding="utf-8", newline=""
    ) as stream:
        moa = {row["drug"]: row["Primary_MOA"] for row in csv.DictReader(stream)}
    classes = {drug: _mechanism_class(drug, moa.get(drug, "")) for drug in drugs}

    records: list[dict[str, Any]] = []
    for drug in drugs:
        numerator = 0.0
        denominator = 0.0
        neighbor_count = 0
        for reference in drugs:
            if reference == drug or classes[reference] == classes[drug]:
                continue
            shared = associations[drug] & associations[reference]
            if not shared:
                continue
            similarity = 0.5 * sum(gene_weights[gene] for gene in shared) * (
                1.0 / sum_weights[drug] + 1.0 / sum_weights[reference]
            )
            numerator += similarity * observed[reference]
            denominator += similarity
            neighbor_count += 1
        if denominator <= 0:
            continue
        records.append(
            {
                "drug": drug,
                "mechanism_class": classes[drug],
                "observed_residual_percentile": observed[drug],
                "lomco_prediction": numerator / denominator,
                "eligible_neighbor_count": neighbor_count,
                "n_associations": len(associations[drug]),
            }
        )

    results = pd.DataFrame.from_records(records)
    statistic = spearmanr(
        results["observed_residual_percentile"],
        results["lomco_prediction"],
    )
    metrics = {
        "analysis_status": "legacy exploratory; not an independent predictive validation",
        "locked_universe_n": len(drugs),
        "evaluable_n": len(results),
        "spearman_rho": float(statistic.correlation),
        "spearman_p": float(statistic.pvalue),
        "leave_out_unit": "mechanism class",
        "similarity_lambda": 1.0,
    }
    return results, metrics


def _save_figure(figure: plt.Figure, stem: Path) -> list[Path]:
    paths = [stem.with_suffix(suffix) for suffix in (".pdf", ".svg", ".png")]
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        kwargs: dict[str, Any] = {}
        if path.suffix == ".png":
            kwargs["dpi"] = 1000
        figure.savefig(path, **kwargs)
    plt.close(figure)
    return paths


def _plot_weight_grid(results: pd.DataFrame, metrics: dict[str, Any]) -> plt.Figure:
    selected = results.loc[results["drug"].isin(SELECTED_DRUGS)].copy()
    order = (
        selected.groupby("drug")["rank"].median().sort_values().index.tolist()
    )
    distributions = [selected.loc[selected["drug"] == drug, "rank"] for drug in order]
    rates = [float(np.mean(values <= 20)) for values in distributions]

    figure, axes = plt.subplots(
        1,
        2,
        figsize=(170 / 25.4, 3.65),
        gridspec_kw={"width_ratios": [1.9, 1.0]},
    )
    box = axes[0].boxplot(
        distributions,
        orientation="horizontal",
        tick_labels=order,
        patch_artist=True,
        showfliers=False,
        widths=0.62,
    )
    for patch in box["boxes"]:
        patch.set_facecolor("#BEE3F8")
        patch.set_edgecolor(PALETTE["blue"])
    for median in box["medians"]:
        median.set_color(PALETTE["red"])
        median.set_linewidth(1.4)
    axes[0].axvline(20.5, color=PALETTE["orange"], linestyle="--", linewidth=1)
    axes[0].set_xlabel("Rank across 203 weight settings")
    axes[0].invert_yaxis()
    axes[0].grid(axis="x", color="#EDF2F7", linewidth=0.7)
    axes[0].set_title("a  Candidate rank distributions", loc="left", weight="bold")

    positions = np.arange(len(order))
    axes[1].barh(positions, rates, color=PALETTE["green"], alpha=0.88)
    axes[1].set_yticks(positions, labels=[])
    axes[1].set_xlim(0, 1)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("Fraction ranked in Top 20")
    axes[1].set_title("b  Top-20 retention", loc="left", weight="bold")
    axes[1].grid(axis="x", color="#EDF2F7", linewidth=0.7)
    for position, value in zip(positions, rates, strict=True):
        axes[1].text(min(value + 0.025, 0.96), position, f"{value:.0%}", va="center", fontsize=7)

    figure.suptitle(
        "Exploratory four-component weight-grid sensitivity",
        x=0.06,
        y=0.975,
        ha="left",
        fontsize=10,
        weight="bold",
    )
    figure.text(
        0.06,
        0.91,
        f"Locked 108-drug universe; {metrics['n_weight_settings']} simplex settings",
        fontsize=8,
        color=PALETTE["gray"],
    )
    figure.text(
        0.06,
        0.01,
        "Components: C_ACC, ACC-relative residual, S_B-neighbor and S_external. "
        "Legacy exploratory analysis; excluded from primary robustness claims.",
        fontsize=7,
        color=PALETTE["gray"],
    )
    figure.tight_layout(rect=(0.04, 0.10, 1, 0.86))
    return figure


def _plot_lomco(results: pd.DataFrame, metrics: dict[str, Any]) -> plt.Figure:
    figure, axis = plt.subplots(figsize=(170 / 25.4, 4.25))
    highlighted = {"CDK4/6", "MEK", "PARP", "Proteasome", "EGFR"}
    color_map = {
        "CDK4/6": PALETTE["red"],
        "MEK": PALETTE["orange"],
        "PARP": PALETTE["purple"],
        "Proteasome": PALETTE["green"],
        "EGFR": PALETTE["blue"],
        "Other classes": "#A0AEC0",
    }
    plot_class = results["mechanism_class"].where(
        results["mechanism_class"].isin(highlighted), "Other classes"
    )
    for label in ["Other classes", "CDK4/6", "MEK", "PARP", "Proteasome", "EGFR"]:
        subset = results.loc[plot_class == label]
        if subset.empty:
            continue
        axis.scatter(
            subset["observed_residual_percentile"],
            subset["lomco_prediction"],
            s=32 if label != "Other classes" else 20,
            alpha=0.9 if label != "Other classes" else 0.55,
            color=color_map[label],
            edgecolor="white",
            linewidth=0.35,
            label=label,
        )
    axis.set_xlabel("Observed ACC-relative residual percentile")
    axis.set_ylabel("LOMCO drug-neighbour prediction")
    axis.set_title(
        "Leave-one-mechanism-class-out assessment of the drug-neighbour score",
        loc="left",
        weight="bold",
    )
    axis.grid(color="#EDF2F7", linewidth=0.7)
    axis.set_axisbelow(True)
    axis.text(
        0.02,
        0.97,
        f"Spearman ρ = {metrics['spearman_rho']:.3f}\n"
        f"P = {metrics['spearman_p']:.3f}; n = {metrics['evaluable_n']}",
        transform=axis.transAxes,
        va="top",
        ha="left",
        fontsize=8.5,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": PALETTE["light"]},
    )
    axis.legend(frameon=False, fontsize=7, ncol=2, loc="lower right")
    figure.text(
        0.08,
        0.015,
        "Mechanism-class members were excluded from each drug's neighbour reference set. "
        "The result does not support independent predictive value.",
        fontsize=7,
        color=PALETTE["gray"],
    )
    figure.tight_layout(rect=(0.04, 0.06, 1, 1))
    return figure


def generate_outputs(project_root: Path, output_root: Path | None = None) -> dict[str, Any]:
    project_root = project_root.resolve()
    destination = output_root.resolve() if output_root else project_root
    figure_data_dir = destination / "figure_data/revision"
    figure_dir = destination / "figures/revision"
    result_dir = destination / "results/supplementary_legacy"
    for directory in (figure_data_dir, figure_dir, result_dir):
        directory.mkdir(parents=True, exist_ok=True)

    grid_results, grid_metrics = compute_weight_grid(project_root)
    lomco_results, lomco_metrics = compute_lomco(project_root)

    grid_data_path = figure_data_dir / "FigS1_four_component_weight_grid.csv"
    lomco_data_path = figure_data_dir / "FigS2_lomco_neighbor_assessment.csv"
    grid_results.to_csv(grid_data_path, index=False)
    lomco_results.to_csv(lomco_data_path, index=False)

    grid_paths = _save_figure(
        _plot_weight_grid(grid_results, grid_metrics),
        figure_dir / "FigS1_four_component_weight_grid",
    )
    lomco_paths = _save_figure(
        _plot_lomco(lomco_results, lomco_metrics),
        figure_dir / "FigS2_lomco_neighbor_assessment",
    )

    metrics_path = result_dir / "supplementary_legacy_metrics.json"
    metrics_path.write_text(
        json.dumps(
            {"figure_s1": grid_metrics, "figure_s2": lomco_metrics},
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    source_paths = [
        project_root / "data/bindex_network/ADRS_final_fullSTRING_ranked.csv",
        project_root / "data/bindex_network/ADRS_v2_ranked.csv",
        project_root / "data/bindex_network/bindex_edges_1304.csv",
        project_root / "data/bindex_network/rACC_399_fullSTRING.csv",
        project_root / "data/bindex_network/Sactivity_124_v1.csv",
        project_root / "data/bindex_network/NCI60_potency_124.csv",
    ]
    manifest_path = result_dir / "run_manifest.json"
    manifest = {
        "analysis": "regenerated exploratory Supplementary Figures S1-S2",
        "status": "descriptive/exploratory; not a primary validation analysis",
        "inputs": {
            str(path.relative_to(project_root)).replace("\\", "/"): _sha256(path)
            for path in source_paths
        },
        "outputs": [
            str(path.relative_to(destination)).replace("\\", "/")
            for path in [grid_data_path, lomco_data_path, *grid_paths, *lomco_paths, metrics_path]
        ],
    }
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return {
        "weight_grid_metrics": grid_metrics,
        "lomco_metrics": lomco_metrics,
        "paths": [
            grid_data_path,
            lomco_data_path,
            *grid_paths,
            *lomco_paths,
            metrics_path,
            manifest_path,
        ],
    }


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument("--output-root", type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    outputs = generate_outputs(args.project_root, args.output_root)
    print(
        "Supplementary legacy figures complete: "
        f"S1 settings={outputs['weight_grid_metrics']['n_weight_settings']}; "
        f"S2 rho={outputs['lomco_metrics']['spearman_rho']:.3f}; "
        f"n={outputs['lomco_metrics']['evaluable_n']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
