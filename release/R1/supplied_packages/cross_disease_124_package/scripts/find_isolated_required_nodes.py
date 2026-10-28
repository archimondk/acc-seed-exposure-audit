"""List required nodes that are not STRING v12 preferred names (run locally).

Prints, for the ACC arm, every gene passed as a required node (399 association
genes plus the 45 disease-only seeds and the 5 therapy-informed genes) that is
absent from 9606.protein.info.v12.0.txt.gz preferred_name, i.e. the genes that
load_string_graph adds as isolated nodes.

Usage (from the ACC repository root):
    python -I find_isolated_required_nodes.py <protein_info.txt.gz>
"""
import csv
import gzip
import json
import sys
from pathlib import Path

info = Path(sys.argv[1])
with gzip.open(info, "rt", encoding="utf-8") as handle:
    header = next(handle).rstrip("\n").split("\t")
    col = header.index("preferred_name")
    names = {line.rstrip("\n").split("\t")[col] for line in handle}
genes = {r["gene"] for r in csv.DictReader(open("data/bindex_network/bindex_genes_399.csv", encoding="utf-8-sig"))}
seeds = set(json.load(open("leakage_audit/observational/seed_set_used.json"))["weights"])
therapy = {"MGMT", "SLFN11", "ABCB1", "SOAT1", "UBA1"}
for label, group in (("399 association genes", genes), ("45 disease-only seeds", seeds), ("therapy-informed genes", therapy)):
    missing = sorted(group - names)
    print(f"{label}: {len(missing)} not in STRING preferred_name -> {missing}")
