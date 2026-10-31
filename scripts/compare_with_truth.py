"""Example data only: compare vafSOMA's L2 counts with the simulated truth.
Usage: python compare_with_truth.py RESULTS_DIR EXAMPLE_DIR"""
import sys, zipfile
from pathlib import Path
import pandas as pd

res, ex = Path(sys.argv[1]), Path(sys.argv[2])
truth = pd.read_csv(ex / "truth_layer_carriers.csv.gz")
with zipfile.ZipFile(ex / "vaf_tables.zip") as z:
    first = sorted(z.namelist(), key=lambda n: int(n.split("_dp_")[1].split("_")[0]))[0]
    sites = set(pd.read_csv(z.open(first)).site_id)           # sites in the input tables
l2 = truth[(truth.layer == "L2") & truth.mutation_id.isin(sites)].copy()
l2["tip"] = l2.organ_id.str.split("_").str[0]
fixed = l2[l2.carrier_fraction >= 0.999]
groups = res / "13_group_variant_calls_window.csv"
if groups.exists():                                            # replicate-group model: one count per tip
    g = pd.read_csv(groups)
    est = g.groupby("replicate_group").mean_P_L2.sum()
    n_leaves = l2.groupby("tip").organ_id.nunique()
    true = fixed.groupby(["tip", "mutation_id"]).organ_id.nunique()
    true = (true == n_leaves.reindex(true.index.get_level_values(0)).to_numpy()).groupby(level=0).sum()
    unit = "tip (fixed in all leaves)"
else:                                                          # per-organ model: one count per leaf
    o = pd.read_csv(res / "08_organ_variant_calls_window.csv")
    est = o.groupby("sample").mean_P_L2.sum()
    true = fixed.groupby("organ_id").size()
    unit = "leaf"
t = pd.DataFrame({"true fixed L2": true, "vafSOMA (posterior-weighted)": est.round(1)}).fillna(0)
t.index.name = unit
print(t.to_string()); print(f"\ntotal: true {t.iloc[:, 0].sum():.0f}, vafSOMA {t.iloc[:, 1].sum():.1f}")
