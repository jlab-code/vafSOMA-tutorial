"""VAF spectrum of one organ at one depth tier, coloured by vafSOMA's call, with the expected layer VAFs.
Usage: python plot_vaf_spectrum.py RESULTS_DIR VAF_TABLES ORGAN TIER OUT.png"""
import sys, zipfile
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

res, tables, organ, tier, out = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], int(sys.argv[4]), sys.argv[5]
with zipfile.ZipFile(tables) as z:
    name = next(n for n in z.namelist() if n.endswith(f"_dp_{tier}_vaf.csv") or n.endswith(f"_dp_{tier}_vaf.txt"))
    t = pd.read_csv(z.open(name), sep=None, engine="python")
t = t[(t["sample"] == organ) & (t["alt_depth"] > 0)]
calls = pd.read_csv(res / "08_organ_variant_calls_window.csv")
calls = calls[calls["sample"] == organ][["pos", "organ_calling_status"]]
t = t.merge(calls, on="pos", how="left").fillna({"organ_calling_status": "not evaluated"})
fit = pd.read_csv(res / "02_fit_parameters_by_tree_depth.csv")
fit = fit[fit["depth_filter"] == tier].iloc[0]

groups = [("called_noise", "background", "#c9c9c9"), ("called_fixed_L2", "fixed in L2", "#2a7de1"),
          ("called_fixed_L1", "fixed in L1", "#e8692d"), ("called_fixed_L3", "fixed in L3", "#8a5ac2")]
fig, ax = plt.subplots(figsize=(6, 3.4))
bins = [i / 50 for i in range(51)]
data = [t.loc[t.organ_calling_status == g, "VAF"] for g, _, _ in groups]
ax.hist(data, bins=bins, stacked=True, color=[c for *_, c in groups],
        label=[f"{lab} ({len(d)})" for (_, lab, _), d in zip(groups, data)])
for layer, col, y in (("L1", "#e8692d", 0.86), ("L2", "#2a7de1", 0.97), ("L3", "#8a5ac2", 0.97)):
    phi = fit[f"mu_{layer}"]
    if pd.notna(phi) and phi > 0:
        ax.axvline(phi, color=col, ls="--", lw=1)
        ax.text(phi, y, f" φ$_{layer[1]}$ = {phi:.3f}", color=col, fontsize=8, va="top",
                transform=ax.get_xaxis_transform())
ax.set_yscale("log"); ax.set_xlabel("VAF"); ax.set_ylabel("positive observations")
ax.set_title(f"{organ}, depth tier ≥ {tier}", fontsize=10, loc="left")
ax.legend(fontsize=8, frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.85))
fig.tight_layout(); fig.savefig(out, dpi=150); print("wrote", out)
