"""Posterior-weighted number of fixed mutations of each layer per depth tier (summed over organs or replicate groups),
with the selected window and the resolvability status. Usage: python plot_layer_counts.py RESULTS_DIR OUT.png"""
import sys
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

res, out = Path(sys.argv[1]), sys.argv[2]
r = pd.read_csv(res / "03_fixed_counts_by_tree_depth.csv")
s = pd.read_csv(res / "05_selected_window_summary.csv").iloc[0]
fig, ax = plt.subplots(figsize=(6, 3.4))
ax.axvspan(s.selected_window_start - 0.3, s.selected_window_end + 0.3, color="#eef3fa", lw=0)
for layer, col in (("L1", "#e8692d"), ("L2", "#2a7de1"), ("L3", "#8a5ac2")):
    status = s.get(f"{layer}_status", "")
    ok = status == "reportable"
    ax.plot(r.depth_filter, r[f"sum_posterior_weighted_{layer}_count"], marker="o", color=col, lw=1.6 if ok else 1,
            ls="-" if ok else ":", mfc=col if ok else "white", label=f"{layer} ({status})")
ax.set_xlabel("depth tier (minimum depth)"); ax.set_ylabel("posterior-weighted fixed mutations")
ax.set_title("shaded: selected window; dotted: layer not resolvable", fontsize=9, loc="left")
ax.legend(fontsize=8, frameon=False); fig.tight_layout(); fig.savefig(out, dpi=150); print("wrote", out)
