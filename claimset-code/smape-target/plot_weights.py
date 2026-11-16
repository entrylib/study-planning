"""Figure 1 of the proofs (H4 in the record, H5 in the paper): the weight y/(y+f)^2 of an actual on the forecast that
minimises expected sMAPE, against the ratio of actual to forecast, relative to its largest value; MAE weighs every actual
the same.

    python3 plot_weights.py --out <folder>
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

S1, S2 = "#2a78d6", "#eb6834"
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#6b6a66", "#e6e5e1", "#fcfcfb"
r = np.exp(np.linspace(np.log(0.01), np.log(100), 2001))
w = 4 * r / (1 + r) ** 2                      # sMAPE weight y/(y+f)^2, relative to its largest value 1/(4f)

fig, ax = plt.subplots(figsize=(7.5, 4.0), facecolor=SURFACE)
ax.set_facecolor(SURFACE)
ax.plot(r, np.ones_like(r), color=S2, lw=2, label="MAE: every actual weighs the same")
ax.plot(r, w, color=S1, lw=2, label="sMAPE: weight $y/(y+f)^2$, relative to its largest value")
for x in (0.1, 10):
    v = 4 * x / (1 + x) ** 2
    ax.plot([x], [v], "o", ms=6, color=S1, mec=SURFACE, mew=2, zorder=3)
    ax.annotate(f"{v:.2f}", (x, v), textcoords="offset points", xytext=(0, 9), ha="center", fontsize=8, color=INK)
ax.annotate("largest at $y = f$", (1, 1), textcoords="offset points", xytext=(38, -26), ha="left", fontsize=8, color=INK,
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
ax.set_xscale("log")
ticks = [0.01, 0.1, 0.25, 0.5, 1, 2, 4, 10, 100]
ax.set_xticks(ticks, ["1/100", "1/10", "1/4", "1/2", "1", "2", "4", "10", "100"])
ax.minorticks_off()
ax.set_ylim(0, 1.12)
ax.set_xlabel("actual / forecast ($y/f$, log scale)", fontsize=9, color=INK)
ax.set_ylabel("weight of the actual on the forecast", fontsize=9, color=INK)
ax.grid(color=GRID, lw=0.6)
ax.tick_params(colors=MUTED, labelsize=8)
for s in ax.spines.values():
    s.set_visible(False)
ax.legend(fontsize=8, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 0.02))
ax.set_title("How much each actual pulls the forecast that minimises the measure (H5)", fontsize=10, color=INK,
             loc="left")
fig.tight_layout()
ap = argparse.ArgumentParser()
ap.add_argument("--out", default=".")
out = os.environ.get("CLAIMSET_RUN_DIR", ap.parse_args().out)
fig.savefig(os.path.join(out, "smape_weights.png"), dpi=200)
