"""Figure of version 4: the asymptotic variance V of each measure's minimiser (H49, H52) against the spread, on log-log
axes, for lognormal, log-Laplace and log-uniform actuals. V is computed by quadrature with the formulas of H49 (the code
of experiment_v4.py); where a moment that H49 needs is infinite (RMSE and MAPE on log-Laplace actuals with scale 1/2 or
more), the line stops. Writes its values to $CLAIMSET_RESULTS and the figure to <run dir>/figures/."""
import json
import math
import os
import sys

import matplotlib
import numpy as np

import experiment_v4 as ex

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

INK, MUTED, GRID, SURFACE = "#0b0b0b", "#6b6a66", "#e6e5e1", "#fcfcfb"
COLOURS = {"RMSE": "#1c5cab", "MAE": "#2a9d55", "MAPE": "#c2571a", "sMAPE": "#8e3aa8"}
FAMILIES = (("lognormal", ex.lognormal, "standard deviation of ln y"),
            ("log-Laplace", ex.log_laplace, "scale of ln y"),
            ("log-uniform", ex.log_uniform, "half-width of ln y"))
SPREADS = np.geomspace(0.1, 8.0, 23)


def main(out):
    res = []
    for family, make, _ in FAMILIES:
        for s in SPREADS:
            vf = make(float(s)).variances()
            for m in ("RMSE", "MAE", "MAPE", "sMAPE"):
                v = vf[m][0]
                if math.isfinite(v):
                    res.append({"variable": "asymptotic variance by formula", "value": float(v),
                                "conditions": {"family": family, "spread": round(float(s), 6), "measure": m}})
        print(family, flush=True)
    path = os.environ.get("CLAIMSET_RESULTS")
    if path:
        with open(path, "w") as fh:
            json.dump(res, fh)
    os.makedirs(os.path.join(out, "figures"), exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.6), sharey=True, facecolor=SURFACE)
    for ax, (family, _, xlabel) in zip(axes, FAMILIES):
        ax.set_facecolor(SURFACE)
        for m, colour in COLOURS.items():
            pts = [(r["conditions"]["spread"], r["value"]) for r in res
                   if r["conditions"]["family"] == family and r["conditions"]["measure"] == m]
            if pts:
                x, y = zip(*pts)
                ax.plot(x, y, color=colour, lw=1.6, label=m)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(1e-3, 1e8)
        ax.set_title(family, fontsize=9, color=INK, loc="left")
        ax.set_xlabel(xlabel, fontsize=8, color=INK)
        ax.grid(color=GRID, lw=0.6)
        ax.tick_params(colors=MUTED, labelsize=7)
        for sp in ax.spines.values():
            sp.set_visible(False)
    axes[0].set_ylabel("asymptotic variance V of ln(minimiser)", fontsize=8, color=INK)
    axes[0].legend(fontsize=7, frameon=False)
    fig.suptitle("How precisely a test set locates each measure's target as the actuals spread out "
                 "(larger V, less precise)", fontsize=9, color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "figures", "precision_growth.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.environ.get("CLAIMSET_RUN_DIR", "."))
