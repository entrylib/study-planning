"""Figure of plan PL2: relative efficiency over test-set size, one panel per setting, with its 99% band and the
asymptotic value (H11) as a dashed line. Above 1 the sMAPE minimiser varies less than the sample median."""
import json
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

LINE, BAND = "#1c5cab", "#3987e5"
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#6b6a66", "#e6e5e1", "#fcfcfb"


def main(out):
    res = json.load(open(os.environ.get("CLAIMSET_RESULTS") or os.path.join(out, "results.json")))
    os.makedirs(os.path.join(out, "figures"), exist_ok=True)
    settings = list(dict.fromkeys(r["conditions"]["setting"] for r in res))
    for family in ("lognormal", "log-Laplace"):
        sub = [s for s in settings if s.startswith(family)]
        if not sub:
            continue
        cols = 4
        rows = (len(sub) + cols - 1) // cols
        fig, axes = plt.subplots(rows, cols, figsize=(3.2 * cols, 2.6 * rows), sharey=True, facecolor=SURFACE)
        for ax, s in zip(np.atleast_1d(axes).ravel(), sub):
            ax.set_facecolor(SURFACE)
            pts = sorted((r["conditions"]["test-set size"], r["value"], r["interval"]["lower"], r["interval"]["upper"])
                         for r in res if r["variable"] == "relative efficiency" and r["conditions"]["setting"] == s)
            n, v, lo, hi = map(np.array, zip(*pts))
            ax.fill_between(n, lo, hi, color=BAND, alpha=0.25, lw=0)
            ax.plot(n, v, color=LINE, lw=1.4, marker="o", ms=2.5)
            asym = next(r["value"] for r in res if r["variable"] == "asymptotic relative efficiency"
                        and r["conditions"]["setting"] == s)
            ax.axhline(asym, color=INK, ls="--", lw=1)
            ax.axhline(1, color=MUTED, ls=":", lw=0.9)
            ax.set_xscale("log")
            ax.set_title(f"{s}: asymptotic {asym:.3g}", fontsize=8, color=INK, loc="left")
            ax.grid(color=GRID, lw=0.6)
            ax.tick_params(colors=MUTED, labelsize=7)
            for sp in ax.spines.values():
                sp.set_visible(False)
        for ax in np.atleast_1d(axes).ravel()[len(sub):]:
            ax.set_visible(False)
        fig.supxlabel("test-set size n", fontsize=8, color=INK)
        fig.supylabel("Var(ln median) / Var(ln sMAPE minimiser)", fontsize=8, color=INK)
        fig.suptitle(f"{family}: relative efficiency of the sMAPE minimiser against the sample median (99% band; "
                     "dashed, large-sample value; above 1, sMAPE varies less)", fontsize=9, color=INK)
        fig.tight_layout()
        fig.savefig(os.path.join(out, "figures", f"efficiency_{family}.png"), dpi=150)
        plt.close(fig)


if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else os.environ.get("CLAIMSET_RUN_DIR", "."))
