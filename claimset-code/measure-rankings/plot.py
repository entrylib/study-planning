"""Figures of the choices: one PNG per outcome distribution, one panel per measure.

In each panel, one violin per test-set size n shows the distribution, over test sets, of the candidate the measure
ranks first: its quantile level for continuous outcomes, its value for counts. A boxplot inside each violin marks the
median, the quartiles and the 5th and 95th percentiles, and a diamond the mean of the choices. A candidate that is not a
quantile (the mean, a minimiser) is placed at its own quantile level, the share of the distribution at or below it. Horizontal lines mark every measure's preferred candidate (the
lowest expected score): the panel's own as a dashed line, the other measures' as dotted lines, each labelled at the
right, so that the gaps between the targets can be read against the spread of the choices."""
import csv
import json
import os
import re

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

FILL, EDGE = "#3987e5", "#1c5cab"                    # one series: one hue (the reference palette's blue)
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#6b6a66", "#e6e5e1", "#fcfcfb"
SHORT = {"RMSE": "RMSE", "MAE": "MAE", "MAPE": "MAPE", "sMAPE": "sMAPE"}


def _samples(rows, scale):
    """The choices over test sets, rebuilt from the shares (each share times the number of test sets is a count)."""
    out = []
    for r in rows:
        k = int(round(float(r["share"]) * int(r["test sets"])))
        out += [float(r[scale])] * k
    return np.array(out)


def main(out):
    rows = list(csv.DictReader(open(os.path.join(out, "selection_shares.csv"))))
    res = json.load(open(os.environ.get("CLAIMSET_RESULTS") or os.path.join(out, "results.json")))
    os.makedirs(os.path.join(out, "figures"), exist_ok=True)
    for setting in dict.fromkeys(r["setting"] for r in rows):
        count = setting.startswith(("Poisson", "intermittent"))
        scale = "value" if count else "quantile level"
        var = "value of the preferred candidate" if count else "quantile level of the preferred candidate"
        targets = {r["conditions"]["measure"]: r["value"] for r in res
                   if r["variable"] == var and r["conditions"]["setting"] == setting}
        sub = [r for r in rows if r["setting"] == setting]
        measures = list(dict.fromkeys(r["measure"] for r in sub))
        sizes = sorted({int(r["test-set size"]) for r in sub})
        fig, axes = plt.subplots(len(measures), 1, figsize=(10, 2.3 * len(measures)), sharex=True, facecolor=SURFACE)
        for ax, m in zip(np.atleast_1d(axes), measures):
            ax.set_facecolor(SURFACE)
            data = [_samples([r for r in sub if r["measure"] == m and int(r["test-set size"]) == n], scale) for n in sizes]
            pos = np.arange(len(sizes))
            spread = [i for i, d in enumerate(data) if d.size and np.ptp(d) > 0]
            if spread:
                vp = ax.violinplot([data[i] for i in spread], positions=pos[spread], widths=0.85,
                                   showextrema=False, bw_method=0.08 if count else 0.3)   # counts: bulges at values
                for b in vp["bodies"]:
                    b.set_facecolor(FILL)
                    b.set_edgecolor(EDGE)
                    b.set_linewidth(0.8)
                    b.set_alpha(0.35)
            ax.boxplot(data, positions=pos, widths=0.18, whis=(5, 95), showfliers=False, patch_artist=True,
                       showmeans=True, meanprops=dict(marker="D", markersize=4, markerfacecolor=INK, markeredgecolor=INK),
                       boxprops=dict(facecolor=SURFACE, edgecolor=INK, linewidth=0.9),
                       medianprops=dict(color=INK, linewidth=1.6), whiskerprops=dict(color=INK, linewidth=0.9),
                       capprops=dict(color=INK, linewidth=0.9))
            lines = []                                    # targets closer than this share one labelled line
            close = 1e-6 if count else 0.005
            for o, v in sorted(targets.items(), key=lambda t: t[1]):
                if lines and v - lines[-1][0] <= close:
                    lines[-1][1].append(o)
                else:
                    lines.append([v, [o]])
            if not count:
                ax.set_ylim(0, 1)
            lo, hi = ax.get_ylim()
            spots = []                                    # label heights in axes fraction, at least GAP apart
            for v, _ in lines:
                y = (v - lo) / (hi - lo)
                spots.append(max(y, spots[-1] + 0.09) if spots else y)
            for (v, names), y in zip(lines, spots):
                own = m in names
                ax.axhline(v, color=INK if own else MUTED, lw=1.1 if own else 0.9, ls="--" if own else ":", zorder=0)
                ax.text(1.005, y, " / ".join(SHORT[o] for o in names), transform=ax.transAxes,
                        fontsize=7, color=INK if own else MUTED, va="center", ha="left")
            ax.set_title(f"{m}: candidate ranked first in each test set", fontsize=9, color=INK, loc="left")
            ax.set_ylabel("value" if count else "quantile level", fontsize=8, color=INK)
            if not count:
                ax.set_ylim(0, 1)
            ax.grid(axis="y", color=GRID, lw=0.6)
            ax.tick_params(colors=MUTED, labelsize=7)
            for s in ax.spines.values():
                s.set_visible(False)
        last = np.atleast_1d(axes)[-1]
        last.set_xticks(np.arange(len(sizes)), [f"{n:,}" for n in sizes])
        last.set_xlabel("test-set size n (outcomes per test set)", fontsize=8, color=INK)
        fig.suptitle(f"{setting}: dashed, the panel's own preferred candidate; dotted, the other measures' "
                     "(lowest expected score)", fontsize=9, color=INK)
        fig.tight_layout(rect=(0, 0, 0.94, 1))
        fig.savefig(os.path.join(out, "figures", re.sub(r"[^A-Za-z0-9.=-]+", "_", setting) + ".png"), dpi=150)
        plt.close(fig)


if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else os.environ.get("CLAIMSET_RUN_DIR", "."))
