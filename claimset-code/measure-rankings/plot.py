"""Figures of the selection shares: one PNG per setting, one panel per measure.

Continuous settings: the share of test sets choosing each candidate, against the candidate's quantile level, one line
per test-set size, with a short tick at the top for the mean quantile level of the choices (their central tendency)
and a dashed line at the measure's preferred candidate. Counts: grouped bars, one group per candidate value."""
import csv
import json
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SIZES = (10, 100, 1000, 100000)                       # the ordinal ramp passes the palette checks for four steps
RAMP = ("#86b6ef", "#3987e5", "#1c5cab", "#0d366b")   # one hue, light to dark: small to large test sets
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#6b6a66", "#e6e5e1", "#fcfcfb"


def _read(out):
    rows = list(csv.DictReader(open(os.path.join(out, "selection_shares.csv"))))
    res = json.load(open(os.environ.get("CLAIMSET_RESULTS") or os.path.join(out, "results.json")))
    return rows, res


def main(out):
    rows, res = _read(out)
    pref = {(r["conditions"]["setting"], r["conditions"]["measure"]): r["value"] for r in res
            if r["variable"] == "quantile level of the preferred candidate"}
    prefv = {(r["conditions"]["setting"], r["conditions"]["measure"]): r["value"] for r in res
             if r["variable"] == "value of the preferred candidate"}
    central = {}
    for r in res:
        c = r["conditions"]
        count = c.get("setting", "").startswith(("Poisson", "intermittent"))
        want = "mean value of the chosen candidates" if count else "mean quantile level of the chosen candidates"
        if r["variable"] == want:
            central[(c["setting"], c["measure"], c["test-set size"])] = r["value"]
    os.makedirs(os.path.join(out, "figures"), exist_ok=True)
    for setting in dict.fromkeys(r["setting"] for r in rows):
        count = setting.startswith(("Poisson", "intermittent"))
        measures = list(dict.fromkeys(r["measure"] for r in rows if r["setting"] == setting))
        fig, axes = plt.subplots(1, len(measures), figsize=(3.4 * len(measures), 3.0), sharey=True, facecolor=SURFACE)
        for ax, m in zip(axes, measures):
            ax.set_facecolor(SURFACE)
            sub = [r for r in rows if r["setting"] == setting and r["measure"] == m]
            if count:
                values = sorted({float(r["value"]) for r in sub})
                width = 0.8 / len(SIZES)
                for j, (n, col) in enumerate(zip(SIZES, RAMP)):
                    share = {float(r["value"]): float(r["share"]) for r in sub if int(r["test-set size"]) == n}
                    ax.bar([i + (j - 1.5) * width for i in range(len(values))], [share.get(v, 0) for v in values],
                           width=width * 0.9, color=col, label=f"n = {n:,}")
                ax.set_xticks(range(len(values)), [f"{v:g}" if v == int(v) else f"{v:.2f}" for v in values],
                              fontsize=7, rotation=0)
                p = prefv.get((setting, m))
                if p is not None and p in values:
                    ax.axvline(values.index(p), color=MUTED, lw=1, ls="--")
                ax.set_xlabel("candidate value", fontsize=8, color=INK)
            else:
                for n, col in zip(SIZES, RAMP):
                    pts = sorted((float(r["quantile level"]), float(r["share"])) for r in sub if int(r["test-set size"]) == n)
                    if not pts:
                        continue
                    xs, ys, acc = [0.0], [0.0], 0.0             # cumulative share: a step at each chosen candidate
                    for x, y in pts:
                        xs += [x, x]
                        ys += [acc, acc + y]
                        acc += y
                    xs.append(1.0)
                    ys.append(acc)
                    ax.plot(xs, ys, color=col, lw=2, label=f"n = {n:,}")
                    c = central.get((setting, m, n))
                    if c is not None:
                        ax.plot([c, c], [1.0, 0.94], color=col, lw=2, transform=ax.get_xaxis_transform(), clip_on=False)
                p = pref.get((setting, m))
                if p is not None:
                    ax.axvline(p, color=MUTED, lw=1, ls="--")
                ax.set_xlim(0, 1)
                ax.set_xlabel("quantile level of the chosen candidate", fontsize=8, color=INK)
            ax.set_title(m, fontsize=9, color=INK)
            ax.grid(axis="y", color=GRID, lw=0.6)
            ax.tick_params(colors=MUTED, labelsize=7)
            for s in ax.spines.values():
                s.set_visible(False)
        axes[0].set_ylabel("share of test sets choosing it" if count else
                           "share of test sets choosing\na candidate at or below this level", fontsize=8, color=INK)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, fontsize=7, frameon=False, loc="center right")
        fig.suptitle(f"{setting}: which candidate each measure ranks first" +
                     ("; dashed: the preferred candidate" if count else
                      "; ticks: mean level of the choices; dashed: the preferred candidate"), fontsize=9, color=INK)
        fig.tight_layout(rect=(0, 0, 0.9, 1))
        fig.savefig(os.path.join(out, "figures", re.sub(r"[^A-Za-z0-9.=-]+", "_", setting) + ".png"), dpi=150)
        plt.close(fig)


if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else os.environ.get("CLAIMSET_RUN_DIR", "."))
