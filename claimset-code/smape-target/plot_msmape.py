"""Figure of plan PL11: the msMAPE target beside the sMAPE target, drawn from the recorded results of its run.

Left: near zero, the target as a share of the median of the positive outcomes, against the share of zeros, for
lognormal outcomes (sigma 1) at median 1 and at median 100, with each measure's zero threshold as a vertical line.
Middle: away from zero, the relative difference of the msMAPE target from the sMAPE target, in per cent, for each
distribution at median about 100. Right: on Poisson counts, the two targets against the mean."""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BLUE, ORANGE = "#1c5cab", "#d4762c"
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#6b6a66", "#e6e5e1", "#fcfcfb"
COLOURS = {"sMAPE": BLUE, "msMAPE": ORANGE}


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=0.6)
    for s in ax.spines.values():
        s.set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=7)


def value(res, variable, **cond):
    hits = [r["value"] for r in res if r["variable"] == variable
            and all(r["conditions"].get(k) == v for k, v in cond.items())]
    return hits[0] if hits else None


def main(path, out):
    res = json.load(open(path))
    os.makedirs(os.path.join(out, "figures"), exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    for ax in axes:
        style(ax)

    ax = axes[0]
    for med, mark in ((1.0, "o"), (100.0, "s")):
        for meas in ("sMAPE", "msMAPE"):
            pts = []
            for p0 in (0.0, 0.2, 0.4, 0.6, 0.8):
                name = f"lognormal median={med:g} sigma=1" + (f" with zeros p0={p0:g}" if p0 else "")
                t = value(res, "target", setting=name, measure=meas)
                if t is not None:
                    pts.append((p0, t / med))
            ax.plot(*zip(*pts), marker=mark, color=COLOURS[meas], linestyle="-" if med == 1 else "--",
                    label=f"{meas}, median {med:g}", markersize=4, linewidth=1)
            z = value(res, "zero threshold", setting=f"lognormal median={med:g} sigma=1", measure=meas)
            if z is not None and med == 1.0:
                ax.axvline(z, color=COLOURS[meas], linewidth=0.8, linestyle=":")
    ax.set_xlabel("share of zero outcomes", fontsize=8, color=INK)
    ax.set_ylabel("target / median of the positive outcomes", fontsize=8, color=INK)
    ax.set_title("Near zero (lognormal, sigma 1): dotted lines are\nthe zero thresholds at median 1", fontsize=8,
                 color=INK, loc="left")
    ax.legend(fontsize=6, frameon=False)

    ax = axes[1]
    names = [r["conditions"]["setting"] for r in res if r["variable"].startswith("relative difference")
             and r["conditions"].get("group") == "away from zero"]
    diffs = [100 * value(res, "relative difference of the msMAPE target from the sMAPE target", setting=n)
             for n in names]
    ax.barh(range(len(names)), diffs, color=ORANGE)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=6)
    ax.axvline(0, color=MUTED, linewidth=0.8)
    ax.set_xlabel("msMAPE target minus sMAPE target, % of the sMAPE target", fontsize=8, color=INK)
    ax.set_title("Away from zero (median about 100)", fontsize=8, color=INK, loc="left")

    ax = axes[2]
    for meas, mark, size in (("sMAPE", "o", 7), ("msMAPE", "x", 6)):
        pts = [(m, value(res, "target", setting=f"Poisson mean={m:g}", measure=meas))
               for m in (0.3, 0.5, 1.0, 1.5, 3.0, 10.0)]
        ax.plot(*zip(*pts), marker=mark, linestyle="none", color=COLOURS[meas], markersize=size, label=meas)
    ax.set_xscale("log")
    ax.set_xlabel("Poisson mean", fontsize=8, color=INK)
    ax.set_ylabel("target", fontsize=8, color=INK)
    ax.set_title("Poisson counts: the two targets coincide", fontsize=8, color=INK, loc="left")
    ax.legend(fontsize=6, frameon=False)

    fig.suptitle("msMAPE beside sMAPE: their targets, computed from the distributions", fontsize=9, color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "figures", "msmape_targets.png"), dpi=150)


if __name__ == "__main__":
    main(sys.argv[1], os.environ.get("CLAIMSET_RUN_DIR") or ".")
