"""Condense the simulations' evidence: summary tables that select and rearrange recorded values (no value is computed
here), the full results of each run in a CSV artefact of its own, and the runs' figures as data plots that name the
observations they show."""
import csv
import math
import hashlib
import json
import re
import shutil
from pathlib import Path


def sig4(v):
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        return v
    if isinstance(v, int) or float(v).is_integer():
        return int(v) if abs(v) < 1e15 else v
    return float(f"{v:.4g}")


def finite(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


def value(res, var, **cond):
    hits = [r for r in res if r["variable"] == var and all(r.get("conditions", {}).get(k) == x for k, x in cond.items())]
    return hits[0]["value"] if hits else None


def table(tid, description, res, row_keys, rows, cols, artefact, conditions=None):
    """rows: list of label tuples (values of row_keys); cols: list of (variable, extra conditions, column spec)."""
    out_rows = []
    for i, lab in enumerate(rows):
        cond = dict(zip(row_keys, lab))
        vals = [sig4(value(res, var, **cond, **extra)) for var, extra, _ in cols]
        if all(v is None for v in vals):
            continue
        out_rows.append({"id": f"r{len(out_rows) + 1}", "label": list(lab), "values": vals})
    t = {"id": tid, "description": description, "row_labels": list(row_keys),
         "columns": [spec for _, _, spec in cols], "rows": out_rows, "artefact": artefact}
    if conditions:
        t["conditions"] = conditions
    return t


def col(var, _as=None, **extra):
    """A column: the record's variable and conditions that pick its values, and how the table names it (`_as`, where
    the paper's name differs from the record's: the Suilin loss was called msMAPE when RN10 ran)."""
    return (var, extra, _as if _as is not None else ({"variable": var, "conditions": extra} if extra else var))


def settings(res, var=None, **cond):
    return list(dict.fromkeys(r["conditions"]["setting"] for r in res if (var is None or r["variable"] == var)
                              and all(r["conditions"].get(k) == x for k, x in cond.items())))


def slug(s):
    return re.sub(r"[^A-Za-z0-9.=-]+", "_", s)


def sha(path):
    return "sha256:" + hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(d, runs_dir, out_dir, master_csv, before_new=None):
    runs_dir, out_dir = Path(runs_dir), Path(out_dir)
    (out_dir / "files" / "figures").mkdir(parents=True, exist_ok=True)
    ev = {tuple(e.get("lab_runs") or []): e for e in d["evidence"]}
    arts = d["artefacts"]
    master = next(a for a in arts if a["uri"] == master_csv)
    with open(out_dir / master_csv, newline="") as fh:
        reader = csv.DictReader(fh)
        header, all_rows = reader.fieldnames, list(reader)

    def new_art(a):
        a["id"] = f"A{len(arts) + 1}"
        arts.append(a)
        return a["id"]

    def run_csv(rn, title, note=""):
        """The run's every result, in a CSV of its own, sorted by setting so that a setting's rows are contiguous."""
        rows = [r for r in all_rows if r["run"] == rn]
        def key(r):
            c = json.loads(r["conditions"] or "{}")
            return (str(c.get("setting", "")), r["quantity"], json.dumps(c, sort_keys=True))
        rows.sort(key=key)
        path = out_dir / "files" / f"smape-target.{rn}.results.csv"
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=header)
            w.writeheader()
            w.writerows(rows)
        aid = new_art({"kind": "data", "uri": f"files/{path.name}", "version": sha(path), "licence": "CC-BY-4.0",
                       "availability": {"level": "public"}, "columns": master.get("columns"),
                       "description": f"Every result of lab run {rn} ({title}) at full precision, one row per result, "
                                      f"sorted by setting: the lab result's ID, its run, the quantity, its conditions "
                                      f"and its value (both as JSON). The evidence's tables give selected values to "
                                      f"four significant digits." + note})
        spans = {}
        for i, r in enumerate(rows, start=2):                  # RFC 7111: the header is row 1
            s = json.loads(r["conditions"] or "{}").get("setting", "")
            spans.setdefault(s, [i, i])[1] = i
        return aid, spans

    def figure(rn, fname, description, shows):
        src = runs_dir / rn / "figures" / fname
        dst = out_dir / "files" / "figures" / re.sub(r"[^A-Za-z0-9._-]", "_", f"{rn}-{fname}").replace("=", "-")
        shutil.copyfile(src, dst)
        return new_art({"kind": "figure", "function": "data_plot", "uri": f"files/figures/{dst.name}",
                        "version": sha(dst), "licence": "CC-BY-4.0", "availability": {"level": "public"},
                        "shows": shows,
                        "description": description + f" Drawn by lab run {rn} from its own results "
                                                     f"(lab/runs/{rn}/figures/{fname})."})

    load = lambda rn: json.load(open(runs_dir / rn / "results.json"))

    # ---------------------------------------------------------------- E13, PL1 (RN1): targets on 34 distributions
    e, res = ev[("RN1",)], load("RN1")
    csv_id, spans = run_csv("RN1", "PL1: the targets of RMSE, MAE, MAPE and sMAPE on 34 distributions")
    eid, allset = e["id"], settings(res)
    measures = ["RMSE", "MAE", "MAPE", "sMAPE"]
    t1 = table("T1", "The sMAPE target and related values derived from each distribution: the lowest local minimum of "
                     "expected sMAPE (rank 1), its quantile level and expected sMAPE, its one-step approximation, the "
                     "pull from above minus the pull from below at the median, the share of zeros at which the target "
                     "switches to 0, the quantile level of the log-mode, and the MAPE target",
               res, ["setting"], [(s,) for s in allset],
               [col("derived local minimum of expected sMAPE", rank=1),
                col("quantile level of a derived local minimum of expected sMAPE", rank=1),
                col("expected sMAPE at a derived local minimum", rank=1),
                col("one-step quantile level of the sMAPE target"),
                col("pull from above minus pull from below at the median"), col("zero threshold of sMAPE"),
                col("quantile level of the log-mode"),
                col("median of the outcome distribution"), col("mean of the outcome distribution")], csv_id)
    t2 = table("T2", "The second local minimum of expected sMAPE, where there is one (rank 2)", res, ["setting"],
               [(s,) for s in settings(res, "derived local minimum of expected sMAPE", rank=2)],
               [col("derived local minimum of expected sMAPE", rank=2),
                col("quantile level of a derived local minimum of expected sMAPE", rank=2),
                col("expected sMAPE at a derived local minimum", rank=2)], csv_id)
    t3 = table("T3", "Each measure's preferred candidate (the candidate forecast with the lowest expected score): its "
                     "quantile level, and the number of candidates whose expected score is infinite", res, ["setting"],
               [(s,) for s in allset],
               [col("quantile level of the preferred candidate", measure=m) for m in measures] +
               [col("candidates whose expected score is infinite", measure=m) for m in ("RMSE", "MAPE")], csv_id)
    pick = ["lognormal sigma=1", "inverse Gaussian shape=10", "bimodal equal modes",
            "lognormal sigma=1 with zeros p0=0.42"]
    t4 = table("T4", "The candidates a measure ranks first over test sets, for four of the distributions at test-set "
                     "sizes 6, 100 and 1000: the mean and the 5th, 50th and 95th percentiles of their quantile levels, "
                     "the share of test sets choosing the preferred candidate, and (MAPE) the test sets in which every "
                     "candidate's score is infinite; every distribution and size is in the run's CSV and its figure",
               res, ["setting", "measure", "test-set size"],
               [(s, m, n) for s in pick for m in measures for n in (6, 100, 1000)],
               [col("mean quantile level of the chosen candidates")] +
               [col("percentile of the chosen candidates' quantile levels", percentile=p) for p in (5, 50, 95)] +
               [col("share of test sets choosing the preferred candidate"),
                col("test sets in which every candidate's score is infinite")], csv_id)
    t5 = table("T5", "The MAPE target derived from each distribution (a median of the actuals weighted by 1/y), where "
                     "it is finite", res, ["setting"],
               [(s,) for s in allset if finite(value(res, "derived MAPE target", setting=s))],
               [col("derived MAPE target")], csv_id)
    t6 = table("T6", "The value of each measure's preferred candidate, for the four distributions of T4",
               res, ["setting", "measure"],
               [(s, m) for s in pick for m in measures
                if finite(value(res, "value of the preferred candidate", setting=s, measure=m))],
               [col("value of the preferred candidate")], csv_id)
    e["tables"], e["results"] = [t1, t2, t3, t4, t5, t6], []
    e.setdefault("artefacts", [])
    fig_ids = []
    for s in allset:
        fname = slug(s) + ".png"
        if not (runs_dir / "RN1" / "figures" / fname).exists():
            continue
        a, b = spans.get(s, [None, None])
        t3row = next((r["id"] for r in t3["rows"] if r["label"] == [s]), None)
        shows = [csv_id] + ([f"{eid}.T3.{t3row}"] if t3row else [])
        if s in pick:
            shows.append(f"{eid}.T4")
        fig_ids.append(figure("RN1", fname,
                              f"{s}: for each measure (RMSE, MAE, MAPE, sMAPE), one violin per test-set size showing "
                              f"how the candidate the measure ranks first is spread over 1000 test sets, as its quantile "
                              f"level, with a boxplot of its median, quartiles and 5th and 95th percentiles and a diamond "
                              f"at its mean; a dashed line marks the measure's preferred candidate and dotted lines the "
                              f"other measures'.", shows))
    e["artefacts"] = list(dict.fromkeys(e["artefacts"] + [csv_id] + fig_ids))
    e13_figs = fig_ids

    # ---------------------------------------------------------------- E14, PL3 (RN2): flatness
    e, res = ev[("RN2",)], load("RN2")
    csv_id, spans = run_csv("RN2", "PL3: flatness of expected sMAPE as the spread of ln y grows")
    eid, allset = e["id"], settings(res, "quantile level of the sMAPE target")
    t1 = table("T1", "The sMAPE target as the spread of ln y grows: its quantile level beside the log-mode's, the "
                     "expected sMAPE there and its large-spread approximation", res, ["setting"],
               [(s,) for s in allset],
               [col("quantile level of the sMAPE target"), col("quantile level of the log-mode"),
                col("expected sMAPE at the target"), col("large-spread approximation of expected sMAPE at the target"),
                col("local minima of expected sMAPE")], csv_id)
    sel = ["lognormal sigma=0.5", "lognormal sigma=4", "log-Gumbel right sigma=4"]
    t2 = table("T2", "Ranking agreement: the share of test sets that give a measure's target a lower score than a "
                     "forecast half a sigma above it, for three of the settings; every setting, offset and size is in "
                     "the run's CSV and its figure",
               res, ["setting", "test-set size"], [(s, n) for s in sel for n in (12, 48, 250, 1000)],
               [col("ranking agreement", measure=m, **{"offset in sigma": 0.5}) for m in measures], csv_id)
    t3 = table("T3", "How much worse a forecast half a sigma or one sigma below or above the target scores, as a "
                     "share of the score at the target, for lognormal outcomes",
               res, ["setting", "offset in sigma"],
               [(s, o) for s in ("lognormal sigma=0.5", "lognormal sigma=2", "lognormal sigma=6")
                for o in (-1.0, -0.5, 0.5, 1.0)],
               [col("excess expected score over the target, relative", measure=m) for m in measures] +
               [col("excess expected sMAPE over the target", measure="sMAPE")], csv_id)
    t4 = table("T4", "Where the sMAPE minimiser of a test set falls: the mean and the 5th, 50th and 95th percentiles of "
                     "its quantile level over test sets",
               res, ["setting", "test-set size"], [(s, n) for s in allset for n in (48, 1000)],
               [col("mean quantile level of the sMAPE minimiser of a test set")] +
               [col("percentile of the quantile level of the sMAPE minimiser of a test set", percentile=p)
                for p in (5, 50, 95)], csv_id)
    e["tables"], e["results"] = [t1, t2, t3, t4], []
    figs = [("flatness_target.png", "The quantile level of the sMAPE target against the spread sigma of ln y, for "
                                    "lognormal and log-Gumbel outcomes, beside the level of the log-mode: the target "
                                    "moves from the median to the mode of ln y as the spread grows.", ["T1"]),
            ("flatness_expected_sMAPE.png", "Expected sMAPE against the quantile level of the forecast, one curve per "
                                            "spread sigma, for each family: the curve flattens towards 200 as sigma "
                                            "grows.", ["T1"]),
            ("flatness_sensitivity.png", "Lognormal outcomes: how much worse a forecast half a sigma or one sigma "
                                         "away from each measure's target scores, as a share of the score at the "
                                         "target, against sigma.", ["T3"]),
            ("flatness_agreement.png", "Ranking agreement: the share of test sets that prefer each measure's target "
                                       "to a forecast away from it, against sigma, for several test-set sizes.", ["T2"]),
            ("flatness_minimiser.png", "Where the sMAPE minimiser of a test set falls, and how widely it spreads: its "
                                       "quantile level over test sets against sigma, for each family.", ["T4"])]
    ids = [figure("RN2", f, desc, [csv_id] + [f"{eid}.{t}" for t in tabs]) for f, desc, tabs in figs]
    e["artefacts"] = list(dict.fromkeys((e.get("artefacts") or []) + [csv_id] + ids))

    # ---------------------------------------------------------------- E15, PL2 (RN7): efficiency
    e, res = ev[("RN7",)], load("RN7")
    csv_id, spans = run_csv("RN7", "PL2: efficiency of the sMAPE minimiser against the sample median")
    eid, allset = e["id"], settings(res, "asymptotic relative efficiency")
    t1 = table("T1", "Relative efficiency of the sMAPE minimiser of a test set against the sample median, at test-set "
                     "sizes 6, 12, 100, 1000 and 10000, beside its large-sample value",
               res, ["setting"], [(s,) for s in allset],
               [col("relative efficiency", **{"test-set size": n}) for n in (6, 12, 100, 1000, 10000)] +
               [col("asymptotic relative efficiency")], csv_id)
    few = [s for s in allset if s.split("=")[-1] in ("0.25", "1", "6")]
    t2 = table("T2", "Bias and variance at test-set sizes 12 and 1000, for spreads 0.25, 1 and 6: the mean of "
                     "ln(sMAPE minimiser) minus ln(median), and the variances of ln(sample median) and ln(sMAPE "
                     "minimiser); every setting and size is in the run's CSV",
               res, ["setting", "test-set size"], [(s, n) for s in few for n in (12, 1000)],
               [col("mean of ln(sMAPE minimiser) minus ln(median)"), col("variance of ln(sample median)"),
                col("variance of ln(sMAPE minimiser)"), col("ratio of mean squared errors about the true median")],
               csv_id)
    e["tables"], e["results"] = [t1, t2], []
    ids = [figure("RN7", f"efficiency_{fam}.png",
                  f"{fam} outcomes: relative efficiency of the sMAPE minimiser of a test set against the sample "
                  f"median, against test-set size, one panel per spread, with its 99% band and the large-sample value "
                  f"as a dashed line; below 1 the sMAPE minimiser varies more than the sample median.",
                  [csv_id, f"{eid}.T1"]) for fam in ("lognormal", "log-Laplace")]
    e["artefacts"] = list(dict.fromkeys((e.get("artefacts") or []) + [csv_id] + ids))

    # ---------------------------------------------------------------- E16, PL11 (RN10): the Suilin loss
    # RN10 ran when the paper called the Suilin loss msMAPE: its results name the measure "msMAPE". The tables name
    # it as the paper now does; the CSV keeps the record's names.
    e, res = ev[("RN10",)], load("RN10")
    csv_id, spans = run_csv("RN10", "PL11: the Suilin and sMAPE targets",
                            " The run named the Suilin loss msMAPE, as version 2 of the paper did: in this file the "
                            "measure \"msMAPE\" and the quantity \"relative difference of the msMAPE target from the "
                            "sMAPE target\" are the Suilin loss's.")
    eid, allset = e["id"], settings(res, "target")
    SU = {"measure": "Suilin"}
    DIFF = "relative difference of the Suilin target from the sMAPE target"
    t1 = table("T1", "The Suilin target beside the sMAPE target, with their quantile levels and the relative "
                     "difference of the Suilin target from the sMAPE target", res, ["setting"], [(s,) for s in allset],
               [col("target", measure="sMAPE"),
                col("target", _as={"variable": "target", "conditions": SU}, measure="msMAPE"),
                col("quantile level of the target", measure="sMAPE"),
                col("quantile level of the target", _as={"variable": "quantile level of the target", "conditions": SU},
                    measure="msMAPE"),
                col("relative difference of the msMAPE target from the sMAPE target", _as=DIFF),
                col("share of zero outcomes")], csv_id)
    t2 = table("T2", "The share of zeros at which each measure's target switches to the forecast 0, for lognormal "
                     "outcomes with sigma 1", res, ["setting"],
               [(s,) for s in settings(res, "zero threshold")],
               [col("zero threshold", measure="sMAPE"),
                col("zero threshold", _as={"variable": "zero threshold", "conditions": SU}, measure="msMAPE")], csv_id)
    e["tables"], e["results"] = [t1, t2], []
    rn11 = runs_dir / "RN11" / "figures" / "msmape_targets.png"
    dst = out_dir / "files" / "figures" / "RN11-msmape_targets.png"
    shutil.copyfile(rn11, dst)
    fid = new_art({"kind": "figure", "function": "data_plot", "uri": f"files/figures/{dst.name}",
                   "version": sha(dst), "licence": "CC-BY-4.0", "availability": {"level": "public"},
                   "shows": [f"{eid}.T1", f"{eid}.T2"],
                   "description": "The Suilin loss beside sMAPE: near zero, the two targets as a share of the median "
                                  "against the share of zeros, for lognormal outcomes at medians 1 and 100, with the "
                                  "zero thresholds at median 1; away from zero, the relative difference of the Suilin "
                                  "target from the sMAPE target for each distribution; on Poisson counts, the two "
                                  "targets against the mean. The figure was drawn when the paper called the Suilin "
                                  "loss msMAPE, and labels it so. Drawn by lab run RN11 from the recorded results of "
                                  "RN10 (lab/runs/RN11/figures/msmape_targets.png)."})
    e["artefacts"] = list(dict.fromkeys((e.get("artefacts") or []) + [csv_id, fid]))

    # version 2's artefacts end here; what version 3 adds takes the numbers after them
    if before_new:
        before_new()

    # ---------------------------------------------------------------- E34, PL14 (RN15): intermittent demand
    e, res = ev[("RN15",)], load("RN15")
    csv_id, spans = run_csv("RN15", "PL14: the targets of sMAPE, msMAPE and the Suilin loss on intermittent demand")
    eid, allset = e["id"], settings(res, "target")
    counts = settings(res, "mode of the outcomes")
    five = ["sMAPE", "msMAPE nu=0.1", "msMAPE nu=0.5", "msMAPE nu=1", "Suilin"]
    sizes = (12, 52, 365)
    t1 = table("T1", "The forecast each measure rewards (its target, computed from the distribution itself) on each "
                     "intermittent demand setting, beside the setting's share of zeros, median, mode and mean; the "
                     "measures are sMAPE, msMAPE with nu = 0.1, 0.5 and 1, and the Suilin loss",
               res, ["setting"], [(s,) for s in allset],
               [col("share of zero outcomes"), col("median of the outcomes"), col("mode of the outcomes"),
                col("mean of the outcomes")] + [col("target", measure=m) for m in five], csv_id)
    t2 = table("T2", "The share of zeros at which each measure's target switches to the forecast 0, for lognormal "
                     "demand sizes (median 1, sigma 1) with zeros", res, ["setting"],
               [(s,) for s in settings(res, "zero threshold")],
               [col("zero threshold", measure=m) for m in five], csv_id)
    t3 = table("T3", "The quantile level of each measure's target and the expected score there (200 times the "
                     "expected term), for each setting", res, ["setting"], [(s,) for s in allset],
               [col("quantile level of the target", measure=m) for m in five] +
               [col("expected score at the target", measure=m) for m in five], csv_id)
    t4 = table("T4", "Test sets: the share of 2000 test sets of n independent actuals in which the forecast with the "
                     "lowest score is 0, for each setting, n and measure", res, ["setting", "n"],
               [(s, n) for s in allset for n in sizes],
               [col("share of test sets whose minimiser is 0", measure=m) for m in five], csv_id)
    t5 = table("T5", "Test sets: the share of 2000 test sets of n independent actuals in which the forecast with the "
                     "lowest score is the measure's target of T1, for each count setting, n and measure",
               res, ["setting", "n"], [(s, n) for s in counts for n in sizes],
               [col("share of test sets whose minimiser is the target", measure=m) for m in five], csv_id)
    t6 = table("T6", "Test sets: the median, over 2000 test sets of n independent actuals, of the forecast with the "
                     "lowest score, for each setting, n and measure", res, ["setting", "n"],
               [(s, n) for s in allset for n in sizes],
               [col("median of the test-set minimisers", measure=m) for m in five], csv_id)
    e["tables"], e["results"] = [t1, t2, t3, t4, t5, t6], []
    fid = figure("RN15", "msmape_intermittent.png",
                 "The forecast with the lowest expected score (the target) against the share of zero actuals, on "
                 "intermittent demand: one panel per family of settings (Poisson counts, negative binomial counts, "
                 "lognormal demand sizes with zeros), one line per measure (sMAPE, msMAPE with nu = 0.1, 0.5 and 1, "
                 "the Suilin loss).", [csv_id, f"{eid}.T1"])
    e["artefacts"] = list(dict.fromkeys((e.get("artefacts") or []) + [csv_id, fid]))

    # columns named by the definition of their variable, as the export names them
    defmap = {x["term"]: x["id"] for x in d["definitions"] if x.get("kind") == "variable"}
    for e in d["evidence"]:
        for t in e.get("tables") or []:
            t["columns"] = [({**c, "variable": defmap.get(c["variable"], c["variable"])} if isinstance(c, dict)
                             else defmap.get(c, c)) for c in t["columns"]]
    return e13_figs
