"""Complete the exported draft of smape-target into version 4 of the claimset (claimset v0.44).

    python3 build4.py SRC V3 OUT

SRC is the draft `claimset lab export` wrote (its files/ beside it), V3 the published version 3 (its files/ beside
it), OUT the claimset to write (its files go to files/ beside it). Version 4 is version 3 with what the record added
since: the claims of H53 (in place of version 3's C6), H54, H49, H52 and the estimate of PL19; their definitions,
proofs and evidence; the exploration counts and the verdicts that give no claim, each with a one-line statement; the
contributors' recent work; new limitations and references. Everything else of version 3 is kept as published.
The study's files are read from claimset/smape-target/ of this repository and its lab's runs from its lab/runs."""
import copy
import csv
import hashlib
import json
import math
import re
import shutil
import sys
from pathlib import Path

SRC, V3, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
REPO = Path(__file__).resolve().parents[3]
STUDY = REPO / "claimset" / "smape-target"
RUNS = STUDY / "lab" / "runs"
X = json.loads(SRC.read_text())
P = json.loads(V3.read_text())
d = copy.deepcopy(P)
OUTDIR = OUT.parent
FILES = OUTDIR / "files"
if FILES.exists():
    shutil.rmtree(FILES)
shutil.copytree(V3.parent / "files", FILES)


def sha(path):
    return "sha256:" + hashlib.sha256(Path(path).read_bytes()).hexdigest()


def walk(x, f):
    if isinstance(x, dict):
        return {k: walk(v, f) for k, v in x.items()}
    if isinstance(x, list):
        return [walk(v, f) for v in x]
    return f(x) if isinstance(x, str) else x


# ------------------------------------------------------------------ the IDs readers have seen
# Claims: matched by the hypothesis (or the estimate's plan) they report; the new ones take the next free numbers.
def key(c):
    if c.get("hypothesis"):
        return ("H", c["hypothesis"]["lab_event"])
    if c.get("estimate"):
        return ("PL", c["estimate"]["lab_plan"])
    return ("A", c["statement"])                      # assumption claims: matched by their text


p_by_key = {key(c): c["id"] for c in P["claims"]}
NEW_CLAIMS = {"C6": "C6", "C17": "C22", "C18": "C23", "C19": "C24", "C25": "C25"}
CMAP = {}
for c in X["claims"]:
    if c["id"] in NEW_CLAIMS:
        CMAP[c["id"]] = NEW_CLAIMS[c["id"]]
    else:
        CMAP[c["id"]] = p_by_key[key(c)]
assert {key(c) for c in X["claims"] if c["id"] not in NEW_CLAIMS} <= set(p_by_key)
assert max(int(c["id"][1:]) for c in P["claims"]) == 21

# Definitions: by term; the new ones after version 3's last (D35).
p_terms = {x["term"]: x["id"] for x in P["definitions"]}
DMAP, nxt = {}, 36
for x in X["definitions"]:
    if x["term"] in p_terms:
        DMAP[x["id"]] = p_terms[x["term"]]
    else:
        DMAP[x["id"]] = f"D{nxt}"
        nxt += 1
NEWDEF = {DMAP[x["id"]]: x for x in X["definitions"] if DMAP[x["id"]] not in p_terms.values()}
assert sorted(NEWDEF, key=lambda s: int(s[1:])) == ["D36", "D37", "D38", "D39", "D40", "D41"], NEWDEF.keys()
TERM_TM = {"TM8": "D8", "TM12": "D12", "TM29": "D2", "TM30": "D3", "TM6": "D6", "TM10": "D10", "TM11": "D11",
           "TM31": "D5", "TM41": "D36", "TM42": "D37", "TM43": "D38", "TM44": "D39", "TM45": "D40", "TM46": "D41"}


def remap(s):
    s = re.sub(r"\[\[(D\d+):", lambda m: f"[[{DMAP.get(m.group(1), m.group(1))}:", s)
    s = re.sub(r"\[\[(TM\d+):", lambda m: f"[[{TERM_TM[m.group(1)]}:", s)
    return s


def dref(x):
    return DMAP.get(x, x) if re.fullmatch(r"D\d+", x) else x


# ------------------------------------------------------------------ definitions
NEWDEF_TEXT = {}
for did, x in NEWDEF.items():
    y = {k: v for k, v in x.items() if k not in ("id", "source")}
    y["id"] = did
    y["definition"] = remap(y["definition"])
    NEWDEF_TEXT[did] = y
NEWDEF_TEXT["D37"]["symbol"] = "$h(0)\\,\\mathrm E[Z^2]/\\mathrm E|Z|$"
NEWDEF_TEXT["D36"]["symbol"] = "$V$"
NEWDEF_TEXT["D38"]["source"] = "B17"
NEWDEF_TEXT["D40"]["symbol"] = "$\\mathrm{ARE}$"
NEWDEF_TEXT["D40"]["definition"] += (" In symbols, $\\mathrm{ARE} = V_{\\mathrm{MAE}}/V_{\\mathrm{sMAPE}}$, with $V$ the "
                                     "[[D36:asymptotic variance of a minimiser]] of each.")
NEWDEF_TEXT["D41"]["source"] = "B19"
NEWDEF_TEXT["D41"]["symbol"] = "$H$"
NEWDEF_TEXT["D41"]["definition"] = (
    "For a random variable $X$ with a finite positive variance $\\sigma^2$: $H = \\mathrm E|X - \\mathrm E X|^{1/2}/"
    "\\sigma^{1/2}$, the half absolute moment scaled by the square root of the standard deviation. It is defined and "
    "studied in entrylib:2026.000016 (D4). It lies in $(0, 1]$. It is lower when the values lie unevenly about the mean, "
    "many near it and a few far out. It is 0.877 for a uniform, 0.822 for a normal and 0.745 for a Laplace distribution. "
    "The Extremity coefficient of the same claimset (D5), $2H^p - 1$ with $p \\approx 5.30$, rises with $H$, so the two "
    "order distributions alike. Here $X = \\ln y$.")
d["definitions"] += [NEWDEF_TEXT[k] for k in sorted(NEWDEF_TEXT, key=lambda s: int(s[1:]))]

# ------------------------------------------------------------------ references
d["references"] += [
    {"id": "B17", "target": "10.1198/jasa.2011.r10138", "id_type": "doi"},
    {"id": "B18", "target": "Rockafellar, R. T. (1970). Convex Analysis. Princeton University Press.",
     "id_type": "reference"},
    {"id": "B19", "target": "entrylib:2026.000016#D4", "id_type": "claimset"},
    {"id": "B20", "target": "entrylib:2026.000016#C19", "id_type": "claimset"},
]

# ------------------------------------------------------------------ artefacts
arts = {a["id"]: a for a in d["artefacts"]}
xarts = {a["id"]: a for a in X["artefacts"]}


def set_file(aid, src, name, **extra):
    dst = FILES / name
    dst.parent.mkdir(parents=True, exist_ok=True)
    if Path(src).resolve() != dst.resolve():
        shutil.copyfile(src, dst)
    arts[aid].update({"uri": f"files/{name}", "version": sha(dst)}, **extra)


# A1, A2: the proofs as the record names them now (AR15, AR16)
for aid, name in (("A1", "proofs.md"), ("A2", "proofs.pdf")):
    rec = [a for a in X["artefacts"] if a["uri"].endswith(name)][-1]
    set_file(aid, STUDY / "files" / name, name)
    assert arts[aid]["version"] == rec["version"], (aid, "the file is not the one the record names")
arts["A1"]["description"] = (
    "Proofs of the formal claims, version 4, in Markdown with LaTeX. Sections H1 to H14 are those of version 3 and "
    "prove C1, C2, C3, C9, C4, C5, C7, C11, C8, C10, C12, C18 and C19, in that order (section H7 proves version 3's C6, "
    "for two mirror-image groups, which version 4 no longer claims). New sections H47 to H50 prove C6, C22, C23 and "
    "C24. Sections H1 to H11 were checked independently by lab actor AC4 (proof-check-report.md) and by the "
    "supervisor, H12 to H14 by AC6 (CK27 to CK31, CK34 to CK37), and H47 to H50 by AC4 (CK38 to CK44), each corrected "
    "as its checker proposed.")
arts["A2"]["description"] = "The same proofs, version 4, rendered as PDF."
# A6: every result of the record, as export wrote it
master = next(a for a in X["artefacts"] if a["uri"].endswith(".results.csv") and "/E" not in a["uri"])
set_file("A6", SRC.parent / master["uri"], "smape-target.results.csv")
# A7: the lab record, at the commit of this export
lab = next(a for a in X["artefacts"] if a["uri"].endswith("lab.jsonl"))
arts["A7"].update({"uri": lab["uri"], "version": lab["version"]})
# A54: the study's code at the commit of the record's last run (RN24), which holds every script of the earlier runs
CODE = json.loads((STUDY / "lab" / "lab.jsonl").read_text().splitlines()[-1]) if False else None
for line in (STUDY / "lab" / "lab.jsonl").read_text().splitlines():
    e = json.loads(line)
    if e.get("type") == "run" and e.get("id") == "RN24":
        CODE = e["body"]["code_version"]
assert CODE and re.fullmatch(r"[0-9a-f]{40}", CODE), CODE
arts["A54"].update({"uri": f"https://github.com/entrylib/study-planning/tree/{CODE}/claimset-code/smape-target",
                    "version": CODE})
arts["A54"]["description"] = arts["A54"]["description"].rstrip(".") + (
    "; and, for version 4, the experiment of the variance formulas (experiment_v4.py) and its figure (plot_v4.py).")
next_a = max(int(a[1:]) for a in arts) + 1


def new_artefact(**a):
    global next_a
    aid = f"A{next_a}"
    next_a += 1
    arts[aid] = {"id": aid, **a}
    return aid


# ------------------------------------------------------------------ evidence of the new runs
EV_RUN = {"RN20": "E42", "RN22": "E43", "RN21": "E44"}
X_EV = {"RN20": "E18", "RN22": "E17", "RN21": "E25"}
xev = {e["id"]: e for e in X["evidence"]}
new_ev = {}
for run, eid in EV_RUN.items():
    e = copy.deepcopy(xev[X_EV[run]])
    old = e["id"]
    e["id"] = eid
    e = walk(e, lambda s: s.replace(f"{old}.", f"{eid}."))
    amap = {}
    for t in e["tables"]:
        t["columns"] = [dref(c) for c in t["columns"]]
        if "extremity" in t["columns"]:
            t["description"] = t["description"].replace("extremity", "extremity (Hypotailedness, D41)", 1)
        t["description"] = re.sub(r"\bD(\d+)\b", lambda m: dref(m.group(0)), t["description"])
        if t.get("rows_file"):
            src = SRC.parent / xarts[t["rows_file"]]["uri"]
            xa = xarts[t["rows_file"]]
            aid = new_artefact(kind="data", uri="", availability={"level": "public"}, licence="CC-BY-4.0",
                               description=xa["description"].replace(f"{old}.", f"{eid}."),
                               columns=[{**col, **({"variable": dref(col["variable"])} if "variable" in col else {})}
                                        for col in xa["columns"]])
            set_file(aid, src, f"smape-target.{eid}.{t['id']}.csv")
            amap[t["rows_file"]] = aid
            t["rows_file"] = aid
        t["artefact"] = "A6"
    e["artefacts"] = ["A6"] + [amap[a] for a in e["artefacts"] if a in amap]
    if "setup" in e:
        e["setup"] = {k: [dref(x) for x in v] for k, v in e["setup"].items()}
    e["description"] = remap(re.sub(r"\(D(\d+)\)", lambda m: f"({dref('D' + m.group(1))})", e["description"]))
    new_ev[eid] = e

# notes the record asks to be said where the values are read
new_ev["E43"]["description"] += (
    " The share of test-set minimisers in the target's group (table T10) is an estimate, not a check: the "
    "hypothesis says nothing about finite test sets (the supervisor's decision after RN20, deviation DV6).")
new_ev["E44"]["description"] += (
    " The values named 'extremity' in this run are Hypotailedness (D41) values; the name was corrected after the run "
    "(deviation DV7), the values are unchanged.")
new_ev["E42"]["description"] += (
    " Three settings converge slowly and are reported as such, not as checks of the formulas (deviation DV5, see L8). "
    "One check against an exact small-sample value failed through a defect of the code, since fixed (deviation DV4).")

# ------------------------------------------------------------------ the two tables of the paper (from RN21)
T4 = {(r["distribution"], r["measure"]): float(r["asymptotic variance by formula"])
      for r in csv.DictReader(open(FILES / "smape-target.E44.T4.csv"))}
FAM = [("lognormal", "lognormal sigma 0.5", "lognormal sigma 2.0", "sd of ln y 0.5", "sd of ln y 2"),
       ("log-Laplace", "log-Laplace scale 0.25", "log-Laplace scale 1.0", "scale 0.25", "scale 1"),
       ("log-uniform", "log-uniform half-width 0.5", "log-uniform half-width 2.0", "half-width 0.5", "half-width 2"),
       ("gamma", "gamma shape 20", "gamma shape 0.5", "shape 20", "shape 0.5"),
       ("inverse Gaussian", "inverse Gaussian mean 1.0 shape 1.0", "inverse Gaussian mean 1.0 shape 0.05",
        "mean 1, shape 1", "mean 1, shape 0.05"),
       ("Weibull", "Weibull shape 2.0", "Weibull shape 0.5", "shape 2", "shape 0.5")]
MEAS = ("RMSE", "MAE", "MAPE", "sMAPE")
rows1 = []
for fam, d1, d2, l1, l2 in FAM:
    for dist, lab_ in ((d1, l1), (d2, l2)):
        vals = [T4.get((dist, m)) for m in MEAS]
        best = min((v, m) for v, m in zip(vals, MEAS) if v is not None)[1]
        rows1.append({"id": f"r{len(rows1) + 1}", "label": [fam, lab_, best],
                      "values": [None if v is None else round(v, 4 if v < 10 else 2) for v in vals]})
new_ev["E44"]["tables"].append({
    "id": "T11", "description": "Table 1. The asymptotic variance V of the log minimiser of each measure (C23), for "
    "six families at a narrow and a wide spread, and the measure with the smallest V. A missing value is an infinite "
    "V: a moment the formula needs is infinite. From table T4 of this run.",
    "row_labels": ["family", "spread", "most precise"], "columns": ["V of RMSE", "V of MAE", "V of MAPE", "V of sMAPE"],
    "rows": rows1})
T5 = {(r["shape"], r["spread"]): r for r in csv.DictReader(open(FILES / "smape-target.E44.T5.csv"))}
e44 = new_ev["E44"]
T6 = {r["label"][0]: r["values"] for r in next(t for t in e44["tables"] if t["id"] == "T6")["rows"]}
assert next(t for t in e44["tables"] if t["id"] == "T6")["columns"] == ["D37", "extremity"]
T1 = {r["label"][0]: r["values"] for r in next(t for t in e44["tables"] if t["id"] == "T1")["rows"]}
EXCESS = {"uniform": -1.2, "normal": 0.0, "logistic": 1.2, "Student t 6": 3.0, "Laplace": 3.0,
          "generalised normal 0.75": round(math.gamma(5 / 0.75) * math.gamma(1 / 0.75) / math.gamma(3 / 0.75) ** 2
                                           - 3, 3)}
ARE1 = {"normal": T1["lognormal sigma 1.0"][1], "Laplace": T1["log-Laplace scale 1.0"][1],
        "uniform": T1["log-uniform half-width 1.0"][1]}
rows2 = []
for shape in ("uniform", "normal", "logistic", "Student t 6", "Laplace", "generalised normal 0.75"):
    small = float(T5[(shape, "0.01")]["(ARE - 1)/(sigma^2/2)"])
    big = T5.get((shape, "100.0"), {}).get("sigma times ARE")
    rows2.append({"id": f"r{len(rows2) + 1}", "label": [shape],
                  "values": [round(T6[shape][0], 3), round(T6[shape][1], 3), EXCESS[shape], round(small, 3) + 0.0,
                             round(ARE1[shape], 3) if shape in ARE1 else None,
                             float(f"{float(big):.3g}") if big else None]})
e44["tables"].append({
    "id": "T12", "description": "Table 2. For six shapes of ln y: the peakedness index (D37) and Hypotailedness (D41). "
    "Then the excess kurtosis, a closed form the run did not compute. Then the small-spread coefficient "
    "(ARE - 1)/(sigma^2/2) at spread 0.01, the ARE at spread 1 and sigma times the ARE at spread 100. With ARE above "
    "1, the sMAPE minimiser is the more precise. Values from tables T1, T5 and T6 of this run; the run did not compute the ARE at spread 1 for "
    "the logistic, Student t and generalised normal shapes, nor sigma ARE for the uniform and generalised normal.",
    "row_labels": ["shape"], "columns": ["D37", "D41", "excess kurtosis", "(ARE - 1)/(sigma^2/2) at spread 0.01",
                                         "ARE at spread 1", "sigma times ARE at spread 100"],
    "rows": rows2})

# the smallest n from which the simulated variance stays within 10% of the formula (from T2 of this run)
sim = {}
for r in csv.DictReader(open(FILES / "smape-target.E44.T2.csv")):
    if r["n times variance over formula"]:
        sim.setdefault((r["distribution"], r["measure"]), {})[int(r["n"])] = float(r["n times variance over formula"])
FROM = {}
for (dist, m), rs in sim.items():
    ns = sorted(rs)
    FROM[(dist, m)] = next((n for i, n in enumerate(ns) if all(abs(rs[k] - 1) <= 0.1 for k in ns[i:])), None)


def by_n(m, n):
    xs = [v for (dd, mm), v in FROM.items() if mm == m]
    return sum(1 for v in xs if v is not None and v <= n), len(xs)


# ------------------------------------------------------------------ the figure (RN24, exploratory)
res24 = json.loads((RUNS / "RN24" / "results.json").read_text())
fig = new_artefact(kind="figure", uri="", availability={"level": "public"}, licence="CC-BY-4.0",
                   description="Figure 1. The asymptotic variance V of the log minimiser of RMSE, MAE, MAPE and sMAPE "
                               "against the spread, on log-log axes, for lognormal, log-Laplace and log-uniform "
                               "actuals (exploratory run RN24, by quadrature with the formulas of C23). Where a moment "
                               "the formula needs is infinite, the line stops.",
                   function="data_plot")
set_file(fig, RUNS / "RN24" / "figures" / "precision_growth.png", "figures/RN24-precision_growth.png")
csv24 = FILES / "smape-target.E45.T1.csv"
with open(csv24, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["family", "spread", "measure", "asymptotic variance by formula"])
    for r in res24:
        c = r["conditions"]
        w.writerow([c["family"], c["spread"], c["measure"], r["value"]])
rows24 = new_artefact(kind="data", uri="", availability={"level": "public"}, licence="CC-BY-4.0",
                      description=f"The {len(res24)} rows of table E45.T1: the values Figure 1 shows (lab run RN24).",
                      columns=[{"name": "family", "dimension": "family"}, {"name": "spread", "dimension": "spread"},
                               {"name": "measure", "dimension": "measure"},
                               {"name": "asymptotic variance by formula",
                                "variable": "asymptotic variance by formula"}])
set_file(rows24, csv24, "smape-target.E45.T1.csv")
arts[fig]["shows"] = [rows24, "E45.T1"]
run24 = next(json.loads(x) for x in (STUDY / "lab" / "lab.jsonl").read_text().splitlines()
             if json.loads(x).get("id") == "RN24")
new_ev["E45"] = {
    "id": "E45", "profile": "quant", "design": "other",
    "description": "Exploratory lab run RN24, after the confirmatory runs. It computes the asymptotic variance V of "
                   "each measure's minimiser by quadrature with the formulas of C23. It does so on 23 spreads from 0.1 "
                   "to 8, for lognormal, log-Laplace and log-uniform actuals, shown as Figure 1. It illustrates C24 and "
                   "decides nothing.",
    "availability": {"level": "public"}, "lab_runs": ["RN24"],
    "setup": {"studied": ["D10", "D11", "D12", "D2"], "cases": ["D39"], "lens": ["D36"]},
    "tables": [{"id": "T1", "description": "V by family, spread and measure (lab run RN24)",
                "row_labels": ["family", "spread", "measure"], "columns": ["asymptotic variance by formula"],
                "rows_file": rows24, "n_rows": len(res24), "artefact": rows24}],
    "artefacts": [fig, rows24],
    "execution": {"code_version": CODE, "software": run24["body"].get("environment", {}).get("software",
                                                                                             {"python": "3.13"})}}

# ------------------------------------------------------------------ proofs -> formal evidence
PAPER = {"H1": "C1", "H2": "C2", "H3": "C3", "H4": "C9", "H5": "C4", "H6": "C5", "H8": "C7", "H9": "C11",
         "H10": "C8", "H11": "C10", "H12": "C12", "H13": "C18", "H14": "C19",
         "H47": "C6", "H48": "C22", "H49": "C23", "H50": "C24"}
PROOF_ID = {"C6": "E38", "C22": "E39", "C23": "E40", "C24": "E41"}
MATH = re.compile(r"(\$\$.+?\$\$|\$[^$]+\$)", re.S)


def sections(text):
    out, cur, name = {}, [], None
    for line in text.splitlines():
        m = re.match(r"^## (H\d+)\. (.*)", line)
        if m or line.startswith("## References"):
            if name:
                out[name] = "\n".join(cur)
            name, cur = (m.group(1) if m else None), [line]
        elif name:
            cur.append(line)
    return out


def tidy(s):
    s = re.sub(r"\s+", " ", s).strip()
    s = re.sub(r"^\*(.+?\.)\*\s", r"\1 ", s)
    parts = MATH.split(s)
    for i in range(0, len(parts), 2):
        parts[i] = re.sub(r"\*\*([^*]+)\*\*", r"\1", parts[i])
        parts[i] = re.sub(r"\*([^*]+)\*", r"\1", parts[i])
        assert "*" not in parts[i], parts[i]
    s = "".join(parts)
    s = re.sub(r"\s*\$\$\s*(.+?)\s*\$\$\s*", lambda m: "\n$$" + m.group(1) + "$$\n", s).strip()
    s = s.replace("This generalises H7, whose groups were mirror images.",
                  "This generalises version 3 of C6, whose groups were mirror images.")
    s = re.sub(r"\bH(\d+)\b", lambda m: PAPER.get(f"H{m.group(1)}", m.group(0)), s)
    return s


RULES = [("dominated convergence", "dominated convergence"), ("bounded convergence", "bounded convergence"),
         ("Taylor", "Taylor's theorem"), ("van der Vaart", "cited theorem"), ("Rockafellar", "cited theorem"),
         ("central limit", "cited theorem"), ("concave", "concavity"), ("By symmetry", "symmetry"),
         ("mirror", "symmetry"), ("continuous", "continuity"), ("increasing", "monotonicity"),
         ("convex", "convexity"), ("Write ", "definition"), ("Take ", "case analysis")]
STEPS = r"steps? (\d+(?:(?:, and |, | and | to )\d+)*)"


def rule_of(text, refs):
    for k, rule in RULES:
        if k in text:
            return rule
    return "earlier result" if refs else ("algebra" if "$" in text else "logical argument")


def step_numbers(listed):
    out = []
    for a, to, b in re.findall(r"(\d+)(?:( to )(\d+))?", listed):
        out += list(range(int(a), int(b) + 1)) if to else [int(a)]
    return out


def steps_of(body):
    lines = body.split("**Proof.**", 1)[1].splitlines()
    intro, steps, closing, cur = [], [], [], None
    for line in lines:
        st = line.strip()
        if not st:
            continue
        if closing or (line.startswith("So ") and not line.startswith(" ")):
            closing.append(st)
            continue
        m = re.match(r"^(\d+)\.\s+(.*)", line)
        if m:
            cur = [f"S{m.group(1)}", m.group(2)]
            steps.append(cur)
        elif cur is None:
            intro.append(st)
        else:
            cur[1] += " " + st
    out = ([["S0", " ".join(intro)]] if intro else []) + steps
    if closing:
        out.append(["SQ", " ".join(closing).replace("∎", "").strip()])
    return out


def formal_evidence(eid, cid, paper_h, body):
    shows = re.search(r"\*What it shows:\*\s*(.*?)\n\n", body, re.S).group(1)
    raw = steps_of(body)
    steps, ids = [], []
    shift = 1 if raw and raw[0][0] == "S0" else 0
    for sid, text in raw:
        t = tidy(text)
        if shift:
            t = re.sub(r"(?<!C\d, )(?<!C\d\d, )\b" + STEPS,
                       lambda m: re.sub(r"\d+", lambda n: str(int(n.group(0)) + 1), m.group(0)), t)
        refs = sorted({PAPER[f"H{x}"] for x in re.findall(r"\bH(\d+)\b", text) if f"H{x}" in PAPER} - {cid})
        local = re.sub(r"H\d+, " + STEPS, "", text)
        cited = [f"S{x}" for listed in re.findall(STEPS, local) for x in step_numbers(listed) if f"S{x}" in ids]
        by = list(dict.fromkeys(cited + refs))
        by += [b for k, b in (("van der Vaart", "B15"), ("Gneiting", "B17"), ("Rockafellar", "B18")) if k in text]
        if sid == "SQ":
            rule, by = "conclusion", [i for i in ids if i != "S0"][-1:] or ids[-1:]
        elif sid == "S0":
            rule = "setup"
        else:
            rule = rule_of(text, refs)
            if not by:
                by = [i for i in ids if i != "S0"][-1:] or [i for i in ids if i == "S0"]
        step = {"id": sid, "statement": t, "rule": rule, "by": by}
        steps.append(step)
        ids.append(sid)
    new = {s["id"]: f"S{i + 1}" for i, s in enumerate(steps)}
    for s in steps:
        s["id"] = new[s["id"]]
        if "by" in s:
            s["by"] = [new.get(b, b) for b in s["by"]]
    return {"id": eid, "profile": "formal",
            "description": f"Proof of {cid} (section {paper_h} of the proof document A1, rendered in A2). "
                           f"{tidy(shows)[:1].upper() + tidy(shows)[1:]}",
            "foundation": "Kolmogorov probability axioms and real analysis",
            "steps": steps, "artefacts": ["A1", "A2"], "availability": {"level": "public"}}


secs = sections((STUDY / "files" / "proofs.md").read_text())
for h, cid in (("H47", "C6"), ("H48", "C22"), ("H49", "C23"), ("H50", "C24")):
    new_ev[PROOF_ID[cid]] = formal_evidence(PROOF_ID[cid], cid, h, secs[h])

# version 3's proof of C6 (E24) proved the mirror-image statement, which version 4 no longer claims
d["evidence"] = [e for e in d["evidence"] if e["id"] != "E24"]
d["evidence"] += [new_ev[k] for k in sorted(new_ev, key=lambda s: int(s[1:]))]

# ------------------------------------------------------------------ claims
xcl = {c["id"]: c for c in X["claims"]}
SUPPORT = {
    "C6": [{"kind": "formal_proof", "evidence": "E38"},
           {"kind": "empirical_quant", "evidence": "E43", "results": ["E43.T9", "E43.T10", "E43.T1"]}],
    "C22": [{"kind": "formal_proof", "evidence": "E39"},
            {"kind": "empirical_quant", "evidence": "E43",
             "results": ["E43.T5", "E43.T6", "E43.T7", "E43.T8", "E43.T1", "E43.T3"]}],
    "C23": [{"kind": "formal_proof", "evidence": "E40"},
            {"kind": "empirical_quant", "evidence": "E42", "results": ["E42.T4", "E42.T2"]},
            {"kind": "cited", "reference": "B15", "locator": "Theorems 3.1, 5.7 and 5.23, Corollary 21.5"},
            {"kind": "cited", "reference": "B17", "locator": "the MAPE target, the (-1)-median"},
            {"kind": "cited", "reference": "B18", "locator": "Theorem 10.8"}],
    "C24": [{"kind": "formal_proof", "evidence": "E41"},
            {"kind": "empirical_quant", "evidence": "E43", "results": ["E43.T5", "E43.T4"]},
            {"kind": "empirical_quant", "evidence": "E45", "results": ["E45.T1"]}],
}
NAMES = {
    "C6": ("Two separated groups of actuals: the target lies in one of them",
           "When the log actuals fall in two separated groups, expected sMAPE is concave across the gap. So the "
           "target lies in one of the groups. With equal weights, every forecast in the gap is a median, and none is "
           "the target."),
    "C22": ("The peakedness index decides which estimate of the median is more precise",
            "Where sMAPE and MAE both target the median, the sMAPE minimiser is the more precise for small spreads "
            "when the peakedness index exceeds 1. The kurtosis does not decide it. For large spreads it loses its "
            "precision. A narrow core of actuals at the median makes it the more precise."),
    "C23": ("How precisely a test set locates each measure's target",
            "For any continuous distribution of positive actuals, the log of each measure's test-set minimiser is "
            "asymptotically normal about the log of the measure's target. Its variance has a closed form for RMSE, "
            "MAE, MAPE and sMAPE."),
    "C24": ("How each measure's precision falls as the actuals spread out",
            "When the log actuals have a fixed shape and a growing spread, MAE's variance grows like the square of "
            "the spread for every shape. RMSE's and MAPE's depend on the tails, sMAPE's on the peak. No measure is "
            "always the slowest or the fastest to lose precision."),
    "C25": ("Each measure's precision, computed and simulated",
            "On 21 distributions, the formulas of C23 and 2000 simulated test sets agree, except in three slow "
            "cases. RMSE is the most precise on most families; MAE on lognormal actuals with a wide spread and on "
            "log-Laplace actuals. Hypotailedness orders the six shapes as the peakedness index does."),
}
c25_statement = (
    "On 21 distributions of positive actuals, the [[D36:asymptotic variance of a minimiser]] $V$ of each of "
    "[[D10:RMSE]], [[D11:MAE]], [[D12:MAPE]] and [[D2:sMAPE]] was computed by the formulas of C23. It was also "
    "estimated from 2000 simulated test sets of $n$ actuals each, for $n$ from 6 to 1000. "
    "The two agree within the plan's tolerance at $n = 1000$ or 4000, except in three slow cases: RMSE on lognormal "
    "actuals with $\\sigma = 2$, and MAE and sMAPE on the log generalised normal shape with power 0.5. "
    f"The simulated variance stays within 10% of the formula from $n = 24$ in {by_n('RMSE', 24)[0]} of "
    f"{by_n('RMSE', 24)[1]} settings for RMSE, {by_n('MAE', 24)[0]} of {by_n('MAE', 24)[1]} for MAE, "
    f"{by_n('MAPE', 24)[0]} of {by_n('MAPE', 24)[1]} for MAPE and {by_n('sMAPE', 24)[0]} of {by_n('sMAPE', 24)[1]} "
    "for sMAPE. Where the density of $\\ln y$ has a corner at its median (log-Laplace actuals), MAE and sMAPE need "
    "250 to 1000 actuals. "
    "In Table 1, RMSE has the smallest $V$ at the narrow spreads and on log-uniform, gamma and Weibull actuals. MAE "
    "has it on log-Laplace actuals and on lognormal actuals with a wide spread. MAPE has it on inverse Gaussian "
    "actuals with shape 0.05. "
    "On the six shapes of Table 2, [[D41:Hypotailedness]] falls as the [[D37:peakedness index]] rises, and the index "
    "crosses 1 at the Laplace shape.")
for xid, cid in NEW_CLAIMS.items():
    x = copy.deepcopy(xcl[xid])
    c = {"id": cid, "name": NAMES[cid][0], "summary": NAMES[cid][1]}
    c["statement"] = c25_statement if cid == "C25" else remap(x["statement"])
    for h, cc in (("H49", "C23"), ("H48", "C22"), ("H50", "C24"), ("H47", "C6")):
        c["statement"] = re.sub(rf"\b{h}\b", cc, c["statement"])
    if "about" in x:
        a = x["about"]
        c["about"] = {"subject": dref(a["subject"]), "relation": a["relation"],
                      "object": [dref(o) for o in a["object"]] if isinstance(a["object"], list) else dref(a["object"]),
                      "comparator": a["comparator"], "conditions": a["conditions"]}
    if cid == "C25":
        c.update({"type": "descriptive", "role": "secondary", "novelty": "no_prior_found"})
    else:
        c.update({"type": "formal", "role": "principal",
                  "novelty": {"C6": "no_prior_found", "C22": "no_prior_found", "C23": "extension",
                              "C24": "no_prior_found"}[cid]})
    c["support"] = SUPPORT.get(cid) or [{"kind": "empirical_quant", "evidence": "E44",
                                         "results": [f"E44.T11.r{i}" for i in range(1, len(rows1) + 1)]
                                         + [f"E44.T12.r{i}" for i in range(1, len(rows2) + 1)]
                                         + ["E44.T2", "E44.T4"]}]
    for k in ("provenance", "hypothesis", "planning", "planned_at", "planned_at_anchor", "test", "estimate", "checks"):
        if k in x:
            c[k] = x[k]
    if "estimate" in c:
        c["estimate"] = walk(c["estimate"], lambda s: dref(s))
    if cid != "C25":
        c["verification"] = "agent_checked"
    d["claims"] = [y for y in d["claims"] if y["id"] != cid] + [c]
claims = {c["id"]: c for c in d["claims"]}
claims["C6"]["depends_on"] = ["C3", "C9", "C4", "C12"]
claims["C22"]["depends_on"] = ["C19", "C5"]
claims["C23"]["depends_on"] = ["C3", "C9", "C12", "C19"]
claims["C24"]["depends_on"] = ["C23", "C19", "C5"]
claims["C25"]["depends_on"] = ["C23"]
claims["C25"]["about"] = {"subject": "D36", "relation": "computed and simulated on 21 distributions for",
                          "object": ["D10", "D11", "D12", "D2"], "comparator": None, "conditions": []}
claims["C25"]["relations"] = [
    {"type": "supports", "target": "C23", "id_type": "local",
     "note": "C23 gives the asymptotic variance of each measure's minimiser in closed form.",
     "results": ["E44.T2", "E44.T4"]}]
claims["C23"]["scope"] = {"conditions": "positive actuals with a continuous distribution; independent actuals; the "
                                         "conditions each part states"}
claims["C24"]["scope"] = {"conditions": "log-location-scale families; limits as the spread grows"}
claims["C22"]["scope"] = {"conditions": "log-symmetric unimodal actuals; small and large spreads, and narrow-core "
                                         "mixtures"}
claims["C6"]["scope"] = {"conditions": "positive actuals whose logs fall in two separated intervals, any shapes "
                                        "and weights"}
claims["C25"]["scope"] = {"conditions": "21 distributions; test sets of 6 to 1000 actuals; 2000 test sets each"}
claims["C22"]["relations"] = [
    {"type": "extends", "reference": "B20", "self": True,
     "note": "C19 of entrylib:2026.000016 shows that the Extremity coefficient and kurtosis order distributions "
             "differently: the Laplace distribution has a smaller excess kurtosis than Student's t with 5 degrees of "
             "freedom but a smaller Extremity."},
    {"type": "refines", "target": "C19", "id_type": "local",
     "note": "C19 gives the asymptotic relative efficiency for log-symmetric unimodal actuals and its limits for "
             "lognormal and log-Laplace actuals.", "results": ["E43.T5"]}]
claims["C23"]["relations"] = [
    {"type": "extends", "reference": "B17",
     "note": "Gneiting (2011) shows that MAPE's target is the (-1)-median, the median of the distribution with density "
             "proportional to that of y divided by y."}]
claims["C6"]["relations"] = [
    {"type": "qualifies", "reference": "B4",
     "note": "C5 of entrylib:2026.000009 found sMAPE's preferred candidate close to the median on the distributions it "
             "studied."},
    {"type": "refines", "target": "C12", "id_type": "local",
     "note": "C12 says when the ordinary median is the target: the weights balance there and enough actuals lie near "
             "it.", "results": ["E43.T9"]}]
claims["C24"]["relations"] = [
    {"type": "extends", "reference": "B15",
     "note": "van der Vaart (1998) gives the asymptotic normality of estimators that minimise a criterion, with their "
             "limiting variance."},
    {"type": "explains", "target": "C10", "id_type": "local",
     "note": "C10 says expected sMAPE flattens as the actuals spread out, so a test set tells forecasts apart less "
             "well.", "results": ["E43.T5"]}]

# ------------------------------------------------------------------ limitations
lims = {l["id"]: l for l in d["limitations"]}
del lims["L5"]                                  # version 3's C6 proved mirror-image groups only; C6 is now general
for c in d["claims"]:
    if "L5" in c.get("limitations", []):
        c["limitations"] = [x for x in c["limitations"] if x != "L5"]
lims["L3"]["affects"] = sorted(set(lims["L3"]["affects"]) | {"C22", "C23", "C24", "C25"}, key=lambda s: int(s[1:]))
lims["L4"]["affects"] = sorted(set(lims["L4"]["affects"]) | {"C25"}, key=lambda s: int(s[1:]))
lims["L8"] = {"id": "L8", "category": "deviation", "affects": ["C23", "C25"], "impact": "reduces_transferability",
              "mitigation": "partial", "mitigation_evidence": "E42",
              "statement": "Lab run RN20 deviated from plan PL18 (deviation DV5). Three settings did not come within the "
                           "plan's tolerance of the formulas of C23 by n = 4000. For RMSE on lognormal actuals with "
                           "sigma 2, the simulated V is 35.2 at n = 1000 and 40.6 at 4000, against the formula's 53.6. "
                           "On the log generalised normal shape with power 0.5, it is 6.56 and 5.36 for MAE (formula "
                           "4), and 4.51 and 3.66 for sMAPE (formula 2.36). Each moves towards the formula as n grows, and the exact small-sample "
                           "values of MAE (6.35 at 1000, 5.49 at 4000) agree with the simulation. As the supervisor "
                           "decided after the run, the formulas are taken as checked where n = 4000 suffices, and "
                           "these three are reported as slow cases. For very heavy tails or a sharp peak, the "
                           "formulas need far larger test sets than are used in practice."}
lims["L9"] = {"id": "L9", "category": "model_specification", "affects": ["C24"], "impact": "could_bias_magnitude",
              "mitigation": "none",
              "statement": "The growth rates of C24 are limits as the spread grows. At moderate spreads they are "
                           "approached slowly: V of sMAPE over sigma cubed is 0.398 at spread 10 and 0.336 at 200 for "
                           "the Laplace shape, against its limit 1/3."}
lims["L10"] = {"id": "L10", "category": "model_specification", "affects": ["C22", "C23", "C24", "C25"],
               "impact": "reduces_transferability", "mitigation": "none",
               "statement": "C22 to C25 are about positive actuals with a continuous distribution. They do not cover "
                            "actuals with zeros or repeated values, such as intermittent demand, which C7, C11 and "
                            "C21 treat."}
lims["L11"] = {"id": "L11", "category": "other", "affects": ["C22"], "impact": "could_bias_magnitude",
               "mitigation": "none",
               "statement": "The peakedness index decides the sign of ARE - 1 only for small spreads (C22 a). At "
                            "moderate spreads the ARE must be computed. Hypotailedness orders the six shapes of "
                            "Table 2 as the index does, but it is not a criterion: across narrow-core mixtures the two "
                            "do not move together."}
d["limitations"] = [lims[k] for k in sorted(lims, key=lambda s: int(s[1:]))]
for lid, l in lims.items():
    for cid in l.get("affects", []):
        c = claims.get(cid)
        if c is not None and lid not in c.get("limitations", []):
            c["limitations"] = sorted(c.get("limitations", []) + [lid], key=lambda s: int(s[1:]))

e1 = next(e for e in d["evidence"] if e["id"] == "E1")
assert "H7 checks C6" in e1["description"]
e1["description"] = e1["description"].replace("H7 checks C6", "H7 checked version 3's C6 (two mirror-image groups)")

# ------------------------------------------------------------------ paper
paper = d["paper"]
expl = copy.deepcopy(X["paper"]["exploration"])
expl["lab_record"] = "A7"
SHORT = {
    ("H11", "PL5"): "A step of the flatness argument, as first worded.",
    ("H16", "PL5"): "An earlier statement of the weighted-median target.",
    ("H9", "PL5"): "An earlier statement of sMAPE near the median.",
    ("H10", "PL5"): "An earlier statement of the asymmetry of sMAPE.",
    ("H27", "PL9"): "Version 3's C6: with two mirror-image groups, the median is a local maximum.",
    ("H18", "PL5"): "An earlier wording of version 3's C6 on two groups of actuals.",
    ("H47", "PL18"): "C6 as first worded: its part (c) said no median is the target, which (d) contradicts for two "
                     "points with p = 1/2.",
    ("H51", "PL20"): "C6 with 'mirror' not yet tied to p = 1/2, as the checker noted.",
    ("H48", "PL18"): "C22 with 'kurtosis 6' not yet said to be excess kurtosis 3; one exact-value check failed through "
                     "a code defect (DV4).",
    ("H50", "PL18"): "C24 as first worded: parts (a) and (c) lacked continuity of h at the median and at z*.",
    ("H39", "PL12"): "C18 as first worded: a condition of part (c) was missing.",
    ("H39", "PL13"): "C18 as first worded, decided again on new runs: the missing condition stood.",
    ("H43", "PL15"): "C18 with its closing sentence too broad.",
    ("H45", "PL16"): "C18 with 'of order' for 'at most of order'.",
    ("H40", "PL12"): "C19 calling the ratio of limiting variances the large-sample relative efficiency.",
}
for n in expl["not_claimed"]:
    if isinstance(n.get("into"), str):
        n["into"] = CMAP.get(n["into"], n["into"])
    k = (n["lab_event"], n["lab_plan"])
    if k in SHORT:
        n["statement"] = SHORT[k]
    else:
        into = n.get("into")
        name = claims[into]["name"] if into in claims else ""
        n["statement"] = f"An earlier wording of {into}: {name[:1].lower() + name[1:]}."
    assert len(n["statement"]) <= 200, n["statement"]
paper["exploration"] = expl
paper["contributions"] = paper["contributions"] + ["C22", "C23", "C24"]
for s in paper["sections"]:
    if s["title"] == "Precision on a test set":
        s["claims"] = ["C23", "C24", "C19", "C22", "C15", "C25"]
paper["summary"] = (
    "Which forecast does sMAPE reward, and how precisely does a test set find it? For positive values sMAPE depends "
    "only on the ratio of actual to forecast: it is the hyperbolic tangent of half the absolute log error. So it scores "
    "a forecast k times too high and one k times too low the same. Being bounded, it has a finite expected value for "
    "every distribution. Its target is a median of the actuals, each weighted by how close it is to the forecast by "
    "ratio. This is the median when the log of the actual is symmetric and unimodal. With two separated groups of "
    "actuals, the target lies in one of them, and the median can be a local maximum. On intermittent demand the target "
    "is 0 once the share of zeros passes a threshold below one half. As the actuals spread out, expected sMAPE "
    "approaches its bound for every forecast, and forecasts become hard to tell apart. For any continuous distribution "
    "of positive actuals, we give the large-sample variance of the test-set minimiser of RMSE, MAE, MAPE and sMAPE. It "
    "grows with the spread like its square for MAE, at a rate set by the tails for RMSE and MAPE, and by the peak for "
    "sMAPE. For small spreads, a peakedness index, not the kurtosis, decides whether the sMAPE minimiser or the sample "
    "median is the more precise.")
paper["highlights"] = ["C1", "C4", "C23", "C22", "C6"]
paper["classification"]["keywords"] = paper["classification"]["keywords"] + ["asymptotic variance",
                                                                             "relative efficiency", "Hypotailedness"]
paper["correction_reason"] = (
    "Version 4 answers a reader's notes on version 3. Corrected: C6 now covers any two separated groups of actuals, "
    "with any shapes and weights; version 3 proved it for mirror-image groups only. Added: how precisely a test set "
    "locates the target of RMSE, MAE, MAPE and sMAPE (C23); how this precision falls as the actuals spread out (C24); "
    "when the sMAPE minimiser is more precise than the sample median (C22); and these computed and simulated on 21 "
    "distributions (C25). Also added: limitations L8 to L11, and the proofs of the new claims.")

# ------------------------------------------------------------------ contributors: the record's run periods and tasks
xc = {c["id"]: c for c in X["contributors"]}
for c in d["contributors"]:
    for k in ("run_period", "models", "ai_tasks"):
        if k in xc[c["id"]]:
            c[k] = xc[c["id"]][k]

d["artefacts"] = [arts[k] for k in sorted(arts, key=lambda s: int(s[1:]))]

# long formulas on a line of their own (warning 44), in every text the validator reads for formulas
sys.path.insert(0, str(REPO.parent / "claimset" / "src"))
from claimset import validate as _v  # noqa: E402


def display_long(text):
    def one(m):
        tex = (m.group(1) or m.group(2) or "").strip()
        if not tex or _v.formula_symbols(tex) <= _v.MAX_INLINE_SYMBOLS:
            return m.group(0)
        return "\n$$" + tex + "$$\n"
    out = _v.INLINE_FORMULA.sub(one, text)
    out = re.sub(r"\$\$\n([.,;:])", r"\1$$\n", out)          # a sentence's punctuation stays with its formula
    out = re.sub(r"\$\$\n\s*([.,;:])", lambda m: m.group(1) + "$$\n", out)
    out = re.sub(r"\$\$([.,;:])\$\$", r"\1", out)
    out = re.sub(r"\$\$\n\s+", "$$\n", out)
    out = re.sub(r"\s+\n\$\$", "\n$$", out)
    return out.strip()


for lst, fields in _v.FORMULA_FIELDS.items():
    for x in d.get(lst) or []:
        for f in fields:
            if isinstance(x.get(f), str):
                x[f] = display_long(x[f])
        if lst == "evidence":
            for st in x.get("steps") or []:
                st["statement"] = display_long(st["statement"])
d["claims"].sort(key=lambda c: int(c["id"][1:]))
d["evidence"].sort(key=lambda e: int(e["id"][1:]))

# definitions in reading order: version 3's order, the new ones after the terms they use
order = [x["id"] for x in P["definitions"]]
pos = order.index("D9") + 1
order[pos:pos] = ["D36", "D39", "D40", "D37", "D41", "D38"]
dd = {x["id"]: x for x in d["definitions"]}
d["definitions"] = [dd[i] for i in order]
OUT.write_text(json.dumps(d, ensure_ascii=False, indent=1) + "\n")
print("wrote", OUT)
