"""Complete the exported draft of smape-target into version 3 of the claimset (claimset v0.44).

    python3 build.py SRC OUT PROOFS

SRC is the draft `claimset lab export` wrote (its results CSV beside it, in files/), OUT the claimset to write (its
files go to files/ beside it), PROOFS the study's proofs.md (proof-suilin.md and the other files beside it, the lab's
runs in ../lab/runs)."""
import hashlib
import json
import re
import shutil
import sys
from pathlib import Path

SRC, OUT, PROOFS = sys.argv[1], sys.argv[2], sys.argv[3]
d = json.loads(Path(SRC).read_text())
FILES = Path(PROOFS).parent
MASTER = "files/smape-target.results.csv"

# ------------------------------------------------------------------ the package's files, beside the claimset
out_files = Path(OUT).parent / "files"
(out_files / "figures").mkdir(parents=True, exist_ok=True)
for name in ("proofs.md", "proofs.pdf", "proof-suilin.md", "smape_weights.png", "proof-check-report.md"):
    shutil.copyfile(FILES / name, out_files / name)
shutil.copyfile(Path(SRC).parent / "files" / (Path(SRC).name.split(".")[0] + ".results.csv"), Path(OUT).parent / MASTER)
LOCAL = {"A1": "files/proofs.md", "A2": "files/proofs.pdf", "A3": "files/smape_weights.png",
         "A4": "files/proof-check-report.md", "A5": "files/proof-suilin.md", "A6": MASTER}
for a in d["artefacts"]:
    if a["id"] in LOCAL:
        assert a["uri"].endswith(LOCAL[a["id"]].split("/")[-1]) or a["id"] == "A6", a["uri"]
        a["uri"] = LOCAL[a["id"]]
        sha = "sha256:" + hashlib.sha256((Path(OUT).parent / a["uri"]).read_bytes()).hexdigest()
        assert a.get("version", sha) == sha, (a["id"], "the file is not the one the record names")
        a["version"] = sha


# ------------------------------------------------------------------ the IDs readers have seen
# Published IDs never change meaning, and a dropped ID is not reused (claimset §9). The record numbers its claims,
# evidence and definitions afresh at each export; here they get the published numbers, and what is new the next free
# ones.
def _walk(x, f):
    if isinstance(x, dict):
        return {f(k) if isinstance(k, str) else k: _walk(v, f) for k, v in x.items()}
    if isinstance(x, list):
        return [_walk(v, f) for v in x]
    return f(x) if isinstance(x, str) else x


# definitions: D31 is the published "formula checks" (added below); the export's D31 to D34 are new terms
DMAP = {"D31": "D32", "D32": "D33", "D33": "D34", "D34": "D35"}
assert [x["term"] for x in d["definitions"] if x["id"] in DMAP] == [
    "signal-to-noise ratio of a comparison", "msMAPE term", "msMAPE target", "intermittent demand settings"]
assert {x["id"]: x["term"] for x in d["definitions"] if x["id"] in ("D28", "D29")} == {"D28": "Suilin term",
                                                                                      "D29": "Suilin target"}


def _dmap(x):
    if x in DMAP:
        return DMAP[x]
    return re.sub(r"\[\[(D\d+):", lambda m: f"[[{DMAP.get(m.group(1), m.group(1))}:", x)


d = _walk(d, _dmap)

# claims, by the hypothesis or plan they report
_lab = {c.get("hypothesis", {}).get("lab_event") or c.get("estimate", {}).get("lab_plan"): c["id"] for c in d["claims"]}
PUBLISHED = {"H21": "C1", "H22": "C2", "H36": "C3", "H37": "C4", "H26": "C5", "H27": "C6", "H28": "C7", "H29": "C8",
             "H24": "C9", "H30": "C10", "H31": "C11", "H38": "C12", "PL1": "C13", "PL4": "C14", "PL2": "C15",
             "PL11": "C16", "H41": "C17", "H46": "C18", "H44": "C19", "H42": "C20", "PL14": "C21"}
assert set(_lab) == set(PUBLISHED), sorted(set(_lab) ^ set(PUBLISHED))
cmap = {_lab[k]: v for k, v in PUBLISHED.items()}

# evidence, by its lab run: E1 the checks of RN8 (one object, below), E11 of RN6, E13 to E16 the four computed
# studies; new: E31 the checks of RN13 (one object), E32 of RN18, E33 of RN16, E34 the study of RN15. E17 (the checks
# of RN9) is gone: the record now decides C17 by the checks of RN13. E18 to E30 and E35 to E37 are the proofs.
_runs = {}
for e in d["evidence"]:
    _runs.setdefault(tuple(e.get("lab_runs") or []), []).append(e["id"])
EV_PUBLISHED = {("RN8",): "E1", ("RN6",): "E11", ("RN1",): "E13", ("RN2",): "E14", ("RN7",): "E15", ("RN10",): "E16",
                ("RN13",): "E31", ("RN18",): "E32", ("RN16",): "E33", ("RN15",): "E34"}
assert set(_runs) == set(EV_PUBLISHED), sorted(_runs)
emap, _spare = {}, 9000
for run, ids in _runs.items():
    emap[ids[0]] = EV_PUBLISHED[run]
    for other in ids[1:]:                       # the same run exported once per claim: merged into the first below
        _spare += 1
        emap[other] = f"E{_spare}"
_idpat = re.compile(r"^([CE]\d+)(\..*)?$")


def _cemap(x):
    m = _idpat.match(x)
    if m and (m.group(1) in cmap or m.group(1) in emap):
        return (cmap.get(m.group(1)) or emap.get(m.group(1))) + (m.group(2) or "")
    return x


d = _walk(d, _cemap)
d["claims"].sort(key=lambda c: int(c["id"][1:]))
d["evidence"].sort(key=lambda e: int(e["id"][1:]))


def _code_artefact():
    """A54, as published, at the commit the last runs ran at."""
    assert len(d["artefacts"]) == 53, len(d["artefacts"])
    d["artefacts"].append({"id": "A54", "kind": "code",
                           "uri": f"https://github.com/entrylib/study-planning/tree/{CODE}/claimset-code/smape-target",
                           "version": CODE, "licence": "CC-BY-4.0", "availability": {"level": "public"},
                           "description": "The study's code: the distributions, the exact quadrature of expected "
                                          "sMAPE, the formula checks, the three simulations (targets.py, "
                                          "efficiency.py, flatness.py), and the targets of the Suilin loss "
                                          "(msmape_targets.py) and of msMAPE on intermittent demand "
                                          "(msmape_intermittent.py)."})


# the code at the commit of the record's last run (RN18), which holds every script of the earlier runs
CODE = next(e["execution"]["code_version"] for e in d["evidence"] if e.get("lab_runs") == ["RN18"])
import condense
condense.run(d, FILES.parent / "lab" / "runs", Path(OUT).parent, MASTER, before_new=_code_artefact)
claims = {c["id"]: c for c in d["claims"]}
ev = {e["id"]: e for e in d["evidence"]}

# paper section of the proofs -> claim of this file
PAPER = {"H1": "C1", "H2": "C2", "H3": "C3", "H4": "C9", "H5": "C4", "H6": "C5", "H7": "C6", "H8": "C7",
         "H9": "C11", "H10": "C8", "H11": "C10", "H12": "C12", "H13": "C18", "H14": "C19"}
# the proofs' evidence objects keep their published numbers; the three new proofs take the next free ones
PROOF_ID = {"C1": "E18", "C2": "E19", "C3": "E20", "C9": "E21", "C4": "E22", "C5": "E23", "C6": "E24", "C7": "E25",
            "C11": "E26", "C8": "E27", "C10": "E28", "C12": "E29", "C17": "E30", "C18": "E35", "C19": "E36",
            "C20": "E37"}


# ------------------------------------------------------------------ proofs -> formal evidence
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


MATH = re.compile(r"(\$\$.+?\$\$|\$[^$]+\$)", re.S)


def tidy(s):
    s = s.split("*Revision.*")[0]
    s = re.sub(r"\s+", " ", s).strip()
    s = re.sub(r"^\*(.+?\.)\*\s", r"\1 ", s)
    # Markdown emphasis, outside mathematics only: in $\bar s(f^*) \le \bar s(f^*)$ the two asterisks are not a pair
    parts = MATH.split(s)
    for i in range(0, len(parts), 2):
        parts[i] = re.sub(r"\*\*([^*]+)\*\*", r"\1", parts[i])
        parts[i] = re.sub(r"\*([^*]+)\*", r"\1", parts[i])
        assert "*" not in parts[i], parts[i]
    s = "".join(parts)
    s = re.sub(r"\s*\$\$\s*(.+?)\s*\$\$\s*", lambda m: "\n$$" + m.group(1) + "$$\n", s).strip()
    s = re.sub(r"\bH(\d+)\b", lambda m: PAPER.get(f"H{m.group(1)}", m.group(0)), s)
    s = s.replace("This is the asymmetry Goodwin and Lawton (1999) describe. ", "")
    s = s.replace(" This is the asymmetry Goodwin and Lawton (1999) describe.", "")
    s = s.replace(" This is an elementary case of Anderson's theorem on symmetric unimodal densities (Anderson, 1955).",
                  "")
    return s


RULES = [("quotient rule", "quotient rule"), ("dominated convergence", "dominated convergence"),
         ("bounded convergence", "bounded convergence"), ("Taylor", "Taylor's theorem"),
         ("Fermi", "standard integral"), ("By symmetry", "symmetry"), ("by symmetry", "symmetry"),
         ("mirror image", "symmetry"), ("Expanding", "algebra"), ("Multiplying", "algebra"),
         ("Dividing", "algebra"), ("Collecting", "algebra"), ("Substituting", "change of variables"),
         ("continuous", "continuity"), ("increasing", "monotonicity"), ("convex", "convexity"),
         ("Write ", "definition"), ("Take ", "case analysis")]


def rule_of(text, refs):
    for key, rule in RULES:
        if key in text:
            return rule
    if refs:
        return "earlier result"
    return "algebra" if "$" in text else "logical argument"


def steps_of(body):
    """[[id, text]] of a proof: an introductory paragraph (S0), numbered steps (S1, Sa1 for part (a)), and the
    closing paragraph (SQ)."""
    lines = body.split("**Proof.**", 1)[1].splitlines()
    intro, steps, closing, part, cur = [], [], [], "", None
    for line in lines:
        st = line.strip()
        if not st:
            continue
        if closing or (line.startswith("So ") and not line.startswith(" ")):
            closing.append(st)
            continue
        if re.fullmatch(r"\(([a-c])\)", st):
            part = st[1]
            continue
        m = re.match(r"^(\d+)\.\s+(.*)", line)
        if m:
            cur = [f"S{part}{m.group(1)}", m.group(2)]
            steps.append(cur)
        elif cur is not None and line.startswith(" "):
            cur[1] += " " + st
        elif cur is None:
            intro.append(st)
        else:
            cur[1] += " " + st
    out = ([["S0", " ".join(intro)]] if intro else []) + steps
    if closing:
        out.append(["SQ", " ".join(closing).replace("∎", "").strip()])
    return out


# "step 3", "steps 2 and 3", "steps 2 to 4": every number of the list is a step
STEPS = r"steps? (\d+(?:(?:, and |, | and | to )\d+)*)"


def step_numbers(listed):
    """The step numbers a list names: "2 to 4" is 2, 3 and 4."""
    out = []
    for a, to, b in re.findall(r"(\d+)(?:( to )(\d+))?", listed):
        out += list(range(int(a), int(b) + 1)) if to else [int(a)]
    return out


def formal_evidence(eid, cid, paper_h, body):
    assert body.count("**Proof.**") == 1 and "*What it shows:*" in body, (cid, paper_h)
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
        prefix = sid[1] if len(sid) > 1 and sid[1] in "abc" else ""
        cited = [f"S{prefix}{x}" for listed in re.findall(STEPS, local) for x in step_numbers(listed)
                 if f"S{prefix}{x}" in ids]
        by = list(dict.fromkeys(cited + refs))
        if sid == "SQ":
            rule, by = "conclusion", [i for i in ids if i != "S0"][-1:] or ids[-1:]
        elif sid == "S0":
            rule = "setup"
        else:
            rule = rule_of(text, refs)
            prev = [i for i in ids if i != "S0" and (i[1:2] in "abc") == bool(prefix) and i[1:2] == (prefix or i[1:2])]
            if not by:
                by = prev[-1:] or [i for i in ids if i == "S0"]
        step = {"id": sid, "statement": t, "rule": rule}
        if by:
            step["by"] = by
        steps.append(step)
        ids.append(sid)
    # renumber S0/Sa1/SQ into S1..Sn, keeping references
    new = {s["id"]: f"S{i + 1}" for i, s in enumerate(steps)}
    for s in steps:
        s["id"] = new[s["id"]]
        if "by" in s:
            s["by"] = [new.get(b, b) for b in s["by"]]
    return {"id": eid, "profile": "formal",
            "description": f"Proof of {cid} (section {paper_h} of the proof document A1, rendered in A2). "
                           f"{tidy(shows)[:1].upper() + tidy(shows)[1:]}",
            "foundation": "Kolmogorov probability axioms and real analysis",
            "steps": steps, "artefacts": ["A1", "A2"] + (["A3"] if cid == "C4" else []),
            "availability": {"level": "public"}}


secs = sections(Path(PROOFS).read_text())
assert list(secs) == [f"H{i}" for i in range(1, 15)], list(secs)
formal_ids = {}
for paper_h, cid in PAPER.items():
    e = formal_evidence(PROOF_ID[cid], cid, paper_h, secs[paper_h])
    d["evidence"].append(e)
    formal_ids[cid] = e["id"]

# C17 and C20: the two proofs of proof-suilin.md (A5), each under a top-level heading
suilin = re.split(r"(?m)^(?=# )", (FILES / "proof-suilin.md").read_text())
suilin = [x for x in suilin if x.strip()]
assert len(suilin) == 2 and suilin[0].startswith("# The Suilin loss beside sMAPE") \
    and suilin[1].startswith("# The msMAPE of Smyl beside sMAPE"), [x[:50] for x in suilin]
for cid, body, which in (("C17", suilin[0], "first"), ("C20", suilin[1], "second")):
    assert not re.search(r"\bH\d+\b", body.split("**Proof.**")[1]), cid     # no section numbers of proofs.md to map
    e = formal_evidence(PROOF_ID[cid], cid, "", body)
    shows = tidy(re.search(r"\*What it shows:\*\s*(.*?)\n\n", body, re.S).group(1))
    e["description"] = f"Proof of {cid} (the {which} proof of the proof document A5). {shows[:1].upper() + shows[1:]}"
    e["artefacts"] = ["A5"]
    d["evidence"].append(e)
    formal_ids[cid] = e["id"]
assert formal_ids == PROOF_ID

# ------------------------------------------------------------------ definitions
defs = {x["id"]: x for x in d["definitions"]}
defs["D5"]["definition"] = re.sub(r"\bH(\d+)\b", lambda m: PAPER.get(f"H{m.group(1)}", m.group(0)),
                                  defs["D5"]["definition"])
d["definitions"].append({"id": "D31", "term": "formula checks", "kind": "method",
                         "definition": "Numerical checks of the formulas of the formal claims "
                                       "(claimset-code/smape-target/check_formulas.py, seed 6; check_jump.py for C11): "
                                       "identities and inequalities on 200000 random pairs, derivatives against finite "
                                       "differences, minimisers, balances and thresholds against quadrature and "
                                       "bisection on every setting of the selection design, and the large-spread "
                                       "expansion against its bound. For the claims new or reworded in version 3 "
                                       "(C3, C4, C12, C17 to C20): check_v3.py (seed 8), check_h39.py and "
                                       "check_corrected.py (seed 8), in the same way."})

# ------------------------------------------------------------------ claims
CLAIM_ROLE = {"C1": "principal", "C2": "supporting", "C3": "principal", "C4": "principal", "C5": "supporting",
              "C6": "principal", "C7": "principal", "C8": "supporting", "C9": "principal", "C10": "principal",
              "C11": "supporting", "C12": "principal", "C13": "supporting", "C14": "supporting", "C15": "principal",
              "C16": "secondary", "C17": "principal", "C18": "principal", "C19": "principal", "C20": "principal",
              "C21": "secondary"}
LENS = {"C17": "D28", "C18": "D32", "C19": "D8", "C20": "D33", "C1": "D1", "C2": "D1", "C3": "D1", "C4": "D6", "C5": "D6", "C6": "D3", "C7": "D6", "C8": "D3",
        "C9": "D3", "C10": "D3", "C11": "D6", "C12": "D6"}
for cid, c in claims.items():
    c["role"] = CLAIM_ROLE[cid]
    if cid in formal_ids:
        c["type"] = "formal"
        c["verification"] = "agent_checked"
        quant = [s for s in c["support"] if s["kind"] == "empirical_quant"]
        c["support"] = [{"kind": "formal_proof", "evidence": formal_ids[cid]}] + quant
        for s in quant:
            e = ev[s["evidence"]]
            # the plans of version 3 record their setup; the earlier check runs get the claim's term
            e["setup"] = ({**e["setup"], "procedure": ["D31"]} if e.get("setup") and "procedure" not in e["setup"]
                          else e.get("setup") if e.get("setup") else {"studied": [LENS[cid]], "procedure": ["D31"]})
            run = ", ".join(e.get("lab_runs") or [])
            e["description"] = (f"Numerical checks of the formulas of {cid} ([[D31:formula checks]]), lab run {run} "
                                f"under the formal plan that decided {cid}; every check is within the plan's "
                                f"tolerances. They check the formulas on concrete numbers; the proof is the evidence "
                                f"for the claim.")

claims["C9"]["statement"] = claims["C9"]["statement"].replace("(unlike RMSE, MAE and MAPE, H12)",
                                                              "(unlike RMSE, MAE and MAPE, C8)")
claims["C5"]["statement"] = claims["C5"]["statement"].replace("[[D6:minimiser of expected sMAPE]]",
                                                              "[[D6:sMAPE target]]")

claims["C13"]["statement"] = (
    "Computed from the distributions themselves, the [[D6:sMAPE target]] on the 34 outcome distributions of plan PL1 "
    "is the median on lognormal outcomes, lies below it when $\\ln Y$ is skewed right and above it when skewed left "
    "(quantile levels from 0.439 to 0.556 on the unimodal positive distributions), is two forecasts, one near each "
    "mode, on equal bimodal outcomes and the larger mode's on unequal ones, lies away from zero on the side of the "
    "mean on normal outcomes straddling zero, and switches to 0 at shares of zero outcomes from 0.343 to 0.468, all "
    "below one half; the one-step approximation of the target's quantile level is within 0.003 of the exact level on "
    "every distribution with one local minimum.")
claims["C13"]["about"] = {"subject": "D6", "relation": "was derived for", "object": "34 outcome distributions",
                          "comparator": None, "conditions": []}
claims["C13"]["scope"] = {"data_type": "outcome distributions given exactly (no data)",
                          "conditions": ["34 distributions: inverse Gaussian, log-skewed, lognormal, bimodal, normal "
                                         "straddling zero, zero-inflated lognormal"]}
claims["C14"]["statement"] = (
    "As the spread $\\sigma$ of $\\ln Y$ grows from 0.1 to 10, the [[D3:expected sMAPE]] at the [[D6:sMAPE target]] "
    "rises towards its bound 200 (from 7.97 to 178.2 on lognormal outcomes), the target's quantile level moves from "
    "the median to that of the mode of $\\ln Y$ (0.372 against the mode's 0.368 at $\\sigma = 10$ on log-Gumbel "
    "outcomes skewed right), and sMAPE separates its target from a forecast half a $\\sigma$ away less reliably than "
    "MAE: at $\\sigma = 4$ with test sets of 12 outcomes it ranks its target first in 0.61 of the test sets, MAE in "
    "0.82.")
claims["C14"]["about"] = {"subject": "D3", "relation": "flattens as the spread grows, by", "object": "D20",
                          "comparator": "D11", "conditions": []}
claims["C14"]["scope"] = {"data_type": "outcome distributions given exactly, and test sets drawn from them",
                          "conditions": ["lognormal and log-Gumbel outcomes, spread of ln y from 0.1 to 10",
                                         "test sets of 12 to 1000 outcomes"]}
claims["C15"]["statement"] = (
    "Where sMAPE and MAE both target the median (log-symmetric actuals), the [[D8:sMAPE minimiser of a test set]] "
    "varies more than the sample median for every spread and test-set size studied: its [[D9:relative efficiency]] "
    "is below 1 and falls as the actuals spread out, at $n = 10000$ from 0.97 to 0.03 on lognormal actuals as "
    "$\\sigma$ grows from 0.25 to 6, and from 1.00 to 0.34 on log-Laplace actuals as $b$ grows from 0.25 to 6, close "
    "to the large-sample values.")
claims["C15"]["about"] = {"subject": "D8", "relation": "is less efficient than", "object": "D9",
                          "comparator": "the sample median", "conditions": ["D26"]}
claims["C15"]["scope"] = {"data_type": "test sets of independent draws from known distributions",
                          "conditions": ["lognormal and log-Laplace outcomes, spread 0.25 to 6",
                                         "test sets of 6 to 10000 outcomes"]}

# intervals in PL2's results: only level, lower, upper (claimset rule 82)
for e in d["evidence"]:
    hit = False
    for r in e.get("results") or []:
        if isinstance(r.get("interval"), dict) and "method" in r["interval"]:
            r["interval"].pop("method")
            hit = True
    for t in e.get("tables") or []:
        for row in t.get("rows") or []:
            for cell in (row.get("cells") or {}).values() if isinstance(row.get("cells"), dict) else []:
                if isinstance(cell, dict) and isinstance(cell.get("interval"), dict):
                    cell["interval"].pop("method", None)
    if hit:
        e["description"] = (e.get("description") or "") + (" The intervals of the relative efficiency are 99% "
                                                            "bootstrap percentile intervals from 2000 resamples.")

# ------------------------------------------------------------------ relations
REL = {
    "C2": [{"type": "qualifies", "target": "entrylib:2026.000008@v2#C33", "id_type": "claimset",
            "note": "C33 says sMAPE penalises a forecast below the outcome more than one above it. That holds for two "
                    "forecasts the same absolute distance away (C2 b); a forecast k times too low and one k times too "
                    "high score the same (C2 a)."}],
    "C12": [{"type": "explains", "target": "entrylib:2026.000009@v1#C5", "id_type": "claimset",
             "note": "C5 found by simulation that sMAPE prefers the median on lognormal and log-Laplace outcomes; C12 "
                     "proves it for every symmetric unimodal ln Y."}],
    "C4": [{"type": "explains", "target": "entrylib:2026.000009@v1#C5", "id_type": "claimset",
            "note": "C5 found sMAPE's preferred candidate close to the median on other distributions; C4 shows why: "
                    "the target is a weighted median, on the side of the median that carries more weight (C12)."}],
    "C7": [{"type": "refines", "target": "entrylib:2026.000009@v1#C6", "id_type": "claimset",
            "note": "C6 found sMAPE preferring zero on counts with more than half zeros; C7 gives the exact threshold, "
                    "which is below one half whenever the positive outcomes vary."}],
    "C15": [{"type": "supports", "target": "entrylib:2026.000009@v1#C11", "id_type": "claimset",
             "note": "C11 found sMAPE's choices in short test sets spreading more than MAE's; C15 measures the "
                     "loss of precision for test sets of 6 to 10000 outcomes."},
            {"type": "supports", "target": "entrylib:2026.000008@v2#C6", "id_type": "claimset",
             "note": "C6 recommends against sMAPE for point forecasts; C15 shows that even where it targets the median "
                     "it estimates it less precisely than the sample median."}],
    "C10": [{"type": "supports", "target": "entrylib:2026.000008@v2#C6", "id_type": "claimset",
             "note": "C6 recommends against sMAPE for point forecasts; C10 shows that on widely spread actuals every "
                     "forecast scores close to the worst value and the target moves to the mode of ln Y."}],
}
# C18 stands to the earlier work as C10 does, which it refines
REL["C18"] = [
    {"type": "supports", "target": "entrylib:2026.000008@v2#C6", "id_type": "claimset",
     "note": "C6 recommends against sMAPE for point forecasts; C18 shows that on widely spread actuals a test set "
             "tells two forecasts apart less and less well under sMAPE."}]
for cid, rels in REL.items():
    claims[cid]["relations"] = rels

# ------------------------------------------------------------------ limitations
nl = len(d["limitations"])
d["limitations"] += [
    {"id": f"L{nl + 1}", "category": "sampling", "affects": ["C14", "C15"], "impact": "reduces_transferability",
     "mitigation": "none",
     "statement": "The test sets are independent draws from a fixed distribution. Real forecast errors are often "
                  "correlated and their distribution can change over a test period, so a real test set carries less "
                  "information than the same number of outcomes here."},
    {"id": f"L{nl + 2}", "category": "model_specification", "affects": ["C13", "C14", "C15"],
     "impact": "reduces_transferability", "mitigation": "none",
     "statement": "Only the distribution families studied were computed; the formal claims (C4, C10) say what holds "
                  "beyond them, the simulations only illustrate it."}]
assert nl == 2 and [(l["id"], l["affects"]) for l in d["limitations"][:2]] == [("L1", ["C13"]), ("L2", ["C15"])]

# ------------------------------------------------------------------ artefacts
arts = {a["id"]: a for a in d["artefacts"]}
for a in d["artefacts"]:
    if a["kind"] != "lab_record":
        a["licence"] = "CC-BY-4.0"
arts["A3"]["function"] = "diagram"
arts["A4"]["description"] = ("The independent checker's reports on the proofs (Claude Fable 5.2, lab actor AC4): "
                             "derivations from the statements alone, then comparisons with the proofs, with its "
                             "counterexamples and the corrections they led to.")
# the record's descriptions of the proof documents, with the claims a reader knows them by
assert arts["A1"]["description"].startswith("Proofs of the formal hypotheses, version 3,")
arts["A1"]["description"] = arts["A1"]["description"].replace("Proofs of the formal hypotheses", "Proofs of the formal claims") + (
    " Sections H1 to H14 prove C1, C2, C3, C9, C4, C5, C6, C7, C11, C8, C10, C12, C18 and C19, in that order. "
    "Sections H1 to H11 were also checked independently by lab actor AC4 (proof-check-report.md) and by the "
    "supervisor.")
assert arts["A5"]["description"].startswith("Proofs of H41 (the Suilin loss beside sMAPE;") \
    and ") and of H42 (the msMAPE of Smyl" in arts["A5"]["description"] \
    and arts["A5"]["description"].endswith("The msMAPE term is written s_nu here; the recorded statement writes p.")
arts["A5"]["description"] = arts["A5"]["description"].replace(
    "Proofs of H41 (the Suilin loss", "Proofs of C17 (record hypothesis H41: the Suilin loss").replace(
    ") and of H42 (the msMAPE of Smyl", ") and of C20 (record hypothesis H42: the msMAPE of Smyl").replace(
    "written s_nu here; the recorded statement writes p.",
    "written s_nu here, as in this claimset; the recorded statement writes p. In version 2 the file was "
    "proof-msmape.md.")

# ------------------------------------------------------------------ paper and contributors
p = d["paper"]
p["title"] = "What does sMAPE reward? An analytical and simulation study"
p["licence"] = "CC-BY-4.0"
p["autonomy"] = "human_supervised"
p["highlights"] = ["C1", "C12", "C4", "C15", "C6"]
p["contributions"] = ["C1", "C12", "C4", "C15", "C6", "C3", "C9", "C10", "C7"]
p["summary"] = (
    "Which forecast does sMAPE reward? The paper proves that sMAPE depends only on the ratio of actual to forecast and "
    "is a hyperbolic tangent of the log error, so its expected value is always finite. Its target is a weighted median "
    "of the actuals, weighting each by how close it is to the forecast by ratio; it is the median when the log of the "
    "actual is symmetric and unimodal, and lies on the side of the median that carries more weight otherwise. With "
    "two separate groups of actuals the median is a local maximum, with enough zeros the target is 0, and across zero "
    "it jumps. As the actuals spread out, expected sMAPE flattens towards its bound and the target moves to the mode "
    "of the log. Simulations show the target on 34 distributions and, where sMAPE and MAE both target the median, "
    "that a test set's sMAPE minimiser is a much noisier estimate of it than the sample median.")
# The guarantor chose the subfield (Statistics and Probability) on 2026-11-17: the classifier had assigned another.
p["classification"] = {"subfield": "openalex:subfields/2613",
                       "keywords": ["forecast evaluation", "error measures", "sMAPE", "point forecasts",
                                    "consistent scoring functions", "proof", "simulation study", "msMAPE",
                                    "Suilin loss", "intermittent demand"]}
p["declarations"] = {"funding": [], "competing_interests": (
    "Ivan Svetunkov runs EntryLib and wrote the review (entrylib:2026.000008) whose claims about sMAPE this study "
    "tests.")}
for c in d["contributors"]:
    if c["id"] == "P1":
        c["tools"] = ["code_execution", "file_access"]
        c["credit_roles"] = ["methodology", "formal_analysis", "software", "investigation", "visualization",
                             "writing_original_draft"]
        c["ai_tasks"] = [{"task": t, "extent": "full"} for t in ("analysis", "code", "writing", "json_preparation")]
    if c["id"] == "P2":
        c["positions"] = ["supervisor", "guarantor"]
        c["credit_roles"] = ["conceptualization", "supervision", "validation", "writing_review_editing"]
    if c["id"] == "P4":                              # the agent of version 3, which also completed this file
        assert c["label"].startswith("version 3") and c["positions"] == ["author"]
        c["tools"] = ["code_execution", "file_access"]
        c["credit_roles"] = ["methodology", "formal_analysis", "software", "investigation", "visualization",
                             "writing_original_draft"]
        c["ai_tasks"] = [{"task": t, "extent": "full"} for t in ("analysis", "code", "writing", "json_preparation")]


# ------------------------------------------------------------------ second pass
# one evidence object for the formula checks of a run, cited by every formal claim it checks
def one_object(run):
    """The first of the evidence objects export wrote for a run (one per claim, alike) stands for all of them."""
    same = [e for e in d["evidence"] if e.get("lab_runs") == [run]]
    keep = same[0]
    assert all(e.get("tables") == keep.get("tables") and e.get("results") == keep.get("results") for e in same), run
    gone = {e["id"] for e in same[1:]}
    d["evidence"] = [e for e in d["evidence"] if e["id"] not in gone]
    for c in d["claims"]:
        for sp in c["support"]:
            if sp.get("evidence") in gone or sp.get("evidence") == keep["id"]:
                sp["evidence"] = keep["id"]
                sp.pop("results", None)
    return keep


keep = one_object("RN8")
assert keep["id"] == "E1"
keep["description"] = ("Numerical checks of the formulas of the formal claims ([[D31:formula checks]]), lab run RN8 "
                       "under the formal plan that decided them; every check is within the plan's tolerances. They "
                       "check the formulas on concrete numbers; the proofs are the evidence for the claims. The rows are labelled by the hypotheses of the record's first plan, "
                       "whose wording kept the factor 200: H1 checks C1; H2 and H10 check C2; H5 checks C5; H7 checks "
                       "C6; H8 checks C7; H12 checks C8; H13 checks C9; H14 checks C10. H3 and H9 checked C3, and H4 "
                       "and H6 checked C4, as version 2 worded them; C3, C4 and C12 as now worded are checked in "
                       "E31. The two results E1.R1 and E1.R2 compare a simulated variance with the large-sample "
                       "variance that C19 now proves.")
keep["setup"] = {"studied": ["D1", "D3", "D4", "D5", "D6"], "procedure": ["D31"]}
# the sparse table of checks (one hypothesis per row, one check per column) as a list of checks
sparse = keep["tables"][0]
keep["results"] = list(keep.get("results") or [])
for row in sparse["rows"]:
    for c, v in zip(sparse["columns"], row["values"]):
        if v is not None:
            var = c["variable"] if isinstance(c, dict) else c
            keep["results"].append({"id": f"R{len(keep['results']) + 1}", "variable": var, "value": v,
                                    "conditions": {"hypothesis": row["label"][0]}})
keep["tables"] = keep["tables"][1:]
for i, t in enumerate(keep["tables"], 1):
    t["id"] = f"T{i}"
assert sorted(c["id"] for c in d["claims"] for sp in c["support"] if sp.get("evidence") == "E1") == sorted(
    ["C1", "C2", "C5", "C6", "C7", "C8", "C9", "C10"])

# RN13: the checks of the claims new or reworded in version 3. Its script labels each check by a hypothesis of the
# record (H3, H4, H34, H35, H39, H40, H42); a table whose rows are such labels, each with its own checks, becomes one
# table per label.
keep = one_object("RN13")
assert keep["id"] == "E31"
tables = []
for t in keep["tables"]:
    if t.get("row_labels") and t["row_labels"][0] == "hypothesis":
        for row in t["rows"]:
            cols = [(c, v) for c, v in zip(t["columns"], row["values"]) if v is not None]
            cond = dict(zip(t["row_labels"], row["label"]))
            names = [c["variable"] if isinstance(c, dict) else c for c, _ in cols]
            tables.append({"description": ", ".join(names) + "; " + ", ".join(f"{k}: {v}" for k, v in cond.items())
                                          + " (lab runs RN13)",
                           "row_labels": [], "columns": [c for c, _ in cols],
                           "rows": [{"id": "r1", "label": [], "values": [v for _, v in cols]}],
                           "conditions": cond, **({"artefact": t["artefact"]} if t.get("artefact") else {})})
    else:
        tables.append(t)
tables.sort(key=lambda t: ["H3", "H4", "H34", "H35", "H42", "H39", "H40"].index(t["conditions"]["hypothesis"]))
for i, t in enumerate(tables, 1):
    t["id"] = f"T{i}"
keep["tables"] = tables
RN13_LABELS = {"H3": "C3", "H4": "C4", "H34": "C12", "H35": "C17", "H42": "C20", "H39": "C18", "H40": "C19"}
assert {t["conditions"]["hypothesis"] for t in tables} == set(RN13_LABELS)
assert sorted(c["id"] for c in d["claims"] for sp in c["support"] if sp.get("evidence") == "E31") == sorted(
    ["C3", "C4", "C12", "C17", "C20"])


def tables_of(e, label):
    ids = [t["id"] for t in e["tables"] if t["conditions"]["hypothesis"] == label]
    assert ids == [f"T{i}" for i in range(int(ids[0][1:]), int(ids[-1][1:]) + 1)], ids
    return ids[0] if len(ids) == 1 else f"{ids[0]} and {ids[1]}" if len(ids) == 2 else f"{ids[0]} to {ids[-1]}"


keep["description"] = (
    "Numerical checks of the formulas of the claims new or reworded in version 3 ([[D31:formula checks]]), lab run "
    "RN13 (check_v3.py, seed 8) under the formal plan that decided C3, C4, C12, C17 and C20; every check is within "
    "the plan's tolerances. They check the formulas on concrete numbers; the proofs are the evidence for the claims. "
    "The tables are labelled by the hypotheses the run's script names: "
    + "; ".join(f"{h} ({tables_of(keep, h)}) checks {c}" for h, c in RN13_LABELS.items() if h not in ("H39", "H40"))
    + f". The tables of H39 ({tables_of(keep, 'H39')}) and H40 ({tables_of(keep, 'H40')}) check the formulas of C18 "
      "and C19, whose wording was decided after later checks (E32, E33).")
e32 = next(e for e in d["evidence"] if e.get("lab_runs") == ["RN18"])
e32["description"] = (
    "Numerical checks of the formulas of C18 ([[D31:formula checks]]), lab run RN18 (check_h39.py) under the formal "
    "plan that decided C18; every check is within the plan's tolerances. They check the formulas on concrete "
    "numbers; the proof is the evidence for the claim. The run's script labels the checks H39, the record's first "
    "wording of C18; the formulas are those of C18.")
e33 = next(e for e in d["evidence"] if e.get("lab_runs") == ["RN16"])
e33["description"] = (
    "Numerical checks of the formulas of C19 ([[D31:formula checks]]), lab run RN16 (check_corrected.py, seed 8) "
    "under the formal plan that decided C19; every check is within the plan's tolerances. They check the formulas "
    "on concrete numbers; the proof is the evidence for the claim. The run's script labels the checks by the "
    f"record's first wordings: H40 ({tables_of(e33, 'H40')} and the result R1) checks C19; H39 "
    f"({tables_of(e33, 'H39')}) checks the formulas of C18, which the same plan tested in an earlier wording.")
for c in d["claims"]:
    if c["id"] in ("C18", "C19"):
        for sp in c["support"]:
            sp.pop("results", None)
ev11 = next(e for e in d["evidence"] if e.get("lab_runs") == ["RN6"])
ev11["setup"] = {"studied": ["D6"], "procedure": ["D31"]}
for c in d["claims"]:
    if c["id"] == "C11":
        for sp in c["support"]:
            sp.pop("results", None)

# proof steps without a premise rest on the definitions they use
for e in d["evidence"]:
    if e.get("profile") != "formal":
        continue
    for st in e["steps"]:
        if not st.get("by"):
            text = st["statement"]
            st["by"] = [k for k, w in (("D10", "RMSE"), ("D11", "MAE"), ("D12", "MAPE")) if w in text and "(a)" not in text] \
                or [{"E30": "D28", "E37": "D33"}.get(e["id"], "D1")]

# dependencies, from the claims each proof cites
fe = {sp["evidence"]: c["id"] for c in d["claims"] for sp in c["support"] if sp["kind"] == "formal_proof"}
for e in d["evidence"]:
    if e.get("profile") == "formal":
        cid = fe[e["id"]]
        deps = sorted({b for st in e["steps"] for b in st.get("by", []) if b.startswith("C") and b != cid},
                      key=lambda x: int(x[1:]))
        if deps:
            claims[cid]["depends_on"] = deps
for cid in ("C2", "C8", "C11", "C13", "C14"):
    claims[cid]["role"] = "secondary"

# intervals: bootstrap and binomial intervals are confidence intervals
def kind(x):
    if isinstance(x, dict):
        if isinstance(x.get("interval"), dict):
            x["interval"].setdefault("kind", "confidence")
        for v in x.values():
            kind(v)
    elif isinstance(x, list):
        for v in x:
            kind(v)
kind(d["evidence"])

# statements: terms marked, long formulas on their own lines
claims["C14"]["statement"] = claims["C14"]["statement"].replace("less reliably than MAE", "less reliably than [[D11:MAE]]") \
    .replace("sMAPE separates its target", "sMAPE's [[D20:ranking agreement]] for its target")
claims["C14"]["statement"] = claims["C14"]["statement"].replace(
    "sMAPE's [[D20:ranking agreement]] for its target from a forecast half a $\\sigma$ away less reliably than [[D11:MAE]]",
    "the [[D20:ranking agreement]] of sMAPE, for its target against a forecast half a $\\sigma$ away, is lower than "
    "that of [[D11:MAE]]")
claims["C15"]["statement"] = claims["C15"]["statement"].replace("(log-symmetric actuals)", "([[D26:log-symmetric outcomes]])")
claims["C2"]["statement"] = (
    "For $y > 0$ and $f > 0$: (a) replacing $r = y/f$ by $1/r$ leaves the [[D1:sMAPE term]] unchanged, so a forecast "
    "$k$ times too high ($f = ky$) and one $k$ times too low ($f = y/k$) score the same; (b) of two forecasts the same "
    "absolute distance $e$ from $y$ ($0 < e < y$), the higher one scores better:\n$$s(y, y+e) = \\frac{e}{2y+e} < "
    "\\frac{e}{2y-e} = s(y, y-e).$$")
def display(cid, pairs):
    """Long formulas of a statement on lines of their own (claimset §8 rule 6): the wording does not change."""
    for a, b in pairs:
        assert claims[cid]["statement"].count(a) == 1, (cid, a)
        claims[cid]["statement"] = claims[cid]["statement"].replace(a, b)


display("C3", [(
    "satisfies $\\frac1n\\sum_t|\\ln y_t - \\ln f^*| - \\frac1n\\sum_t|\\ln y_t - \\ln m| \\le \\frac1{12}\\cdot\\frac1n\\sum_t|\\ln y_t - \\ln f^*|^3$: its mean",
    "satisfies\n$$\\frac1n\\sum_t|\\ln y_t - \\ln f^*| - \\frac1n\\sum_t|\\ln y_t - \\ln m| \\le \\frac1{12}\\cdot\\frac1n\\sum_t|\\ln y_t - \\ln f^*|^3:$$\nits mean")])
display("C4", [(
    ": $\\mathrm E[w(Y,f)\\mathbf 1\\{Y<f\\}] = \\mathrm E[w(Y,f)\\mathbf 1\\{Y>f\\}]$.",
    ":\n$$\\mathrm E[w(Y,f)\\,\\mathbf 1\\{Y<f\\}] = \\mathrm E[w(Y,f)\\,\\mathbf 1\\{Y>f\\}].$$\n")])
display("C12", [(
    "in addition, $2g(\\ln m) > \\mathrm E[\\operatorname{sech}^2((\\ln Y - \\ln m)/2)\\tanh(|\\ln Y - \\ln m|/2)]$, where",
    "in addition,\n$$2g(\\ln m) > \\mathrm E\\Big[\\operatorname{sech}^2\\Big(\\frac{\\ln Y - \\ln m}{2}\\Big)\\tanh\\Big(\\frac{|\\ln Y - \\ln m|}{2}\\Big)\\Big],$$\nwhere")])
display("C18", [
    ("as $\\sigma$ grows, $\\mathrm E[D] = 4\\ln 2\\,(h(c_1) - h(c_2))/\\sigma + O(\\sigma^{-3})$ and $\\operatorname{Var}(D) = 4(2\\ln 2 - 1)(h(c_1)+h(c_2))/\\sigma + O(\\sigma^{-2})$. So",
     "as $\\sigma$ grows,\n$$\\mathrm E[D] = \\frac{4\\ln 2\\,(h(c_1) - h(c_2))}{\\sigma} + O(\\sigma^{-3}), \\qquad \\operatorname{Var}(D) = \\frac{4(2\\ln 2 - 1)(h(c_1)+h(c_2))}{\\sigma} + O(\\sigma^{-2}).$$\nSo"),
    ("the ratio is $\\frac{2\\ln 2}{\\sqrt{2\\ln 2 - 1}}\\,\\frac{|h(c_1)-h(c_2)|}{\\sqrt{h(c_1)+h(c_2)}}\\,\\sigma^{-1/2}\\,(1 + O(1/\\sigma))$: it falls",
     "the ratio is\n$$\\frac{2\\ln 2}{\\sqrt{2\\ln 2 - 1}}\\,\\frac{|h(c_1)-h(c_2)|}{\\sqrt{h(c_1)+h(c_2)}}\\,\\sigma^{-1/2}\\,(1 + O(1/\\sigma)):$$\nit falls"),
    ("the ratio tends to $|2H(c_2) - 1|/(2\\sqrt{H(c_2)(1 - H(c_2))})$ when",
     "the ratio tends to\n$$\\frac{|2H(c_2) - 1|}{2\\sqrt{H(c_2)(1 - H(c_2))}}$$\nwhen")])
display("C19", [(
    "is $\\mathrm{ARE} = (2g(\\ln m) - A)^2/(4g(\\ln m)^2 B)$; that the",
    "is\n$$\\mathrm{ARE} = \\frac{(2g(\\ln m) - A)^2}{4\\,g(\\ln m)^2\\,B};$$\nthat the")])

# a disputed check answered by a correction that the checker then confirmed
for c in d["claims"]:
    for ck in c.get("checks") or []:
        if ck["outcome"] == "disputed" and "resolution" not in ck:
            ck["resolution"] = ("The statement was corrected, and the checker confirmed the corrected wording in a "
                                "later check listed here.")

defs = {x["id"]: x for x in d["definitions"]}
defs["D6"]["definition"] = ("The value or values of $f$ that minimise expected sMAPE for a given distribution of the "
                            "actual $Y$, $\\arg\\min_f \\mathrm E[s(Y,f)]$: the functional sMAPE rewards. It may be a "
                            "single value or several.")
for a in d["artefacts"]:
    a.setdefault("licence", "CC-BY-4.0")
d["paper"]["summary"] = "Which forecast does sMAPE reward? sMAPE has been criticised for its asymmetry and its behaviour near zero, but the forecast that minimises it had been found only in special cases. We show that for positive values sMAPE depends only on the ratio of actual to forecast, not on their difference: it is the hyperbolic tangent of half the absolute log error. So it scores a forecast k times too high and one k times too low the same, and is asymmetric only for equal absolute errors. Being bounded, it has a finite expected value for every distribution. Its target is a median of the actuals, each weighted by how close it is to the forecast by ratio. This is the median when the log of the actual is symmetric and unimodal; otherwise it can move to the side of the median that carries more weight. With two well-separated groups of actuals, the median becomes a local maximum. On intermittent demand the target is 0 once the share of zeros passes a threshold below one half. Otherwise sMAPE aims at a weighted median of the non-zero actuals. As the actuals spread out, expected sMAPE approaches its bound for every forecast, so forecasts become hard to tell apart, and the target moves from the median to the mode of the log. Where sMAPE and MAE both target the median, a test set's sMAPE minimiser estimates it less precisely than the sample median, increasingly so as the actuals spread out."  # the guarantor's text, 2026-11-17

# C9's proof uses only the definition for zero and negative forecasts (no cycle with C7)
for e in d["evidence"]:
    if e.get("profile") == "formal" and fe.get(e["id"]) == "C9":
        for st in e["steps"]:
            if "C7" in st.get("by", []):
                st["by"] = [b for b in st["by"] if b != "C7"] + ["D1"]
                st["statement"] = st["statement"].replace(" (C7)", "")
claims["C9"]["depends_on"] = [x for x in claims["C9"].get("depends_on", []) if x != "C7"]
if not claims["C9"]["depends_on"]:
    del claims["C9"]["depends_on"]
# planned estimates cite their tables, not every value
for cid in ("C13", "C14", "C15"):
    for sp in claims[cid]["support"]:
        e = next(x for x in d["evidence"] if x["id"] == sp["evidence"])
        tabs = [f"{e['id']}.{t['id']}" for t in e.get("tables") or []]
        if tabs:
            sp["results"] = tabs

claims["C8"]["statement"] = claims["C8"]["statement"].replace(
    "every difference $\\mathrm E[(Y-f)^2 - (Y-g)^2] = (f-g)(f+g-2\\,\\mathrm E[Y])$ is finite",
    "every difference\n$$\\mathrm E[(Y-f)^2 - (Y-g)^2] = (f-g)(f+g-2\\,\\mathrm E[Y])$$\nis finite")
claims["C13"]["relations"] = [{"type": "supports", "target": "C4", "id_type": "local",
                               "note": "Computed on 34 distributions, the target lies on the side of the median that "
                                       "carries more weight, as C4 and C12 prove."}]
claims["C14"]["relations"] = [{"type": "supports", "target": "C10", "id_type": "local",
                               "note": "Computed for spreads of ln Y from 0.1 to 10, expected sMAPE flattens and the "
                                       "target moves to the mode of ln Y, as C10 proves for large spreads."}]

for e in d["evidence"]:
    if e.get("profile") == "formal" and fe.get(e["id"]) == "C8":
        heads = {"S1": ("(a) RMSE. ", "D10"), "S5": ("(b) MAE. ", "D11"), "S8": ("(c) MAPE. ", "D12")}
        for st in e["steps"]:
            if st["id"] in heads:
                st["statement"] = heads[st["id"]][0] + st["statement"]
                st["by"] = [heads[st["id"]][1]]

d["references"] = [
    {"id": "B1", "target": "10.1016/S0169-2070(99)00007-2", "id_type": "doi"},
    {"id": "B2", "target": "Anderson, T. W. (1955). The integral of a symmetric unimodal function over a symmetric "
                           "convex set and some probability inequalities. Proceedings of the American Mathematical "
                           "Society, 6(2), 170-176.", "id_type": "reference"}]
claims["C2"].setdefault("relations", []).append(
    {"type": "refines", "reference": "B1",
     "note": "Goodwin and Lawton (1999) showed that sMAPE is asymmetric: for the same absolute error it scores a forecast "
             "above the actual better than one below. C2 proves this (b) and shows that by ratio sMAPE is symmetric (a)."})
claims["C5"]["novelty"] = "replication"
claims["C5"].setdefault("relations", []).append(
    {"type": "replicates", "reference": "B2",
     "note": "C5 is an elementary case of Anderson's theorem on symmetric unimodal densities, proved here directly for "
             "the expected sMAPE."})

# the figure of the weights is a diagram, not evidence
for e in d["evidence"]:
    if e.get("artefacts") and "A3" in e["artefacts"] and e.get("profile") == "formal":
        e["artefacts"] = [a for a in e["artefacts"] if a != "A3"]

# planned estimates cite the rows and results holding the numbers their statements give
def obs(e):
    for r in e.get("results") or []:
        if isinstance(r.get("value"), (int, float)):
            yield f"{e['id']}.{r['id']}", [r["value"]], r.get("conditions") or {}
    for t in e.get("tables") or []:
        for row in t.get("rows") or []:
            yield f"{e['id']}.{t['id']}.{row['id']}", [v for v in row.get("values") or [] if isinstance(v, (int, float))
                                                       and not isinstance(v, bool)], {"label": row.get("label")}

def rows_for(cid, numbers, keep=lambda cond: True):
    e = next(x for x in d["evidence"] if x["id"] == claims[cid]["support"][0]["evidence"])
    found = []
    for num in numbers:
        dec = len(num.split(".")[1]) if "." in num else 0
        for ref, vals, cond in obs(e):
            if keep(cond) and any(round(v, dec) == float(num) for v in vals) and ref not in found:
                found.append(ref)
                break
    return found

claims["C13"]["statement"] = (
    "Computed from the distributions themselves, the [[D6:sMAPE target]] on the 34 outcome distributions of plan PL1 "
    "is the median on lognormal outcomes, lies below it when $\\ln Y$ is skewed right and above it when skewed left "
    "(quantile levels from 0.439 to 0.556 on the unimodal positive distributions), is two forecasts, one near each "
    "mode, on equal bimodal outcomes and the larger mode's on unequal ones, lies away from zero on the side of the "
    "mean on normal outcomes straddling zero, and switches to 0 at shares of zero outcomes from 0.343 to 0.468, all "
    "below one half. On test sets the forecast sMAPE ranks first spreads widely when they are small: on lognormal "
    "outcomes with $\\sigma = 1$, its quantile level runs from 0.16 to 0.82 (5th to 95th percentile over test sets) "
    "with 6 outcomes and from 0.47 to 0.53 with 1000.")
claims["C13"]["support"][0]["results"] = rows_for("C13", ["0.439", "0.556", "0.343", "0.468", "0.16", "0.82", "0.47",
                                                          "0.53"])
claims["C14"]["support"][0]["results"] = rows_for("C14", ["7.97", "178.2", "0.372", "0.368", "0.61", "0.82"])
claims["C15"]["support"][0]["results"] = rows_for("C15", ["0.97", "0.03", "1.00", "0.34"])
print("cited", {k: claims[k]["support"][0]["results"] for k in ("C13", "C14", "C15")})

# every cited work once in references, cited by its ID (claimset §3.13): the EntryLib claims too
d["references"][1] = {"id": "B2", "target": "10.1090/S0002-9939-1955-0069229-1", "id_type": "doi"}
known = {r["target"]: r["id"] for r in d["references"]}
for c in d["claims"]:
    for rel in c.get("relations") or []:
        if rel.get("id_type") == "claimset":
            t = rel.pop("target"); rel.pop("id_type")
            if t not in known:
                known[t] = f"B{len(d['references']) + 1}"
                d["references"].append({"id": known[t], "target": t, "id_type": "claimset"})
            rel["reference"] = known[t]
            rel.pop("note", None) if False else None

# published work on sMAPE, each checked against Crossref and, where open, its text
LIT = [("B8", "10.1016/0169-2070(93)90079-3"), ("B9", "10.1016/j.ijforecast.2006.03.001"),
       ("B10", "10.1007/s10618-022-00894-5"), ("B11", "10.1371/journal.pone.0174202"),
       ("B12", "10.1057/jors.2014.103"), ("B13", "10.3390/forecast8020024")]
for rid, doi in LIT:
    assert rid not in {r["id"] for r in d["references"]}
    d["references"].append({"id": rid, "target": doi, "id_type": "doi"})
def rel(cid, typ, rid, note):
    claims[cid].setdefault("relations", []).append({"type": typ, "reference": rid, "note": note})
rel("C2", "refines", "B9", "Hyndman and Koehler (2006) note that for the same actual sMAPE penalises low forecasts more "
    "heavily than high ones; C2 shows this holds for equal absolute errors and not for equal ratios.")
rel("C2", "refines", "B10", "Hewamalage et al. (2023) note that sMAPE penalises underestimates more than overestimates "
    "for the same actual; C2 shows this holds for equal absolute errors and not for equal ratios.")
rel("C2", "refines", "B11", "Chen et al. (2017) show a heavier penalty on a forecast 10% under the actual than on one 10% "
    "over; C2 shows that forecasts the same factor too high and too low score the same.")
rel("C3", "extends", "B12", "Tofallis (2015) relates sMAPE to the log accuracy ratio among relative measures; C3 shows "
    "the sMAPE term is exactly the hyperbolic tangent of half the absolute log ratio.")
claims["C3"]["novelty"] = "extension"
rel("C4", "addresses", "B12", "Tofallis (2015) states that minimising sMAPE does not produce an established measure of "
    "location; C4 identifies the target as a median weighted by y/(y+f)^2.")
rel("C4", "qualifies", "B10", "Hewamalage et al. (2023) say sMAPE may tend towards selecting a slightly overestimating "
    "model; C4 shows the target is a weighted median, which lies on the side of the median that carries more weight, "
    "above or below it (C12).")
rel("C7", "extends", "B13", "Kolassa (2026) shows that sMAPE elicits a zero forecast for sufficiently intermittent "
    "actuals; C7 gives the share of zeros at which it does, which is below one half.")
claims["C7"]["novelty"] = "extension"
rel("C9", "extends", "B10", "Hewamalage et al. (2023) note that sMAPE is bounded by 200; C9 adds that its expected value "
    "is therefore finite for every distribution and that a best forecast exists.")
claims["C9"]["novelty"] = "extension"
rel("C15", "supports", "B13", "Kolassa (2026) concludes that sMAPE should not be used; C15 shows that even where it "
    "targets the median it estimates it less precisely than the sample median.")
rel("C10", "supports", "B13", "Kolassa (2026) concludes that sMAPE should not be used; C10 shows that on widely spread "
    "actuals every forecast scores close to the worst value.")
rel("C18", "supports", "B13", "Kolassa (2026) concludes that sMAPE should not be used; C18 shows that on widely spread "
    "actuals a test set tells two forecasts apart less and less well under sMAPE.")
next(x for x in d["definitions"] if x["id"] == "D2")["source"] = "B8"

rel("C1", "extends", "B12", "Tofallis (2015) lists the sMAPE term, (f-g)/[(f+g)/2], among Tornqvist et al.'s measures "
    "of relative change; C1 shows it depends on the actual and the forecast only through their ratio.")
claims["C1"]["novelty"] = "extension"

D2 = next(x for x in d["definitions"] if x["id"] == "D2")
D2["definition"] = D2["definition"].rstrip() + (" Introduced as the symmetric MAPE by Makridakis (1993), whose "
                                               "denominator had no absolute values; for positive actuals and "
                                               "forecasts the two agree.")

# ------------------------------------------------------------------ the Suilin loss (C16, C17)
# Version 2 called this measure msMAPE; that name is now Smyl's measure only (C20, C21).
d["references"].append({"id": "B14", "target": "10.1016/j.ijforecast.2020.06.008", "id_type": "doi"})
next(x for x in d["definitions"] if x["id"] == "D28")["source"] = "B14"
rel("C17", "extends", "B14", "Hewamalage, Bergmeir and Bandara (2021) evaluate forecasts with the Suilin loss, Suilin's "
    "modification of sMAPE for series with values near zero; C17 shows that it keeps sMAPE's behaviour away from zero "
    "and becomes an absolute error near zero.")
claims["C17"]["novelty"] = "extension"
claims["C16"]["statement"] = (
    "Computed from the distributions themselves, the [[D29:Suilin target]] differs little from the [[D6:sMAPE target]] "
    "away from zero and on counts, and differs near zero: it is within 0.3% of the sMAPE target on the six "
    "distributions with median about 100, and the same as the sMAPE target on every [[D30:Poisson counts]] distribution "
    "studied; on lognormal actuals with median 1 and $\\sigma = 1$ it is 0.985 instead of 1, and the share of zeros at "
    "which it switches to the forecast 0 is 0.470 instead of sMAPE's 0.394, while at median 100 the two shares are "
    "0.3937 and 0.3940.")
claims["C16"]["about"] = {"subject": "D29", "relation": "differs from", "object": "D6", "comparator": None,
                          "conditions": ["D30"]}
claims["C16"]["scope"] = {"data_type": "outcome distributions given exactly (no data)",
                          "conditions": ["positive distributions at median about 100, lognormal with zeros at medians 1 "
                                         "and 100, Poisson counts with means 0.3 to 10"]}
claims["C16"]["support"][0]["results"] = rows_for("C16", ["0.985", "0.470", "0.394", "0.3937", "0.3940"])
# the estimand names its variable as the table now does (the record: "... of the msMAPE target ...")
est = claims["C16"]["estimate"]["estimands"]["difference"]
assert est["variable"] == "relative difference of the msMAPE target from the sMAPE target"
est["variable"] = "relative difference of the Suilin target from the sMAPE target"
claims["C16"]["relations"] = [{"type": "supports", "target": "C17", "id_type": "local",
                               "note": "Computed on actual distributions, the Suilin target stays at the sMAPE target "
                                       "away from zero and moves near zero, as C17's two cases imply."}]
for l in d["limitations"]:
    if "C13" in l["affects"] and "C16" not in l["affects"] and l["category"] == "model_specification":
        l["affects"].append("C16")
d["paper"]["contributions"].append("C17")
e16 = next(e for e in d["evidence"] if e.get("lab_runs") == ["RN10"])
assert e16["description"].count("msMAPE") == 2 and "msmape_targets.py" in e16["description"]
e16["description"] = e16["description"].replace("the expected msMAPE term", "the expected Suilin term").replace(
    "relative difference of the msMAPE target", "relative difference of the Suilin target") + (
    " The plan, the run's script and its results file name the Suilin loss msMAPE, as version 2 of this paper did.")
claims["C17"]["about"] = {"subject": "D28", "relation": "behaves like, except near zero,", "object": "D1",
                          "comparator": None, "conditions": []}

# ------------------------------------------------------------------ version 3: C18 to C21
# C18, C19 and C20 are proved, like C17; C21 is computed, like C16.
d["paper"]["contributions"] += ["C18", "C19", "C20"]
assert all(claims[c]["type"] == "formal" and claims[c]["verification"] == "agent_checked" for c in ("C18", "C19", "C20"))

# works cited for the first time
d["references"] += [
    {"id": "B15", "target": "10.1017/CBO9780511802256", "id_type": "doi"},        # van der Vaart (1998)
    {"id": "B16", "target": "Smyl, S. (2025). Sparse-Proof sMAPE. Foresight: The International Journal of Applied "
                            "Forecasting, 78, 29–30. https://ideas.repec.org/a/for/ijafaa/y2025i78p29-30.html",
     "id_type": "reference"}]
assert [r["id"] for r in d["references"]] == [f"B{i}" for i in range(1, 17)]

# C19 applies two published theorems: cited support (claimset §5, Verification), named by the steps that use them
e36 = next(e for e in d["evidence"] if e["id"] == formal_ids["C19"])
used = {st["id"]: st for st in e36["steps"] for key in ("van der Vaart (1998, Theorem 5.7)",
                                                       "van der Vaart, 1998, Theorem 5.23", "The same theorem applies")
        if key in st["statement"]}
assert len(used) == 3, sorted(used)
for st in used.values():
    st["by"] = st["by"] + ["B15"]
    if "Theorem 5.7" not in st["statement"]:        # the step on Theorem 5.7 gives the argument itself
        st["rule"] = "cited theorem"
claims["C19"]["support"].append({"kind": "cited", "reference": "B15", "locator": "Theorems 5.7 and 5.23"})

# C20: the msMAPE of Smyl, as Kolassa (2026, Definition 4) states it
next(x for x in d["definitions"] if x["id"] == "D33")["source"] = "B13"
claims["C20"]["about"] = {"subject": "D33", "relation": "equals, except for actuals below its constant,",
                          "object": "D1", "comparator": None, "conditions": []}
rel("C20", "extends", "B13", "Kolassa (2026, Proposition 3) shows that the minimisers of expected msMAPE on counts are "
    "integers; C20 shows where the msMAPE term equals the sMAPE term and that it removes the jump at the zero "
    "forecast.")
rel("C20", "extends", "B16", "Smyl (2025) proposed msMAPE, the sparse-proof sMAPE, as Kolassa (2026) reports; Smyl's "
    "article was not read. C20 states properties of its term as Kolassa (2026, Definition 4) defines it.")
claims["C20"]["novelty"] = "extension"

# C21: the planned estimate of PL14 (run RN15), stated as C16 is: what its tables show
E34 = next(e for e in d["evidence"] if e.get("lab_runs") == ["RN15"])
assert E34["id"] == "E34" and E34["description"].endswith("a panel per measure); no test is made and no threshold is set.")
E34["description"] += " The figure as drawn (A56) has a panel per family of settings and a line per measure."


def cell(table, setting, measure, n=None):
    t = next(t for t in E34["tables"] if t["id"] == table)
    row = next(r for r in t["rows"] if r["label"] == ([setting] if n is None else [setting, n]))
    col = next(i for i, c in enumerate(t["columns"]) if isinstance(c, dict) and c["conditions"].get("measure") == measure)
    return row["values"][col]


_counts = [r["label"][0] for r in next(t for t in E34["tables"] if t["id"] == "T5")["rows"] if r["label"][1] == 12]
_all = [r["label"][0] for r in next(t for t in E34["tables"] if t["id"] == "T1")["rows"]]
_nu = ["msMAPE nu=0.1", "msMAPE nu=0.5", "msMAPE nu=1"]
assert len(_all) == 14 and len(_counts) == 10
# (a) what the statement says of the targets on counts
assert all(cell("T1", s_, "msMAPE nu=0.1") == cell("T1", s_, "Suilin") == cell("T1", s_, "sMAPE") for s_ in _counts)
_differs = lambda m: {s_: (cell("T1", s_, "sMAPE"), cell("T1", s_, m)) for s_ in _counts
                      if cell("T1", s_, m) != cell("T1", s_, "sMAPE")}
assert _differs("msMAPE nu=0.5") == {"Poisson mean=0.7": (0, 1), "negative binomial mean=2 size=0.5": (0, 1),
                                     "Poisson mean=1.5": (2, 1)}
assert _differs("msMAPE nu=1") == {**_differs("msMAPE nu=0.5"), "Poisson mean=0.5": (0, 1),
                                   "negative binomial mean=1 size=0.5": (0, 1)}
# (b) the zero thresholds
_thr = {m: cell("T2", "lognormal sizes with zeros", m) for m in ["sMAPE", "Suilin"] + _nu}
assert [f"{_thr[m]:.3f}" for m in ("sMAPE", "msMAPE nu=0.1", "Suilin", "msMAPE nu=0.5", "msMAPE nu=1")] == [
    "0.394", "0.419", "0.470", "0.519", "0.565"], _thr
# (c) test sets: 0 no more often under msMAPE than under sMAPE, nor under a larger nu
for s_ in _all:
    for n_ in (12, 52, 365):
        v = [cell("T4", s_, m, n_) for m in ["sMAPE"] + _nu]
        assert v[0] >= v[1] >= v[2] >= v[3], (s_, n_, v)
assert [cell("T4", "Poisson mean=0.7", m, 52) for m in ["sMAPE"] + _nu] == [0.6575, 0.523, 0.1525, 0.0265]
claims["C21"]["statement"] = (
    "The forecasts that sMAPE, msMAPE with $\\nu = 0.1$, 0.5 and 1, and the Suilin loss reward were computed from the "
    "distributions themselves on the 14 [[D35:intermittent demand settings]]. (a) On each of the 10 count "
    "distributions, the [[D34:msMAPE target]] with $\\nu = 0.1$ and the [[D29:Suilin target]] equal the "
    "[[D6:sMAPE target]]. With $\\nu = 0.5$ the msMAPE target is 1 instead of 0 on Poisson counts with mean 0.7 and on "
    "negative binomial counts with mean 2, and 1 instead of 2 on Poisson counts with mean 1.5. With $\\nu = 1$ it is "
    "also 1 instead of 0 on Poisson counts with mean 0.5 and on negative binomial counts with mean 1. (b) For "
    "lognormal demand sizes with zeros, the share of zeros at which the target switches to the forecast 0 is 0.394 "
    "for sMAPE, 0.419 for msMAPE with $\\nu = 0.1$, 0.470 for the Suilin loss, 0.519 for msMAPE with $\\nu = 0.5$ and "
    "0.565 with $\\nu = 1$. (c) In test sets of 12, 52 and 365 actuals, the forecast with the lowest score is 0 no "
    "more often under msMAPE than under sMAPE, and no more often under a larger $\\nu$. With 52 actuals from Poisson "
    "counts with mean 0.7, it is 0 in 0.6575 of the test sets under sMAPE, 0.523 with $\\nu = 0.1$, 0.1525 with "
    "$\\nu = 0.5$ and 0.0265 with $\\nu = 1$.")
C21_SUMMARY = (
    "On counts, msMAPE with ν = 0.1 and the Suilin loss reward the same forecast as sMAPE. With ν = 0.5 or 1 it is 0 "
    "less often. On non-integer demand the zero forecast needs a larger share of zeros: sMAPE 0.394; msMAPE 0.419, "
    "0.519, 0.565; Suilin 0.470.")
claims["C21"]["about"] = {"subject": "D34", "relation": "differs from", "object": "D6", "comparator": "D29",
                          "conditions": ["D35"]}
claims["C21"]["scope"] = {
    "data_type": "outcome distributions given exactly, and test sets drawn from them",
    "conditions": ["14 settings: 10 count distributions (Poisson with means 0.1 to 3, negative binomial with size 0.5 "
                   "and means 0.5 to 2) and lognormal demand sizes (median 1, sigma 1) with shares of zeros 0.2 to 0.8",
                   "msMAPE with nu = 0.1, 0.5 and 1",
                   "test sets of 12, 52 and 365 independent actuals, 2000 test sets per setting and size"]}
# the observations the statement reports: (a) the targets on the 10 count distributions, (b) the zero thresholds,
# (c) the shares of test sets in every setting (the whole table) and the row its example quotes
_row = lambda table, label: next(f"E34.{table}.{r['id']}" for t in E34["tables"] if t["id"] == table
                                 for r in t["rows"] if r["label"] == label)
claims["C21"]["support"][0]["results"] = ([_row("T1", [s_]) for s_ in _counts]
                                          + [_row("T2", ["lognormal sizes with zeros"]), "E34.T4",
                                             _row("T4", ["Poisson mean=0.7", 52])])

# limitations: version 2's L3 and L4 also hold for C21; three new ones
for l in d["limitations"]:
    if l["id"] in ("L3", "L4"):
        assert "C21" not in l["affects"]
        l["affects"].append("C21")
assert len(d["limitations"]) == 4
d["limitations"] += [
    {"id": "L5", "category": "model_specification", "affects": ["C6"], "impact": "reduces_transferability",
     "mitigation": "none",
     "statement": "C6 is proved for two groups that mirror each other about the median on the log scale. It is not "
                  "proved for groups of different shape or spread."},
    {"id": "L6", "category": "measurement", "affects": ["C20", "C21"], "impact": "reduces_transferability",
     "mitigation": "none",
     "statement": "The definition of msMAPE is taken from Kolassa (2026, Definition 4). Smyl's article was not read. "
                  "Smyl gives no value for the constant; ν = 0.1, 0.5 and 1 are this paper's choices."},
    {"id": "L7", "category": "other", "affects": ["C19"], "impact": "reduces_transferability",
     "mitigation": "partial", "mitigation_evidence": "E15",
     "statement": "That the relative efficiency of finite test sets tends to the asymptotic value is not shown; the "
                  "simulation of C15 agrees with it at n = 10000."}]

# ------------------------------------------------------------------ what a reader finds their way by
# (claimset v0.44 as amended on 2026-11-17: names, summaries, symbols, the results of a relation). Version 2's, with
# the names and summaries the guarantor asked for in version 3 (C2 to C5, C12, C16, C17) and those of C18 to C21.
NAMES = {
    "C1": "sMAPE depends only on the ratio of actual to forecast",
    "C2": "Too high and too low by the same factor score the same",
    "C3": "sMAPE is a hyperbolic tangent of the log error",
    "C4": "sMAPE targets a weighted median",
    "C5": "Symmetric, unimodal logs: the target is the median",
    "C6": "Two separate groups: the median is a local maximum",
    "C7": "Zero actuals: when the target is zero",
    "C8": "Where RMSE, MAE and MAPE fail, and what they then target",
    "C9": "Expected sMAPE is always finite and has a minimiser",
    "C10": "Expected sMAPE flattens as the actuals spread out",
    "C11": "Actuals of both signs: the target jumps across zero",
    "C12": "When sMAPE targets the ordinary median",
    "C13": "The sMAPE target on 34 distributions",
    "C14": "Flattening, computed",
    "C15": "sMAPE finds the median less precisely than the sample median",
    "C16": "The Suilin target beside the sMAPE target",
    "C17": "The Suilin loss behaves like sMAPE except near zero",
    "C18": "sMAPE loses its power to tell forecasts apart as the actuals spread out",
    "C19": "The sMAPE minimiser in large samples",
    "C20": "msMAPE is sMAPE except for small actuals",
    "C21": "What msMAPE rewards on intermittent demand",
}
SUMMARIES = {
    "C1": "Only the ratio of the actual to the forecast matters to sMAPE, not their size.",
    # the earlier part shortened to keep the summary within 50 words with the guarantor's addition
    "C2": "A forecast twice too high and one half too low score the same. For equal absolute errors the higher one "
          "scores better. This compares forecasts for one known actual. Which forecast is best when the actual is "
          "uncertain is claim C4: sMAPE does not in general reward forecasting too high.",
    "C3": "sMAPE measures the error on the log scale. For small errors it is half the absolute log error; for large "
          "errors it levels off. For small errors its best forecast is close to the median, but not in general the "
          "median.",
    "C4": "For any distribution, the forecast sMAPE rewards balances the actuals below it and above it. Each actual "
          "is weighted by how close it is to the forecast by ratio.",
    "C5": "When the logarithm of the actual is symmetric with a single peak, the weights balance at the median, and "
          "the median is the forecast sMAPE rewards. Lognormal and log-Laplace actuals are examples.",
    "C6": "When the actuals form two separate groups, the median falls in the gap between them, far from both. So "
          "sMAPE scores it worse than forecasts near either group.",
    "C7": "The worst possible score goes to forecasts of the wrong sign and to zero forecasts. With enough zeros "
          "among the actuals, sMAPE prefers the forecast 0. The share of zeros needed is below one half.",
    "C8": "RMSE, MAE and MAPE each have an infinite expected value, for every forecast, when a moment they need is "
          "infinite. What each then rewards depends on the measure.",
    "C9": "The expected sMAPE of a forecast is always finite, even when the actuals have a heavy tail that makes "
          "the expected RMSE, MAE or MAPE infinite. And there is always a best forecast.",
    "C10": "For widely spread actuals, every forecast scores close to the worst value. The differences between "
           "forecasts shrink, and the forecast sMAPE rewards moves from the median to the most likely value on the "
           "log scale.",
    "C11": "When actuals can be negative or positive, the forecast sMAPE rewards is never near zero. It is on the "
           "same side of zero as the mean. So a small shift of the mean across zero makes it jump from a negative "
           "value to a positive one.",
    "C12": "The forecast sMAPE rewards lies above the ordinary median when the actuals above it carry more weight, "
           "and below it when those below do. When the weights balance and enough actuals lie near the median, the "
           "median is rewarded.",
    "C13": "Computed on 34 distributions, the forecast sMAPE rewards is the median on lognormal data. It moves away "
           "from the median on skewed data, sits near a mode on two-group data, and becomes zero once enough of the "
           "data are zeros, always fewer than half.",
    "C14": "As the data spread out, the best possible sMAPE rises towards its worst value, the forecast it rewards "
           "moves from the median towards the most likely value on the log scale, and sMAPE ranks its own best "
           "forecast first less often than MAE does.",
    "C15": "Even where sMAPE and MAE both reward the median, the forecast that minimises sMAPE on a test set varies "
           "more than the sample median. The loss of precision grows as the data spread out.",
    "C16": "Computed on actual distributions, the forecast the Suilin loss rewards is almost the same as sMAPE's "
           "away from zero and on counts. It differs near zero, where the Suilin loss needs a larger share of zeros "
           "before it rewards a zero forecast.",
    "C17": "The Suilin loss behaves like sMAPE except near zero, where it becomes an absolute error.",
    "C18": "For widely spread actuals a test set tells two forecasts apart less and less well under sMAPE: the "
           "signal-to-noise ratio of their score difference falls like one over the square root of the spread. The "
           "absolute log error and MAE keep their power.",
    "C19": "Where the log of the actual is symmetric with one peak, the forecast that minimises sMAPE on a test set "
           "is asymptotically normal about the median. It is asymptotically as precise as the sample median for "
           "narrow distributions and loses its precision as the actuals spread out.",
    "C20": "Smyl's msMAPE equals sMAPE for actuals of at least its constant ν. Below ν its denominator stops "
           "shrinking. A zero actual then no longer gives every non-zero forecast the worst score, and the expected "
           "score has no jump at the zero forecast.",
    "C21": C21_SUMMARY,
}
for cid, c in claims.items():
    c["name"], c["summary"] = NAMES[cid], SUMMARIES[cid]
    # in reading order of the fields: id, name, summary, statement, the rest
    keep = dict(c)
    c.clear()
    c.update({k: keep[k] for k in ("id", "name", "summary", "statement") if k in keep})
    c.update({k: v for k, v in keep.items() if k not in c})
SYMBOLS = {"D1": "$s(y,f)$", "D4": "$z$", "D5": "$w(y,f)$", "D8": "$\\hat f_n$", "D9": "$\\mathrm{RE}(n)$",
           "D28": "$\\tilde s(y,f)$", "D33": "$s_\\nu(y,f)$"}
for x in d["definitions"]:
    if x["id"] in SYMBOLS:
        x["symbol"] = SYMBOLS[x["id"]]


def rows(eid, tables, settings):
    """The rows of the named tables of an evidence object whose first label is one of the settings."""
    e_ = next(x for x in d["evidence"] if x["id"] == eid)
    out = [f"{eid}.{t['id']}.{r['id']}" for t in e_["tables"] if t["id"] in tables for r in t["rows"]
           if r.get("label") and r["label"][0] in settings]
    assert out, (eid, tables, settings)
    return out


def local(cid, typ, target, note, results):
    rels = claims[cid].setdefault("relations", [])
    have = next((r for r in rels if r.get("id_type") == "local" and r["target"] == target), None)
    if have is None:
        have = {"type": typ, "target": target, "id_type": "local", "note": note}
        rels.append(have)
    if results:
        have["results"] = results


ZEROS = [f"lognormal sigma={sg} with zeros p0={p0}" for sg, ps in
         (("0.5", ("0.25", "0.42", "0.47", "0.6")), ("1", ("0.25", "0.37", "0.42", "0.6"))) for p0 in ps]
E13 = next(e["id"] for e in d["evidence"] if e.get("lab_runs") == ["RN1"])
E14 = next(e["id"] for e in d["evidence"] if e.get("lab_runs") == ["RN2"])
E15 = next(e["id"] for e in d["evidence"] if e.get("lab_runs") == ["RN7"])
E16 = next(e["id"] for e in d["evidence"] if e.get("lab_runs") == ["RN10"])
# the simulations beside the proved claims: T1 holds the targets, T2 a second minimum, T3 the rows the figures show
local("C13", "supports", "C4", None,
      rows(E13, ("T1", "T3"), ["log-skewed right sd=1 skew=+1", "log-skewed left sd=1 skew=-1",
                               "inverse Gaussian shape=10"]))
local("C13", "supports", "C5", "Computed on lognormal outcomes, the target is the median, as C5 proves.",
      rows(E13, ("T1", "T3"), ["lognormal sigma=0.5", "lognormal sigma=1"]))
local("C13", "supports", "C6", "Computed on bimodal outcomes, the target is two forecasts, one near each mode, when "
      "the modes are equal, and the larger mode's when they are not, as C6 proves.",
      rows(E13, ("T1", "T2", "T3"), ["bimodal equal modes", "bimodal left mode higher", "bimodal right mode higher"]))
local("C13", "supports", "C7", "Computed on lognormal outcomes with zeros, the target switches to 0 at a share of "
      "zeros below one half, as C7 proves.", rows(E13, ("T1", "T3"), ZEROS))
local("C13", "supports", "C11", "Computed on normal outcomes straddling zero, the target lies away from zero, on the "
      "side of the mean, as C11 proves.",
      rows(E13, ("T1", "T2", "T3"), [f"normal straddling zero mean={m}" for m in ("0", "0.5", "1")]))
local("C14", "supports", "C10", None, [f"{E14}.T1"])
local("C15", "qualifies", "C12", "Where C12 makes the median the target, the forecast that minimises sMAPE on a test "
      "set finds it less precisely than the sample median.", [f"{E15}.T1"])
local("C16", "supports", "C17", None, [f"{E16}.T1", f"{E16}.T2"])

# ------------------------------------------------------------------ version 3: relations between the paper's claims
E33 = next(e for e in d["evidence"] if e.get("lab_runs") == ["RN16"])
local("C4", "qualifies", "C2", "C2 compares two forecasts of one known actual. For an uncertain actual, the forecast "
      "sMAPE rewards is the weighted median of C4, not a forecast that is too high.", None)
local("C18", "refines", "C10", "C10 shows that the scores of forecasts at fixed quantile levels move together like "
      "1/σ. C18 adds the noise and shows that the signal-to-noise ratio falls like σ^(-1/2), while that of the "
      "absolute log error does not change.", None)
# C19's own evidence that bears on C15: the formula at the 16 distributions of C15, and its distance from the
# large-sample values C15's run (RN7) recorded
_t = {t["columns"][0] if isinstance(t["columns"][0], str) else t["columns"][0]["variable"]: t["id"] for t in E33["tables"]}
local("C19", "explains", "C15", "C15 found by simulation that the sMAPE minimiser of a test set varies more than the "
      "sample median, the more so the wider the spread. C19 gives the limiting variances of both and the rates at "
      "which their ratio falls.",
      [f"E33.{_t['large-sample relative efficiency, formula']}",
       f"E33.{_t['distributions where A is not below 2 g at the median']}"])
local("C21", "supports", "C20", "Computed on intermittent demand, msMAPE needs a larger share of zeros than sMAPE "
      "before its target is 0, in line with C20 (b): a zero actual no longer gives every non-zero forecast the worst "
      "score.", ["E34.T1", "E34.T2"])

# ------------------------------------------------------------------ sections, and the definitions in reading order
d["paper"]["sections"] = [
    {"title": "What sMAPE measures", "claims": ["C1", "C2", "C3"]},
    {"title": "What sMAPE rewards for positive actuals", "claims": ["C9", "C4", "C5", "C12", "C6"]},
    {"title": "Zeros and negative actuals", "claims": ["C7", "C11"]},
    {"title": "Widely spread actuals", "claims": ["C10", "C18"]},
    {"title": "Other measures", "claims": ["C8"]},
    {"title": "sMAPE modifications", "claims": ["C17", "C20", "C16", "C21"]},
    {"title": "Precision on a test set", "claims": ["C15", "C19"]},
    {"title": "Numerical checks", "claims": ["C13", "C14"]}]
assert sorted(c for sec in d["paper"]["sections"] for c in sec["claims"]) == sorted(claims)

# sMAPE first and in full, then its term; a term after those it names; the paper's objects before designs and settings
defs = {x["id"]: x for x in d["definitions"]}
assert defs["D2"]["definition"].startswith(
    "Symmetric mean absolute percentage error of a point forecast $f$ over a test set of actuals $y_1, \\dots, y_n$: "
    "$\\frac{200}{n}\\sum_t s(y_t, f)$, 200 times the mean of the sMAPE terms. Introduced as the symmetric MAPE")
defs["D2"]["definition"] = (
    "Symmetric mean absolute percentage error of a point forecast $f$ over a test set of actuals $y_1, \\dots, y_n$: "
    "$\\frac{200}{n}\\sum_t |y_t - f|/(|y_t|+|f|)$, where a term with $y_t = f = 0$ counts as 0. Introduced as the "
    "symmetric MAPE by Makridakis (1993), whose denominator had no absolute values; for positive actuals and "
    "forecasts the two agree.")
assert defs["D1"]["definition"].endswith("(including $y = f = 0$, where the fraction is $0/0$). sMAPE is 200 times its mean.")
defs["D1"]["definition"] = (
    "The contribution of one actual $y$ and one point forecast $f$ to [[D2:sMAPE]]: $s(y,f) = |y-f|/(|y|+|f|)$, "
    "between 0 and 1, with $s = 0$ when $y = f$ (including $y = f = 0$, where the fraction is $0/0$).")
assert defs["D8"]["definition"] == ("The forecast that minimises the sMAPE of one test set, found numerically (the "
                                    "sample counterpart of the sMAPE target).")
defs["D8"]["definition"] = (
    "The forecast that minimises the sMAPE of one test set $y_1, \\dots, y_n$, "
    "$\\hat f_n = \\arg\\min_f \\frac1n\\sum_t s(y_t, f)$, found numerically (the sample counterpart of the "
    "[[D6:sMAPE target]]).")
ORDER = ["D2", "D1", "D3", "D6", "D4", "D5", "D7",                 # sMAPE, its term, its expectation and its target
         "D10", "D11", "D12", "D19",                               # the other measures
         "D28", "D29", "D33", "D34",                               # the modifications of sMAPE
         "D8", "D9", "D20", "D32",                                 # on a test set
         "D13", "D14", "D15", "D16", "D17", "D18", "D31",          # how the values were obtained
         "D21", "D22", "D23", "D24", "D25", "D26", "D27", "D30", "D35"]     # the distributions studied
assert sorted(ORDER) == sorted(defs), sorted(set(ORDER) ^ set(defs))
d["definitions"] = [defs[k] for k in ORDER]

# version 3 of the published paper: its identifier, and what changed and why (claimset §9)
d["paper"]["id"] = "entrylib:2026.000013"
d["paper"]["correction_reason"] = (
    "Version 3 answers a reader's comments. Corrected: C3 no longer says sMAPE is minimised by the median for small "
    "errors, only near it, with a bound; C4 keeps the weighted median, and C12 now holds what was said about the "
    "ordinary median; the measure called msMAPE in version 2 is the Suilin loss. Added: how well a test set tells two "
    "forecasts apart as the actuals spread out (C18); the large-sample distribution of the sMAPE minimiser (C19); "
    "Smyl's msMAPE (C20) and what it rewards on intermittent demand (C21). Also added: sections, the definition of "
    "sMAPE first, and a limitation on C6.")

import notation
import notation
d = notation.apply(d)
for _b in notation.leftovers(d):
    print("NOTATION LEFT:", _b)
Path(OUT).write_text(json.dumps(d, ensure_ascii=False, indent=1) + "\n")
print("written", OUT, len(d["claims"]), "claims", len(d["evidence"]), "evidence")
