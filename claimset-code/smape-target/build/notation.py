"""One notation across the claimset (agreed 2026-11-17, extended for version 3 on 2026-11-18): y is the actual, as a
value and as a random quantity, with E the expectation over it; m is the median of y; g is the density of ln y; s is
the sMAPE term, s~ the Suilin term (called the msMAPE term in version 2) and s_nu the msMAPE term of Smyl, which the
record writes p; epsilon is the random error of a stated model; the constant of the Suilin term is delta and that of
the msMAPE term nu. Applied last, to the finished document."""
import re

SYM_Y = re.compile(r"(?<![A-Za-z\\])Y(?![A-Za-z])")


def _rep(text, pairs, where):
    for a, b in pairs:
        n = text.count(a)
        assert n >= 1, (where, a)
        text = text.replace(a, b)
    return text


CLAIMS = {
    "C2": [("absolute distance $e$ from $y$ ($0 < e < y$)", "absolute distance $d$ from $y$ ($0 < d < y$)"),
           ("$$s(y, y+e) = \\frac{e}{2y+e} < \\frac{e}{2y-e} = s(y, y-e).$$",
            "$$s(y, y+d) = \\frac{d}{2y+d} < \\frac{d}{2y-d} = s(y, y-d).$$")],
    "C4": [("Let $Y$ be positive and take no single value", "Let the actual $y$ be positive and take no single value")],
    "C5": [("Let $\\ln Y$ have a density $g$ symmetric about $m$, $g(m+v) = g(m-v)$, and unimodal: $g(m+u) \\ge g(m+v)$ "
            "whenever $|u| \\le |v|$, strictly on a set of positive probability. Then the median $e^m$ of $Y$ is the only",
            "Let $\\ln y$ have a density $g$ symmetric about $\\ln m$, $g(\\ln m+v) = g(\\ln m-v)$, and unimodal: "
            "$g(\\ln m+v_1) \\ge g(\\ln m+v_2)$ whenever $|v_1| \\le |v_2|$, strictly on a set of positive probability. "
            "Then the median $m$ of $y$ is the only")],
    "C6": [("Let $\\ln Y$ have two modes, at $m - a$ and $m + a$ ($a > 0$), with the actuals concentrated near them: a "
            "mixture of a distribution near $m - a$ and its mirror image about $m$,",
            "Let $\\ln y$ have two modes, at $\\ln m - a$ and $\\ln m + a$ ($a > 0$), with the actuals concentrated near "
            "them: a mixture of a distribution near $\\ln m - a$ and its mirror image about $\\ln m$,"),
           ("the median $e^m$ lies in the gap", "the median $m$ lies in the gap"),
           ("these approach $e^{m-a}$ and $e^{m+a}$ as the groups narrow", "these approach $m\\,e^{-a}$ and $m\\,e^{a}$ as the groups narrow")],
    "C7": [("Let $Y$ be 0 with probability $p_0$ and otherwise positive, and $m_+$ the lowest",
            "Let the actual $y$ be 0 with probability $p_0$ and otherwise positive, and $s_+$ the lowest"),
           ("$p_0 > (1 - m_+)/(2 - m_+)$", "$p_0 > (1 - s_+)/(2 - s_+)$"), ("whenever $m_+ > 0$", "whenever $s_+ > 0$")],
    "C8": [("For a positive actual $Y$:", "For a positive random actual $y$:"),
           ("every difference\n$$\\mathrm E[(Y-f)^2 - (Y-g)^2] = (f-g)(f+g-2\\,\\mathrm E[Y])$$\nis finite",
            "every difference between two forecasts $f_1$ and $f_2$,\n"
            "$$\\mathrm E[(y-f_1)^2 - (y-f_2)^2] = (f_1-f_2)(f_1+f_2-2\\,\\mathrm E[y]),$$\nis finite")],
    "C9": [("For a positive random actual $Y$", "For a positive random actual $y$")],
    "C10": [("has a fixed density $g$ with", "has a fixed density $h$ with"),
            ("4\\ln 2\\,g(c)/\\sigma", "4\\ln 2\\,h(c)/\\sigma"), ("When $g$ has a single highest point", "When $h$ has a single highest point")],
    "C12": [("Let $Y$ be positive and take no single value", "Let the actual $y$ be positive and take no single value")],
    "C17": [("Where $|y|+|f| \\ge 0.5$, $\\tilde s(y,f) = s(y,f)\\,\\frac{|y|+|f|}{|y|+|f|+\\delta}$: the sMAPE term",
             "Where $|y|+|f| \\ge 0.5$,\n$$\\tilde s(y,f) = s(y,f)\\,\\frac{|y|+|f|}{|y|+|f|+\\delta}:$$\nthe sMAPE term")],
    # the msMAPE term of Smyl: s_nu in the paper and in its proof (proof-suilin.md), p in the recorded statement
    "C20": [("[[D33:msMAPE term]] $p$ is the", "[[D33:msMAPE term]] $s_\\nu$ is the"),
            ("$p(y,f) = s(y,f)$", "$s_\\nu(y,f) = s(y,f)$"),
            ("$p(y,f) = |f-y|/(f+\\nu)$", "$s_\\nu(y,f) = |f-y|/(f+\\nu)$"),
            ("$0 \\le p \\le 1$, with $p < 1$ where", "$0 \\le s_\\nu \\le 1$, with $s_\\nu < 1$ where"),
            ("$p(y,f)$ is continuous", "$s_\\nu(y,f)$ is continuous"),
            ("$p$ differs from $s$", "$s_\\nu$ differs from $s$")],
}

DEFS = {
    "D3": "The expected sMAPE of a fixed forecast $f$ under the distribution of the actual $y$: "
          "$200\\,\\mathrm{E}[s(y, f)]$, where $\\mathrm E$ is the expectation over the random actual $y$. Computed from "
          "the distribution (by integration or exact sums), never from a sample.",
    "D4": "For a positive actual $y$ and a positive forecast $f$, $z = \\ln y - \\ln f$.",
    "D7": "For an outcome distribution with quantile function $Q$ and a forecast $f$, the pull of the outcome at each "
          "quantile level $u$, $u \\mapsto w(Q(u), f)$ on $(0,1)$; the minimiser's quantile level balances its integrals "
          "below and above.",
    "D9": "For a test-set size $n$ and an outcome distribution where sMAPE and MAE both target the median: "
          "$\\mathrm{RE}(n) = \\operatorname{Var}(\\ln \\hat m_n) / \\operatorname{Var}(\\ln \\hat f_n)$ over test sets, "
          "where $\\hat m_n$ is the sample median (the MAE minimiser) and $\\hat f_n$ the sMAPE minimiser of the same test "
          "set. Above 1, sMAPE's minimiser varies less; below 1, MAE's does.",
    "D10": "Root mean squared error of a point forecast $f$ over a test set of actuals $y_1, \\dots, y_n$: "
           "$\\sqrt{\\frac1n\\sum_t (y_t-f)^2}$.",
    "D11": "Mean absolute error of a point forecast $f$ over a test set of actuals $y_1, \\dots, y_n$: "
           "$\\frac1n\\sum_t |y_t-f|$.",
    "D12": "Mean absolute percentage error of a point forecast $f$ over a test set of actuals $y_1, \\dots, y_n$: "
           "$\\frac{100}{n}\\sum_t |y_t-f|/|y_t|$; infinite when an actual is zero.",
    "D15": "A value computed from the outcome distribution itself by numerical integration of its density "
           "(Gauss-Legendre quadrature and a convolution on a grid of 524288 points of $\\ln y$), never from a sample: "
           "the expected scores, the local minima of expected sMAPE, the formulas of the hypotheses.",
    "D16": "For each outcome distribution: derive the candidate forecasts and their expected scores; then, for each "
           "test-set size $n$, draw test sets of $n$ independent outcomes, score every candidate by every measure on each "
           "test set, and count the candidate with the lowest score as the measure's choice (a tie of $j$ candidates "
           "counts $1/j$ each; a test set in which every candidate scores infinite chooses none).",
    "D17": "For each outcome distribution and test-set size $n$: draw test sets of $n$ independent outcomes; in each, "
           "compute the sample median and the exact sMAPE minimiser (an outcome of the test set, since the test set's "
           "sMAPE is concave in $\\ln f$ between consecutive outcomes); over test sets, compare the variances of their "
           "logarithms.",
    "D18": "For families whose spread $\\sigma$ of $\\ln y$ grows from 0.1 to 10: derive the sMAPE target, the expected "
           "sMAPE of forecasts at fixed quantile levels and the large-spread approximation of H14; on test sets, count how "
           "often each measure scores the target lower than a forecast away from it, and where the test set's sMAPE "
           "minimiser falls.",
    "D20": "For one distribution, one measure, a test-set size $n$, the measure's target and a forecast at another "
           "quantile level: the probability that a test set of $n$ outcomes gives the target the strictly lower score "
           "(ties count as disagreement); 0.5 is a coin toss.",
    "D21": "Inverse Gaussian outcomes with mean 100 and shape 1000, 200, 100, 50, 20 or 10 (coefficient of variation 0.32 "
           "to 3.16): positive, skewed right, every moment and $\\mathrm E[1/y]$ finite.",
    "D22": "$y = 100\\,e^x$ with $x = \\pm\\sigma(G-\\alpha)/\\sqrt\\alpha$, $G \\sim \\mathrm{Gamma}(\\alpha, 1)$: $\\ln y$ with "
           "standard deviation $\\sigma$ (0.5 or 1; sd in the settings' names) and skewness $\\pm 2/\\sqrt\\alpha$ (-2, -1, "
           "-0.5, +0.5, +1, +2), and the lognormal of the same $\\sigma$ (skewness 0). Some have infinite moments: "
           "$\\mathrm E[y]$ infinite ($\\sigma = 1$, skewness +2), $\\mathrm E[y^2]$ infinite ($\\sigma = 0.5$, skewness +2; "
           "$\\sigma = 1$, skewness +1), $\\mathrm E[1/y]$ infinite ($\\sigma = 1$, skewness -2).",
    "D23": "$y = 100\\,e^x$, $x$ normal with mean $-1$ (probability $p$) or $+1$, standard deviation 0.25; $p = 0.5$ "
           "(equal modes), 0.65 (left mode higher), 0.35 (right mode higher).",
    "D24": "Normal outcomes with standard deviation 1 and mean 0, 0.5 or 1: 50%, 31% and 16% negative; "
           "$\\mathrm E[1/|y|]$ infinite.",
    "D25": "Lognormal outcomes (median 100, $\\sigma$ 0.5 or 1) replaced by 0 with probability $p_0$: 0.25, 0.42, 0.47, 0.6 "
           "for $\\sigma = 0.5$ and 0.25, 0.37, 0.42, 0.6 for $\\sigma = 1$, on both sides of the zero threshold of H8 "
           "(0.4469 and 0.3940).",
    "D26": "Lognormal outcomes with $\\sigma$, and log-Laplace outcomes with scale $b$, each 0.25, 0.5, 1, 1.5, 2, 3, 4 or 6, "
           "median 100: $\\ln y$ symmetric and strictly unimodal, so that sMAPE and MAE both target the median (H5).",
    "D27": "$\\ln y = \\ln 100 + \\sigma\\varepsilon$ with $\\varepsilon$ standard normal (lognormal), standard Gumbel "
           "(log-Gumbel right, $\\ln y$ skewed right) or minus standard Gumbel (log-Gumbel left), for $\\sigma$ = 0.1, 0.25, "
           "0.5, 1, 1.5, 2, 3, 4, 6 and 10.",
}

# definitions whose text the record changed or added for version 3: the record's wording, with the paper's notation
DEF_EDITS = {
    "D28": [("to the Suilin loss: $\\tilde s(y,f) = |y-f|/\\max(|y|+|f|+\\delta,\\ 0.5+\\delta)$ with $\\delta = 0.1$; "
             "the Suilin loss is",
             "to the Suilin loss:\n$$\\tilde s(y,f) = \\frac{|y-f|}{\\max(|y|+|f|+\\delta,\\ 0.5+\\delta)}, \\qquad "
             "\\delta = 0.1;$$\nthe Suilin loss is"),
            ("calls it the Suilin loss.", "calls it the Suilin loss. Version 2 of this paper called it msMAPE.")],
    "D33": [("(Kolassa writes it $m$, which is the median here): $p(y,f) = |f-y|/(f+\\max\\{y, \\nu\\})$, which is 0",
             "(Kolassa writes it $m$, which is the median here):\n$$s_\\nu(y,f) = \\frac{|f-y|}{f+\\max\\{y, \\nu\\}},$$\n"
             "which is 0")],
}

# the record's variable names, as the tables and results show them (the CSV files keep the record's own names)
LABELS = [("(f-g)(f+g-2EY)", "(f1-f2)(f1+f2-2E[y])"), ("800 ln2 g)", "800 ln2 h)"),
          ("between m(0, f) and", "between s~(0, f) and"), ("between m and", "between s~ and"), ("+eps", "+delta"),
          ("(0.5+eps)", "(0.5+delta)")]
NOTES = [("lists the sMAPE term, (f-g)/[(f+g)/2], among", "lists the sMAPE term, |y-f|/[(|y|+|f|)/2], among")]


def _labels(text):
    for a, b in LABELS:
        text = text.replace(a, b)
    return text


XREF = {"D14": [("rank first (H12)", "rank first (C8)")],
        "D18": [("large-spread approximation of H14", "large-spread approximation of C10")],
        "D19": [("small-error limit (H9)", "small-error limit (C3)")],
        "D25": [("zero threshold of H8", "zero threshold of C7")],
        "D26": [("both target the median (H5)", "both target the median (C5)")]}


# words the record's checkers wrote stay as they wrote them
VERBATIM = {"checks", "exploration"}


def apply(d):
    claims = {c["id"]: c for c in d["claims"]}
    for cid, pairs in CLAIMS.items():
        claims[cid]["statement"] = _rep(claims[cid]["statement"], pairs, cid)
    defs = {x["id"]: x for x in d["definitions"]}
    for did, text in DEFS.items():
        defs[did]["definition"] = text
    for did, pairs in DEF_EDITS.items():
        defs[did]["definition"] = _rep(defs[did]["definition"], pairs, did)
    # terms the record marks by its own term IDs ([[TM1:sMAPE term]]) are marked by their definitions
    by_term = {x["term"]: x["id"] for x in d["definitions"]}
    for x in d["definitions"]:
        x["definition"] = re.sub(r"\[\[TM\d+:([^\]]+)\]\]", lambda m: f"[[{by_term[m.group(1)]}:{m.group(1)}]]",
                                 x["definition"])
    # the record's hypothesis numbers, which a reader of the paper cannot look up, as the claims they became
    for did, pairs in XREF.items():
        defs[did]["definition"] = _rep(defs[did]["definition"], pairs, did)
    for e in d["evidence"]:
        if e.get("description", "").find("SMAPE") >= 0:
            e["description"] = e["description"].replace("SMAPE", "sMAPE")

    def walk(x, path):
        if isinstance(x, dict):
            return {k: (v if k in VERBATIM else walk(v, path + [k])) for k, v in x.items()}
        if isinstance(x, list):
            return [walk(v, path) for v in x]
        if isinstance(x, str):
            return SYM_Y.sub("y", x) if ("$" in x or "ln Y" in x or "E[" in x) else x
        return x
    for e in d["evidence"]:
        for r in e.get("results") or []:
            r["variable"] = _labels(r["variable"])
        for t in e.get("tables") or []:
            t["columns"] = [_labels(c) if isinstance(c, str) else ({**c, "variable": _labels(c["variable"])} if isinstance(c, dict) and isinstance(c.get("variable"), str) else c) for c in t.get("columns") or []]
            if t.get("description"):
                t["description"] = _labels(t["description"])
    for c in d["claims"]:
        for r in c.get("relations") or []:
            for a, b in NOTES:
                if a in r.get("note", ""):
                    r["note"] = r["note"].replace(a, b)
    for key in ("claims", "definitions", "evidence", "limitations", "warrants", "open_questions"):
        if key in d:
            d[key] = walk(d[key], [key])
    return d


def leftovers(d):
    """Symbols of the old notation that remain outside the checkers' own words."""
    bad = []
    pats = [SYM_Y, re.compile(r"m_\+"), re.compile(r"e\^\{?m"), re.compile(r"\$m\(y"), re.compile(r"\[\[TM\d+"),
            re.compile(r"(?<![A-Za-z\\_])p\(y, ?f\)")]

    def walk(x, path):
        if isinstance(x, dict):
            for k, v in x.items():
                if k not in VERBATIM:
                    walk(v, path + "." + k)
        elif isinstance(x, list):
            for i, v in enumerate(x):
                walk(v, f"{path}[{i}]")
        elif isinstance(x, str):
            for p in pats:
                if p.search(x) and not (p is SYM_Y and not ("$" in x or "ln Y" in x or "E[" in x)):
                    bad.append((path, p.pattern, x[max(0, p.search(x).start() - 40):p.search(x).end() + 40]))
            if "\\epsilon" in x and "write the constant" not in x:
                bad.append((path, "epsilon", x[:80]))
    walk(d, "")
    # the msMAPE term is s_nu: no bare p where the paper speaks of it (p_0, the share of zeros, stays)
    bare_p = re.compile(r"(?<![A-Za-z\\_{^])p(?![A-Za-z_0-9{])")
    math = re.compile(r"\$\$.+?\$\$|\$[^$]+\$", re.S)
    claims = {c["id"]: c for c in d["claims"]}
    proofs = {sp["evidence"] for sp in claims["C20"]["support"] if sp["kind"] == "formal_proof"}
    texts = [("claims.C20.statement", claims["C20"]["statement"])]
    texts += [(f"definitions.{x['id']}", x["definition"]) for x in d["definitions"]
              if x["term"] in ("msMAPE term", "msMAPE target")]
    texts += [(f"evidence.{e['id']}.{st['id']}", st["statement"]) for e in d["evidence"] if e["id"] in proofs
              for st in e["steps"]]
    for where, text in texts:
        for m in math.finditer(text):
            if bare_p.search(m.group(0)):
                bad.append((where, "bare p", m.group(0)[:80]))
    return bad
