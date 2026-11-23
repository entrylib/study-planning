"""Write every division a/b inside a formula ($...$ or $$...$$) as \\frac{a}{b}, as a paper sets it.

The numerator is the whole term before the slash, back to an operator or relation at the same depth (a product divided
by something is the same either way); the denominator is the one factor after it: a group, a number, or a letter or
command with its arguments, sub- and superscripts and an attached argument list, as in h(0) or M(2\\sigma)^2. Outer
parentheses that only delimited an operand are dropped. Text outside formulas is left as it is."""
import re

OPEN, CLOSE = "([{", ")]}"
PAIR = {")": "(", "]": "[", "}": "{"}
STOP_CHARS = set("+-=<>,;&:")
STOP_CMDS = {"le", "ge", "leq", "geq", "ne", "neq", "approx", "to", "in", "sim", "pm", "mp", "mid", "quad", "qquad",
             "ll", "gg", "le", "propto", "equiv", "rightarrow", "Rightarrow", "iff", "implies", "times", "text",
             "colon", "cdots", "ldots", "dots"}
SPACE_CMDS = {",", ";", "!", ":", " "}
ARG_CMDS = {"mathrm", "operatorname", "mathbf", "mathit", "mathcal", "mathbb", "boldsymbol", "hat", "bar", "tilde",
            "widehat", "overline", "dot", "vec"}


def _back_group(s, j):
    """Index of the opening of the group that closes at s[j] (a bracket or an absolute-value bar)."""
    c = s[j]
    if c == "|":
        k = j - 1
        depth = 0
        while k >= 0:
            if s[k] in CLOSE:
                depth += 1
            elif s[k] in OPEN:
                depth -= 1
            elif s[k] == "|" and depth == 0:
                return k
            k -= 1
        return None
    depth = 0
    k = j
    while k >= 0:
        if s[k] in CLOSE and (k == 0 or s[k - 1] != "\\"):
            depth += 1
        elif s[k] in OPEN and (k == 0 or s[k - 1] != "\\"):
            depth -= 1
            if depth == 0:
                return k
        k -= 1
    return None


def _fwd_group(s, i):
    c = s[i]
    if c == "|":
        k = i + 1
        depth = 0
        while k < len(s):
            if s[k] in OPEN:
                depth += 1
            elif s[k] in CLOSE:
                depth -= 1
            elif s[k] == "|" and depth == 0:
                return k
            k += 1
        return None
    close = {"(": ")", "[": "]", "{": "}"}[c]
    depth = 0
    k = i
    while k < len(s):
        if s[k] == c and (k == 0 or s[k - 1] != "\\"):
            depth += 1
        elif s[k] == close and s[k - 1] != "\\":
            depth -= 1
            if depth == 0:
                return k
        k += 1
    return None


def _numerator_start(s, i):
    """Start of the term that ends just before s[i] (the slash)."""
    j = i - 1
    start = i
    while j >= 0:
        c = s[j]
        if c in CLOSE or (c == "|" and _back_group(s, j) is not None and _is_closing_bar(s, j)):
            k = _back_group(s, j)
            if k is None:
                return None
            start, j = k, k - 1
            continue
        if c in OPEN:                       # the enclosing group opens here
            break
        if c in STOP_CHARS:
            break
        if c.isspace():
            # a space continues the term only after a command (\mathrm E) or before more of it
            k = j
            while k >= 0 and s[k].isspace():
                k -= 1
            m = re.search(r"\\([A-Za-z]+)$", s[:k + 1])
            if m and m.group(1) in STOP_CMDS:
                break
            if k < 0 or s[k] in STOP_CHARS or s[k] in OPEN:
                break
            j = k
            continue
        if c.isalpha():
            m = re.search(r"\\([A-Za-z]+)$", s[:j + 1])
            if m:
                if m.group(1) in STOP_CMDS:
                    break
                start, j = m.start(), m.start() - 1
                continue
        if c == "\\":                        # a spacing command such as \, (the comma is handled above)
            start, j = j, j - 1
            continue
        if j > 0 and s[j - 1] == "\\" and c in SPACE_CMDS:
            start, j = j - 1, j - 2
            continue
        start, j = j, j - 1
    return start if start < i else None


def _is_closing_bar(s, j):
    """Whether the bar at s[j] closes an absolute value (an even number of bars before it at this depth)."""
    depth = 0
    bars = 0
    k = j - 1
    while k >= 0:
        if s[k] in CLOSE:
            depth += 1
        elif s[k] in OPEN:
            if depth == 0:
                break
            depth -= 1
        elif s[k] == "|" and depth == 0:
            bars += 1
        k -= 1
    return bars % 2 == 1


def _denominator_end(s, i):
    """End (exclusive) of the one factor that starts at s[i] (after the slash and any spaces)."""
    n = len(s)
    k = i
    if k >= n:
        return None
    if s[k] in "([{" or s[k] == "|":
        e = _fwd_group(s, k)
        if e is None:
            return None
        k = e + 1
    elif s[k].isdigit() or s[k] == ".":
        while k < n and (s[k].isdigit() or s[k] == "."):
            k += 1
        return k                                     # a number takes nothing attached: 1/2x is (1/2)x
    elif s[k] == "\\":
        m = re.match(r"\\([A-Za-z]+)", s[k:])
        if not m:
            return None
        name = m.group(1)
        k += m.end()
        if name in ARG_CMDS:
            if k < n and s[k] == "{":
                k = _fwd_group(s, k) + 1
            elif k < n and s[k] == " " and k + 1 < n and s[k + 1].isalpha():
                k += 2
        elif name == "frac":
            for _ in range(2):
                if k < n and s[k] == "{":
                    k = _fwd_group(s, k) + 1
                elif k < n:
                    k += 1
        elif name in ("ln", "log", "exp", "sqrt", "sin", "cos", "tanh", "cosh", "sinh"):
            if k < n and s[k] == "{":
                k = _fwd_group(s, k) + 1
            elif k < n and s[k] == " " and k + 1 < n and (s[k + 1].isalnum()):
                k += 1
                while k < n and s[k].isalnum():
                    k += 1
                return k
    elif s[k].isalpha():
        k += 1
    else:
        return None
    while k < n:                                     # sub/superscripts, primes, an attached argument list
        if s[k] in "^_":
            k += 1
            if k < n and s[k] == "{":
                k = _fwd_group(s, k) + 1
            elif k < n and s[k] == "\\":
                m = re.match(r"\\[A-Za-z]+", s[k:])
                k += m.end() if m else 1
            else:
                k += 1
        elif s[k] == "'":
            k += 1
        elif s[k] in "([" or (s[k] == "|" and k > i):
            e = _fwd_group(s, k)
            if e is None:
                break
            k = e + 1
        else:
            break
    return k


def _strip(x):
    x = x.strip()
    for lo, hi in (("\\big(", "\\big)"), ("\\Big(", "\\Big)"), ("\\left(", "\\right)")):
        if x.startswith(lo) and x.endswith(hi) and x.count(lo) == 1:
            return x[len(lo):-len(hi)].strip()
    if len(x) >= 2 and x[0] == "(" and _fwd_group(x, 0) == len(x) - 1:
        return x[1:-1].strip()
    return x


def _in_script(m, i):
    """Whether position i is inside a ^{...} or _{...} group, where a fraction would be unreadable."""
    depth_starts = []
    for k, c in enumerate(m[:i]):
        if c == "{" and (k == 0 or m[k - 1] != "\\"):
            depth_starts.append(k)
        elif c == "}" and (k == 0 or m[k - 1] != "\\") and depth_starts:
            depth_starts.pop()
    return any(k > 0 and m[k - 1] in "^_" for k in depth_starts)


def convert_math(m):
    """One formula's body, every slash outside sub- and superscripts turned into \\frac; None if one cannot be
    parsed."""
    guard = 0
    while guard < 50:
        guard += 1
        slashes = [k for k, c in enumerate(m) if c == "/" and not _in_script(m, k)]
        if not slashes:
            return m
        i = slashes[0]
        a = _numerator_start(m, i)
        j = i + 1
        while j < len(m) and m[j] == " ":
            j += 1
        b = _denominator_end(m, j)
        if a is None or b is None:
            return None
        num, den = _strip(m[a:i]), _strip(m[j:b])
        if not num or not den:
            return None
        m = m[:a] + "\\frac{" + num + "}{" + den + "}" + m[b:]
    return None


MATH = re.compile(r"(\$\$.+?\$\$|\$[^$]+?\$)", re.S)


def convert_text(text, failures=None):
    def one(mo):
        whole = mo.group(0)
        d = 2 if whole.startswith("$$") else 1
        body = whole[d:-d]
        if "/" not in body:
            return whole
        out = convert_math(body)
        if out is None:
            if failures is not None:
                failures.append(whole)
            return whole
        return whole[:d] + out + whole[-d:]
    return MATH.sub(one, text)


def convert_doc(doc, failures=None, skip_keys=("uri", "version", "id", "target", "doi", "url")):
    if isinstance(doc, dict):
        return {k: (v if k in skip_keys else convert_doc(v, failures, skip_keys)) for k, v in doc.items()}
    if isinstance(doc, list):
        return [convert_doc(v, failures, skip_keys) for v in doc]
    if isinstance(doc, str) and "$" in doc and "/" in doc:
        return convert_text(doc, failures)
    return doc
