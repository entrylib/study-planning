"""The lab records only numbers as result values. Booleans become 1 (true) or 0 (false); a value that does not exist
("none": a measure that targets nothing, an expected score that is infinite) becomes a result named
"<variable> is undefined" with value 1, keeping its conditions (and the reason, where one is given)."""


def numeric(results):
    out = []
    for r in results:
        v = r["value"]
        if isinstance(v, bool):
            out.append({**r, "value": int(v)})
        elif v == "none":
            out.append({**r, "variable": f"{r['variable']} is undefined", "value": 1})
        else:
            out.append(r)
    return out
