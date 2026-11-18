"""The deciding run for the corrected wording of the signal-to-noise and large-sample efficiency hypotheses: their
formulas did not change, so it repeats check_h39.py and the efficiency checks (H40) of check_v3.py.

    python3 check_corrected.py --seed 8

Writes its results to $CLAIMSET_RESULTS when the lab runs it."""
import argparse
import json
import os
import sys
import tempfile

import check_h39
import check_v3


def collect(module, argv):
    target = os.environ.get("CLAIMSET_RESULTS")
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "results.json")
        os.environ["CLAIMSET_RESULTS"] = path
        old = sys.argv
        sys.argv = [module.__name__ + ".py", *argv]
        try:
            module.main()
        finally:
            sys.argv = old
            if target is None:
                del os.environ["CLAIMSET_RESULTS"]
            else:
                os.environ["CLAIMSET_RESULTS"] = target
        with open(path) as fh:
            return json.load(fh)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=8)
    a = ap.parse_args()
    res = collect(check_h39, [])
    res += [r for r in collect(check_v3, ["--seed", str(a.seed)]) if r["conditions"].get("hypothesis") == "H40"]
    path = os.environ.get("CLAIMSET_RESULTS")
    if path:
        with open(path, "w") as fh:
            json.dump(res, fh)


if __name__ == "__main__":
    main()
