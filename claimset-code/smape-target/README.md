# Code of the study smape-target

Which point forecast minimises expected sMAPE, and why; where sMAPE and MAE both target the median, which one's minimiser
varies less; and how expected sMAPE flattens as the spread of the outcomes grows. The study's lab record is on the branch
`claimset/smape-target`; its hypotheses are formal results (proved in the study's proof document) whose predictions these
scripts compute exactly and check on simulated test sets.

- `dists.py`: the outcome distributions, each with a sampler, the density of ln y and which moments are finite.
- `theory.py`: what the derivations predict, computed from each distribution by quadrature (never from a sample):
  expected scores, every local minimum of expected sMAPE, the one-step formula, the zero threshold, the large-spread
  approximation, the large-sample variances.
- `measures.py`: the measures, scoring many fixed forecasts on one test set at once.
- `minimiser.py`: the exact sMAPE minimiser of a test set (always one of its outcomes).
- `targets.py` (plan PL1), `efficiency.py` (PL2), `flatness.py` (PL3): the simulations; `plot_*.py`: their figures.
- `experiment_v4.py`: the numerical experiment of version 4: the formulas of H47 to H50 by quadrature and on simulated test sets of 6 to 1000 actuals, exact small-sample variances where they exist, and the Hypotailedness (entrylib:2026.000016) of each shape.
- `test_smape.py`: checks against exact results (`python3 -m unittest test_smape`).

Requirements: `pip install -r requirements.txt` (Python 3.10+).
