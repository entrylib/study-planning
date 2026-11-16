# Code of the study measure-rankings

Which fixed candidate forecast does each error measure (RMSE, MAE, MAPE, sMAPE) rank first, by its expected score and
in test sets of a given size? The candidates are fixed from the true outcome distribution before any test set is drawn
(its quantiles at levels 0.01 to 0.99, its mean, and the minimisers of the expected MAPE and sMAPE); no forecast is
estimated from data. The study's lab record is on the branch `claimset/measure-rankings`.

- `settings.py`: the outcome distributions (positive continuous ones, with and without outliers; counts).
- `reference.py`: the candidates and their expected scores (one large reference sample, or exact sums for counts).
- `measures.py`: the four measures, scoring many fixed candidates on one test set at once.
- `simulate.py`: the simulation; writes `results.json`, `candidates.csv`, `selection_shares.csv` and figures.
- `plot.py`: the figures. `test_rankings.py`: checks against exact results.

Requirements: `pip install -r requirements.txt` (Python 3.10+).
