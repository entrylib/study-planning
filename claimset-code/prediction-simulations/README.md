# Code of the study prediction-simulations

Simulations testing what point forecasts RMSE, MAE, MAPE and sMAPE favour, for the claims of the guide
`prediction.md` (claimset `entrylib:2026.000008`). The study's lab record is on the branch
`claimset/prediction-simulations`; this code is run from the repository's root under that record
(`claimset lab run --plan …`).

- `measures.py`: the four measures, scoring many candidate forecasts at once.
- `oracle.py`: which candidate each measure selects when outcomes come from a known distribution.
- `models.py`: ETS(A,N,N) on logs (statsforecast), its median and mean forecasts scored by RMSE and MAE.
- `aggregation.py`: summed RMSE and mean RMSSE over series on different scales.

Requirements: `pip install -r requirements.txt` (Python 3.10+).
