# trend_calculator.py

import numpy as np

def calculate_trend(args):
    x, y, annual_mean, time_values_init = args
    chla_values = annual_mean['chla'].loc[dict(x=x, y=y)]
    mask = np.isfinite(chla_values)
    chla_values = chla_values[mask]
    time_values = time_values_init[mask]

    if chla_values.size == 0:
        return x, y, np.nan, np.nan, np.nan

    trend = np.polyfit(time_values, chla_values, 1)
    predicted_values = np.polyval(trend, time_values)
    rmse = np.sqrt(np.mean((chla_values - predicted_values)**2))

    return x, y, trend[0], trend[1], rmse