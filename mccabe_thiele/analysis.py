"""Sensitivity studies: stages and cost against reflux, and feed-stage placement."""

from dataclasses import replace
import math
import numpy as np
import pandas as pd

from .column import design, forced_feed_stage_count, PinchError
from .economics import costs


def reflux_sweep(vle, spec, factors=None, econ=None, param_min=None):
    """Design the column over a range of R/R_min (boil-up/minimum for strippers).

    Returns a DataFrame with stage counts, duties and, when component
    properties are available, the cost breakdown.
    """
    if factors is None:
        factors = np.round(np.concatenate([[1.005, 1.01], np.linspace(1.02, 1.5, 17), np.linspace(1.6, 3.0, 8)]), 3)
    stripping = spec.configuration == "stripping"
    rows = []
    for f in factors:
        s = replace(spec, reflux_ratio=None, boilup_ratio=None,
                    **({"boilup_factor": f} if stripping else {"reflux_factor": f}))
        try:
            d = design(vle, s, param_min)
        except PinchError:
            continue
        param_min = d.param_min
        row = {"factor": f, "ratio": d.param, "N_theoretical": d.N,
               "theoretical_trays": d.theoretical_trays, "actual_trays": d.actual_trays,
               "Q_C_kW": d.Q_C, "Q_R_kW": d.Q_R}
        if vle.has_properties:
            row.update(costs(d, econ))
        rows.append(row)
    return pd.DataFrame(rows)


def optimum(sweep):
    """Row of a reflux sweep with the lowest total annual cost."""
    if "total_annual_cost" not in sweep or not len(sweep):
        return None
    i = sweep["total_annual_cost"].idxmin()
    best = sweep.loc[i].copy()
    best["at_bound"] = i in (sweep.index[0], sweep.index[-1])
    return best


def feed_stage_sensitivity(d, extra=6):
    """Theoretical stages needed when the single feed is moved away from the optimum."""
    if len(d.spec.feeds) != 1 or d.spec.side_draws or d.spec.configuration != "full":
        return None
    n_max = int(math.ceil(d.N)) + extra
    rows = [{"feed_stage": k, "N_theoretical": forced_feed_stage_count(d.vle, d.spec, d.param, k)}
            for k in range(1, n_max + 1)]
    return pd.DataFrame(rows)
