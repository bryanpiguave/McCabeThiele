"""Rating an existing column: given the number of stages, the feed stages, the
reflux ratio and the distillate rate, find the product compositions."""

from dataclasses import dataclass, replace
import numpy as np
from scipy.optimize import brentq

from .column import ColumnModel, step_stages


@dataclass
class Rating:
    x_D: float
    x_W: float
    D: float
    W: float
    stepping: object
    model: object


def rate(vle, spec, n_stages, feed_stages, R, D):
    """Product compositions of a full column with a partial reboiler.

    n_stages     equilibrium stages, including the reboiler (and a partial condenser)
    feed_stages  stage of each feed, counted from the top (same order as spec.feeds)
    R            reflux ratio L/D
    D            distillate rate, kmol/h
    """
    if spec.configuration != "full" or spec.reboiler != "partial" or spec.side_draws:
        raise ValueError("Rating supports full columns with a partial reboiler and no side draws")
    feeds = sorted(spec.feeds, key=lambda f: -f.z)
    order = [spec.feeds.index(f) for f in feeds]
    forced = [feed_stages[i] for i in order]
    F = sum(f.flow for f in feeds)
    Fz = sum(f.flow * f.z for f in feeds)
    W = F - D
    if not 0 < D < F:
        raise ValueError("Need 0 < D < total feed")

    def run(x_D):
        x_W = (Fz - D * x_D) / W
        s = replace(spec, feeds=feeds, x_D=x_D, x_W=x_W)
        model = ColumnModel(vle, s, R)
        st = step_stages(vle, model, spec.murphree, forced=forced, n_fixed=n_stages,
                         ideal_last=True, ideal_first=spec.condenser == "partial")
        return st.stages[-1][1] - x_W, x_W, st, model

    lo = max(Fz / D - W / D, 0.0) + 1e-9      # x_W <= 1
    hi = min(Fz / D, 1.0) - 1e-9              # x_W >= 0
    grid = np.linspace(lo, hi, 400)
    res = [run(x)[0] for x in grid]
    root = None
    for a, b, ra, rb in zip(grid, grid[1:], res, res[1:]):
        if ra * rb <= 0:
            root = brentq(lambda x: run(x)[0], a, b, xtol=1e-12)
    if root is None:
        raise ValueError("No consistent product compositions found for this column")
    _, x_W, st, model = run(root)
    return Rating(root, x_W, D, W, st, model)
