"""Ponchon-Savarit method: stage stepping on the enthalpy-composition diagram.

Unlike McCabe-Thiele, it does not assume constant molar overflow, so the
difference in latent heats and the sensible heat are accounted for.

Enthalpy model (reference: pure liquids at 25 degC, heat of mixing neglected)
  saturated liquid  h(x) = cp_L(x) (T_bubble - T_ref)
  saturated vapor   H(y) = cp_L(y) (T_dew - T_ref) + sum y_i dHvap_i(T_dew)

Supports one feed, a total condenser and a partial reboiler.
"""

from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq

from .column import KW_PER_KMOLH_KJMOL, PinchError

T_REF = 298.15


@dataclass
class PonchonResult:
    R_min: float
    R: float
    N: float
    feed_stage: int
    Q_C: float                 # kW
    Q_R: float                 # kW
    delta_D: tuple             # (x, h) difference points, kJ/mol
    delta_W: tuple
    feed_point: tuple
    stages: list               # (stage, x_n, y_n)
    liquid_curve: tuple        # (x, h)
    vapor_curve: tuple         # (y, H)
    L_over_V_top: float
    L_over_V_bottom: float


class EnthalpyDiagram:
    def __init__(self, vle):
        if not (vle.has_temperature and vle.has_properties):
            raise ValueError("Ponchon-Savarit needs a thermodynamic VLE model with component properties")
        self.vle = vle
        x, y, T = vle.x, vle.y, vle.T
        self.xg, self.hg = x, vle.cp_liq(x) * (T - T_REF)
        Hg = vle.cp_liq(y) * (T - T_REF) + vle.latent_heat(y, T)
        mono = slice(0, len(vle._yi))
        self.yg, self.Hg = y[mono], Hg[mono]

    def h(self, x):
        return np.interp(x, self.xg, self.hg)

    def H(self, y):
        return np.interp(y, self.yg, self.Hg)

    def tie_extension(self, x, X):
        """Enthalpy where the tie line from liquid x, extended, reaches composition X."""
        y = self.vle.y_eq(x)
        return self.h(x) + (self.H(y) - self.h(x)) * (X - x) / (y - x)


def _line(p, q):
    (x1, h1), (x2, h2) = p, q
    return lambda x: h1 + (h2 - h1) * (x - x1) / (x2 - x1)


def ponchon_savarit(vle, x_D, x_W, feed, R=None, reflux_factor=1.5, max_stages=300):
    ed = EnthalpyDiagram(vle)
    F, z, q = feed.flow, feed.z, feed.q
    D = F * (z - x_W) / (x_D - x_W)
    W = F - D
    h_D, h_W = float(ed.h(x_D)), float(ed.h(x_W))
    H_1 = float(ed.H(x_D))
    h_F = float(ed.H(z) - q * (ed.H(z) - ed.h(z)))
    Fp = (z, h_F)

    # Minimum reflux: highest difference point implied by any extended tie line
    x_F_tie = brentq(lambda x: ed.tie_extension(x, z) - h_F, x_W + 1e-9, z - 1e-9) \
        if (ed.tie_extension(x_W + 1e-9, z) - h_F) * (ed.tie_extension(z - 1e-9, z) - h_F) < 0 else z
    xs = np.linspace(x_W + 1e-6, x_D - 1e-6, 2000)
    rect = xs[xs >= x_F_tie]
    strip = xs[xs < x_F_tie]
    cand = list(ed.tie_extension(rect, x_D))
    h_dW = ed.tie_extension(strip, x_W)
    cand += list(h_F + (h_F - h_dW) * (x_D - z) / (z - x_W))
    h_dD_min = max(cand)
    R_min = (h_dD_min - h_D) / (H_1 - h_D) - 1
    if R is None:
        R = reflux_factor * R_min

    # Difference points for the operating reflux
    Qc_per_D = (R + 1) * (H_1 - h_D)
    dD = (x_D, h_D + Qc_per_D)
    dW = (x_W, _line(dD, Fp)(x_W))
    Q_C = D * Qc_per_D * KW_PER_KMOLH_KJMOL
    Q_R = W * (h_W - dW[1]) * KW_PER_KMOLH_KJMOL
    op = _line(dD, dW)
    x_switch = brentq(lambda x: ed.h(x) - op(x), x_W, x_D)

    def next_vapor(x, delta):
        line = _line(delta, (x, float(ed.h(x))))
        y_hi = x_D if delta is dD else 1.0
        ys = np.linspace(x + 1e-9, y_hi, 400)
        g = ed.H(ys) - line(ys)
        idx = np.nonzero(np.sign(g[:-1]) != np.sign(g[1:]))[0]
        if len(idx) == 0:
            raise PinchError("Ponchon-Savarit stepping pinched at x = %.4f" % x)
        i = idx[0]
        return brentq(lambda y: ed.H(y) - line(y), ys[i], ys[i + 1])

    stages, y, x_prev, delta, feed_stage = [], x_D, x_D, dD, None
    for n in range(1, max_stages + 1):
        x = float(vle.x_eq(y))
        stages.append((n, x, y))
        if x <= x_W:
            N = n - 1 + (x_prev - x_W) / (x_prev - x)
            break
        if x >= x_prev - 1e-12:
            raise PinchError("Ponchon-Savarit stepping pinched at x = %.4f" % x)
        if feed_stage is None and x < x_switch:
            feed_stage, delta = n, dW
        y = next_vapor(x, delta)
        x_prev = x
    else:
        raise PinchError("x_W not reached in %d stages" % max_stages)

    # Internal L/V at the top and bottom of the column (not constant here):
    # L1/V2 = (x_D - y2)/(x_D - x1) and L_m/V_(m+1) = (y_(m+1) - x_W)/(x_m - x_W)
    if len(stages) > 1:
        LV_top = (x_D - stages[1][2]) / (x_D - stages[0][1])
        LV_bot = (stages[-1][2] - x_W) / (stages[-2][1] - x_W)
    else:
        LV_top = LV_bot = float("nan")
    return PonchonResult(R_min, R, N, feed_stage, Q_C, Q_R, dD, dW, Fp, stages,
                         (ed.xg, ed.hg), (ed.yg, ed.Hg), LV_top, LV_bot)
