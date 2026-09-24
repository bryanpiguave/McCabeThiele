"""Vapor-liquid equilibrium for binary systems.

Each VLE object stores the equilibrium curve on a composition grid (light
component mole fractions) and interpolates on it. It can come from
  * a thermodynamic model: modified Raoult's law, y_i P = x_i gamma_i Psat_i(T),
    with an ideal, Margules, Van Laar, Wilson or NRTL activity model,
  * a constant relative volatility, or
  * a table (CSV with columns X, Y and optionally T in degC).
"""

import os
import numpy as np
import pandas as pd

from .components import ATM_KPA, get_component, binary_parameters

# Grid clustered near x = 0 and x = 1 so high-purity specs are resolved
GRID = 0.5 * (1 - np.cos(np.pi * np.linspace(0, 1, 1601)))


# ---------------------------------------------------------------- activity models
def gamma_ideal(x1, T, p):
    one = np.ones_like(x1)
    return one, one


def gamma_margules(x1, T, p):
    A12, A21 = p["A12"], p["A21"]
    x2 = 1 - x1
    return (np.exp(x2**2 * (A12 + 2 * (A21 - A12) * x1)),
            np.exp(x1**2 * (A21 + 2 * (A12 - A21) * x2)))


def gamma_van_laar(x1, T, p):
    A12, A21 = p["A12"], p["A21"]
    x2 = 1 - x1
    den = A12 * x1 + A21 * x2
    return (np.exp(A12 * (A21 * x2 / den) ** 2),
            np.exp(A21 * (A12 * x1 / den) ** 2))


def gamma_wilson(x1, T, p):
    L12, L21 = p["Lambda12"], p["Lambda21"]
    x2 = 1 - x1
    t = L12 / (x1 + L12 * x2) - L21 / (x2 + L21 * x1)
    return (np.exp(-np.log(x1 + L12 * x2) + x2 * t),
            np.exp(-np.log(x2 + L21 * x1) - x1 * t))


def gamma_nrtl(x1, T, p):
    """tau_ij = a_ij + b_ij / T, with b_ij in K."""
    x2 = 1 - x1
    t12 = p.get("a12", 0.0) + p.get("b12", 0.0) / T
    t21 = p.get("a21", 0.0) + p.get("b21", 0.0) / T
    alpha = p.get("alpha", 0.3)
    G12, G21 = np.exp(-alpha * t12), np.exp(-alpha * t21)
    ln1 = x2**2 * (t21 * (G21 / (x1 + x2 * G21)) ** 2 + t12 * G12 / (x2 + x1 * G12) ** 2)
    ln2 = x1**2 * (t12 * (G12 / (x2 + x1 * G12)) ** 2 + t21 * G21 / (x1 + x2 * G21) ** 2)
    return np.exp(ln1), np.exp(ln2)


ACTIVITY_MODELS = {"raoult": gamma_ideal, "ideal": gamma_ideal, "margules": gamma_margules,
                   "van_laar": gamma_van_laar, "wilson": gamma_wilson, "nrtl": gamma_nrtl}


# ---------------------------------------------------------------- VLE container
class VLE:
    """Binary equilibrium curve y*(x), optionally with bubble temperatures."""

    def __init__(self, x, y, T=None, components=None, P=ATM_KPA, label=""):
        self.x = np.asarray(x, float)
        self.y = np.asarray(y, float)
        self.T = None if T is None else np.asarray(T, float)   # K
        self.components = components                          # (light, heavy) Component or None
        self.P = P                                             # kPa
        self.label = label
        # y must increase with x for the inverse x*(y); stop at an azeotrope
        dy = np.diff(self.y)
        stop = np.argmax(dy <= 0) if np.any(dy <= 0) else len(dy)
        self._xi, self._yi = self.x[:stop + 1], self.y[:stop + 1]

    @property
    def has_temperature(self):
        return self.T is not None

    @property
    def has_properties(self):
        return self.components is not None

    def y_eq(self, x):
        return np.interp(x, self.x, self.y)

    def x_eq(self, y):
        return np.interp(y, self._yi, self._xi)

    def T_bubble(self, x):
        """Bubble-point temperature in K."""
        if self.T is None:
            raise ValueError("This VLE has no temperature data")
        return np.interp(x, self.x, self.T)

    def T_dew(self, y):
        """Dew-point temperature in K."""
        return self.T_bubble(self.x_eq(y))

    def alpha(self, x):
        x = np.clip(x, 1e-9, 1 - 1e-9)
        y = np.clip(self.y_eq(x), 1e-12, 1 - 1e-12)
        return y * (1 - x) / (x * (1 - y))

    def azeotropes(self):
        """Interior points where y = x, as a list of (x, T in K or None)."""
        d = self.y - self.x
        found = []
        for i in range(1, len(d) - 2):
            if d[i] == 0 or d[i] * d[i + 1] < 0:
                xa = self.x[i] - d[i] * (self.x[i + 1] - self.x[i]) / (d[i + 1] - d[i])
                found.append((xa, self.T_bubble(xa) if self.has_temperature else None))
        return found

    # ------------------------------------------------ mixture properties
    def _mix(self, x, prop):
        c1, c2 = self._require_properties()
        return x * prop(c1) + (1 - x) * prop(c2)

    def _require_properties(self):
        if self.components is None:
            raise ValueError("Component properties are needed; set system.components in the config")
        return self.components

    def T_or_estimate(self, x):
        """Bubble temperature (K), or a mole-fraction average of the boiling points."""
        if self.has_temperature:
            return self.T_bubble(x)
        c1, c2 = self._require_properties()
        return x * c1.tsat(self.P) + (1 - x) * c2.tsat(self.P)

    def latent_heat(self, x, T=None):
        """Mixture latent heat in kJ/mol."""
        c1, c2 = self._require_properties()
        T = self.T_or_estimate(x) if T is None else T
        return x * c1.dHvap(T) + (1 - x) * c2.dHvap(T)

    def mw(self, x):
        return self._mix(x, lambda c: c.mw)

    def cp_liq(self, x):
        """kJ/(mol K)"""
        return self._mix(x, lambda c: c.cp_liq) / 1000

    def cp_vap(self, x):
        """kJ/(mol K)"""
        return self._mix(x, lambda c: c.cp_vap) / 1000

    def rho_liq(self, x):
        """kg/m3 (ideal mixing by volume)"""
        c1, c2 = self._require_properties()
        w1 = x * c1.mw / (x * c1.mw + (1 - x) * c2.mw)
        return 1 / (w1 / c1.rho_liq + (1 - w1) / c2.rho_liq)

    def rho_vap(self, y, T):
        """kg/m3 (ideal gas)"""
        return self.P * self.mw(y) / (8.314 * T)

    def viscosity(self, x, T):
        """Liquid viscosity in cP (log-mole-fraction mixing)."""
        c1, c2 = self._require_properties()
        return np.exp(x * np.log(c1.viscosity(T)) + (1 - x) * np.log(c2.viscosity(T)))


# ---------------------------------------------------------------- constructors
def model_vle(light, heavy, P=ATM_KPA, model="raoult", params=None):
    """Bubble-point curve from modified Raoult's law, solved by vectorized bisection."""
    c1, c2 = get_component(light), get_component(heavy)
    model = model.lower()
    if model not in ACTIVITY_MODELS:
        raise KeyError("Unknown VLE model '%s'. Available: %s" % (model, ", ".join(ACTIVITY_MODELS)))
    gamma = ACTIVITY_MODELS[model]
    if model not in ("raoult", "ideal") and not params:
        params = binary_parameters(light, heavy, model)

    x1 = GRID.copy()
    lo = np.full_like(x1, min(c1.tsat(P), c2.tsat(P)) - 80)
    hi = np.full_like(x1, max(c1.tsat(P), c2.tsat(P)) + 80)

    def excess(T):
        g1, g2 = gamma(x1, T, params)
        return x1 * g1 * c1.psat(T) + (1 - x1) * g2 * c2.psat(T) - P

    for _ in range(60):
        mid = 0.5 * (lo + hi)
        above = excess(mid) > 0
        hi = np.where(above, mid, hi)
        lo = np.where(above, lo, mid)
    T = 0.5 * (lo + hi)
    g1, _ = gamma(x1, T, params)
    y1 = x1 * g1 * c1.psat(T) / P
    label = "%s-%s, %s, %.1f kPa" % (c1.name, c2.name, model, P)
    return VLE(x1, y1, T, (c1, c2), P, label)


def constant_alpha_vle(alpha, components=None, P=ATM_KPA):
    x = GRID.copy()
    y = alpha * x / (1 + (alpha - 1) * x)
    comps = tuple(get_component(c) for c in components) if components else None
    return VLE(x, y, None, comps, P, "constant alpha = %.3g" % alpha)


def table_vle(path, components=None, P=ATM_KPA):
    data = pd.read_csv(path).sort_values("X")
    T = data["T"].values + 273.15 if "T" in data.columns else None
    comps = tuple(get_component(c) for c in components) if components else None
    return VLE(data["X"].values, data["Y"].values, T, comps, P, "table %s" % os.path.basename(path))


def make_vle(system):
    """Build a VLE from the `system` section of a config dict."""
    comps = system.get("components")
    P = float(system.get("pressure_kPa", ATM_KPA))
    vle = system.get("vle", {}) or {}
    model = vle.get("model", "raoult").lower()
    if model == "table":
        return table_vle(vle["file"], comps, P)
    if model == "constant_alpha":
        return constant_alpha_vle(float(vle["alpha"]), comps, P)
    if not comps or len(comps) != 2:
        raise ValueError("system.components must list two components (light first) for model '%s'" % model)
    return model_vle(comps[0], comps[1], P, model, vle.get("params"))
