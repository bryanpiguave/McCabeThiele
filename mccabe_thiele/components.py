"""Pure-component property data and binary interaction parameters.

Values are typical handbook numbers (Perry's, Reid-Prausnitz-Poling, Seader &
Henley) and are meant for preliminary design, not for final equipment sizing.

Antoine equation: log10(P / mmHg) = A - B / (C + T / degC)
"""

from dataclasses import dataclass
import numpy as np

MMHG_TO_KPA = 0.133322
ATM_KPA = 101.325


@dataclass(frozen=True)
class Component:
    name: str
    mw: float            # kg/kmol
    antoine: tuple       # (A, B, C), mmHg and degC
    Tb: float            # normal boiling point, K
    Tc: float            # critical temperature, K
    dHvap_b: float       # latent heat at Tb, kJ/mol
    cp_liq: float        # liquid heat capacity, J/(mol K)
    cp_vap: float        # vapor heat capacity near Tb, J/(mol K)
    rho_liq: float       # liquid density, kg/m3
    mu_25: float         # liquid viscosity at 25 degC, cP
    mu_b: float          # liquid viscosity at Tb, cP

    def psat(self, T):
        """Vapor pressure in kPa at T in K."""
        A, B, C = self.antoine
        return MMHG_TO_KPA * 10 ** (A - B / (C + T - 273.15))

    def tsat(self, P):
        """Saturation temperature in K at P in kPa."""
        A, B, C = self.antoine
        return B / (A - np.log10(P / MMHG_TO_KPA)) - C + 273.15

    def dHvap(self, T):
        """Latent heat in kJ/mol at T in K (Watson correlation)."""
        T = np.minimum(T, self.Tc - 1e-6)
        return self.dHvap_b * ((self.Tc - T) / (self.Tc - self.Tb)) ** 0.38

    def viscosity(self, T):
        """Liquid viscosity in cP at T in K (Andrade fit through two points)."""
        T1, T2 = 298.15, self.Tb
        B = np.log(self.mu_25 / self.mu_b) / (1 / T1 - 1 / T2)
        return self.mu_25 * np.exp(B * (1 / T - 1 / T1))


COMPONENTS = {c.name: c for c in [
    Component("benzene",  78.11, (6.90565, 1211.033, 220.790), 353.2, 562.0, 30.72, 136.0,  96.0, 876.0, 0.604, 0.316),
    Component("toluene",  92.14, (6.95464, 1344.800, 219.482), 383.8, 591.8, 33.18, 157.0, 125.0, 867.0, 0.560, 0.248),
    Component("ethanol",  46.07, (8.20417, 1642.890, 230.300), 351.4, 513.9, 38.56, 112.0,  78.0, 789.0, 1.074, 0.440),
    Component("water",    18.015, (8.07131, 1730.630, 233.426), 373.15, 647.1, 40.65, 75.3,  34.0, 997.0, 0.890, 0.282),
    Component("methanol", 32.04, (8.08097, 1582.271, 239.726), 337.7, 512.6, 35.21,  81.0,  50.0, 792.0, 0.544, 0.330),
    Component("n-hexane", 86.18, (6.87601, 1171.170, 224.408), 341.9, 507.6, 28.85, 195.0, 160.0, 655.0, 0.300, 0.200),
    Component("n-heptane", 100.20, (6.89677, 1264.900, 216.544), 371.6, 540.2, 31.77, 224.0, 185.0, 684.0, 0.387, 0.200),
]}

# Binary interaction parameters, keyed by (light, heavy) component names.
# Van Laar and Margules constants are for ln(gamma), from Seader & Henley.
BINARY_PARAMETERS = {
    ("ethanol", "water"): {"van_laar": {"A12": 1.6798, "A21": 0.9227},
                           "margules": {"A12": 1.6022, "A21": 0.7947}},
    ("methanol", "water"): {"van_laar": {"A12": 0.9011, "A21": 0.5486},
                            "margules": {"A12": 0.7923, "A21": 0.5434}},
}


def get_component(name):
    try:
        return COMPONENTS[name.lower()]
    except KeyError:
        raise KeyError("Unknown component '%s'. Available: %s" % (name, ", ".join(COMPONENTS))) from None


def binary_parameters(light, heavy, model):
    params = BINARY_PARAMETERS.get((light.lower(), heavy.lower()), {}).get(model)
    if params is None:
        raise KeyError("No built-in %s parameters for %s-%s; give them in the config under vle.params"
                       % (model, light, heavy))
    return dict(params)
