"""Preliminary column sizing and cost estimate.

Diameter from the Souders-Brown flooding velocity, height from the tray count,
installed costs from Guthrie's correlations as given by Douglas (1988),
scaled with the Marshall & Swift index. Accuracy is about +/-30 %, which is
enough to compare design options and locate the optimum reflux ratio.
"""

from dataclasses import dataclass
import math

M_TO_FT = 3.2808
M2_TO_FT2 = 10.764


@dataclass
class Economics:
    ms_index: float = 2250.0           # Marshall & Swift index (M&S 2011 escalated to ~2024 with CEPCI)
    tray_spacing_m: float = 0.61
    extra_height_m: float = 4.0        # disengagement space and sump
    flooding_fraction: float = 0.8
    souders_brown_C: float = 0.09      # m/s, for about 0.6 m tray spacing
    downcomer_fraction: float = 0.12
    material_factor: float = 1.0       # carbon steel = 1
    U_condenser: float = 0.8           # kW/(m2 K)
    U_reboiler: float = 1.0            # kW/(m2 K)
    steam_T_C: float = 160.0
    cw_in_C: float = 30.0
    cw_out_C: float = 45.0
    steam_cost_per_GJ: float = 12.0
    cw_cost_per_GJ: float = 0.5
    hours_per_year: float = 8000.0
    capital_charge: float = 1 / 3      # 1/yr: annualized fraction of installed cost


def _diameter(vle, V, x_liq, y_vap, T, e):
    rho_V = vle.rho_vap(y_vap, T)
    rho_L = vle.rho_liq(x_liq)
    u_flood = e.souders_brown_C * math.sqrt((rho_L - rho_V) / rho_V)
    Q = V / 3600 * vle.mw(y_vap) / rho_V                     # m3/s
    area = Q / (e.flooding_fraction * u_flood) / (1 - e.downcomer_fraction)
    return math.sqrt(4 * area / math.pi)


def _exchanger_cost(area_m2, e):
    return e.ms_index / 280 * 101.3 * (area_m2 * M2_TO_FT2) ** 0.65 * (2.29 + e.material_factor)


def costs(design, e=None):
    """Sizes and costs of a designed column. Returns a dict (USD, m, kW)."""
    e = e or Economics()
    vle, m = design.vle, design.model
    if not vle.has_properties:
        raise ValueError("Costing needs component properties; set system.components in the config")
    T_top, T_bot = float(vle.T_or_estimate(m.x_D)), float(vle.T_or_estimate(m.x_W))
    diam = max(_diameter(vle, m.sections[0].V, m.x_D, m.x_D, T_top, e),
               _diameter(vle, m.sections[-1].V, m.x_W, float(vle.y_eq(m.x_W)), T_bot, e))
    trays = design.actual_trays
    height = trays * e.tray_spacing_m + e.extra_height_m
    D_ft, H_ft = diam * M_TO_FT, height * M_TO_FT
    shell = e.ms_index / 280 * 101.9 * D_ft ** 1.066 * H_ft ** 0.802 * (2.18 + e.material_factor)
    internals = e.ms_index / 280 * 4.7 * D_ft ** 1.55 * H_ft * e.material_factor

    Q_C, Q_R = design.Q_C or 0.0, design.Q_R or 0.0
    cond = reb = 0.0
    A_C = A_R = 0.0
    if Q_C > 0:
        dT1, dT2 = T_top - 273.15 - e.cw_in_C, T_top - 273.15 - e.cw_out_C
        if dT2 <= 0:
            raise ValueError("Condenser temperature %.1f degC is below the cooling-water outlet"
                             % (T_top - 273.15))
        lmtd = (dT1 - dT2) / math.log(dT1 / dT2)
        A_C = Q_C / (e.U_condenser * lmtd)
        cond = _exchanger_cost(A_C, e)
    if Q_R > 0:
        dT = e.steam_T_C - (T_bot - 273.15)
        if dT <= 0:
            raise ValueError("Steam at %.0f degC can't boil the bottoms at %.1f degC"
                             % (e.steam_T_C, T_bot - 273.15))
        A_R = Q_R / (e.U_reboiler * dT)
        reb = _exchanger_cost(A_R, e)

    # Direct (live) steam is charged at its latent heat
    steam_kW = Q_R + m.steam * 40.65 / 3.6
    GJ_per_kW_year = 3600 * e.hours_per_year / 1e6
    steam_cost = steam_kW * GJ_per_kW_year * e.steam_cost_per_GJ
    cw_cost = Q_C * GJ_per_kW_year * e.cw_cost_per_GJ
    capital = shell + internals + cond + reb
    return {
        "diameter_m": diam, "height_m": height, "actual_trays": trays,
        "condenser_area_m2": A_C, "reboiler_area_m2": A_R,
        "shell_cost": shell, "tray_cost": internals, "condenser_cost": cond, "reboiler_cost": reb,
        "capital_cost": capital, "annualized_capital": capital * e.capital_charge,
        "steam_cost": steam_cost, "cooling_water_cost": cw_cost,
        "operating_cost": steam_cost + cw_cost,
        "total_annual_cost": capital * e.capital_charge + steam_cost + cw_cost,
    }
