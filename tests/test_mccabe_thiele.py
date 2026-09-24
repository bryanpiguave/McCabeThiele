import math
import os
from dataclasses import replace

import numpy as np
import pytest

from mccabe_thiele import (ColumnSpec, Feed, SideDraw, design, model_vle, constant_alpha_vle, table_vle,
                           rate, ponchon_savarit, reflux_sweep, optimum, costs, VLE)
from mccabe_thiele.column import feed_q, ColumnModel
from mccabe_thiele.shortcut import fenske, underwood
from mccabe_thiele import config as cfgmod
from mccabe_thiele.cli import run

HERE = os.path.dirname(__file__)
EXAMPLES = os.path.join(HERE, "..", "examples")


@pytest.fixture(scope="module")
def bt():
    return model_vle("benzene", "toluene")


@pytest.fixture(scope="module")
def ew():
    return model_vle("ethanol", "water", model="van_laar")


def spec(**kw):
    base = dict(feeds=[Feed(100, 0.6, 1.0)], x_D=0.95, x_W=0.05, reflux_factor=1.5)
    base.update(kw)
    return ColumnSpec(**base)


# ------------------------------------------------------------- thermodynamics
def test_pure_component_boiling_points(bt):
    assert bt.T_bubble(1.0) - 273.15 == pytest.approx(80.1, abs=0.2)
    assert bt.T_bubble(0.0) - 273.15 == pytest.approx(110.6, abs=0.2)


def test_ethanol_water_azeotrope(ew):
    (x, T), = ew.azeotropes()
    assert x == pytest.approx(0.90, abs=0.03)          # experimental 0.894
    assert T - 273.15 == pytest.approx(78.2, abs=0.5)  # experimental 78.15 degC


def test_dew_and_bubble_are_consistent(bt):
    x = 0.4
    assert bt.T_dew(bt.y_eq(x)) == pytest.approx(bt.T_bubble(x), abs=1e-3)


def test_azeotrope_between_products_is_rejected(ew):
    with pytest.raises(ValueError, match="Azeotrope"):
        design(ew, spec(feeds=[Feed(100, 0.5, 1)], x_D=0.95, x_W=0.02))


# ------------------------------------------------------------- regression on the original script
def test_original_script_results():
    vle = table_vle(os.path.join(EXAMPLES, "benzene_toluene_constant_alpha.csv"))
    d = design(vle, ColumnSpec(feeds=[Feed(100, 0.6, 0.45)], x_D=0.95, x_W=0.10, reflux_factor=1.6))
    assert d.param_min == pytest.approx(1.20, abs=0.005)
    assert d.param == pytest.approx(1.92, abs=0.005)
    assert d.N == pytest.approx(9.58, abs=0.01)
    assert d.stepping.event_stages == [5]
    assert d.N_min == pytest.approx(5.82, abs=0.01)


# ------------------------------------------------------------- consistency checks
def test_mass_balance_closes(bt):
    s = spec(feeds=[Feed(60, 0.7, 1), Feed(40, 0.35, 0.5)], side_draws=[SideDraw(10, 0.5)], x_D=0.97, x_W=0.03)
    m = design(bt, s).model
    assert m.D + m.W + 10 == pytest.approx(100)
    assert m.D * 0.97 + m.W * 0.03 + 10 * 0.5 == pytest.approx(60 * 0.7 + 40 * 0.35)
    # The bottom operating line passes through (x_W, x_W)
    assert m.sections[-1].y(0.03) == pytest.approx(0.03)


def test_underwood_matches_graphical_rmin_for_constant_alpha():
    vle = constant_alpha_vle(2.5)
    for q in (1.0, 0.5, 1.2):
        d = design(vle, spec(feeds=[Feed(100, 0.5, q)]))
        assert d.param_min == pytest.approx(underwood(0.95, 0.5, q, 2.5), rel=2e-3)


def test_fenske_matches_total_reflux_stepping():
    vle = constant_alpha_vle(2.5)
    d = design(vle, spec())
    assert d.N_min == pytest.approx(fenske(0.95, 0.05, 2.5), abs=0.3)


def test_tangent_pinch_detected(ew):
    d = design(ew, spec(feeds=[Feed(100, 0.1, 1)], x_D=0.85, x_W=0.005))
    assert d.pinch_kind == "tangent pinch"
    # The naive feed-pinch estimate is far too low
    y_q = float(ew.y_eq(0.1))
    slope = (0.85 - y_q) / (0.85 - 0.1)
    assert d.param_min > 1.3 * slope / (1 - slope)


def test_stages_decrease_with_reflux(bt):
    sweep = reflux_sweep(bt, spec(), [1.1, 1.5, 2.0, 3.0])
    assert list(sweep["N_theoretical"]) == sorted(sweep["N_theoretical"], reverse=True)


def test_q_from_temperature(bt):
    Tb = float(bt.T_bubble(0.6)) - 273.15
    Td = float(bt.T_dew(0.6)) - 273.15
    assert feed_q(bt, 0.6, Tb - 1e-6) == pytest.approx(1.0, abs=1e-4)
    assert feed_q(bt, 0.6, Td + 1e-6) == pytest.approx(0.0, abs=1e-4)
    assert feed_q(bt, 0.6, 40) > 1
    assert 0 < feed_q(bt, 0.6, 0.5 * (Tb + Td)) < 1
    assert feed_q(bt, 0.6, 130) < 0


def test_partial_condenser_counts_as_stage(bt):
    total = design(bt, spec())
    partial = design(bt, spec(condenser="partial"))
    assert partial.N == pytest.approx(total.N)
    assert partial.theoretical_trays == pytest.approx(total.theoretical_trays - 1)


def test_murphree_efficiency(bt):
    ideal = design(bt, spec())
    assert design(bt, spec(murphree=1.0)).actual_trays == math.ceil(ideal.theoretical_trays)
    real = design(bt, spec(murphree=0.6))
    assert real.actual_trays > ideal.theoretical_trays / 0.8


def test_oconnell_efficiency_is_typical(bt):
    d = design(bt, spec(overall_efficiency="oconnell"))
    assert 0.45 < d.overall_efficiency < 0.7


def test_direct_steam(ew):
    d = design(ew, spec(feeds=[Feed(100, 0.1, 1)], x_D=0.75, x_W=0.002, reboiler="steam"))
    m = d.model
    assert m.sections[-1].y(0.002) == pytest.approx(0.0, abs=1e-12)   # line through (x_W, 0)
    assert m.sections[-1].V == pytest.approx(m.steam)
    assert m.D + m.W == pytest.approx(100 + m.steam)
    assert d.Q_R is None


def test_stripping_column(bt):
    d = design(bt, ColumnSpec(feeds=[Feed(100, 0.3, 1)], x_W=0.02, configuration="stripping",
                              boilup_factor=1.5))
    m = d.model
    assert m.D * m.x_D + m.W * m.x_W == pytest.approx(30)
    assert d.Q_C is None and d.Q_R > 0
    assert d.param == pytest.approx(1.5 * d.param_min)


def test_rectifying_column(bt):
    d = design(bt, ColumnSpec(feeds=[Feed(100, 0.5, 0.0)], x_D=0.9, configuration="rectifying",
                              reflux_factor=1.5))
    m = d.model
    assert m.D * 0.9 + m.W * m.x_W == pytest.approx(50)
    assert d.Q_R is None and d.Q_C > 0


# ------------------------------------------------------------- rating
def test_rating_round_trip(bt):
    s = spec()
    d = design(bt, s)
    n = math.ceil(d.N)
    r = rate(bt, s, n, d.stepping.event_stages, d.param, d.model.D)
    assert r.x_D >= 0.95 - 1e-6 and r.x_W <= 0.05 + 1e-6
    # Designing for the purities the column reaches needs exactly its stage count
    d2 = design(bt, replace(s, x_D=r.x_D, x_W=r.x_W, reflux_ratio=d.param))
    assert d2.N == pytest.approx(n, abs=0.05)


# ------------------------------------------------------------- Ponchon-Savarit
def test_ponchon_reduces_to_mccabe_thiele_with_constant_molar_overflow(bt):
    c1, c2 = bt.components
    c1 = replace(c1, cp_liq=0.0, Tc=1e9, dHvap_b=30.0)
    c2 = replace(c2, cp_liq=0.0, Tc=1e9, dHvap_b=30.0)
    cmo = VLE(bt.x, bt.y, bt.T, (c1, c2))
    f = Feed(100, 0.6, 0.45)
    d = design(cmo, spec(feeds=[f], x_W=0.1, reflux_factor=1.6))
    p = ponchon_savarit(cmo, 0.95, 0.1, f, R=d.param)
    assert p.R_min == pytest.approx(d.param_min, rel=1e-3)
    assert p.N == pytest.approx(d.N, abs=1e-3)
    assert p.Q_R == pytest.approx(d.Q_R, rel=1e-6)


def test_ponchon_energy_balance(bt):
    f = Feed(100, 0.6, 1.0)
    p = ponchon_savarit(bt, 0.95, 0.05, f, R=1.2)
    D = 100 * (0.6 - 0.05) / 0.9
    from mccabe_thiele.ponchon import EnthalpyDiagram
    ed = EnthalpyDiagram(bt)
    h_in = 100 * p.feed_point[1] + p.Q_R * 3.6
    h_out = D * ed.h(0.95) + (100 - D) * ed.h(0.05) + p.Q_C * 3.6
    assert h_in == pytest.approx(h_out, rel=1e-6)


# ------------------------------------------------------------- economics
def test_cost_optimum_is_interior(bt):
    sweep = reflux_sweep(bt, spec(overall_efficiency=0.6))
    best = optimum(sweep)
    assert not best["at_bound"]
    assert 1.02 < best["factor"] < 2.0


def test_column_diameter_scales_with_throughput(bt):
    small = costs(design(bt, spec()))
    big = costs(design(bt, spec(feeds=[Feed(400, 0.6, 1.0)])))
    assert big["diameter_m"] == pytest.approx(2 * small["diameter_m"], rel=1e-6)


# ------------------------------------------------------------- config and CLI
@pytest.mark.parametrize("name", ["benzene_toluene", "original_constant_alpha", "ethanol_water",
                                  "ethanol_water_direct_steam", "two_feeds_side_draw", "stripper",
                                  "enricher"])
def test_examples_run(name, tmp_path):
    cfg = cfgmod.load(os.path.join(EXAMPLES, name + ".yaml"))
    cfg.setdefault("output", {})["directory"] = str(tmp_path)
    res = run(cfg, plots=False, verbose=False)
    assert res["design"].N > 0
    assert all(os.path.exists(f) for f in res["files"])


def test_unknown_config_key_is_reported():
    with pytest.raises(KeyError, match="reflux_facter"):
        cfgmod.build({"system": {"components": ["benzene", "toluene"]},
                      "feeds": [{"flow": 1, "z": 0.5, "q": 1}],
                      "specs": {"x_D": 0.9, "x_W": 0.1, "reflux_facter": 1.2}})


def test_invalid_specs_are_rejected(bt):
    with pytest.raises(ValueError):
        design(bt, spec(x_D=0.5, x_W=0.6))
    with pytest.raises(ValueError, match="below the minimum"):
        design(bt, spec(reflux_ratio=0.1, reflux_factor=None))
