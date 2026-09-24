"""Interactive what-if explorer.

    streamlit run app.py
"""

import io
import tempfile

import pandas as pd
import streamlit as st

from mccabe_thiele import ColumnSpec, Feed, design, make_vle, PinchError, ponchon_savarit
from mccabe_thiele.analysis import reflux_sweep, optimum, feed_stage_sensitivity
from mccabe_thiele.components import COMPONENTS, BINARY_PARAMETERS
from mccabe_thiele.economics import Economics, costs
from mccabe_thiele.plotting import (plot_mccabe_thiele, plot_reflux_sweep, plot_temperature_profile,
                                    plot_feed_stage_sensitivity, plot_ponchon)
from mccabe_thiele.report import summary, stage_table

st.set_page_config(page_title="McCabe-Thiele designer", layout="wide")
st.title("McCabe–Thiele column designer")

# ------------------------------------------------------------------ inputs
sb = st.sidebar
sb.header("System")
names = list(COMPONENTS)
light = sb.selectbox("Light component", names, index=names.index("benzene"))
heavy = sb.selectbox("Heavy component", names, index=names.index("toluene"))
models = ["raoult", "constant_alpha"] + sorted({m for p in BINARY_PARAMETERS.values() for m in p})
model = sb.selectbox("VLE model", models)
P = sb.number_input("Pressure (kPa)", 10.0, 1000.0, 101.325)
vle_cfg = {"model": model}
if model == "constant_alpha":
    vle_cfg["alpha"] = sb.number_input("Relative volatility", 1.05, 20.0, 2.45)

sb.header("Column")
configuration = sb.selectbox("Configuration", ["full", "stripping", "rectifying"])
condenser = sb.selectbox("Condenser", ["total", "partial"])
reboiler = sb.selectbox("Reboiler", ["partial", "steam"])

sb.header("Feed")
F = sb.number_input("Flow (kmol/h)", 1.0, 1e5, 100.0)
z = sb.slider("Light-component fraction z", 0.01, 0.99, 0.60, 0.01)
by_T = sb.checkbox("Specify feed temperature instead of q")
if by_T:
    T_F = sb.number_input("Feed temperature (°C)", -50.0, 400.0, 60.0)
    feed = Feed(F, z, T_C=T_F)
else:
    feed = Feed(F, z, sb.slider("q (liquid fraction)", -0.5, 1.5, 1.0 if configuration != "rectifying" else 0.0, 0.05))

sb.header("Specifications")
x_D = sb.slider("Distillate x_D", 0.50, 0.999, 0.95, 0.001, format="%.3f")
x_W = sb.slider("Bottoms x_W", 0.001, 0.50, 0.05, 0.001, format="%.3f")
factor = sb.slider("R/R_min (or boil-up / minimum)", 1.01, 4.0, 1.3, 0.01)
eff_mode = sb.selectbox("Tray efficiency", ["none", "Murphree", "overall (number)", "O'Connell"])
murphree = overall = None
if eff_mode == "Murphree":
    murphree = sb.slider("E_MV", 0.2, 1.0, 0.7, 0.05)
elif eff_mode == "overall (number)":
    overall = sb.slider("E_o", 0.2, 1.0, 0.6, 0.05)
elif eff_mode == "O'Connell":
    overall = "oconnell"

sb.header("Economics")
econ = Economics(steam_cost_per_GJ=sb.number_input("Steam ($/GJ)", 0.0, 100.0, 12.0),
                 cw_cost_per_GJ=sb.number_input("Cooling water ($/GJ)", 0.0, 10.0, 0.5),
                 capital_charge=sb.number_input("Capital charge (1/yr)", 0.05, 1.0, 0.333))

# ------------------------------------------------------------------ calculation
try:
    vle = make_vle({"components": [light, heavy], "pressure_kPa": P, "vle": vle_cfg})
    spec = ColumnSpec(feeds=[feed], x_D=x_D if configuration != "stripping" else None,
                      x_W=x_W if configuration != "rectifying" else None,
                      reflux_factor=factor, boilup_factor=factor, configuration=configuration,
                      condenser=condenser, reboiler=reboiler, murphree=murphree, overall_efficiency=overall)
    d = design(vle, spec)
except (ValueError, KeyError, PinchError) as e:
    st.error(str(e))
    st.stop()

cost = costs(d, econ) if vle.has_properties else None
sweep = reflux_sweep(vle, spec, econ=econ, param_min=d.param_min)
best = optimum(sweep) if configuration != "rectifying" else None
p = None
if (configuration == "full" and reboiler == "partial" and condenser == "total" and vle.has_temperature):
    try:
        p = ponchon_savarit(vle, x_D, x_W, feed, R=d.param)
    except (ValueError, PinchError):
        p = None

c = st.columns(5)
c[0].metric("Minimum ratio", "%.3f" % d.param_min, d.pinch_kind, delta_color="off")
c[1].metric("Operating ratio", "%.3f" % d.param)
c[2].metric("Theoretical stages", "%.2f" % d.N)
c[3].metric("Actual trays", "%d" % d.actual_trays)
if cost:
    c[4].metric("Total annual cost", "$%.0fk/yr" % (cost["total_annual_cost"] / 1e3))
for w in d.warnings:
    st.warning(w)

tmp = tempfile.mkdtemp()
tabs = st.tabs(["Diagram", "Stages", "Reflux & cost", "Feed location", "Ponchon–Savarit", "Summary"])
with tabs[0]:
    st.image(plot_mccabe_thiele(d, tmp)[1])
with tabs[1]:
    st.dataframe(stage_table(d), use_container_width=True)
    files = plot_temperature_profile(d, tmp)
    if files:
        st.image(files[1], width=500)
with tabs[2]:
    for f in plot_reflux_sweep(sweep, d, tmp)[1::3]:
        st.image(f)
    st.dataframe(sweep, use_container_width=True)
with tabs[3]:
    table = feed_stage_sensitivity(d)
    if table is None:
        st.info("Available for a full column with one feed and no side draws.")
    else:
        st.image(plot_feed_stage_sensitivity(table, d, tmp)[1])
with tabs[4]:
    if p is None:
        st.info("Needs a full column with a total condenser, partial reboiler and a temperature-aware VLE model.")
    else:
        st.write("R_min = %.4f, stages at R = %.3f: %.2f (McCabe–Thiele: %.2f)" % (p.R_min, p.R, p.N, d.N))
        st.image(plot_ponchon(p, vle, tmp)[1])
with tabs[5]:
    rows = summary(d, cost, best, p)
    df = pd.DataFrame(rows, columns=["Item", "Value", "Unit"])
    st.table(df)
    buf = io.BytesIO()
    with pd.ExcelWriter(buf) as xw:
        df.to_excel(xw, sheet_name="Summary", index=False)
        stage_table(d).to_excel(xw, sheet_name="Stages", index=False)
        sweep.to_excel(xw, sheet_name="Reflux sweep", index=False)
    st.download_button("Download Excel report", buf.getvalue(), "distillation_report.xlsx")
