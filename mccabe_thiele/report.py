"""Text report, stage-by-stage table and Excel/CSV export."""

import math
import pandas as pd

from .column import Feed


def stage_table(d):
    """One row per equilibrium stage: compositions, temperature and internal flows."""
    vle, model, spec = d.vle, d.model, d.spec
    n_last = len(d.stepping.stages)
    names = {}
    for ev, st in zip(model.events, d.stepping.event_stages):
        if st is not None:
            names.setdefault(st, []).append(ev.name)
    rows = []
    for i, (n, x, y, s) in enumerate(d.stepping.stages):
        sec = model.sections[s]
        kind = "Tray"
        if n == 1 and spec.configuration != "stripping" and spec.condenser == "partial":
            kind = "Partial condenser"
        if n == n_last and spec.configuration in ("full", "stripping") and spec.reboiler == "partial":
            kind = "Reboiler"
        rows.append({
            "stage": n, "type": kind, "section": sec.name,
            "feed/draw": ", ".join(names.get(n, [])),
            "x": x, "y": y, "y_eq(x)": float(vle.y_eq(x)),
            "T_C": d.temperatures[i] - 273.15 if d.temperatures else None,
            "L_kmol_h": sec.L, "V_kmol_h": sec.V,
        })
    return pd.DataFrame(rows)


def _included(spec):
    parts = []
    if spec.configuration in ("full", "stripping") and spec.reboiler == "partial":
        parts.append("reboiler")
    if spec.configuration in ("full", "rectifying") and spec.condenser == "partial":
        parts.append("partial condenser")
    return "incl. " + " and ".join(parts) if parts else ""


def summary(d, cost=None, sweep_opt=None, ponchon=None):
    """Key results as (label, value, unit) rows."""
    spec, m = d.spec, d.model
    strip = spec.configuration == "stripping"
    rows = [("VLE", d.vle.label, ""), ("Configuration", spec.configuration, "")]
    for f in spec.feeds:
        rows.append(("%s: flow, z, q" % f.name, "%.2f, %.4f, %.3f" % (f.flow, f.z, f.q), "kmol/h, -, -"))
    for s in spec.side_draws:
        rows.append(("%s (%s)" % (s.name, s.phase), "%.2f at x = %.4f" % (s.flow, s.x), "kmol/h"))
    rows += [
        ("Distillate D, x_D", "%.3f, %.4f" % (m.D, m.x_D), "kmol/h, -"),
        ("Bottoms W, x_W", "%.3f, %.4f" % (m.W, m.x_W), "kmol/h, -"),
    ]
    if m.steam:
        rows.append(("Direct steam", "%.3f" % m.steam, "kmol/h"))
    name = "boil-up ratio V/W" if strip else "reflux ratio R"
    rows += [
        ("Minimum " + name, "%.4f" % d.param_min, ""),
        ("Pinch", "%s at x = %.4f" % (d.pinch_kind, d.pinch[0]) if d.pinch else "-", ""),
        ("Operating " + name, "%.4f (%.2f x min)" % (d.param, d.param / d.param_min), ""),
    ]
    for sec in m.sections:
        rows.append(("%s section L, V" % sec.name, "%.2f, %.2f" % (sec.L, sec.V), "kmol/h"))
    if not strip and spec.configuration == "full" and spec.reboiler == "partial":
        rows.append(("Boil-up ratio V/W", "%.3f" % (m.sections[-1].V / m.W), ""))
    rows += [
        ("Theoretical stages", "%.2f" % d.N, _included(spec)),
        ("Theoretical trays", "%.2f" % d.theoretical_trays, ""),
        ("Optimal feed/draw stage(s)", ", ".join("%s: %s" % (e.name, s) for e, s in
                                                 zip(m.events, d.stepping.event_stages)) or "-", "from top"),
    ]
    if d.N_min is not None:
        rows.append(("Minimum stages N_min (total reflux)", "%.2f" % d.N_min, ""))
    if d.murphree_stepping is not None:
        rows.append(("Murphree efficiency E_MV", "%.2f" % spec.murphree, ""))
    if d.overall_efficiency is not None:
        rows.append(("Overall efficiency E_o", "%.3f" % d.overall_efficiency,
                     "O'Connell" if str(spec.overall_efficiency).lower() == "oconnell" else ""))
    no_eff = d.murphree_stepping is None and d.overall_efficiency is None
    rows.append(("Actual trays", "%d" % d.actual_trays, "assuming 100 % efficiency" if no_eff else ""))
    if d.murphree_stepping is not None or d.overall_efficiency is not None:
        rows.append(("Actual feed/draw tray(s)", ", ".join(str(t) for t in d.feed_trays), "from top"))
    if d.temperatures:
        rows.append(("Top / bottom temperature", "%.1f / %.1f" % (d.temperatures[0] - 273.15,
                                                             d.temperatures[-1] - 273.15), "degC"))
    if d.Q_C is not None:
        rows.append(("Condenser duty Q_C", "%.1f" % d.Q_C, "kW"))
    if d.Q_R is not None:
        rows.append(("Reboiler duty Q_R", "%.1f" % d.Q_R, "kW"))
    if d.shortcut:
        s = d.shortcut
        rows += [("Shortcut (FUG): alpha_avg", "%.3f" % s["alpha"], ""),
                 ("Shortcut: N_min (Fenske)", "%.2f" % s["N_min"], ""),
                 ("Shortcut: R_min (Underwood)", "%.4f" % s["R_min"], ""),
                 ("Shortcut: N (Gilliland)", "%.2f" % s["N"] if math.isfinite(s["N"])
                  else "n/a: R is below the constant-alpha R_min", ""),
                 ("Shortcut: feed stage (Kirkbride)", "%.1f" % s["feed_stage"] if math.isfinite(s["N"])
                  else "n/a", "")]
        a_top, a_bot = float(d.vle.alpha(spec.x_D)), float(d.vle.alpha(spec.x_W))
        if max(a_top, a_bot) / min(a_top, a_bot) > 1.3:
            rows.append(("Warning", "alpha varies from %.2f to %.2f; the shortcut method is unreliable"
                         % (min(a_top, a_bot), max(a_top, a_bot)), ""))
    if ponchon is not None:
        rows += [("Ponchon-Savarit: R_min", "%.4f" % ponchon.R_min, ""),
                 ("Ponchon-Savarit: stages at same R", "%.2f" % ponchon.N, ""),
                 ("Ponchon-Savarit: feed stage", "%d" % ponchon.feed_stage, ""),
                 ("Ponchon-Savarit: Q_C, Q_R", "%.1f, %.1f" % (ponchon.Q_C, ponchon.Q_R), "kW"),
                 ("Ponchon-Savarit: L/V top, bottom", "%.3f, %.3f" % (ponchon.L_over_V_top,
                                                                    ponchon.L_over_V_bottom), "")]
    if cost:
        rows += [("Column diameter", "%.2f" % cost["diameter_m"], "m"),
                 ("Column height", "%.1f" % cost["height_m"], "m"),
                 ("Condenser / reboiler area", "%.1f / %.1f" % (cost["condenser_area_m2"],
                                                                cost["reboiler_area_m2"]), "m2"),
                 ("Installed capital cost", "%.0f" % cost["capital_cost"], "USD"),
                 ("Utility cost", "%.0f" % cost["operating_cost"], "USD/yr"),
                 ("Total annual cost", "%.0f" % cost["total_annual_cost"], "USD/yr")]
    if sweep_opt is not None:
        label = "Cost-optimal boil-up / minimum" if strip else "Cost-optimal R/R_min"
        rows.append((label, "%.3f (%s = %.3f)" % (sweep_opt["factor"], "V/W" if strip else "R",
                                                 sweep_opt["ratio"]), ""))
        if sweep_opt.get("at_bound"):
            rows.append(("Warning", "the cost optimum is at the edge of the swept range; widen "
                                    "analysis.reflux_factors", ""))
    for w in d.warnings:
        rows.append(("Warning", w, ""))
    return rows


def format_summary(rows):
    width = max(len(r[0]) for r in rows) + 2
    return "\n".join("%-*s %s %s" % (width, label, value, unit) for label, value, unit in rows)


def export(d, output_dir, basename, rows, sweep=None, feed_table=None, rating=None, excel=True):
    """Write the stage table (CSV) and a workbook with every table."""
    import os
    os.makedirs(output_dir, exist_ok=True)
    stages = stage_table(d)
    paths = [os.path.join(output_dir, basename + "_stages.csv")]
    stages.to_csv(paths[0], index=False)
    if excel:
        path = os.path.join(output_dir, basename + "_report.xlsx")
        with pd.ExcelWriter(path) as xw:
            pd.DataFrame(rows, columns=["Item", "Value", "Unit"]).to_excel(xw, sheet_name="Summary", index=False)
            stages.to_excel(xw, sheet_name="Stages", index=False)
            if sweep is not None and len(sweep):
                sweep.to_excel(xw, sheet_name="Reflux sweep", index=False)
            if feed_table is not None:
                feed_table.replace(math.inf, None).to_excel(xw, sheet_name="Feed stage", index=False)
            if rating is not None:
                pd.DataFrame(rating, columns=["Item", "Value", "Unit"]).to_excel(xw, sheet_name="Rating", index=False)
        paths.append(path)
    return paths
