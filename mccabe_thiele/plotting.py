"""Figures in publication style: Helvetica (or the nearest fallback), bold text
of at least 16 pt, no titles (the file name describes the figure), and every
figure saved as PDF, PNG (300 dpi) and SVG."""

import os
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .column import ColumnModel, Feed

COLORS = {"equilibrium": (0.05, 0.35, 0.75), "q": (0.05, 0.80, 0.60),
          "rectifying": (0.90, 0.30, 0.05), "stripping": (0.90, 0.50, 0.70),
          "middle": (0.90, 0.60, 0.05), "stages": (0.05, 0.05, 0.05),
          "murphree": (0.35, 0.70, 0.90), "vapor": (0.90, 0.30, 0.05)}
MIDDLE_COLORS = [(0.90, 0.60, 0.05), (0.55, 0.35, 0.75), (0.60, 0.45, 0.20), (0.35, 0.60, 0.25)]
LEGEND = {"weight": "bold", "size": 16}


def setup_style(fontsize=18):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "Nimbus Sans", "DejaVu Sans"],
        "font.weight": "bold", "axes.labelweight": "bold", "mathtext.default": "bf",
        "font.size": fontsize, "xtick.labelsize": fontsize, "ytick.labelsize": fontsize,
        "pdf.fonttype": 42, "svg.fonttype": "none",
    })


def save_all_formats(fig, output_dir, basename):
    os.makedirs(output_dir, exist_ok=True)
    paths = []
    for ext in ("pdf", "png", "svg"):
        path = os.path.join(output_dir, "%s.%s" % (basename, ext))
        fig.savefig(path, dpi=300, bbox_inches="tight")
        paths.append(path)
    plt.close(fig)
    return paths


def _axes(xlabel, ylabel, figsize=(10, 8)):
    setup_style()
    fig, ax = plt.subplots(figsize=figsize)
    ax.grid(True, alpha=0.5)
    ax.minorticks_on()
    ax.grid(visible=True, which="minor", color="#999999", linestyle="-", alpha=0.2)
    ax.set_xlabel(xlabel, fontsize=22)
    ax.set_ylabel(ylabel, fontsize=22)
    return fig, ax


def system_slug(vle):
    if vle.components:
        return "_".join(c.name.replace("-", "") for c in vle.components)
    return "binary"


def _light_name(vle):
    return vle.components[0].name.capitalize() if vle.components else "Light component"


def _section_color(k, n):
    if k == 0 and n > 1:
        return COLORS["rectifying"]
    if k == n - 1:
        return COLORS["stripping"] if n > 1 else COLORS["rectifying"]
    return MIDDLE_COLORS[(k - 1) % len(MIDDLE_COLORS)]


def _plot_sections(ax, model, style="-", alpha=1.0, label=True):
    n = len(model.sections)
    bounds = [model.x_start] + list(model.switches) + [model.x_end]
    for k, sec in enumerate(model.sections):
        x = np.array([bounds[k + 1], bounds[k]])
        ax.plot(x, sec.y(x), style, color=_section_color(k, n), linewidth=2.5 if style == "-" else 2,
                alpha=alpha, label=("%s line" % sec.name if label else None))


def plot_mccabe_thiele(d, output_dir, basename=None):
    vle, model, spec = d.vle, d.model, d.spec
    light = _light_name(vle)
    fig, ax = _axes("%s mole fraction in liquid phase" % light,
                    "%s mole fraction in vapor phase" % light, (10, 9))
    ax.plot([0, 1], [0, 1], color="gray", linewidth=1.5)
    ax.plot(vle.x, vle.y, color=COLORS["equilibrium"], linewidth=2.5, label="Equilibrium curve")

    # Minimum-reflux operating lines
    try:
        _plot_sections(ax, ColumnModel(vle, spec, d.param_min), "--", 0.6, label=False)
        ax.plot([], [], "--", color="gray", linewidth=2,
                label="Minimum %s" % ("boil-up" if spec.configuration == "stripping" else "reflux"))
    except ValueError:
        pass
    if d.pinch:
        ax.plot(*d.pinch, "s", color="k", markersize=8, markerfacecolor="none", markeredgewidth=2,
                label="Pinch (%s)" % d.pinch_kind)

    # q-lines to the point where the operating lines meet
    feeds = [e for e in model.events if isinstance(e, Feed)] or spec.feeds
    for i, f in enumerate(feeds):
        if spec.configuration == "full":
            k = model.events.index(f)
            xe = model.switches[k]
            ye = model.sections[k].y(xe)
        elif spec.configuration == "stripping":
            xe, ye = model.x_start, model.y_start
        else:
            xe = model.x_end
            ye = model.sections[0].y(xe)
        ax.plot([f.z, xe], [f.z, ye], color=COLORS["q"], linewidth=2.5, label="q-line" if i == 0 else None)
        ax.plot(f.z, f.z, "o", color="k", markersize=9, zorder=5)
        ax.annotate("$z_{F%s}$" % ("" if len(feeds) == 1 else i + 1), (f.z, f.z), textcoords="offset points",
                    xytext=(8, -22), fontsize=20)

    # Side draws: the operating lines meet at the draw composition
    for i, sd in enumerate(spec.side_draws):
        k = model.events.index(sd)
        xe = model.switches[k]
        ye = model.sections[k].y(xe)
        ax.plot([xe, xe], [xe, ye], ":", color="gray", linewidth=2.5, label="Side draw" if i == 0 else None)
        ax.plot(xe, ye, "D", color="gray", markersize=9, zorder=5)

    _plot_sections(ax, model)
    st = d.stepping
    ax.plot(st.X, st.Y, color=COLORS["stages"], linewidth=1.5, label="Ideal stages (%.1f)" % st.N)
    if d.murphree_stepping is not None:
        ms = d.murphree_stepping
        ax.plot(ms.X, ms.Y, color=COLORS["murphree"], linewidth=1.5,
                label="Actual stages, $E_{MV}$ = %.2f (%d)" % (spec.murphree, ms.N))

    for label, x in (("$x_D$", model.x_D), ("$x_W$", model.x_W)):
        if 0 < x < 1:
            ax.plot(x, x, "o", color="k", markersize=9, zorder=5)
            offset = (8, -22) if x > 0.06 else (14, 4)
            ax.annotate(label, (x, x), textcoords="offset points", xytext=offset, fontsize=20)

    ax.legend(loc="lower right", frameon=False, prop=LEGEND)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    basename = basename or "mccabe_thiele_%s_theoretical_stages" % system_slug(vle)
    return save_all_formats(fig, output_dir, basename)


def plot_rating(r, vle, output_dir, basename=None):
    light = _light_name(vle)
    fig, ax = _axes("%s mole fraction in liquid phase" % light,
                    "%s mole fraction in vapor phase" % light, (10, 9))
    ax.plot([0, 1], [0, 1], color="gray", linewidth=1.5)
    ax.plot(vle.x, vle.y, color=COLORS["equilibrium"], linewidth=2.5, label="Equilibrium curve")
    model = r.model
    _plot_sections(ax, model)
    ax.plot(r.stepping.X, r.stepping.Y, color=COLORS["stages"], linewidth=1.5,
            label="Existing column (%d stages)" % r.stepping.N)
    for label, x in (("$x_D$", r.x_D), ("$x_W$", r.x_W)):
        ax.plot(x, x, "o", color="k", markersize=9, zorder=5)
        ax.annotate(label, (x, x), textcoords="offset points", xytext=(8, -22), fontsize=20)
    ax.legend(loc="lower right", frameon=False, prop=LEGEND)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    basename = basename or "mccabe_thiele_%s_rating_existing_column" % system_slug(vle)
    return save_all_formats(fig, output_dir, basename)


def plot_reflux_sweep(sweep, d, output_dir):
    stripping = d.spec.configuration == "stripping"
    xlabel = "Boil-up ratio / minimum" if stripping else "Reflux ratio, $R/R_{min}$"
    slug = system_slug(d.vle)
    paths = []
    fig, ax = _axes(xlabel, "Number of stages")
    ax.plot(sweep["factor"], sweep["N_theoretical"], "o-", color=COLORS["equilibrium"],
            linewidth=2.5, markersize=7, label="Theoretical stages")
    ax.plot(sweep["factor"], sweep["actual_trays"], "s-", color=COLORS["rectifying"],
            linewidth=2.5, markersize=7, label="Actual trays")
    ax.legend(loc="upper right", frameon=False, prop=LEGEND)
    paths += save_all_formats(fig, output_dir, "number_of_stages_vs_reflux_ratio_%s" % slug)

    if "total_annual_cost" in sweep and d.spec.configuration != "rectifying":
        fig, ax = _axes(xlabel, "Annual cost (thousand USD/yr)")
        k = 1e-3
        ax.plot(sweep["factor"], sweep["annualized_capital"] * k, "o-", color=COLORS["equilibrium"],
                linewidth=2.5, markersize=7, label="Annualized capital")
        ax.plot(sweep["factor"], sweep["operating_cost"] * k, "s-", color=COLORS["rectifying"],
                linewidth=2.5, markersize=7, label="Utilities")
        ax.plot(sweep["factor"], sweep["total_annual_cost"] * k, "D-", color=COLORS["stages"],
                linewidth=2.5, markersize=7, label="Total")
        best = sweep.loc[sweep["total_annual_cost"].idxmin()]
        ax.plot(best["factor"], best["total_annual_cost"] * k, "*", color=COLORS["q"], markersize=22,
                zorder=5, label="Optimum (%.2f)" % best["factor"])
        ax.legend(loc="best", frameon=False, prop=LEGEND)
        paths += save_all_formats(fig, output_dir, "total_annual_cost_vs_reflux_ratio_%s" % slug)
    return paths


def plot_temperature_profile(d, output_dir):
    if d.temperatures is None:
        return []
    n = [s[0] for s in d.stepping.stages]
    fig, ax = _axes("Temperature (°C)", "Stage (from top)", (8, 9))
    ax.plot(np.array(d.temperatures) - 273.15, n, "o-", color=COLORS["rectifying"], linewidth=2.5,
            markersize=8)
    for fs in d.stepping.event_stages:
        if fs is not None:
            ax.axhline(fs, color=COLORS["q"], linestyle="--", linewidth=2)
    ax.invert_yaxis()
    return save_all_formats(fig, output_dir, "stage_temperature_profile_%s" % system_slug(d.vle))


def plot_feed_stage_sensitivity(table, d, output_dir):
    """Ball-and-stick plot of stages needed against feed location."""
    if table is None:
        return []
    ok = table[np.isfinite(table["N_theoretical"])]
    fig, ax = _axes("Feed stage (from top)", "Theoretical stages needed")
    ax.grid(False, axis="x")
    best = ok.loc[ok["N_theoretical"].idxmin()]
    ax.vlines(ok["feed_stage"], 0, ok["N_theoretical"], color=COLORS["equilibrium"], linewidth=2.5)
    ax.plot(ok["feed_stage"], ok["N_theoretical"], "o", color=COLORS["equilibrium"], markersize=12)
    ax.plot(best["feed_stage"], best["N_theoretical"], "o", color=COLORS["rectifying"], markersize=14,
            label="Optimal feed stage")
    bad = table[~np.isfinite(table["N_theoretical"])]
    if len(bad):
        ax.plot(bad["feed_stage"], np.zeros(len(bad)), "x", color="k", markersize=12, markeredgewidth=2.5,
                label="Infeasible (pinch)")
    ax.set_ylim(-0.6, ok["N_theoretical"].max() * 1.1)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.legend(loc="upper left", frameon=False, prop=LEGEND)
    return save_all_formats(fig, output_dir, "theoretical_stages_vs_feed_stage_%s" % system_slug(d.vle))


def plot_ponchon(p, vle, output_dir):
    light = _light_name(vle)
    fig, ax = _axes("%s mole fraction" % light, "Enthalpy (kJ/mol)", (10, 9))
    ax.plot(*p.liquid_curve, color=COLORS["equilibrium"], linewidth=2.5, label="Saturated liquid")
    ax.plot(*p.vapor_curve, color=COLORS["vapor"], linewidth=2.5, label="Saturated vapor")
    xl, hl = p.liquid_curve
    yv, hv = p.vapor_curve
    for i, (_, x, y) in enumerate(p.stages):
        ax.plot([x, y], [np.interp(x, xl, hl), np.interp(y, yv, hv)], color=COLORS["stages"],
                linewidth=1.2, label="Tie lines (%.1f stages)" % p.N if i == 0 else None)
    ax.plot([p.delta_D[0], p.delta_W[0]], [p.delta_D[1], p.delta_W[1]], "--", color=COLORS["q"],
            linewidth=2, label="$\\Delta_D$ - F - $\\Delta_W$")
    for label, pt in (("$\\Delta_D$", p.delta_D), ("$\\Delta_W$", p.delta_W), ("F", p.feed_point)):
        ax.plot(*pt, "o", color="k", markersize=9, zorder=5)
        ax.annotate(label, pt, textcoords="offset points", xytext=(10, -5), fontsize=20)
    ax.set_xlim(0, 1)
    ax.legend(loc="lower right", frameon=False, prop=LEGEND)
    return save_all_formats(fig, output_dir, "ponchon_savarit_enthalpy_composition_%s" % system_slug(vle))
