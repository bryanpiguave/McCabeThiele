"""Command-line interface.

    python -m mccabe_thiele examples/benzene_toluene.yaml
    python -m mccabe_thiele examples/benzene_toluene.yaml --xD 0.98 --R-factor 1.3 --out results
    python -m mccabe_thiele case.yaml --set efficiency.murphree=0.65 --set feeds.0.T_C=40
"""

import argparse
import sys

from . import config as cfgmod
from .analysis import reflux_sweep, optimum, feed_stage_sensitivity
from .column import design
from .components import COMPONENTS
from .economics import costs
from .plotting import (plot_mccabe_thiele, plot_reflux_sweep, plot_temperature_profile,
                       plot_feed_stage_sensitivity, plot_ponchon, plot_rating, system_slug)
from .ponchon import ponchon_savarit
from .rating import rate
from .report import summary, format_summary, export

SHORTCUTS = {"xD": "specs.x_D", "xW": "specs.x_W", "xF": "feeds.0.z", "q": "feeds.0.q",
             "F": "feeds.0.flow", "TF": "feeds.0.T_C", "R": "specs.reflux_ratio",
             "R_factor": "specs.reflux_factor", "P": "system.pressure_kPa",
             "murphree": "efficiency.murphree", "out": "output.directory"}


def run(cfg, plots=True, verbose=True):
    """Run a full case. Returns a dict with the design and every analysis."""
    vle, spec, econ = cfgmod.build(cfg)
    analysis = cfg.get("analysis", {}) or {}
    out = cfg.get("output", {}) or {}
    out_dir = out.get("directory", "results")
    slug = system_slug(vle)
    res = {"design": design(vle, spec)}
    d = res["design"]

    if vle.has_properties and analysis.get("economics", True):
        res["costs"] = costs(d, econ)
    if analysis.get("reflux_sweep", True):
        res["sweep"] = reflux_sweep(vle, spec, analysis.get("reflux_factors"), econ, d.param_min)
        # In a rectifying column R also changes the bottoms, so costs aren't comparable
        if spec.configuration != "rectifying":
            res["optimum"] = optimum(res["sweep"])
    if analysis.get("feed_stage_sensitivity", True):
        res["feed_table"] = feed_stage_sensitivity(d)
    if (analysis.get("ponchon_savarit", True) and spec.configuration == "full" and len(spec.feeds) == 1
            and not spec.side_draws and spec.reboiler == "partial" and spec.condenser == "total"
            and vle.has_temperature and vle.has_properties):
        res["ponchon"] = ponchon_savarit(vle, spec.x_D, spec.x_W, spec.feeds[0], R=d.param)

    rows = summary(d, res.get("costs"), res.get("optimum"), res.get("ponchon"))
    rating_rows = None
    if cfg.get("rating"):
        r = cfg["rating"]
        rt = rate(vle, spec, int(r["n_stages"]), list(r["feed_stages"]), float(r["reflux_ratio"]),
                  float(r["distillate_flow"]))
        res["rating"] = rt
        rating_rows = [("Stages (incl. reboiler)", "%d" % r["n_stages"], ""),
                       ("Feed stage(s)", ", ".join(map(str, r["feed_stages"])), "from top"),
                       ("Reflux ratio", "%.4f" % r["reflux_ratio"], ""),
                       ("Distillate rate", "%.3f" % rt.D, "kmol/h"),
                       ("Achieved x_D", "%.5f" % rt.x_D, ""),
                       ("Achieved x_W", "%.5f" % rt.x_W, "")]
    res["summary"] = rows

    if verbose:
        print(format_summary(rows))
        if rating_rows:
            print("\nRating of the existing column")
            print(format_summary(rating_rows))

    files = []
    if plots:
        files += plot_mccabe_thiele(d, out_dir)
        files += plot_temperature_profile(d, out_dir)
        if res.get("sweep") is not None and len(res["sweep"]):
            files += plot_reflux_sweep(res["sweep"], d, out_dir)
        files += plot_feed_stage_sensitivity(res.get("feed_table"), d, out_dir)
        if "ponchon" in res:
            files += plot_ponchon(res["ponchon"], vle, out_dir)
        if "rating" in res:
            files += plot_rating(res["rating"], vle, out_dir)
    files += export(d, out_dir, "distillation_%s" % slug, rows, res.get("sweep"), res.get("feed_table"),
                    rating_rows, excel=out.get("excel", True))
    res["files"] = files
    if verbose:
        print("\nWrote %d files to %s/" % (len(files), out_dir))
    return res


def main(argv=None):
    p = argparse.ArgumentParser(prog="mccabe_thiele", description="McCabe-Thiele distillation design")
    p.add_argument("config", nargs="?", help="YAML or JSON design case")
    p.add_argument("--list-components", action="store_true", help="show the built-in components")
    for key, path in SHORTCUTS.items():
        p.add_argument("--" + key.replace("_", "-"), dest=key, help="override %s" % path)
    p.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                   help="override any config entry, e.g. efficiency.overall=oconnell")
    p.add_argument("--no-plots", action="store_true")
    args = p.parse_args(argv)

    if args.list_components:
        for c in COMPONENTS.values():
            print("%-10s MW %6.2f  Tb %6.1f degC" % (c.name, c.mw, c.Tb - 273.15))
        return 0
    if not args.config:
        p.error("a config file is required")
    cfg = cfgmod.load(args.config)
    for key, path in SHORTCUTS.items():
        if getattr(args, key) is not None:
            cfgmod.set_path(cfg, path, getattr(args, key))
    for item in args.set:
        k, _, v = item.partition("=")
        cfgmod.set_path(cfg, k.strip(), v.strip())
    if args.R is not None:
        cfg["specs"].pop("reflux_factor", None)
    try:
        run(cfg, plots=not args.no_plots)
    except (ValueError, KeyError) as e:
        print("Error: %s" % e, file=sys.stderr)
        return 1
    return 0
