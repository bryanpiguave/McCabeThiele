"""Load a design case from YAML or JSON and turn it into model objects."""

import json
import os
from dataclasses import fields

import yaml

from .column import ColumnSpec, Feed, SideDraw
from .economics import Economics
from .thermo import make_vle

SECTIONS = {"system", "feeds", "side_draws", "specs", "column", "efficiency", "economics",
            "analysis", "rating", "output"}


def load(path):
    with open(path) as fh:
        cfg = json.load(fh) if path.lower().endswith(".json") else yaml.safe_load(fh)
    base = os.path.dirname(os.path.abspath(path))
    vle = cfg.get("system", {}).get("vle", {}) or {}
    if vle.get("file") and not os.path.isabs(vle["file"]) and not os.path.exists(vle["file"]):
        vle["file"] = os.path.join(base, vle["file"])
    return cfg


def set_path(cfg, dotted, value):
    """Set cfg['a']['b'] from 'a.b'; list items by index, e.g. feeds.0.z."""
    keys = dotted.split(".")
    node = cfg
    for k in keys[:-1]:
        node = node[int(k)] if isinstance(node, list) else node.setdefault(k, {})
    last = keys[-1]
    value = yaml.safe_load(value) if isinstance(value, str) else value
    if isinstance(node, list):
        node[int(last)] = value
    else:
        node[last] = value


def _pick(cls, data, where):
    data = dict(data or {})
    allowed = {f.name for f in fields(cls)}
    unknown = set(data) - allowed
    if unknown:
        raise KeyError("Unknown key(s) in %s: %s. Allowed: %s"
                       % (where, ", ".join(sorted(unknown)), ", ".join(sorted(allowed))))
    return cls(**data)


def build(cfg):
    """Return (vle, spec, economics) from a config dict."""
    unknown = set(cfg) - SECTIONS
    if unknown:
        raise KeyError("Unknown config section(s): %s" % ", ".join(sorted(unknown)))
    vle = make_vle(cfg.get("system", {}))
    feeds = [_pick(Feed, f, "feeds[%d]" % i) for i, f in enumerate(cfg.get("feeds", []))]
    if not feeds:
        raise ValueError("At least one feed is required")
    for i, f in enumerate(feeds):
        if f.name == "Feed" and len(feeds) > 1:
            f.name = "Feed %d" % (i + 1)
    draws = [_pick(SideDraw, s, "side_draws[%d]" % i) for i, s in enumerate(cfg.get("side_draws") or [])]
    col = dict(cfg.get("specs", {}))
    col.update(cfg.get("column", {}) or {})
    eff = cfg.get("efficiency", {}) or {}
    if "murphree" in eff:
        col["murphree"] = eff["murphree"]
    if "overall" in eff:
        col["overall_efficiency"] = eff["overall"]
    spec = _pick(ColumnSpec, dict(col, feeds=feeds, side_draws=draws), "specs/column")
    econ = _pick(Economics, cfg.get("economics"), "economics")
    return vle, spec, econ
