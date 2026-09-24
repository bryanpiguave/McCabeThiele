"""McCabe-Thiele design of binary distillation columns.

Assumes constant molar overflow. Flows are in kmol/h, duties in kW.
Compositions are light-component mole fractions.

Supported configurations
  full        condenser + column + reboiler (or direct steam), any number of
              feeds and liquid/vapor side draws
  stripping   feed enters the top stage, reboiler at the bottom, no reflux
  rectifying  vapor feed enters below the bottom stage, condenser at the top

Every section of the column has an operating line V y = L x + d, where d is
the net light-component flow leaving through the top of that section.
"""

from dataclasses import dataclass, field
import math
import numpy as np

KW_PER_KMOLH_KJMOL = 1 / 3.6    # (kmol/h)(kJ/mol) -> kW


class PinchError(RuntimeError):
    """The stepping can't reach the target: the operating line touches the equilibrium curve."""


# ---------------------------------------------------------------- inputs
@dataclass
class Feed:
    flow: float                 # kmol/h
    z: float                    # light-component mole fraction
    q: float = None             # liquid fraction; computed from T_C when None
    T_C: float = None           # feed temperature, degC
    name: str = "Feed"


@dataclass
class SideDraw:
    flow: float                 # kmol/h
    x: float                    # composition of the withdrawn stream
    phase: str = "liquid"       # liquid | vapor
    name: str = "Side draw"


@dataclass
class ColumnSpec:
    feeds: list
    x_D: float = None
    x_W: float = None
    reflux_ratio: float = None           # L/D; overrides reflux_factor
    reflux_factor: float = 1.5           # R / R_min
    configuration: str = "full"          # full | stripping | rectifying
    condenser: str = "total"             # total | partial
    reboiler: str = "partial"            # partial | steam
    side_draws: list = field(default_factory=list)
    boilup_ratio: float = None           # V/W (stripping columns); overrides boilup_factor
    boilup_factor: float = 1.5           # boil-up / minimum boil-up (stripping columns)
    murphree: float = None               # Murphree vapor efficiency of the trays
    overall_efficiency: object = None    # number, or "oconnell"
    max_stages: int = 300


# ---------------------------------------------------------------- operating lines
@dataclass
class Section:
    L: float
    V: float
    d: float                    # V y = L x + d
    name: str = ""

    def y(self, x):
        return (self.L * x + self.d) / self.V

    def x(self, y):
        return (self.V * y - self.d) / self.L


def _event_key(ev, vle):
    if isinstance(ev, Feed):
        return ev.z
    return ev.x if ev.phase == "liquid" else vle.x_eq(ev.x)


def _apply_event(L, V, d, ev):
    if isinstance(ev, Feed):
        return L + ev.q * ev.flow, V - (1 - ev.q) * ev.flow, d - ev.flow * ev.z
    if ev.phase == "liquid":
        return L - ev.flow, V, d + ev.flow * ev.x
    return L, V + ev.flow, d + ev.flow * ev.x


def _intersect(a, b):
    """x where two operating lines cross."""
    num = b.d / b.V - a.d / a.V
    den = a.L / a.V - b.L / b.V
    return num / den if abs(den) > 1e-14 else float("nan")


def q_line_intersection(section, feed):
    """Point where an operating line meets the feed q-line."""
    if math.isclose(feed.q, 1.0):
        x = feed.z
    else:
        m = feed.q / (feed.q - 1)
        # section: y = s x + c ; q-line: y = m (x - z) + z
        s, c = section.L / section.V, section.d / section.V
        x = (c - feed.z + m * feed.z) / (m - s)
    return x, section.y(x)


# ---------------------------------------------------------------- feed condition
def feed_q(vle, z, T_C):
    """Liquid fraction q of a feed at temperature T_C."""
    T = T_C + 273.15
    Tb, Td = float(vle.T_bubble(z)), float(vle.T_dew(z))
    if T <= Tb:
        return 1 + vle.cp_liq(z) * (Tb - T) / vle.latent_heat(z, Tb)
    if T >= Td:
        return -vle.cp_vap(z) * (T - Td) / vle.latent_heat(z, Td)
    # Two-phase: liquid on the bubble curve at T, vapor in equilibrium with it
    x = np.interp(T, vle.T[::-1], vle.x[::-1])
    y = vle.y_eq(x)
    return 1 - (z - x) / (y - x)


# ---------------------------------------------------------------- column model
class ColumnModel:
    """Operating lines of a column for a given reflux (or boil-up) ratio."""

    def __init__(self, vle, spec, param):
        self.vle, self.spec, self.param = vle, spec, param
        cfg = spec.configuration
        builder = {"full": self._full, "stripping": self._stripping,
                   "rectifying": self._rectifying}.get(cfg)
        if builder is None:
            raise ValueError("Unknown configuration '%s'" % cfg)
        builder(param)
        self.valid = all(s.L > 0 and s.V > 0 for s in self.sections) and self.D > 0 and self.W > 0

    # -- full column ------------------------------------------------
    def _full(self, R):
        s, vle = self.spec, self.vle
        events = sorted(s.feeds + s.side_draws, key=lambda e: -_event_key(e, vle))
        F = sum(f.flow for f in s.feeds)
        Fz = sum(f.flow * f.z for f in s.feeds)
        Sd = sum(d.flow for d in s.side_draws)
        Sdx = sum(d.flow * d.x for d in s.side_draws)
        self.steam = 0.0
        if s.reboiler == "steam":
            # unknowns D, W, S: totals, light component, bottom vapor = steam
            vap_in = sum((1 - f.q) * f.flow for f in s.feeds)
            vap_out = sum(d.flow for d in s.side_draws if d.phase == "vapor")
            A = np.array([[1, 1, -1], [s.x_D, s.x_W, 0], [R + 1, 0, -1]])
            b = np.array([F - Sd, Fz - Sdx, vap_in - vap_out])
            self.D, self.W, self.steam = np.linalg.solve(A, b)
        else:
            A = np.array([[1, 1], [s.x_D, s.x_W]])
            self.D, self.W = np.linalg.solve(A, np.array([F - Sd, Fz - Sdx]))
        L, V, d = R * self.D, (R + 1) * self.D, self.D * s.x_D
        self.sections = [Section(L, V, d, "Rectifying")]
        self.events = events
        for i, ev in enumerate(events):
            L, V, d = _apply_event(L, V, d, ev)
            last = i == len(events) - 1
            self.sections.append(Section(L, V, d, "Stripping" if last else "Middle %d" % (i + 1)))
        self.switches = [_intersect(a, b) for a, b in zip(self.sections, self.sections[1:])]
        self.x_start, self.y_start = s.x_D, s.x_D
        self.x_end = s.x_W
        self.x_D, self.x_W = s.x_D, s.x_W

    # -- stripping column -------------------------------------------
    def _stripping(self, boilup):
        s = self.spec
        if len(s.feeds) != 1 or s.side_draws:
            raise ValueError("A stripping column takes exactly one feed and no side draws")
        f = s.feeds[0]
        self.W = f.q * f.flow / (1 + boilup)
        self.D = f.flow - self.W
        self.steam = 0.0
        sec = Section(self.W * (1 + boilup), self.W * boilup, -self.W * s.x_W, "Stripping")
        self.sections, self.switches, self.events = [sec], [], []
        self.x_start, self.y_start = q_line_intersection(sec, f)
        self.x_end = self.x_W = s.x_W
        self.x_D = (f.flow * f.z - self.W * s.x_W) / self.D

    # -- rectifying column ------------------------------------------
    def _rectifying(self, R):
        s = self.spec
        if len(s.feeds) != 1 or s.side_draws:
            raise ValueError("A rectifying column takes exactly one feed and no side draws")
        f = s.feeds[0]
        if f.q >= 1:
            raise ValueError("A rectifying column needs a feed with vapor (q < 1)")
        self.D = (1 - f.q) * f.flow / (R + 1)
        self.W = f.flow - self.D
        self.steam = 0.0
        sec = Section(R * self.D, (R + 1) * self.D, self.D * s.x_D, "Rectifying")
        self.sections, self.switches, self.events = [sec], [], []
        self.x_start = self.y_start = self.x_D = s.x_D
        self.x_end, _ = q_line_intersection(sec, f)
        self.x_W = (f.flow * f.z - self.D * s.x_D) / self.W

    # -- helpers ----------------------------------------------------
    def section_index(self, x):
        k = 0
        while k < len(self.switches) and x < self.switches[k]:
            k += 1
        return k

    def y_op(self, x):
        return self.sections[self.section_index(x)].y(x)

    def min_gap(self):
        """Smallest distance y*(x) - y_op(x) over the column, and where it occurs."""
        if not self.valid or self.x_end <= 0 or self.x_end >= self.x_start:
            return -1.0, None
        sw = [x for x in self.switches if not math.isnan(x)]
        if len(sw) != len(self.switches) or sw != sorted(sw, reverse=True):
            return -1.0, None
        x = np.linspace(self.x_end, self.x_start, 4001)
        idx = np.searchsorted(-np.array(sw), -x, side="right") if sw else np.zeros(len(x), int)
        y = np.empty_like(x)
        for k, sec in enumerate(self.sections):
            m = idx == k
            y[m] = sec.y(x[m])
        gap = self.vle.y_eq(x) - y
        i = int(np.argmin(gap))
        return gap[i], (x[i], y[i])

    def feasible(self):
        return self.min_gap()[0] > 0


def find_minimum(vle, spec):
    """Minimum reflux (or boil-up) ratio by bisection on feasibility.

    The pinch can be at a feed (operating lines meet on the equilibrium curve)
    or a tangent pinch inside a section; both are found the same way.
    """
    hi = 1.0
    while not ColumnModel(vle, spec, hi).feasible():
        hi *= 2
        if hi > 1e5:
            raise ValueError("The separation is infeasible at any reflux: check for an azeotrope "
                             "between the product compositions")
    lo = 0.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if ColumnModel(vle, spec, mid).feasible():
            hi = mid
        else:
            lo = mid
    return (hi,) + pinch_info(vle, spec, hi)


def pinch_info(vle, spec, param_min):
    """Pinch point at the minimum reflux (or boil-up) and whether it sits at a feed."""
    model = ColumnModel(vle, spec, param_min)
    _, pinch = model.min_gap()
    kind = "tangent pinch"
    for x_sw in model.switches + ([model.x_start] if spec.configuration == "stripping" else []) \
            + ([model.x_end] if spec.configuration == "rectifying" else []):
        if pinch is not None and abs(pinch[0] - x_sw) < 2e-3:
            kind = "feed pinch"
    return pinch, kind


# ---------------------------------------------------------------- stage stepping
@dataclass
class Stepping:
    X: np.ndarray               # staircase x
    Y: np.ndarray               # staircase y
    stages: list                # (stage, x_n, y_n, section index)
    N: float                    # stage count (fractional when all stages are ideal)
    event_stages: list          # stage where each feed / side draw is located


def step_stages(vle, model, murphree=None, forced=None, n_fixed=None, ideal_first=False,
                ideal_last=False, max_stages=300):
    """Step off stages from the top down.

    murphree    Murphree vapor efficiency of the trays (None = ideal stages)
    forced      stage number of each feed / side draw (None = optimal location)
    n_fixed     step exactly this many stages (rating); otherwise stop at x_end
    ideal_first first stage is an equilibrium stage (partial condenser)
    ideal_last  last stage is an equilibrium stage (partial reboiler)
    """
    E = 1.0 if murphree is None else murphree
    x_prev, y = model.x_start, model.y_start
    X, Y = [x_prev], [y]
    stages, events = [], [None] * len(model.switches)
    s = 0
    n_limit = n_fixed if n_fixed is not None else max_stages
    for n in range(1, n_limit + 1):
        y = min(max(y, 0.0), 1.0)
        x_ideal = float(vle.x_eq(y))
        ideal = E >= 1 or (ideal_first and n == 1) or (ideal_last and (
            (n_fixed is None and x_ideal <= model.x_end) or n == n_fixed))
        if ideal:
            x = x_ideal
        else:
            sec = model.sections[s]
            yp = sec.y(vle.x) + E * (vle.y - sec.y(vle.x))
            yp = np.maximum.accumulate(yp)
            x = float(np.interp(y, yp, vle.x))
        X.append(x)
        Y.append(y)
        stages.append((n, x, y, s))
        if n_fixed is None and x <= model.x_end:
            N = n - 1 + (x_prev - model.x_end) / (x_prev - x) if E >= 1 else n
            X.append(x)
            Y.append(model.sections[s].y(x))
            return Stepping(np.array(X), np.array(Y), stages, N, events)
        if n_fixed is None and x >= x_prev - 1e-12:
            raise PinchError("Pinch: the operating line touches the equilibrium curve at x = %.4f" % x)
        while s < len(model.switches) and (
                (forced is None and x < model.switches[s]) or (forced is not None and n >= forced[s])):
            events[s] = n
            s += 1
        y = model.sections[s].y(x)
        X.append(x)
        Y.append(y)
        x_prev = x
    if n_fixed is not None:
        return Stepping(np.array(X), np.array(Y), stages, n_fixed, events)
    raise PinchError("x_end not reached in %d stages; the reflux ratio is too close to the minimum"
                     % max_stages)


# ---------------------------------------------------------------- design
@dataclass
class Design:
    spec: ColumnSpec
    vle: object
    model: ColumnModel
    param_min: float            # R_min (or minimum boil-up ratio)
    param: float                # R (or boil-up ratio)
    pinch: tuple
    pinch_kind: str
    stepping: Stepping
    N_min: float = None
    murphree_stepping: Stepping = None
    overall_efficiency: float = None
    temperatures: list = None   # K, per theoretical stage
    Q_C: float = None           # kW
    Q_R: float = None           # kW
    shortcut: dict = None
    warnings: list = field(default_factory=list)

    @property
    def N(self):
        return self.stepping.N

    @property
    def extra_stages(self):
        """Equilibrium stages that are not trays (reboiler, partial condenser)."""
        s = self.spec
        n = 1 if (s.configuration in ("full", "stripping") and s.reboiler == "partial") else 0
        if s.configuration in ("full", "rectifying") and s.condenser == "partial":
            n += 1
        return n

    @property
    def theoretical_trays(self):
        return self.N - self.extra_stages

    @property
    def actual_trays(self):
        if self.murphree_stepping is not None:
            return self.murphree_stepping.N - self.extra_stages
        if self.overall_efficiency is not None:
            return math.ceil(self.theoretical_trays / self.overall_efficiency - 1e-9)
        return math.ceil(self.theoretical_trays - 1e-9)

    @property
    def feed_trays(self):
        """Actual tray number (from the top) of each feed / side draw."""
        st = self.murphree_stepping or self.stepping
        top = 1 if (self.spec.configuration != "stripping" and self.spec.condenser == "partial") else 0
        trays = [None if e is None else e - top for e in st.event_stages]
        if self.murphree_stepping is None and self.overall_efficiency is not None:
            trays = [None if t is None else math.ceil(t / self.overall_efficiency - 1e-9) for t in trays]
        return trays


def oconnell_efficiency(vle, x_D, x_W):
    """Overall tray efficiency from the O'Connell correlation."""
    x_avg = 0.5 * (x_D + x_W)
    T_avg = 0.5 * (vle.T_or_estimate(x_D) + vle.T_or_estimate(x_W))
    alpha = math.sqrt(vle.alpha(x_D) * vle.alpha(x_W))
    mu = vle.viscosity(x_avg, T_avg)
    return float(min(1.0, 0.492 * (mu * alpha) ** -0.245))


def validate(vle, spec):
    warnings = []
    if spec.configuration == "stripping":
        f = spec.feeds[0]
        if not 0 < spec.x_W < f.z < 1:
            raise ValueError("Need 0 < x_W < z_F < 1")
    elif spec.configuration == "rectifying":
        f = spec.feeds[0]
        if not 0 < f.z < spec.x_D < 1:
            raise ValueError("Need 0 < z_F < x_D < 1")
    else:
        if not 0 < spec.x_W < spec.x_D < 1:
            raise ValueError("Need 0 < x_W < x_D < 1")
        for f in spec.feeds:
            if not spec.x_W < f.z < spec.x_D:
                raise ValueError("%s: need x_W < z < x_D (z = %.3f)" % (f.name, f.z))
    for f in spec.feeds:
        if f.flow <= 0:
            raise ValueError("%s: flow must be positive" % f.name)
    if spec.murphree is not None and not 0 < spec.murphree <= 1:
        raise ValueError("Murphree efficiency must be in (0, 1]")
    x_lo = spec.x_W if spec.x_W is not None else 0
    x_hi = spec.x_D if spec.x_D is not None else 1
    for x_az, T_az in vle.azeotropes():
        msg = "Azeotrope at x = %.4f" % x_az + ("" if T_az is None else ", %.1f degC" % (T_az - 273.15))
        if x_lo < x_az < x_hi:
            raise ValueError(msg + " lies between the product compositions; the separation is impossible "
                             "by ordinary distillation")
        warnings.append(msg)
    return warnings


def design(vle, spec, param_min=None):
    """Design a column. Pass param_min (R_min or minimum boil-up) to skip its search."""
    for f in spec.feeds:
        if f.q is None:
            if f.T_C is None:
                raise ValueError("%s: give q or T_C" % f.name)
            f.q = float(feed_q(vle, f.z, f.T_C))
    warnings = validate(vle, spec)
    cfg = spec.configuration

    if param_min is None:
        param_min, pinch, kind = find_minimum(vle, spec)
    else:
        pinch, kind = pinch_info(vle, spec, param_min)
    if cfg == "stripping":
        param = spec.boilup_ratio if spec.boilup_ratio is not None else spec.boilup_factor * param_min
    else:
        param = spec.reflux_ratio if spec.reflux_ratio is not None else spec.reflux_factor * param_min
    if param <= param_min:
        raise ValueError("%s = %.4g is below the minimum %.4g"
                         % ("Boil-up ratio" if cfg == "stripping" else "Reflux ratio", param, param_min))
    model = ColumnModel(vle, spec, param)

    partial_cond = cfg != "stripping" and spec.condenser == "partial"
    reboiler = cfg in ("full", "stripping") and spec.reboiler == "partial"
    stepping = step_stages(vle, model, max_stages=spec.max_stages)
    d = Design(spec, vle, model, param_min, param, pinch, kind, stepping, warnings=warnings)

    if cfg == "full":
        total = ColumnModel.__new__(ColumnModel)
        total.sections, total.switches = [Section(1.0, 1.0, 0.0, "Total reflux")], []
        total.x_start = total.y_start = spec.x_D
        total.x_end = spec.x_W
        d.N_min = step_stages(vle, total, max_stages=spec.max_stages).N
        if len(spec.feeds) == 1 and not spec.side_draws and spec.reboiler == "partial":
            from .shortcut import fug
            d.shortcut = fug(vle, spec.x_D, spec.x_W, spec.feeds[0], param, model.D, model.W)

    if spec.murphree is not None and spec.murphree < 1:
        d.murphree_stepping = step_stages(vle, model, spec.murphree, ideal_first=partial_cond,
                                          ideal_last=reboiler, max_stages=spec.max_stages)
    elif spec.overall_efficiency is not None:
        if str(spec.overall_efficiency).lower() == "oconnell":
            d.overall_efficiency = oconnell_efficiency(vle, model.x_D, model.x_W)
        else:
            d.overall_efficiency = float(spec.overall_efficiency)

    if vle.has_temperature:
        d.temperatures = [float(vle.T_bubble(x)) for _, x, _, _ in stepping.stages]
    if vle.has_properties:
        top, bottom = model.sections[0], model.sections[-1]
        if cfg != "stripping":
            condensed = top.L if partial_cond else top.V
            d.Q_C = condensed * vle.latent_heat(model.x_D) * KW_PER_KMOLH_KJMOL
        if reboiler:
            d.Q_R = bottom.V * vle.latent_heat(model.x_W) * KW_PER_KMOLH_KJMOL
    return d


def forced_feed_stage_count(vle, spec, param, feed_stage):
    """Theoretical stages when the (single) feed is placed on a given stage; inf if it pinches."""
    model = ColumnModel(vle, spec, param)
    try:
        st = step_stages(vle, model, forced=[feed_stage], max_stages=spec.max_stages)
    except PinchError:
        return float("inf")
    # A feed below the last stage never enters the column
    return st.N if st.event_stages[0] is not None else float("inf")
