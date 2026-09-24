"""The original benzene-toluene example, now running on the mccabe_thiele package.

For other cases use the command-line tool, e.g.
    python -m mccabe_thiele examples/benzene_toluene.yaml
"""

from mccabe_thiele import ColumnSpec, Feed, design, table_vle
from mccabe_thiele.plotting import plot_mccabe_thiele

x_benz_feed = 0.6
x_benz_dist = 0.95
x_benz_bott = 0.1
q = 0.45
reflux_factor = 1.6

vle = table_vle("examples/benzene_toluene_constant_alpha.csv", components=["benzene", "toluene"])
d = design(vle, ColumnSpec(feeds=[Feed(100.0, x_benz_feed, q)], x_D=x_benz_dist, x_W=x_benz_bott,
                           reflux_factor=reflux_factor))

print("Minimum reflux ratio         %.2f" % d.param_min)
print("Operating reflux (%.1f Rmin)  %.2f" % (reflux_factor, d.param))
print("Theoretical stages           %.2f" % d.N)
print("Feed stage                   %d" % d.stepping.event_stages[0])
print("Minimum stages (Nmin)        %.2f" % d.N_min)

plot_mccabe_thiele(d, ".", "mccabe_thiele_benzene_toluene_theoretical_stages")
