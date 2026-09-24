import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import LineString, Point

#Figure style: Helvetica (or Arial if not installed), bold, >= 16 pt
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "Nimbus Sans", "DejaVu Sans"]
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["mathtext.default"] = "bf"
plt.rcParams["font.size"] = 18
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"

#Color palette
COLORS = {"equilibrium": (0.05, 0.35, 0.75), "q": (0.05, 0.80, 0.60),
          "rectifying": (0.90, 0.30, 0.05), "stripping": (0.90, 0.50, 0.70),
          "stages": (0.05, 0.05, 0.05)}


def intersection(x1_array, y1_array, x2_array, y2_array):
    line1 = LineString([i for i in zip(x1_array, y1_array)])
    line2 = LineString([i for i in zip(x2_array, y2_array)])
    inter = line1.intersection(line2)
    if not isinstance(inter, Point) or inter.is_empty:
        raise ValueError("Expected a single intersection point, got: " + inter.geom_type)
    return inter.x, inter.y


def step_stages(operating_line, x_D, x_W, x_switch=None, max_stages=100):
    """Step off stages from x_D down to x_W between the equilibrium curve and
    the operating line. Returns the staircase points, the number of
    theoretical stages (including the fraction of the last one) and the
    feed stage."""
    X_steps, Y_steps = [x_D], [x_D]
    x, y = x_D, x_D
    feed_stage = None
    for stage in range(1, max_stages + 1):
        #Horizontal to the equilibrium curve
        x_new = np.interp(y, Data_y, Data_x)
        X_steps.append(x_new)
        Y_steps.append(y)
        if x_new <= x_W:
            fraction = (x - x_W) / (x - x_new)
            X_steps.append(x_new)
            Y_steps.append(x_new)
            return np.array(X_steps), np.array(Y_steps), stage - 1 + fraction, feed_stage
        if x_switch is not None and feed_stage is None and x_new < x_switch:
            feed_stage = stage
        #Vertical to the operating line
        x, y = x_new, operating_line(x_new)
        X_steps.append(x)
        Y_steps.append(y)
    raise RuntimeError("x_W not reached in %d stages; the reflux ratio is too close to the minimum" % max_stages)


data = pd.read_csv("distillation.csv")
Data_x=data["X"].values
Data_y=data["Y"].values
n =len(Data_x)

Figure, Graph1 = plt.subplots(figsize=(10,9))
Graph1.grid(True, alpha=0.5)
Graph1.minorticks_on()
Graph1.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)

Graph1.plot([0,1],[0,1],color="gray",linewidth=1.5)
Graph1.plot(Data_x,Data_y,color=COLORS["equilibrium"],linewidth=2.5,label="Equilibrium curve")


q = 0.45
reflux_factor = 1.6

x_benz_feed = 0.6
x_benz_dist = 0.95
x_benz_bott = 0.1

#Plot of the compositions
Graph1.plot([x_benz_feed,x_benz_dist,x_benz_bott],
            [x_benz_feed,x_benz_dist,x_benz_bott],"o",color="k",markersize=9,zorder=5)
Graph1.annotate("$X_F$",(x_benz_feed,x_benz_feed),textcoords="offset points",xytext=(8,-22),fontsize=20)
Graph1.annotate("$X_w$",(x_benz_bott,x_benz_bott),textcoords="offset points",xytext=(8,-22),fontsize=20)
Graph1.annotate("$X_D$",(x_benz_dist,x_benz_dist),textcoords="offset points",xytext=(8,-22),fontsize=20)

#q-line (vertical if q = 1)
if np.isclose(q,1):
    X_aux = np.array([x_benz_feed,x_benz_feed])
    Y_aux = np.array([0.0,1.0])
else:
    X_aux = np.array([0.0,1.0])
    Y_aux = (q/(q-1))*(X_aux - x_benz_feed)+x_benz_feed


#Intersection of the q-line with the equilibrium curve (pinch point)
x_inter,y_inter = intersection(X_aux,Y_aux,Data_x,Data_y)


#Plot of the q-line
Graph1.plot([x_benz_feed,x_inter],[x_benz_feed,y_inter],color=COLORS["q"],linewidth=2.5,label="q-line")

#Minimum reflux
m_ROL_min=(y_inter- x_benz_dist)/(x_inter- x_benz_dist)
R_D_min = m_ROL_min/(1-m_ROL_min)
R_D_operating = reflux_factor*R_D_min

Graph1.plot([x_inter,x_benz_dist],[y_inter,x_benz_dist],"--",color=COLORS["rectifying"],linewidth=2,alpha=0.6,
            label="Rectifying ($R_{min}$)")
Graph1.plot([x_benz_bott,x_inter],[x_benz_bott,y_inter],"--",color=COLORS["stripping"],linewidth=2,alpha=0.6,
            label="Stripping ($R_{min}$)")


#Actual rectifying line
m_ROL = R_D_operating/(R_D_operating+1)
def rectifying_line(x):
    return m_ROL*(x-x_benz_dist) + x_benz_dist

#Intersection of the rectifying line with the q-line
if np.isclose(q,1):
    x_inter = x_benz_feed
else:
    m_q = q/(q-1)
    x_inter = ((m_ROL-1)*x_benz_dist + (1-m_q)*x_benz_feed)/(m_ROL-m_q)
y_inter = rectifying_line(x_inter)

#Actual stripping line
m_SOL = (y_inter-x_benz_bott)/(x_inter-x_benz_bott)
def stripping_line(x):
    return m_SOL*(x-x_benz_bott) + x_benz_bott

def operating_line(x):
    return rectifying_line(x) if x >= x_inter else stripping_line(x)


Graph1.plot([x_inter,x_benz_dist],[y_inter,x_benz_dist],color=COLORS["rectifying"],linewidth=2.5,label="Rectifying line")
Graph1.plot([x_benz_bott,x_inter],[x_benz_bott,y_inter],color=COLORS["stripping"],linewidth=2.5,label="Stripping line")


#Stage stepping
X_steps,Y_steps,stages,feed_stage = step_stages(operating_line,x_benz_dist,x_benz_bott,x_switch=x_inter)
Graph1.plot(X_steps,Y_steps,color=COLORS["stages"],linewidth=1.5,label="Stages (%.1f)" % stages)

#Minimum stages at total reflux
_,_,min_stages,_ = step_stages(lambda x: x,x_benz_dist,x_benz_bott)


print("Minimum reflux ratio         %.2f" % R_D_min)
print("Operating reflux (%.1f Rmin)  %.2f" % (reflux_factor,R_D_operating))
print("Theoretical stages           %.2f" % stages)
print("Feed stage                   %d" % feed_stage)
print("Minimum stages (Nmin)        %.2f" % min_stages)


Graph1.legend(loc="lower right",frameon=False,prop={"weight":"bold","size":16})
Graph1.set_xlabel("Benzene mole fraction in liquid phase",fontsize=22)
Graph1.set_ylabel("Benzene mole fraction in vapor phase",fontsize=22)
Graph1.set_xlim(0,1)
Graph1.set_ylim(0,1.02)

#No title on the figure: the file name describes the plot
filename = "mccabe_thiele_benzene_toluene_theoretical_stages"
for ext in ("pdf","png","svg"):
    Figure.savefig(filename + "." + ext, dpi=300, bbox_inches="tight")
