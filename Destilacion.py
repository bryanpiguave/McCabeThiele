import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import LineString, Point

#Estilo de figura: Helvetica (o Arial si no está instalada), negrita, >= 16 pt
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "Nimbus Sans", "DejaVu Sans"]
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["mathtext.default"] = "bf"
plt.rcParams["font.size"] = 18
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"

#Paleta de colores
COLORES = {"equilibrio": (0.05, 0.35, 0.75), "q": (0.05, 0.80, 0.60),
           "enriquecimiento": (0.90, 0.30, 0.05), "agotamiento": (0.90, 0.50, 0.70),
           "etapas": (0.05, 0.05, 0.05)}


def intersection(x1_array, y1_array, x2_array, y2_array):
    line1 = LineString([i for i in zip(x1_array, y1_array)])
    line2 = LineString([i for i in zip(x2_array, y2_array)])
    inter = line1.intersection(line2)
    if not isinstance(inter, Point) or inter.is_empty:
        raise ValueError("Se esperaba una única intersección, se obtuvo: " + inter.geom_type)
    return inter.x, inter.y


def escalonar(linea_operacion, x_D, x_W, x_cambio=None, max_etapas=100):
    """Escalona etapas desde x_D hasta x_W entre la curva de equilibrio y la
    línea de operación. Devuelve los puntos de la escalera, el número de
    etapas teóricas (con fracción de la última) y la etapa de alimentación."""
    X_esc, Y_esc = [x_D], [x_D]
    x, y = x_D, x_D
    etapa_alim = None
    for etapa in range(1, max_etapas + 1):
        #Horizontal hasta la curva de equilibrio
        x_nuevo = np.interp(y, Data_y, Data_x)
        X_esc.append(x_nuevo)
        Y_esc.append(y)
        if x_nuevo <= x_W:
            fraccion = (x - x_W) / (x - x_nuevo)
            X_esc.append(x_nuevo)
            Y_esc.append(x_nuevo)
            return np.array(X_esc), np.array(Y_esc), etapa - 1 + fraccion, etapa_alim
        if x_cambio is not None and etapa_alim is None and x_nuevo < x_cambio:
            etapa_alim = etapa
        #Vertical hasta la línea de operación
        x, y = x_nuevo, linea_operacion(x_nuevo)
        X_esc.append(x)
        Y_esc.append(y)
    raise RuntimeError("No se alcanzó x_W en %d etapas; el reflujo está muy cerca del mínimo" % max_etapas)


archivo = pd.read_csv("Destilacion.csv")
Data_x=archivo["X"].values
Data_y=archivo["Y"].values
n =len(Data_x)

Figure, Graph1 = plt.subplots(figsize=(10,9))
Graph1.grid(True, alpha=0.5)
Graph1.minorticks_on()
Graph1.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)

Graph1.plot([0,1],[0,1],color="gray",linewidth=1.5)
Graph1.plot(Data_x,Data_y,color=COLORES["equilibrio"],linewidth=2.5,label="Curva de equilibrio")


q = 0.45
factor_reflujo = 1.6

x_ben_alim = 0.6
x_ben_deriv= 0.95
x_ben_res  = 0.1

#Gráfico de fracciones
Graph1.plot([x_ben_alim,x_ben_deriv,x_ben_res],
            [x_ben_alim,x_ben_deriv,x_ben_res],"o",color="k",markersize=9,zorder=5)
Graph1.annotate("$X_F$",(x_ben_alim,x_ben_alim),textcoords="offset points",xytext=(8,-22),fontsize=20)
Graph1.annotate("$X_w$",(x_ben_res,x_ben_res),textcoords="offset points",xytext=(8,-22),fontsize=20)
Graph1.annotate("$X_D$",(x_ben_deriv,x_ben_deriv),textcoords="offset points",xytext=(8,-22),fontsize=20)

#Recta q (vertical si q = 1)
if np.isclose(q,1):
    X_aux = np.array([x_ben_alim,x_ben_alim])
    Y_aux = np.array([0.0,1.0])
else:
    X_aux = np.array([0.0,1.0])
    Y_aux = (q/(q-1))*(X_aux - x_ben_alim)+x_ben_alim


#Intersección de la recta q con la curva de equilibrio (punto de pinch)
x_inter,y_inter = intersection(X_aux,Y_aux,Data_x,Data_y)


#Gráfica de recta q
Graph1.plot([x_ben_alim,x_inter],[x_ben_alim,y_inter],color=COLORES["q"],linewidth=2.5,label="Recta q")

#Reflujo mínimo
m_LOE_min=(y_inter- x_ben_deriv)/(x_inter- x_ben_deriv)
R_D_min = m_LOE_min/(1-m_LOE_min)
R_D_trabajo = factor_reflujo*R_D_min

Graph1.plot([x_inter,x_ben_deriv],[y_inter,x_ben_deriv],"--",color=COLORES["enriquecimiento"],linewidth=2,alpha=0.6,
            label="Enriquecimiento ($R_{min}$)")
Graph1.plot([x_ben_res,x_inter],[x_ben_res,y_inter],"--",color=COLORES["agotamiento"],linewidth=2,alpha=0.6,
            label="Agotamiento ($R_{min}$)")


#Linea de enriquecimiento real
m_LOE = R_D_trabajo/(R_D_trabajo+1)
def recta_enriquecimiento(x):
    return m_LOE*(x-x_ben_deriv) + x_ben_deriv

#Intersección de la recta de enriquecimiento con la recta q
if np.isclose(q,1):
    x_inter = x_ben_alim
else:
    m_q = q/(q-1)
    x_inter = ((m_LOE-1)*x_ben_deriv + (1-m_q)*x_ben_alim)/(m_LOE-m_q)
y_inter = recta_enriquecimiento(x_inter)

#Recta de agotamiento real
m_LOA = (y_inter-x_ben_res)/(x_inter-x_ben_res)
def recta_agotamiento(x):
    return m_LOA*(x-x_ben_res) + x_ben_res

def recta_operacion(x):
    return recta_enriquecimiento(x) if x >= x_inter else recta_agotamiento(x)


Graph1.plot([x_inter,x_ben_deriv],[y_inter,x_ben_deriv],color=COLORES["enriquecimiento"],linewidth=2.5,label="Recta de enriquecimiento")
Graph1.plot([x_ben_res,x_inter],[x_ben_res,y_inter],color=COLORES["agotamiento"],linewidth=2.5,label="Recta de agotamiento")


#Escalonamiento de etapas
X_esc,Y_esc,etapas,etapa_alim = escalonar(recta_operacion,x_ben_deriv,x_ben_res,x_cambio=x_inter)
Graph1.plot(X_esc,Y_esc,color=COLORES["etapas"],linewidth=1.5,label="Etapas (%.1f)" % etapas)

#Etapas mínimas a reflujo total
_,_,etapas_min,_ = escalonar(lambda x: x,x_ben_deriv,x_ben_res)


print("Reflujo mínimo           %.2f" % R_D_min)
print("Reflujo real (%.1f Rmin)  %.2f" % (factor_reflujo,R_D_trabajo))
print("Etapas teóricas          %.2f" % etapas)
print("Etapa de alimentación    %d" % etapa_alim)
print("Etapas mínimas (Nmin)    %.2f" % etapas_min)


Graph1.legend(loc="lower right",frameon=False,prop={"weight":"bold","size":16})
Graph1.set_xlabel("Fracción molar de benceno en fase líquida",fontsize=22)
Graph1.set_ylabel("Fracción molar de benceno en fase de vapor",fontsize=22)
Graph1.set_xlim(0,1)
Graph1.set_ylim(0,1.02)

#Sin título en la figura: el nombre del archivo describe el gráfico
nombre = "mccabe_thiele_benceno_tolueno_etapas_teoricas"
for ext in ("pdf","png","svg"):
    Figure.savefig(nombre + "." + ext, dpi=300, bbox_inches="tight")
