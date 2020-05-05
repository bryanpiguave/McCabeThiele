import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import LineString


def intersection(x1_array, y1_array, x2_array, y2_array):
    line1 = LineString([i for i in zip(x1_array, y1_array)])
    line2 = LineString([i for i in zip(x2_array, y2_array)])
    x_inter = line1.intersection(line2).x
    y_inter = line1.intersection(line2).y
    return x_inter, y_inter


archivo = pd.read_excel("Destilacion.xlsx")
Data_x=archivo["X"].values
Data_y=archivo["Y"].values
n =len(Data_x)

Figure =plt.figure(dpi =200)
Figure.suptitle("Destilación benceno tolueno")
Grid =Figure.add_gridspec(1,1)
Graph1 = Figure.add_subplot(Grid[0,0])
Graph1.grid(True)
Graph1.minorticks_on()
Graph1.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)

x = np.linspace(0,max(Data_x),n)
Graph1.plot(x,x,Data_x,Data_y)
Graph1.scatter(Data_x,Data_y,alpha=0.5)


q = 0.45

x_ben_alim = 0.6
x_ben_deriv= 0.95
x_ben_res  = 0.1

#Gráfico de fracciones 
Graph1.plot(x_ben_alim,x_ben_alim,"-bo",
            x_ben_deriv,x_ben_deriv,"-bo",
            x_ben_res,x_ben_res,"-bo")
Graph1.text(0.62,0.55,"$X_F$")
Graph1.text(0.1,0.05,"$X_w$")
Graph1.text(0.95,0.87,"$X_D$")

#Recta q
X_aux = np.linspace(0.3,x_ben_alim,n,endpoint=True)
Y_aux = (q/(q-1))*(X_aux - x_ben_alim)+x_ben_alim


#Intersección
x_inter,y_inter = intersection(X_aux,Y_aux,Data_x,Data_y)


#Gráfica de recta q 
X_q = np.linspace(x_inter,x_ben_alim,n,endpoint=True)
Y_q =(q/(q-1))*(X_q - x_ben_alim)+x_ben_alim
Graph1.plot(X_q,Y_q,label="Recta q")
#Linea de enriquecimiento

m_LOE_min=(y_inter- x_ben_deriv)/(x_inter- x_ben_deriv)
X_LOE_min=np.linspace(x_inter,x_ben_deriv,n)
Y_LOE_min=m_LOE_min*(X_LOE_min-x_ben_deriv)+x_ben_deriv
Graph1.plot(X_LOE_min,Y_LOE_min,alpha=0.8,label="Recta de enriquecimiento")


#Recta de agotamiento
m_LOA_min = (y_inter-x_ben_res)/(x_inter-x_ben_res)
X_LOA_min = np.linspace(x_ben_res,x_inter,n)
Y_LOA_min = m_LOA_min*(X_LOA_min-x_ben_res)+x_ben_res
Graph1.plot(X_LOA_min,Y_LOA_min,alpha=0.9,label="Recta de agotamiento")
Graph1.legend()
Graph1.set_xlabel("Fracción molar de benceno en fase líquida")
Graph1.set_ylabel("Fracción molar de benceno en fase de vapor")


R_D_min = round(-m_LOE_min/(m_LOE_min-1),2)
print("Reflujo mínimo ",R_D_min)
R_D_trabajo = round(1.6*R_D_min,2)
print("Reflujo real ", R_D_trabajo)

#Linea de enriquecimiento real 
x_trabajo = np.linspace(0,x_ben_deriv,n)
y_trabajo = (R_D_trabajo/(R_D_trabajo+1))*(x_trabajo-x_ben_deriv) + x_ben_deriv


x_inter,y_inter =intersection(x_trabajo,y_trabajo,X_aux,Y_aux)
m_LOE = (y_inter-x_ben_deriv)/(x_inter - x_ben_deriv) 
X_LOE = np.linspace(x_inter,x_ben_deriv,n)
Y_LOE = m_LOE*(X_LOE-x_ben_deriv)+x_ben_deriv

m_LOA = (y_inter-x_ben_res)/(x_inter-x_ben_res)
X_LOA = np.linspace(x_ben_res,x_inter,n)
Y_LOA = m_LOA*(X_LOA-x_ben_res)+x_ben_res


Graph1.plot(X_LOE,Y_LOE)
Graph1.plot(X_LOA,Y_LOA)

Figure.savefig("Destilacion.png")