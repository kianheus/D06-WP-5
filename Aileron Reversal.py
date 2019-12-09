from math import *
import numpy as np
from matplotlib import pyplot as plt
# Given by the previous iterations
C_r = 3.695625413           # Root chord            [m]
C_t = 1.107118055           # Tip chord             [m]
b = 24.01371734             # Wing span             [m]
Taper= C_t / C_r            # Taper ratio
S =3.620008069175999        # Aileron effective area [m^2]
G = 27E9
################################################
# Inputs
t_TopPlate = 0.01
t_BottomPlate= 0.01
t_SidePlates = 0.01
# Assuming hat  stringers
L_Stringer_TopBar =0.004
L_Stringer_L_Shape = 0.004
h_Stringer=0.004
t_Stringer =0.002
#Stringer properties
h_L_Shape= h_Stringer - t_Stringer
Area_L_Shape = (h_L_Shape*t_Stringer) + ((L_Stringer_L_Shape-t_Stringer)*t_Stringer)
Area_TopBar = (t_Stringer)*L_Stringer_TopBar
Area_Stringer = (2*Area_L_Shape)+Area_TopBar
Centroid_Stringer_2 = (L_Stringer_L_Shape * t_Stringer * (L_Stringer_L_Shape/2))+((h_L_Shape-t_Stringer)*t_Stringer*(L_Stringer_L_Shape-t_Stringer+(t_Stringer/2))) + \
                    (Area_TopBar*((L_Stringer_TopBar/2)+(L_Stringer_L_Shape-t_Stringer))) +\
                    ((h_L_Shape-t_Stringer)*t_Stringer*(L_Stringer_TopBar+(t_Stringer/2)+(L_Stringer_L_Shape-t_Stringer))) + \
                    (L_Stringer_L_Shape*t_Stringer * (L_Stringer_TopBar+(L_Stringer_L_Shape-t_Stringer)+((L_Stringer_L_Shape/2)-t_Stringer)))
Centroid_Stringer = Centroid_Stringer_2/Area_Stringer
#Equation inputs
rho = 0.37956556562264265
dCl_dEpsilon = 6.6          #per rad 
dCm_dEpsilon = -3          #per rad
dCl_dAlpha =6.66                   #per rad
Epsilon= 3
##############################################################3
y_range= np.arange(0.75,0.95,0.01)           # Calculating step range
C_y_range = C_r-(C_r * (1-Taper) *y_range)      # Calculating chord lengths
L_TopPlate = ( 0.6 * C_y_range) - (0.2 * C_y_range) # Calculating the length of the top plate
Area_TopPlate = L_TopPlate * t_TopPlate 
L_BottomPlate= ( 0.6 * C_y_range) - (0.2 * C_y_range)
Area_BottomPlate= L_BottomPlate * t_BottomPlate
H_SidePlate_Left = 0.09065 * C_y_range
H_SidePlate_Right = 0.08116 * C_y_range
h_SidePlate_Average = (H_SidePlate_Left + H_SidePlate_Right)/2
Area_SidePlate_Left = H_SidePlate_Left * t_SidePlates
Area_SidePlate_Right = H_SidePlate_Right * t_SidePlates
Area_SidePlate_Average = (Area_SidePlate_Right +  Area_SidePlate_Left)/2
Total_Area = (4*Area_Stringer) + Area_TopPlate + Area_BottomPlate + Area_SidePlate_Right + Area_SidePlate_Left
# Datum is the side closest to LE
x_cord_HorizontalPlates = L_TopPlate/2 
x_cord_VerticalPlate_Left = t_SidePlates/2
x_cord_VerticalPlate_Right =L_TopPlate - (t_SidePlates/2)
# Assuming hat stringers consisting of 2 L shapes and one Bar
Centroid_WingBox = ((Area_TopPlate * x_cord_HorizontalPlates) + (Area_BottomPlate * x_cord_HorizontalPlates) +(Area_SidePlate_Left * x_cord_VerticalPlate_Left) \
                   + (Area_SidePlate_Right* x_cord_VerticalPlate_Right) + (2 * Area_Stringer * t_SidePlates) + ( 2 * Area_Stringer * ( L_TopPlate - t_SidePlates))\
                   )/ Total_Area
e_Value = (((0.6*C_y_range) + Centroid_WingBox)/C_y_range) -0.25
ce = e_Value* C_y_range

Ixx_TopPlate = ((1/12)* L_TopPlate*(t_TopPlate**3)) + (Area_TopPlate*(h_SidePlate_Average**2))
Ixx_SidePlate = ((1/12)*t_SidePlates* (h_SidePlate_Average**3))
Izz_TopPlate = ((1/12)*t_TopPlate*(L_TopPlate**3))
Izz_SidePlate= ((1/12)*h_SidePlate_Average *(t_SidePlates**3))+(Area_SidePlate_Average*((L_TopPlate/2)**2))
Ixx_Total = (2*Ixx_TopPlate)+(2*Ixx_SidePlate)
Izz_Total = (2*Izz_TopPlate) + ( 2*Ixx_SidePlate)
J = Izz_Total + Ixx_Total
K = G * J
####################################################################
Ratio = 0.5*rho*S
Vr = np.sqrt((-K*dCl_dEpsilon)/(Ratio*C_y_range*dCm_dEpsilon*dCl_dAlpha))                                       #Reversal speed
Delta_L = Ratio*Vr*Epsilon*(((Ratio*C_y_range*dCm_dEpsilon*dCl_dAlpha)+(K*dCl_dEpsilon))/(K-Ratio*Vr*ce*dCl_dAlpha))
plt.plot(y_range,Delta_L)
plt.show()

# Va sea level = 138.6938
# Va cruise = 249.1620
# G = 27 GPa
# Reynolds number = 9415768
# rho cruise = 0.37956556562264265

