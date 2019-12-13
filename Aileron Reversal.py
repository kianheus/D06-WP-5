from math import *
import numpy as np
from matplotlib import pyplot as plt
# Given by the previous iterations
C_r = 3.695625413           # Root chord            [m]
C_t = 1.107118055           # Tip chord             [m]
b = 24.01371734             # Wing span             [m]
Taper= C_t / C_r            # Taper ratio
S =57        # Aileron effective area [m^2]
G = 27E9
x=0
################################################
# Inputs
t_TopPlate = 0.004
t_BottomPlate= 0.004
t_SidePlates = 0.004
# Assuming hat  stringers
L_Stringer_TopBar =0.010
L_Stringer_L_Shape = 0.010
h_Stringer=0.010
t_Stringer =0.004
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
dCl_dEpsilon = 2.813795732          #per rad -8.92
dCm_dEpsilon = -0.425707          #per rad 2.254
dCl_dAlpha =6.65802                   #per rad 7.07
Epsilon= 1
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

#Fun stuff
a = h_SidePlate_Average*1000
b=L_TopPlate*1000
t=t_SidePlates*1000
J_Ratio_1=2*(t**2)*((b-2)**2)*((a-t)**2)
J_Ratio_2 = (a*t)+(b*t)-(2*(t**2))
J_1= (J_Ratio_1/J_Ratio_2)*1E-12
J_2= 99819557.3955E-12
K = G*J_1

####################################################################
Ratio = 0.5*rho*S
##Ratio_2=(Ratio*(V**2)*C_y_range*dCm_dEpsilon*dCl_dAlpha)+(K*dCl_dEpsilon)
##Ratio_3= K-(Ratio*ce*(V**2)*dCl_dAlpha)
Vr = np.sqrt((-K*dCl_dEpsilon)/(Ratio*C_y_range*dCm_dEpsilon*dCl_dAlpha))   #Reversal speed
print(Vr[-1])
#Delta_L = (Ratio*(V**2)*Epsilon)*(Ratio_2/Ratio_3)
#Delta_Lr = dCl_dEpsilon *Epsilon *Ratio * (V**2)

##Effectiveness = Delta_L/Delta_Lr
##plt.plot(y_range,Vr)
##plt.show()

################################################################
V = np.arange(0,300,1)
Wierd_Ratio= rho*S*ce*dCl_dAlpha
Vd = np.sqrt((2*K)/Wierd_Ratio)
print(Vd[-1])
Speed_Ratio = (V/(Vr[-1]))**2
Speed_Ratio_2 = (V/(Vd[-1]))**2
Effectiveness = (1-Speed_Ratio)/(1-Speed_Ratio_2)
print(max(Effectiveness))
##print(Effectiveness)
while x <len(Effectiveness):
    if Effectiveness[x]==max(Effectiveness) or Effectiveness[x] == min(Effectiveness):
        print(V[x],Effectiveness[x])
    x=x+1

plt.plot(V,Effectiveness)
plt.show()
# Va sea level = 138.6938
# Va cruise = 249.1620
# G = 27 GPa
# Reynolds number = 9415768
# rho cruise = 0.37956556562264265

