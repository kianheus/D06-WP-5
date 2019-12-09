from math import *
#--------------------Wing parameters--------------------

#Geometry
b   = 24.9                      #Span [m]
S   = 69.3                      #Surface [m^2]
C_root =  4.28                  #Root chord [m]
C_tip  =  1.28                  #Tip chord [m]
taper  = C_tip/C_root           #Taper ratio

#Sweep angle formulas:

#Quarter chord sweep [enter in degrees]
sweep_quarter_chord = 28.8
sweep_quarter_chord = sweep_quarter_chord * pi/180

#Leading edge sweep [radians]
sweep_LE = atan(tan(sweep_quarter_chord) + 0.25*2*C_root*(1-taper)/b)

#Half chord sweep [radians]
sweep_half_chord = atan(tan(sweep_LE) - 0.5*2*C_root*(1-taper)/b)             

#---------------------C_Lalpha---------------------------

#Constants
Cla = 0.1121                    #Airfoil Cl-alpha slope [1/degrees]
M   = 0.77                      #Mach number
n   = 0.95                      #Airfoil efficiency factor

#Formulas
AR = b*b/S                      #Aspect ratio
B = sqrt(1-M*M)                 #Mach number correction

#Wing lift gradient
CLa = round((Cla*AR)/(2+sqrt(4+(AR*B/n)**2*1+tan(sweep_half_chord)**2/(B*B))),4)

#---------------------CL max-----------------------------

#Constants
CLClratio   = 0.8               #CLmax/Clmax ratio
Clmax       = 1.397             #Clmax airfoil

#Wing CLmax
CLmax = round(CLClratio * Clmax,4)

#---------------------Stall angle------------------------

#Constants
a0L = -4.49                     #Zero lift angle of attack [degrees]
alpha_CLclean = 2               #Increase in alpha at CLmax [degrees]

#Stall angle of attack
alpha_stall = round(CLmax/CLa + a0L + alpha_CLclean,4)

#---------------------Mach drag divergence --------------

#Constants
k_a = 0.935                                 #Factor for supercritical airfoils
t_c_stream = 0.1*cos(sweep_quarter_chord)   #Thickness to chord ratio in stream direction
CLdes = 0.46                                #CLdesign

#Drag divergence Mach number
M_dd = round(k_a/cos(sweep_LE) - t_c_stream/(cos(sweep_LE))**2 - CLdes/(10*(cos(sweep_LE))**3),4)

#--------------------------------------------------------

print("CLa: ",CLa,", CLmax: ",CLmax,", Stall AoA: ",alpha_stall,", Drag divergence Mach: ",M_dd)
print("Ready")
