from math import*
import numpy as np
import matplotlib.pyplot as plt

""""
This script will try to calculate the critical gust loads.
Plot the (delta-n , t) graph and the corresponding gust loading diagram.

You are now entering the realm where inconsistent units feast upon aspiring engineers, 
here dreams are smashed into pieces within the blink of an eye. No Second changes. 
If you decide to proceed may the force be with you
 
STRUCTURE
-constants
 *nature
 *aircraft parameter

-independent calculations
-dependent calculations/loops
-selection of lists
-printing results
-plotting

Note all speeds are TAS unless otherwise specified

"""
#-----------------------------------------------------------------------------------------------------------------------
R     =  8.314510           # gasconstant [J/(mol*K)]

M_air =  0.0289645          # Molair mass air [kg/mol]

g     =  9.80665            # gravitational acceleration [m/s^2]

rho_0 = 1.225               # density at sealevel [kg/m^3]

gamma = 1.4                 # adiabatic index of air [-]

M_C = 0.77                  # Design cruise mach number[-]

altitude_cruise = 10668     # [m]

S = 57.7                    # surface area wing [m^2]

MAC = 2.63                  # Mean aerodynamic chord [m]

MTOW = 30502 * g            # [N]

MLW = 29142.1 * g           # ASSUMED TO BE MTOW - 0.7*FUELWEIGHT!![N]

MZFW = MTOW - 4533 * g

OEW = 21963 * g             # Operational empty weight [N]

CL_clean = 1.1              # CL clean config [-]

dCL = 1.08

CL_a = 0.0762213            # CL_a at M=0 from excel calc with prandtl[-]

Zmo  = 12000                # Assumed Ceiling [m]
#----------------------------------------------------------------------------------------------------------

# INDEPENDENT COMPUTATIONS
R_air = R / M_air           # shorthand convention [J/kg*K]
V_C   = 228.31 * sqrt(0.3795655/rho_0)

# Computation of flight profile alleviation factor
R1 = MLW / MTOW
R2 = MZFW / MTOW
F_gz = 1 - (Zmo / 76200)
F_gm = sqrt(R2 * tan(pi * R1 / 4))
F_g = 0.5 * (F_gz + F_gm)
#------------------------------------------------------------------------------------------------------------
#list of weights considered
#Weights = [OEW,MTOW,MZFW]

"""If list of only weights is made it says Kt = 0 at the second indexed weight....
No idea why this is """
Weights = []
Weights.append(MTOW)
Weights.append(OEW)
Weights.append(MZFW)
#Gust gradients needed to be considered [m] in increments of 1 m/s
gustgrad = np.arange(9, max(107, 12.5 * MAC) + 1, 1)
# list to find iteration with biggest dn
maxdnvb = []
maxdnvc = []
maxdnvd = []
for h in (1, 2, 3):

    if h == 1:  # Sealevel conditions

        T = 288.15  # Temperature in [K]

        rho = 1.225  # density [kg/m^3]

        altitude = 0  # altitude in [m]

    elif h == 2:  # intermediate altitute 2000ft GENERICLY CHOSEN!!

        T = 284.1876

        rho = 1.1550945

        altitude = 609.6

    elif h == 3:  # cruise altitude conditions 35000ft

        T = 218.808

        rho = 0.3795655

        altitude = 10668  # altitude in [m]

    a = sqrt(gamma * R_air * T) # speed of sound

    V_Ctas = V_C * sqrt(rho_0/rho)

    V_D1 = V_C / 0.8
    V_D2 = (M_C / 0.8) * a* sqrt(rho / rho_0)
    V_D3 = (M_C + 0.05) * a * sqrt(rho / rho_0)
    if max(V_D1, V_D2 )*sqrt(rho_0/rho) >= 0.95 * a:
        V_D = V_D3
    else:                              # diving velocity if not bound by compress effects [m/s]
        V_D = max(V_D2, V_D1)
    V_Dtas = V_D * sqrt(rho_0/rho)
    velocity = np.arange(0.01,round(V_Dtas)+1,1)

    for W in Weights:

        # computation stall speed [m/s]
        V_S1 = sqrt(2 * W / (rho * (CL_clean) * S))

        for V in velocity:

            # Prandt-Glauert correction
            M = V / a
            CL_a = CL_a / sqrt(1 - M ** 2)

            # Calculation of design speed of maximum gust intensity
            mu = (2 * W / S) / (rho * MAC * CL_a * g)
            K_G = (0.88 * mu) / (5.3 + mu)

            if h < 4572:
                Uref = 17.07 - (h * 3.66 / 4572)
            elif h >= 4572 and h < 18288:
                Uref = 13.41 - ((h - 4572) * 7.05 / 13.716)
            else:
                Uref = 6.36

            if round(V) == round(V_Dtas):
                Uref = 0.5 * Uref

            V_B = V_S1 * sqrt(1 + (K_G * rho_0 * Uref * V_C * CL_a) / (2 * W / S))

            if V >= V_B:
                for  H in gustgrad:
                    #Computation of design velocity
                    Uds = Uref * F_g * (H/107)**(1/6)

                    #Computation of load factor change
                    omega = pi * V/H
                    Kt =  2 * W/(S * CL_a * rho * V * g)
                    if Kt != 0 :

                        t = np.linspace(0,(2 * pi/omega),100)                          # 100 Time increments over the interval
                        #dnp = (Uds/(2 * g)) * (omega * np.sin( omega * t) + ((np.exp(-t/Kt)/Kt) - (np.cos(omega *t)/Kt)
                                                                         #- omega * np.sin(omega * t))/(1 + (omega * Kt)**(-2)))
                        dnp = Uds / (2 * g) * (omega * np.sin(omega * t) + (
                                    np.exp(-t / Kt) / Kt - np.cos(omega * t) / Kt - omega * np.sin(omega * t)) / (
                                                             1 + (omega * Kt) ** -2))

                        if  round(V) == round(V_B):
                            maxdnvb.append([max(abs(dnp)),altitude,W,H,V])
                        elif round(V) == round(V_Ctas):
                            maxdnvc.append([max(abs(dnp)),altitude,W,H,V])
                        elif round(V) == round(V_Dtas):
                            maxdnvd.append([max(abs(dnp)),altitude,W,H,V])

                    else:
                        print(Kt,W,Weights)
#-------------------------------------------------------------------------------------------------------------------------
#plots to be produced

V_Blist = []
for i in range(0,len(maxdnvb)):
    V_Blist.append(maxdnvb[i][0])
indexvb = V_Blist.index(max(V_Blist))
V_Bm = maxdnvb[indexvb]
print("max dn at VB\\ altitude = ",V_Bm[1],"\\ Weight = ",V_Bm[2],"\\ H = ",V_Bm[3],"\\V_B = ",V_Bm[4])

V_Clist = []
for i in range(0,len(maxdnvc)):
    V_Clist.append(maxdnvc[i][0])
indexvc = V_Clist.index(max(V_Clist))
V_Cm = maxdnvc[indexvc]
print("max dn at VC\\ altitude = ",V_Cm[1],"\\ Weight = ",V_Cm[2],"\\ H = ",V_Cm[3],"\\V_C = ",V_Cm[4])

V_Dlist = []
print(len(maxdnvd))
for i in range(0,len(maxdnvd)):
    V_Dlist.append(maxdnvd[i][0])
indexvd = V_Dlist.index(max(V_Dlist))
V_Dm = maxdnvd[indexvd]
print("max dn at VD\\ altitude = ",V_Dm[1],"\\ Weight = ",V_Dm[2],"\\ H = ",V_Dm[3],"\\V_D = ",V_Dm[4])