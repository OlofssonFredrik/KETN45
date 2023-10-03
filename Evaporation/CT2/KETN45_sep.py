import math
import numpy as np
from scipy.optimize import root, fsolve
from H_steam import H_steam
from hL_black_liquor import hL_black_liquor
from BPE import BPE
from k_black_liquor import k_black_liquor

def known():
    F=2.5 # kg/s                     Feed flux
    xF=0.150   # kg dry matter/kg total   Feed dry matter content 
    xL1=0.200 # kg dry matter/kg total   Dry matter content from last evaporator
    Tf=95     # degrees C                Feed temperature
    A1=30     #m2
    A2= 45      #m2
    A3=50       #m2
    k1=2500      #W/m2/K                Apparell  Overall heat transfer coeff
    k2=1500    #W/m2/K               Apparell   Overall heat transfer coeff
    k3=1200    # W/m2/K               Apparell   Overall heat transfer coeff
    Ps = 320000 #Pressure of fresh steam in Pa

    return F, xF, Tf, Ps, xL1, k1, k2, k3, A1, A2, A3

def Hv(T):
    H_v = 2496.4+2.26*T-7.34808*10**-3*T**2+3.38602*10**-5*T**3-8.40678*10**-8*T**4
    return H_v
    
def hL(T):
    h_L=4.19*T
    return h_L

def Tsat(Psat):
    Tsat = 3816.44/(18.3036-np.log(Psat/133.32))-227.03
    return Tsat
 

def evaporator(X):
    [S,V1,V2,V3,L1,L2,L3,T1,T2,T3,xL2, xL3]=X

    [F, xF, Tf, Ps, xL1, k1, k2, k3, A1, A2, A3]=known()

    Ts = Tsat(Ps)

    Hs=Hv(Ts)
    HV1=Hv(T1) #V1
    HV2=Hv(T2) #V2
    HV3=Hv(T3) #V3

    hf = hL(Tf) #feed
    hL1 = hL(T1) #L1
    hL2 = hL(T2) #L2
    hL3 = hL(T3) #L3
    
    hk1 = hL(T1) #k1
    hk2 = hL(T2) #k2
    hk3 = hL(T3) #k3


    
    Y=X*0
    
    #evaporator 1
    Y[0]=V1+L1-F #tot MB
    Y[1]=L1*xL1 - F*xF #MB solids
    Y[2]=V1*HV1 + L1*hL1 - F*hf - S*(Hs-hk1)    #EB för evaporator 1
    Y[3]=k1*A1*(Ts-T1) - S*(Hs-hk1)             #EB för heat ex 1
    
    
    #evaporator 2
    Y[4]=V2+L2-L1 #tot MB
    Y[5]=L2*xL2 - L1*xL1 #MB solids
    Y[6]=V2*HV2 + L2*hL2 - L1*hL1 - V1*(HV1-hk2) #EB för evaporator 2
    Y[7]=k2*A2*(T1-T2) - V1*(HV1-hk2)      #EB för heat ex 2
    

    #evaporator 3
    Y[8]=V3+L3-L2 #tot MB
    Y[9]=L3*xL3 - L2*xL2 #MB solids
    Y[10]=V3*HV3 + L3*hL3 - L2*hL2 - V2*(HV2-hk3)  #EB för evaporator 3
    Y[11]=k3*A3*(T2-T3) - V2*(HV2-hk3)       #EB för heat ex 3
    
    return Y


[F, xF, Tf, Ps, xL1, k1, k2, k3, A1, A2, A3]=known()

# guess=np.array([7,7,7,7,7,7,7,200,90,90,0.425,0.425,0.425]) 
guess=np.array([2.5,1,1,1,1,1,1,90,90,90,0.425,0.425]) 

#sol = root(evaporator, guess, method='hybr')
sol = fsolve(evaporator, guess)
print(f"sol: {sol}")


[S,V1,V2,V3,L1,L2,L3,T1,T2,T3,xL2, xL3]=sol
  
print('Steam flux', f'{S:.2f}', 'kg/s')
print('Vapor flux for 1', f'{V1:.2f}', 'kg/s')
print('Vapor flux for 2', f'{V2:.2f}', 'kg/s')
print('Vapor flux for 3', f'{V3:.2f}', 'kg/s')
print('Liquid flux for 1', f'{L1:.2f}', 'kg/s')
print('Liquid flux for 2', f'{L2:.2f}', 'kg/s')
print('Liquid flux for 3', f'{L3:.2f}', 'kg/s')
print('Temperature for 1', f'{T1:.2f}', 'Celsius')
print('Temperature for 2', f'{T2:.2f}', 'Celsius')
print('Temperature for 3', f'{T3:.2f}', 'Celsius')
print('Molar fraction for L2', f'{xL2:.2f}', 'kg/kg')
print('Molar fraction for L3', f'{xL3:.2f}', 'kg/kg')
print('S/Vtot',S/(V1+V2+V3))
print('Vtot/S',(V1+V2+V3)/S)