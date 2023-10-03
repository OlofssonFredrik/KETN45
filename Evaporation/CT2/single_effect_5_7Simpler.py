
import math
import numpy as np
from scipy.optimize import root
from H_steam import H_steam
from hL_black_liquor import hL_black_liquor
from BPE import BPE
from k_black_liquor import k_black_liquor


def known():

    F=15 # kg/s                     Feed flux
    xF=0.150   # kg dry matter/kg total   Feed dry matter content 
    xL3=0.700 # kg dry matter/kg total   Dry matter content from last evaporator
    Tf=40     # degrees C                Feed temperature
    TL=45.81  #  degrees C      Temperature for vapor flow from last evaporator
    k=1.274    #  kW/m2/K                  Overall heat transfer coeff
    Ts=158.83  # degrees C      Temp of fresh steam
    T3=TL
    return F,xF,Tf,T3,xL3,k,Ts

def evaporator(X):
    #Calculates residuals of mass and energy balances for an evaporator
  

    # ================ Initialisation ========================= 
    #

    [S,V1,V2,V3,L1,L2,L3,A,T1,T2,xL1,xL2]=X
    # 
    [F,xF,Tf,T3,xL3,k,Ts]=known()

    # =============== Calculations start here ==================
    #


    [Hs,dummy1,dummy2]=H_steam(Ts,-1) 
    [HV1,dummy1,dummy2]=H_steam(T1,-1)
    [HV2,dummy1,dummy2]=H_steam(T2,-1)
    [HV3,dummy1,dummy2]=H_steam(T3,-1)


    hf=hL_black_liquor(xF,Tf)
    hL1=hL_black_liquor(xL1,T1)
    hL2=hL_black_liquor(xL2,T2)
    hL3=hL_black_liquor(xL3,T3)
    hk1=hL_black_liquor(0,Ts)
    hk2=hL_black_liquor(0,T1)
    hk3=hL_black_liquor(0,T2)
        

    # ================ Calculating residuals: ==================
    
    Y=X*0
    
    #evaporator 1
    Y[0]=V1+L1-F #tot MB
    Y[1]=L1*xL1 - F*xF #MB solids
    Y[2]=V1*HV1 + L1*hL1 - F*hf - S*(Hs-hk1)    #EB för evaporator 1
    Y[3]=k*A*(Ts-T1) - S*(Hs-hk1)             #EB för heat ex 1
    
    
    #evaporator 2
    Y[4]=V2+L2-L1 #tot MB
    Y[5]=L2*xL2 - L1*xL1 #MB solids
    Y[6]=V2*HV2 + L2*hL2 - L1*hL1 - V1*(HV1-hk2) #EB för evaporator 2
    Y[7]=k*A*(T1-T2) - V1*(HV1-hk2)      #EB för heat ex 2
    

    #evaporator 3
    Y[8]=V3+L3-L2 #tot MB
    Y[9]=L3*xL3 - L2*xL2 #MB solids
    Y[10]=V3*HV3 + L3*hL3 - L2*hL2 - V2*(HV2-hk3)  #EB för evaporator 3
    Y[11]=k*A*(T2-T3) - V2*(HV2-hk3)       #EB för heat ex 3
    
    return Y


[F,xF,Tf,T3,xL3,k,Ts]=known()

guess=np.array([7,7,7,7,7,7,7,200,90,90,0.425,0.425]) 

sol = root(evaporator, guess, method='hybr')

if not sol.success:
    print('Iteration not successful:',sol.message)
else:
    print('Iteration successful:',sol.message)
   
    [S,V1,V2,V3,L1,L2,L3,A,T1,T2,xL1,xL2]=sol.x[:]
   
    print('Steam flux',S,'kg/s' )
    print('Vapor flux for 1',V1,'kg/s' )
    print('Vapor flux for 2',V2,'kg/s' )
    print('Vapor flux for 3',V3,'kg/s' )
    print('Liquid flux for 1',L1,'kg/s' )
    print('Liquid flux for 2',L2,'kg/s' )
    print('Liquid flux for 3',L3,'kg/s' )
    print('Temperature for 1',T1,'Celsius')
    print('Temperature for 2',T2,'Celsius')
    print('Temperature for 3',T3,'Celsius')
    print('Molar fraction for L1',xL1,'kg/kg')
    print('Molar fraction for L2',xL2,'kg/kg')
    print('Molar fraction for L3',xL3,'kg/kg')
    print('Area',A,'m2' )
    print('S/Vtot',S/(V1+V2+V3))
    print('Vtot/S',(V1+V2+V3)/S)

    [Hs,P,dummy]=H_steam(Ts,-1)
    [HV3,P,dummy]=H_steam(T3,-1)
    print('Ts',Ts,'degree C')
    print('Steam enthalpy Hs',Hs,'kJ/kg')
    
    with open("output.txt", "w") as f:
      print(S,V3,L3,A,Ts,Hs,HV3, file=f)