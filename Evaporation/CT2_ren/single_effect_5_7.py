# ? Comments that begin with a question mark are hints to students using
# ? this code as a starting point for calculating something different
import math
import numpy as np
from scipy.optimize import root
from hL_gelatine import hL_gelatine #? Delete this line when working with CT2
from H_steam import H_steam
#? from hL_black_liquor import hL_black_liquor
#? from BPE import BPE
#? from k_black_liquor import k_black_liquor


def Ts_from_Ps_water(Ps):
    # ? Note: In task 5.7 we were given the steam pressure
    # ? You wont be needing Antoines equation in the compulsory task
    # Calculation of vapor temperature (saturated steam)
    Ts=1730.630/(-math.log10(Ps)+10.19625)-233.426 # Antoine equation, see handbook
    return Ts


def evaporator(X,F,xF,Tf,Ps,TL,xL,k_constant,use_BPE,vary_k):
    #Calculates residuals of mass and energy balances for an evaporator
    #
    #    Input: X contains guesses of operating conditions in this order
    #       condition      unit     explanation
    #       S,             kg/s     Steam flux
    #       V,             kg/s     Vapor flux (pure water)
    #       L,             kg/s     Liquid flux (pure water)
    #       A               m2      Heat exchanger area
    #
    #    Other input (known parameters)
    #       F              kg/s     Feed flow
    #       xF             kg/kg    Solids fraction in feed
    #       Tf             Celsius  Feed temperature
    #       Ps             Pa       Fresh steam pressure
    #       TL             Celsius  Product stream temperature
    #       xL             kg/kg    Solids fraction in product stream
    #       k_constant     J/m2,K   Overall heat transfer coefficient
    #       use_BPE        boolean  Take boiling point elevation into account?
    #       vary_k         boolean  Let k depend on e.g. T, x?
    #
    #    Output: residuals that should be zero if all data in X are
    #    internally consistent
    #
    #    Assumptions:
    #       * xV=0, i.e. vapor flow from evaporator is pure water
    #       * heat capacity for water is constant within the temperature range
    #       * involatile substance is a liquid with a constant heat capacity
    #       * no mixing enthalpy, i.e. enthalpy of flows containing water and
    #         the involatile substance can be calculated using heat capacities
    #         for pure water and pure involatile substance



    # ================ Initialisation ========================= 
    #
    # Copy X-values to variables that are easier to understand
    [S,V,L,A]=X

    # =============== Calculations start here ==================
    #
    if use_BPE:
        raise Exception('Boiling point elevation not implemented')
        # ? Delete line above and write some code here to implement
        # ? handling of boiling point elevation

    # ? Note: In this task we were given the steam pressure
    # ? You wont be needing Antoines equation in the compulsory task
    # Calculation of vapor temperature (saturated steam)
    Ts=Ts_from_Ps_water(Ps)

    # Calculation of water vapor enthalpies 
    # ? NOTE: To calculate H for overheated steam, 2nd argument cannot be -1
    [Hs,dummy1,dummy2]=H_steam(Ts,-1) 
    [HV,dummy1,dummy2]=H_steam(TL,-1)

    #Calculation of other enthalpies
    # ? Note that hLGelatine calculates the liquid enthalpy for gelatine
    # ? solutions. To do calculations for black liquor, you will need to call
    # ? hLBlackLiquor instead
    hf=hL_gelatine(xF,Tf)
    hL=hL_gelatine(xL,TL)
    hk=hL_gelatine(0,Ts)
        
    if vary_k:
        raise Exception('Varying k-values not implemented')
        # ? Delete line above and write some code here to implement
        # ? handling of boiling point elevation
    else:
        k=k_constant

    # ================ Calculating residuals: ==================
    #
    # 
    # evaporator 1
    Y=X*0
    Y[0]=V+L-F # MB total
    Y[1]=F*xF-L*xL # MB solids
    Y[2]=S*(Hs-hk)+F*hf-V*HV-L*hL # EB over evaporator 
    Y[3]=S*(Hs-hk)-k*A*(Ts-TL) # EB over heat exchanger
    # Note: Which equation is used for which position in Y does not matter
    # You may choose an order you think makes the code the easiest to read
    return Y

# ? Current data: excercise 5.7 in course KETF10
F=1e3/3600# kg/s                     Feed flux
xF=0.05   # kg dry matter/kg total   Feed dry matter content 
xL_out=0.25# kg dry matter/kg total   Dry matter content from last evaporator
Tf=80     # degrees C                Feed temperature
Ps=200e3  # Pa                       Pressure of steam used
T=99.61  #  degrees C      Temperature for vapor flow from last evaporator
k=1.6    #  kW/m2/K                  Overall heat transfer coeff
use_BPE=False # Should boiling point elevation be taken into account?
vary_k=False # Should k-value depend on temperature etc.?

# Note: In Python a so called tuple can contain different data types
# thus we can have a boolean in the same tuple as some floats
known=(F,xF,Tf,Ps,T,xL_out,k,use_BPE,vary_k) 
guess=np.array([F, F, F, 100]) #Set initial guesses

sol = root(evaporator, guess, args=known, method='hybr')

if not sol.success:
    print('Iteration not successful:',sol.message)
else:
    print('Iteration successful:',sol.message)
# Translate results to variable names easier to understand
[S, V, L, A]=sol.x[:]

# You can use the command "print" to print results
print('Steam flux',S,'kg/s' )
print('Vapor flux',V,'kg/s' )
print('Liquid flux',L,'kg/s' )
print('Area',A,'m2' )
print('S/Vtot',S/V)
print('Vtot/S',V/S)

# Recalculate other results
# ? Note: A common beginners method to retrieve these kind of values is
# ? to print out values inside the routine called by root
# ? There are two problems why that should be avoided:
# ? 1. There is in general no guarante that the solver, in its very
# ?    last call of the function, will use the correct solution
# ? 2. Printing out things to the screen every time a function is called
# ?    by the solver slows down the calculation considerably
Ts=Ts_from_Ps_water(Ps)
[Hs,P,dummy]=H_steam(Ts,-1)
[HV,P,dummy]=H_steam(T,-1)
print('Ts',Ts,'degree C')
print('Steam enthalpy Hs',Hs,'kJ/kg')

#? Variables in Python can easily be printed to files
#? Note: the "w" argument WRITEs over the file if it already exists
#? If you instead want to APPEND your new results to end of file, change to "a"
with open("output.txt", "w") as f:
  print(S,V,L,A,Ts,Hs,HV, file=f)
