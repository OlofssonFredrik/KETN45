import math
from H_steam import H_steam

def hL_black_liquor(xL,TL):
    # NOTE to students: You are NOT supposed to change anything in this file
    # Enthalpy hL (kJ/kg) of black liquor
    # with mass fraction x at temperature TL (Celcius)
    #
    # This routine can handle vectors as input, in which case x and TL
    # must have the same dimensions. The output will have the same dimensions
    # as the input.

    # To increase stability when routine is used in iterations,
    # take x as zero if x<0, take x as 1 if x>1
    xx=xL
    if xL<0:
        xx=0
    elif xL>1:
        xx=1

    Cp_solids=1.926
    h_solids=Cp_solids*xx*(TL-0)

    h_water=0

    [dummy1,dummy2,h_water]=H_steam(TL,-1)
    h_water=max(0,h_water)
    hL=h_water*(1-xx)+h_solids
    return hL


