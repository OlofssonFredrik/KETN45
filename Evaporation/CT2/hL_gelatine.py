def hL_gelatine(x,T):
    # Calculates enthalpy of two-component liquids assuming ideal mixture
    # Can handle x and T as vectors of same length
    # x      Solids content (weight fraction)
    # T      Temperature
    # CHANGE Cp-VALUES BELOW TO SOLVE DIFFERENT PROBLEM
    #
    Cp_water=4.19  # kJ/kg,K  heat capacity water 
    Cp_solids=4.19  # kJ/kg,K heat capacity for non-volatile substance
    hL=(x*Cp_solids+(1-x)*Cp_water)*T # kJ/kg Enthalpy of liquid
    return hL
