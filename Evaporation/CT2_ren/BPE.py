import math

def BPE(xL):
    # Boiling point elevation for black liquor with mass fraction x
    # TAPPI equation by Ray, Rao, Bansal, & Mohanty 1992
    #    Gray literature, quoted by Bhargava (see below)
    # Curve fit by Bhargava, Khanam, Mohanty & Ray 2008
    #    doi:10.1016/j.compchemeng.2007.10.012
    #
    # This routine can handle vectors as input. 
    # The output will have the same dimensions as the input.

    # To increase stability when routine is used in iterations,
    # take x as zero if x<0 and as 1 if x>1
    xx=xL
    if xL<0:
        xx=0
    elif xL>1:
        xx=1
    
    beta=20*(0.2+xx)**2
    return beta
