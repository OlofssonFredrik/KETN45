def k_black_liquor(xF,xL,T,F,L):
    ## This routine is developed for educational purposes only
    # NOTE to students: You are NOT supposed to change anything in this file
    #
    # The routine is tailored to be used for Compulsory Task 2 (CT2) in the
    # course KETF10 Separation processes.
    #
    # The use of this routine outside the scope of CT2 is _not_ advised.
    #
    # xF, xL      Solids fraction (weight fraction)
    # T           Temperature (Celcius)
    # F, L        Liquid flow rates in and out respectively (kg/s)
    # k           Apparent overall heat transfer coefficient (kW/m2,K)
    # mu          Dynamic viscosity (Pa s)
    #
    # xF, xL, T, F and L may be scalars or (equally long) vectors
    # k, mu has the same dimension as the input parameters
    #
    # You are recommended to use this routine as a black box, supplying you
    # with overall heat transfer coefficient k to use in CT2.
    # The relationship used for k is based on research articles, but it needs
    # good viscosity values as input.
    # The second output, mu, is NOT based on research articles
    #
    
    #
    # The routine takes the apparent overall heat transfer coefficient to be
    # a function of the mass flow rate per unit width (specific mass flow rate)
    # gamma and the viscosity of the solution as
    #
    # k=a*gamma**b*viscosity**c
    #
    # In Karlsson et al 2013 a=558, b=0.073, c=-0.21
    # In Johansson et al 2007 a=201, b=0.26, c=-0.41  (cited in Karlsson…)
    # Here we will use other values of a, b and c, selected to give
    # reasonable results for CT2
    #
    # A critical issue with this approach is that viscosity needs to be
    # measured. Viscosity for black liquor is dependent not only on
    # the temperature (T) and the solids content (i.e. xF and xL) but also
    # on the specific cooking conditions (as it affects the concentration
    # of inorganics like Cl and SO4, molecular weight of lignin, etc.)
    # Compare Zaman and Fricke 1996, Ind. Chem. Eng. Res. 590-597
    #
    # Many use an expression based on a combination of two theories:
    # absolute reaction rates and free volume concepts:
    #
    # viscosity = A*T**.5*exp(B/T+C/T0 /(T-T0))
    #
    # as described by Zaman and Fricke 1994. If the reference tempareture T0
    # is set to be close to the glass temperature, the values of the other
    # parameters (A, B and C) often get similar values for different liquids.
    #
    # A different equation is suggested by Zaman
    # and Fricke for high solids content liquids:
    #
    # ln (viscosity) = G0+G1*Z+G2*Z**2 where
    # Z=x/(x+1)*(1/T)
    #
    # Here, however, we will use a good curve fit for pure water and then
    # fake a similar curve fit for the solids and add the two together.
    #
    # DOI for references
    # 
    # Karlsson et al 2013 DOI:10.1016/J.IJHEATMASSTRANSFER.2013.07.003
    # Johansson et al 2007. Grey publication, referenced by Karlsson 2013:
    #    M. Johansson, L. Vamling, L. Olausson, An implementation-oriented heat
    #    transfer model for black liquor evaporation, in: Proceedings of the
    #    International Chemical Recovery Conference 2007, May 29–June 1, 2007,
    #    Quebec City, Canada, 2007, pp. 155–158.
    # Zaman and Fricke 1996, Ind. Chem. Eng. Res. 590-597
    #    DOI 10.1021/ie950202h
    #
    # Students are recommended to stop reading the code here :o)
    #

    TK=T+273.2

    # Let use the average flow rate and the average concentration as our
    # key parameters in the equations
    # To increase stability when routine is used in iterations,
    # take x as zero if x<0 and flow as zero if <0
    xxF=xF
    xxL=xL
    FF=F
    LL=L
    if xF<0:
        xxF=0

    if xL<0:
        xxL=0

    if F<0:
        FF=0

    if L<0:
        LL=0

    x_avg=(xxF+xxL)/2
    F_avg=(FF+LL)/2

    # Data for black liquor ABAFX011,12 (Zaman and Fricke 1994):
    #G0=6.313
    #G1=-3.668E4
    #G2=2.726E7
    # ABAFX013,14
    #G0=22.125
    #G1=-7.4027E4
    #G2=4.961E7
    #Z=x_avg./(x_avg+1).*(1./TK)
    #mu=exp(G0+G1.*Z+G2.*Z.**2)
    # Equation above may work for high x, but for a large x range
    # it does not seem to produce reasonable results

    # Try instead using good curve fit for water:
    # Al-Shemmeri, Tarik (2012). Engineering Fluid Mechnanics. 
    # Ventus Publishing ApS. pp. 17–18. ISBN 978-87-403-0114-4.
    #
    # and fake a similar curve fit for the solids part
    # and then add the two together assuming viscosity is additive
    A=0.00002414
    B=247.8
    C=140
    AA=A*1e7
    BB=309.75
    CC=C
    n=7
    mu=(1-x_avg)*A*10**(B/(TK-C))+x_avg*AA*10**(BB/(TK-CC))*(x_avg/(x_avg+1)*300/TK)**n

    # Curve fit for apparent overall heat transfer coeff:
    # Gamma=mass flow rate per unit width (specific mass flow rate) and thus
    # dependent on the specific design of the evaporator
    # Johansson et al
    # h=201 * Gamma**0.26* my**-0.41
    # 0.3<Gamma < 1.1 kg/ms
    # 0.001 <my < 0.024 Pa s
    #
    # Karlsson et al 2013
    # h=558 * Gamma**0.073 * my**-021
    a=501
    b=0.26
    c=-0.21
    d=5 # Just a value to get reasonable values for CT2 changed 2020-07-03
    gamma=F_avg/d
    k=a*gamma**b*mu**c/1e3 # kW/m2,K
    return [k,mu]
