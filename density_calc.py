#Density calculation based on the stream salinity and temperature

#input: 
 # T- Temperayure in C
 # S- salinity in g/kg 

#output: density in kg/m3 

def density_calc(T, S):
    A= 8.24493E-1 - 4.0899E-3*T + 7.6438E-5*T**2 -8.2467E-7*T**3 + 5.3675E-9*T**4
    B = -5.724E-3 + 1.0227E-4*T - 1.6546E-6*T**2
    C = 4.8314E-4
    rho = 1000*(1 -(T+288.9414)/(508929.2*(T+68.12963))*(T-3.9863)**2)
    rhos = rho + A*S + B*S**(3/2) + C*S**2
    return rhos