def density_calc(T, S):
    """
        Calculate the density of water based on its temperature and salinity.
    
        This function uses empirical equations to calculate the density of seawater 
        based on temperature (T) and salinity (S) using the UNESCO 1983 (EOS 80) formula.
    
        Parameters:
        -----------
        T : float
            Temperature of the water in degrees Celsius (°C).
        S : float
            Salinity of the water in grams per kilogram (g/kg).
    
        Returns:
        --------
        float
            The calculated density of the seawater in kilograms per cubic meter (kg/m³).
        
        Notes:
        ------
        The calculation is based on the following empirical formulas:
        - A: Polynomial coefficients for temperature dependence
        - B: Polynomial coefficients for salinity dependence
        - C: Second-order polynomial term for salinity dependence
        - rho: Pure water density at the given temperature
        - rhos: Final seawater density considering both temperature and salinity
    """
    A= 8.24493E-1 - 4.0899E-3*T + 7.6438E-5*T**2 -8.2467E-7*T**3 + 5.3675E-9*T**4
    B = -5.724E-3 + 1.0227E-4*T - 1.6546E-6*T**2
    C = 4.8314E-4
    rho = 1000*(1 -(T+288.9414)/(508929.2*(T+68.12963))*(T-3.9863)**2)
    rhos = rho + A*S + B*S**(3/2) + C*S**2
    return rhos