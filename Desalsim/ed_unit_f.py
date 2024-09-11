import numpy as np
from . import scaleup
from . import constants 
import math
from .density_calc import density_calc 

#%% Calculations 
"""
A class used to represent Mass Balance for Electrodialysis

...

Attributes
----------
    MWs : float
        Molecular weight of NaCl (g/mol)
    MWw : float
        Molecular weight of water (g/mol)
    R : float
        Resistance of rinse stream (ohm)
    Rp : float
        Resistance of polarization (ohm)
    A : float
        Active area of cell-pair (m^2)
    F : float
        Faraday constant (C/mol)
    T : float
        Temperature in Kelvin
    dp : float
        Parameter dp
    npump : float
        Pump efficiency
    rho_w : float
        Density of water (kg/m^3)
    D : float
        Diffusion coefficient (m^2/s)
    tcu : float
        tcu parameter
    veloc : float
        Velocity (m^2/s)
    h : float
        Height (mm)
    Sh : float
        Sh parameter
    Mem_eff : float
        Membrane efficiency
    Ncp : integer
        Number of cell-pairs
    Qed_in : float
        Inlet flow rate (l/h)
    Qed_in_c : float
        Concentrate inlet flow rate (l/h)
    Qed_in_d : float
        Dilute inlet flow rate (l/h)


Methods
-------
    Ts_cp(S): 
        Transport number for salt in concentrate compartment
    Tw_cp(Sc, Sd): 
        Transport number for water in concentrate compartment
    Ls_cp(Sc, Sd): 
        Permeability for salt in concentrate compartment
    Lw_cp(S): 
        Permeability for water in concentrate compartment
    p_osmo(S): 
        Osmotic pressure in concentrate compartment
    p_osmo2(S): 
        Alternative osmotic pressure calculation
    dC(Ts_cp): 
        Change in concentration
"""
class ElectrodialysisCalc:
    def __init__(self):
        pass

    def Ts_cp(S):
        """Transport number for salt in concentrate compartment"""
        Ts_cp=-4e-6*S**2 + 4e-5*S + 0.96
        return Ts_cp
    
    def Tw_cp(Sc,Sd):
        """Transport number for water in concentrate compartment"""
        return -4e-5*Sc**2 - 1.9e-2*Sd + 11.2
    
    def Ls_cp(Sc, Sd):
        """Permeability for salt in concentrate compartment"""
        return min(2e-12*Sd**2 - 3e-10*Sd + 6e-8, 2e-12*Sc**2 - 3e-10*Sc + 6e-8)
    
    def Lw_cp(S):
        """Permeability for water in concentrate compartment"""
        return 5*S**(-0.416)
    
    
    def p_osmo(S, T, MWs): 
        """calculation for Osmotic pressure in concentrate compartment"""
        C1=S/MWs*constants.MW_Na/constants.MW_Na
        C2=S/MWs*constants.MW_cl/constants.MW_cl
        sum_Ci=sum([C1, C2])
        p_osmo=0.0831446261815324*sum_Ci*T
        return p_osmo
    
    def dC(Ts_cp, tcu, D, Ij,  h, Sh ):
        """Change in concentration"""
        F = 96485.3329  # Faraday constant (C/mol)
        Tcu=(Ts_cp+1)/2
        dC=-(Tcu-tcu)/D*(Ij/F)*(2*h/1000/Sh)
        return dC

