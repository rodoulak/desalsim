import math
import numpy as np
from desalsim.density_calc import density_calc 
from desalsim import constants 
import math
#
    #molecular weight
MW_Na=constants.MW_Na
MW_Cl=constants.MW_cl
MW_SO4=constants.MW_so4
MW_K=constants.MW_K
MW_Ca=constants.MW_Ca
MW_Mg=constants.MW_Mg
MW_HCO3=constants.MW_HCO3
MW_values = [MW_Na, MW_Cl, MW_K, MW_Mg, MW_Ca, MW_SO4]

#%%
class NFMass:
    """
    A class used to represent Mass Balance for Nanofiltration Unit

    ...

    Attributes
    ----------
    comp : str
        Component name
    Cfeedi : float
        Ion concentration in the feed (g/L)
    rjr : float
        Rejection of the ion by the membrane
    Wrec : float
        Water recovery in the first pass
    Qf : float
        Feed flow rate (kg/h)
    Qperm : float
        Permeate flow rate (kg/h)
    Qconc : float
        Concentrate flow rate (kg/h)
    Cpermi : float
        Ion concentration in the permeate (g/L)
    Cconci : float
        Ion concentration in the concentrate (g/L)

    Methods
    -------
    calculate_perm():
        Calculates the permeate and concentrate flow rates, and their ion concentrations
    """

    def __init__(self, comp, Cfeedi, rjr, Wrec, Qf):
        self.comp = comp
        self.Cfeedi = Cfeedi
        self.rjr = rjr
        self.Wrec = Wrec
        self.Qf = Qf
        self.calculate_perm()

    def calculate_perm(self):
        self.Qperm = self.Wrec * self.Qf
        self.Qconc = self.Qf - self.Qperm
        self.Cpermi = (1 - self.rjr) * self.Cfeedi
        self.Cconci = (self.Qf * self.Cfeedi - self.Qperm * self.Cpermi) / self.Qconc

#%% osmotic pressure 
class OsmoticPressure:
    """
        Calculate osmotic pressure for a solution using Gibbs equation based on water activity.

         Attributes
         ----------
            Ci_in : float
                Concentration of ions in the solution (mol/L).
            MW_values : float
                Ions molar mass in g/mol.
        Methods
        -------    
        calculate_molalities():
            Calculate molality of each component
        calculate_moles_of_water():
            Calculate moles of water
        calculate_total_moles_of_solute():
            Calculate the total moles of solute
        calculate_osmotic_pressure():
            Calculates the osmotic pressure of a solution.
            
        Returns
        ------- 
            osmotic_pressure_bar : float
                Osmotic pressure of a solution (bar). 
    """
    def __init__(self, Ci_in, MW_values, T=298):
        """
        Initialize the calculator with concentrations and temperature.
        
        Parameters:
        - Ci_in (list): Concentrations of each component in g/L.
        - T (float): Temperature in Kelvin. Default is 298 K.
        """
        self.Ci_in = Ci_in  # Concentrations in g/L
        self.T = T  # Temperature in Kelvin
        self.R = 0.0821  # Ideal gas constant in L·atm/(mol·K)
        self.V = 0.018015  # Molar volume of water in L/mol
        self.M_w = 0.018015  # Molar mass of water in kg/mol
        self.molar_masses = MW_values   # Molar masses in g/mol for Na, Cl, K, Mg, Ca, SO4

    def calculate_molalities(self):
        """
        Convert the concentration from g/L to mol/kg of water (assuming solution density ~ 1 kg/L).
        
        Returns:
        - molalities (list): Molality of each component in mol/kg.
        """
        molalities = [Ci / mol_mass for Ci, mol_mass in zip(self.Ci_in, self.molar_masses)]
        return molalities

    def calculate_moles_of_water(self):
        """
        Calculate moles of water (assuming 1 kg of water).
        
        Returns:
        - moles_of_water (float): Moles of water.
        """
        return 1 / self.M_w

    def calculate_total_moles_of_solute(self, molalities):
        """
        Calculate the total moles of solute.
        
        Parameters:
        - molalities (list): Molalities of each component.
        
        Returns:
        - total_moles_of_solute (float): Total moles of solute.
        """
        return sum(molalities)

    def calculate_mole_fraction_of_water(self, molalities):
        """
        Calculate the mole fraction of water.
        
        Parameters:
        - molalities (list): Molalities of each component.
        
        Returns:
        - X_water (float): Mole fraction of water.
        """
        moles_of_water = self.calculate_moles_of_water()
        total_moles_of_solute = self.calculate_total_moles_of_solute(molalities)
        return moles_of_water / (moles_of_water + total_moles_of_solute)

    def calculate_osmotic_pressure(self):
        """
        Calculate the osmotic pressure using the Gibbs equation.
        
        Returns:
        - osmotic_pressure_bar (float): Osmotic pressure in bar.
        """
        molalities = self.calculate_molalities()
        X_water = self.calculate_mole_fraction_of_water(molalities)
        a_w = X_water  # Approximating water activity as mole fraction of water
        osmotic_pressure = -(self.R * self.T / self.V) * math.log(a_w)  # Osmotic pressure in atm
        osmotic_pressure_bar = osmotic_pressure * 1.01325  # Convert atm to bar
        return osmotic_pressure_bar
            
#%%Energy consumption 
class NfEnergy:   
    """
        Calculate energy consumption for nanofiltration.
        
            Attributes
            ----------
            P_osmo_c : float
                Osmotic pressure of concentrate stream (bar).
            P_osmo_f : float
                Osmotic pressure of feed stream (bar).
            P_osmo_p : float
                Osmotic pressure of permeate stream (bar).
            dp : float
                Pressure drop (bar).
            d_p:  float
                Permeate stream density (kg/m3).
            Qperm : float
                Permeate flow rate(kg/h).
            Qf : float
                Concentrate flow rate (kg/h).
            d_in : float
                Feed stream density(kg/m3).
            n : float
                Pump efficiency(-).
            Methods
            -------
             calculate_energy_consumption:  
                Calculates the Energy consumption of the nanofiltration unit and the specific energy consumption  
            Returns
            -------
                Papplied : float
                    Applied pressure (bar) 
                Ppump : float
                    Power for applied pressure (W)
                E_el_nf : float
                    Electricity consumption KW
                Spec : float
                    Specific Energy Consumption (Kwh/m3 of permeate)
                SEC_el_feed : float
                    Specific Energy Consumption (Kwh/m3 of feed)
        """
    def __init__(self, P_osmo_c, P_osmo_f, P_osmo_p, dp, d_p, Qperm, Qf, d_in, n):
        self.P_osmo_c = P_osmo_c         
        self.P_osmo_f = P_osmo_f 
        self.P_osmo_p = P_osmo_p 
        self.dp = dp 
        self.d_p = d_p
        self.Qperm = Qperm 
        self.Qf = Qf 
        self.d_in = d_in
        self.n =n
        self.calculate_energy_consumption()

    def calculate_energy_consumption(self):
        """
        Calculate the energy consumption of the nanofiltration unit.

        Returns:
            dict: A dictionary containing the calculated values.
        """
        # Calculate the applied pressure (bar)
        Papplied = (self.P_osmo_c + self.P_osmo_f) / 2 - self.P_osmo_p + self.dp 

        # Calculate the power for the pump
        Ppump = Papplied * self.Qperm / self.d_p * 1e5 / 3600 

        # Calculate the electrical energy consumption (KWh)
        self.E_el_nf = (Ppump / 1000 / self.n ) 

        # Calculate the specific energy consumption for the permeate (KWh/m3 of permeate)
        Spec = self.E_el_nf / (self.Qperm / self.d_p) 

        # Calculate the specific energy consumption for the feed (KWh/m3 of feed)
        SEC_el_feed = self.E_el_nf / (self.Qf / self.d_in) 

        return {
            "Applied pressure (Bar)": round(Papplied,2),
            "Power for pump (KW)": round(Ppump/1000,2) ,
            "E_el_nf (KW)": round(self.E_el_nf,2),
            "Specific Energy Consumption (KWh/m3 of permeate)": round(Spec,2),
            "Specific Energy Consumption (KWh/m3 of feed)": round(SEC_el_feed,2)
        }
#%%
#Example usage
#Feed concentration
components = ['Na', 'Cl', 'K', 'Mg', 'Ca', 'SO4']
Ci_in = [12.33, 21.67, 0.45, 1.39, 0.45, 3.28]
z_values = [1, -1, 1, 2, 2, -2]


#Constants
R=8.314 #gas constant (units: J / mol·K)
T=20+273 #Operating temperature (units: K)

    #molecular weight
MW_Na=constants.MW_Na
MW_Cl=constants.MW_cl
MW_SO4=constants.MW_so4
MW_K=constants.MW_K
MW_Ca=constants.MW_Ca
MW_Mg=constants.MW_Mg
MW_HCO3=constants.MW_HCO3
MW_values = [MW_Na, MW_Cl, MW_K, MW_Mg, MW_Ca, MW_SO4]
mg_in = sum(Ci_in)
#Feed flow density 
d_in = density_calc(T-273, mg_in)  # kg/m3

#Feed flowrate
Qsw = 3000 / 24 * d_in #m3/d
Qf = Qsw  # kg/hr

#Asuumptions  
rjr_values = [0.16, 0.29, 0.21, 0.98, 0.95, 0.98] #Ions rejection rates based on membrane characteristics (units: -)
Wrec = 0.7 # Water recovery based on membrane characteristics (units: -)
n=0.8 #pump efficiency (units: -)
dp=2 # pressure drop (units: bar)


# Function to create NFMass objects for different components
def create_nfmass_objects(components, C_in, rjr_values, Wrec, Qf):
    return [NFMass(comp, Ci, rjr, Wrec, Qf) for comp, Ci, rjr in zip(components, C_in, rjr_values)]

# Create NFMass objects for different components
nfmass_objects = create_nfmass_objects(components, Ci_in, rjr_values, Wrec, Qf)

Cconc = [nf_mass.Cconci for nf_mass in nfmass_objects]
Cperm = [nf_mass.Cpermi for nf_mass in nfmass_objects]
Qperm = nfmass_objects[0].Qperm  # kg/hr
Qconc = nfmass_objects[0].Qconc  # kg/hr
print("Permeate stream flow rate is "+str(round(Qperm,2))+"kg/hr")
print("Permeate stream total concentration is "+str(round(sum(Cperm),2))+"g/l")
print("-----------------------------------------")
print("Concentrate stream flow rate is "+str(round(Qconc,2))+"kg/hr")
print("Concentrate stream total concentration is "+str(round(sum(Cconc),2))+"g/l")
print("-----------------------------------------")

# Calculate Osmotic Pressure
P_osmo_f = OsmoticPressure(Ci_in, MW_values).calculate_osmotic_pressure()
P_osmo_p = OsmoticPressure(Cperm, MW_values).calculate_osmotic_pressure()
P_osmo_c = OsmoticPressure(Cconc, MW_values).calculate_osmotic_pressure()

d_p=density_calc(T-273, sum(Cperm))

#Calculate Energy consumption 
elec=NfEnergy(P_osmo_c, P_osmo_f, P_osmo_p, dp, d_p, Qperm, Qf, d_in,n)
result=elec.calculate_energy_consumption()
for key, value in result.items():
        print(f"{key}: {value}")