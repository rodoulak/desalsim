import math
import numpy as np
from density_calc import density_calc
import constants 
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
        Calculate osmotic pressure for a solution.

         Attributes
         ----------
            C1, C2 C3, C4, C5, C6 : float
                Concentration of ions in the solution (mol/L).
            z1, z2,z3, z4, z5, z6 : int
                Charge of ions in the solution.
        Methods
        -------    
        osmotic_pressure_calculation():
            Calculates the osmotic pressure of a solution.
        Returns
        ------- 
            p_osmo : float
                Osmotic pressure of a solution (bar). 
    """
    def __init__(self, C_values, z_values, T):
        self.Ci = C_values
        self.z = z_values
        self.T = T
        self.sum_Ci=sum(self.Ci)

    def osmotic_pressure_calculation(self):
            """
            Calculates the osmotic pressure of a solution.
            """
            mi=[self.Ci[0]*1000/(MW_Na*1000*((1e+6-self.sum_Ci*1000)/1e+6)), self.Ci[1]*1000/(MW_Cl*1000*((1e+6-self.sum_Ci*1000)/1e+6)), self.Ci[2]*1000/(MW_K*1000*((1e+6-self.sum_Ci*1000)/1e+6)), self.Ci[3]*1000/(MW_Mg*1000*((1e+6-self.sum_Ci*1000)/1e+6)), self.Ci[4]*1000/(MW_Ca*1000*((1e+6-self.sum_Ci*1000)/1e+6)), self.Ci[5]*1000/(MW_SO4*1000*((1e+6-self.sum_Ci*1000)/1e+6))]
            self.mizi_2=[]
            for i in range(6) :
                self.mizi_2.append(mi[i]*self.z[i]**2)
            SI = sum(self.mizi_2) / 2
            B = -348.662 / self.T + 6.72817 - 0.971307 * math.log(self.T)
            C = 40.5016 / self.T - 0.721404 + 0.103915 * math.log(self.T)
            D = 5321 / self.T + 233.76 - 0.9297 * self.T + 0.001417 * self.T ** 2 - 0.0000008292 * self.T ** 3
            S = 1.17202 * (sum(self.mizi_2) / sum(self.Ci)) * 0.9982 ** 0.5 * (23375.556 / (D * self.T)) ** 1.5
            fi=1-S/(3.375*SI)*((1+1.5*SI**0.5)-2*math.log(1+1.5*SI**0.5)-1/(1+1.5*SI**0.5))+B*sum(mi)/2+C*(sum(mi)/2)**2

            #Calculate Osmotic pressure and convert units from psi to bar            
            p_osmo= 1.205 * fi * self.T * sum(mi)/ 14.3 
            return p_osmo
            
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
