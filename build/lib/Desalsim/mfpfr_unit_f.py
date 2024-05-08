import math
from Desalsim import scaleup
from Desalsim.density_calc import density_calc 
import math
from Desalsim import constants 
#%%
    #molecular weight
MW_Na=constants.MW_Na
MW_Cl=constants.MW_cl
MW_SO4=constants.MW_so4
MW_K=constants.MW_K
MW_Ca=constants.MW_Ca
MW_Mg=constants.MW_Mg
MW_HCO3=constants.MW_HCO3
MW_values = [MW_Na, MW_Cl, MW_K, MW_Mg, MW_Ca, MW_SO4]
MW_MgOH=constants.MW_MgOH
MW_CaOH=constants.MW_CaOH

#%%
class MFPFRCALC:
    """
    A class used to represent Mass Balance for percipitation step with the addition of alkaline solution 

    Attributes
    ----------
    MW : float
        Molecular weight of the solute (g/mol)
    Ci_in : float
        Initial concentration of the solute (g/L)
    conv_1: float
        Conversion rate for Magnisium percipitation in the first step 
    conv_2: float
        Conversion rate for Calcium percipitation in the second step    
    C_NaOH_1: float
        Concentration of NaOH solution for first step (mol/L)
    Qin
    QMg_in 
    QNaOH_1 : float 
        Volumetric flow rate of sodium hydroxide in step 1 in L/h
    M_MgOH2_1: float 
        Outlet mass flow rate of magnesium hydroxide produced in kg/h
    Qtot_out_1 : float 
        Outlet volumetric flow rate in L/h
    Mtot_out_1 : float 
        Step 1 outlet mass flow rate in kg/h
    magma_d_1 : float 
        magma density: the quantity of solids produced per volume of slurry kg/l
    ph_1 : float 
        ph of solution during first step 
    kps_MgOH : float 
        Product solubility of Mg(OH)2 
    Ci_out_1 : float 
        The outlet ion concentration from step 1 in mol/L
    QNaOH_2_st : float 
        The stoichiometric volumetric flow rate of sodium hydroxide for the second step L/hr
    QNaOH_2_add : float
        The added volumetric flow rate of sodium hydroxide needed to reach a pH = 13  L/h
    M_CaOH2_2: float 
        The outlet mass flow rate of calcium hydroxide produced during the 2nd step kg/hr 
    M_MgOH2_2 : float 
        The outlet mass flow rate of magnesium hydroxide produced during the 2nd step kg/hr 
    magma_d_2 : float    
        Magma density: the quantity of solids produced per volume of slurry kg/l
    Qtot_out_2 : float 
        Total outlet volumetric flow rate for 2nd step in L/h
    Ci_out_2 : float 
        The outlet ion concentration from step 2 in mol/L
      
    Methods
    -------
    calc_step1():
        Calculates the flowrates, the concentration of the streams and the requirements of alkaline solution in the 1st step 
    calc_step2():
        Calculates the flowrates, the concentration of the streams and the requirements of alkaline solution in the 2nd step 
        """
    def __init__(self, Qin, Cin_mfpfr, C_NaOH_1, C_NaOH_2, conv_1, conv_2):
        self.Qin=Qin
        Cc1, Cc2, Cc3, Cc4, Cc5, Cc6=Cin_mfpfr
        self.CNa_in=Cc1/MW_Na #Na concentration (g/l) to (mol/l)
        self.CCl_in=Cc2/MW_Cl #Cl concentration (g/l) to (mol/l)
        self.CK_in=Cc3/MW_K #K concentration (g/l) to (mol/l)
        self.CMg_in=Cc4/MW_Mg #Mg concentration (g/l) to (mol/l)
        self.CCa_in=Cc5/MW_Ca #Ca concentration (g/l) to (mol/l)
        self.CSO4_in=Cc6/MW_SO4 #SO4 concentration (g/l) to (mol/l)
        self.C_NaOH_1=C_NaOH_1 #Molar concentration of NaOH (mol/l) used in the first precipitation step 
        self.C_NaOH_2=C_NaOH_2 #Molar concentration of NaOH (mol/l) used in the first precipitation step 
        self.conv_1=conv_1 #Conversion rate of Mg the first precipitation step 
        self.conv_2=conv_2 #Conversion rate of Mg the first precipitation step 
        self.d_in=density_calc(25, sum(Cin_mfpfr))
    
    def calc_step1(self, kps_MgOH, d_mgoh_2):
        #Calculate the molar flow rate of magnesium in the reactor during the 1° stepin mol/h
        self.QMg_in=self.Qin*self.CMg_in 
        
        #Calculate the volumetric flow rate of sodium hydroxide in L/h
        self.QNaOH_1=(self.Qin*self.CMg_in*(self.conv_1/100)*2)/self.C_NaOH_1 
        
        #Calculate outlet mass flow rate of magnesium hydroxide produced in kg/h
        self.M_MgOH2_1=(self.Qin*self.CMg_in*(self.conv_1/100)*MW_MgOH)/1000 
        
        #Calculate outlet volumetric flow rate in L/h
        self.Qtot_out_1=self.Qin+self.QNaOH_1
        
        #Calculate outlet mass flow rate in kg/h
        self.Mtot_out_1=self.Qin*self.d_in+self.QNaOH_1*1.04-self.M_MgOH2_1  
        
        #Calculate magma density: the quantity of solids produced per volume of slurry kg/l
        self.magma_d_1=self.M_MgOH2_1/self.Qtot_out_1 
        
        #Calculate ph of solution during the first step 
        self.ph_1=14+math.log10(2*(kps_MgOH/4)**(1/3))
        
        #Calculate total outlet volumetric flow rate in L/h
        self.Qtot_out_1=self.Qin+self.QNaOH_1-self.M_MgOH2_1/d_mgoh_2 
        
        #Calculates the outlet ion concentration in mol/L
        self.CMg_out_1=(self.QMg_in*(1-self.conv_1/100))/self.Qtot_out_1 
        self.CNa_out_1=(self.Qin*self.CNa_in+self.QNaOH_1*self.C_NaOH_1)/self.Qtot_out_1
        self.CCl_out_1=(self.Qin*self.CCl_in)/self.Qtot_out_1
        self.CK_out_1=(self.Qin*self.CK_in)/self.Qtot_out_1
        self.CCa_out_1=(self.Qin*self.CCa_in)/self.Qtot_out_1
        self.CSO4_out_1=(self.Qin*self.CSO4_in)/self.Qtot_out_1
        
        return self.QMg_in, self.QNaOH_1, self.M_MgOH2_1, self.Qtot_out_1, self.Mtot_out_1, self.magma_d_1, self.ph_1, self.CMg_out_1, self.CNa_out_1, self.CCl_out_1, self.CK_out_1, self.CCa_out_1, self.CSO4_out_1
    
    def calc_step2(self, d_mgoh_2, d_caoh_2 ):
       #Calculate the molar flow rate of calcium in the reactor during the 2nd stepin mol/h
       self.QCa_in_2=self.Qtot_out_1*self.CCa_out_1 
       
       #Calculate the stoichiometric volumetric flow rate of sodium hydroxide for the second step L/hr
       self.QNaOH_2_st=(self.Qtot_out_1*(self.CCa_out_1*(self.conv_2/100)+self.CMg_out_1*1)*2)/self.C_NaOH_2 
       
       #Calculate concentration of the hydroxide ion in mol/L for a ph=13 solution  
       self.COH_ph13=10**(-(14-13)) 
       self.COH_st=10**(-(14-self.ph_1))

       #Calculate the added volumetric flow rate of sodium hydroxide needed to reach a pH = 13  L/h
       self.QNaOH_2_add=((0.0216-self.COH_ph13)*(self.QNaOH_2_st+self.Qtot_out_1))/(self.COH_ph13-self.C_NaOH_2)
       
       #Calculate the total outlet volumetric flow rate from 2nd step L/h
       self.Qtot_out_2=self.Qtot_out_1+self.QNaOH_2_st+self.QNaOH_2_add
      
       #Calculate the outlet mass flow rate of calcium hydroxide produced during the 2nd step kg/hr 
       self.M_CaOH2_2=(self.Qtot_out_1*self.CCa_out_1*(self.conv_2/100)*MW_CaOH)/1000 #kg/hr
      
       #Calculate the outlet mass flow rate of magnesium hydroxide produced during the 2nd step kg/hr 
       self.M_MgOH2_2=self.Qtot_out_1*self.CMg_out_1*MW_MgOH/1000 
       
       #Calculate magma density: the quantity of solids produced per volume of slurry kg/l
       self.magma_d_2=(self.M_MgOH2_2+self.M_CaOH2_2)/self.Qtot_out_2 
      
       #Calculate total outlet volumetric flow rate for 2nd step in L/h
       self.Qout_2=self.Qtot_out_2-self.M_CaOH2_2/d_caoh_2-self.M_MgOH2_2/d_mgoh_2
      
       
       #Calculates the outlet ion concentration in mol/L
       self.CNa_out_2=(self.Qtot_out_1*self.CNa_out_1+(self.QNaOH_2_st+self.QNaOH_2_add)*self.C_NaOH_2)/self.Qtot_out_2       
       self.CCa_out_2=(self.QCa_in_2*(1-self.conv_2/100))/self.Qtot_out_2
       self.CCl_out_2=(self.Qtot_out_1*self.CCl_out_1)/self.Qtot_out_2
       self.CK_out_2=(self.Qtot_out_1*self.CK_out_1)/self.Qtot_out_2
       self.CMg_out_2=0
       self.CSO4_out_2=(self.Qtot_out_1*self.CSO4_out_1)/self.Qtot_out_2
       
       #Calculate ph of solution during the second step 
       self.ph_2=14+math.log10(0.1)
    
       return self.QCa_in_2, self.QNaOH_2_st, self.COH_ph13, self.QNaOH_2_add, self.Qtot_out_2, self.M_CaOH2_2, self.M_MgOH2_2, self.magma_d_2, self.Qout_2, self.CNa_out_2, self.CCa_out_2, self.CCl_out_2, self.CK_out_2, self.CMg_out_2, self.CSO4_out_2,  self.ph_2

class HClAddition:
    def __init__(self, Qout_2, Cout_all_m, MW_cl, ph_2, HCl_conc):
        self.Qout_2 = Qout_2
        self.Cout_all_m = Cout_all_m
        self.MW_cl = MW_cl
        self.ph_2 = ph_2
        self.HCl_conc=HCl_conc

    def calculate_HCl_addition(self, Cout_mfpfr_g):
        OH_initial = math.pow(10, -14) / math.pow(10, -self.ph_2)
        OH_final = math.pow(10, -14) / math.pow(10, -7)
        
        #Calculate the volume of HCl 1M in liters, Qout_f
        QHCl = (self.Qout_2 * OH_initial - self.Qout_2 * OH_final) / (OH_final + self.HCl_conc)  
        self.Qout_f = self.Qout_2+QHCl

        # Calculate the outlet concentration of chloride ions
        C_cl_out = (self.Cout_all_m[1] * self.Qout_2 + QHCl * self.HCl_conc) / self.Qout_f  # mol/l
        self.Cout_all_m[1] = C_cl_out
        C_cl_out = C_cl_out * self.MW_cl  # g/l
        Cout_mfpfr_g[1] = C_cl_out

        for i in range(2, 6):
            Cout_mfpfr_g[i] = Cout_mfpfr_g[i] * self.Qout_2 / self.Qout_f

        Cout_mfpfr_g[0] = Cout_mfpfr_g[0] * self.Qout_2 / self.Qout_f
        # Return the volume of HCl added and the outlet concentration of chloride ions
        return QHCl, Cout_mfpfr_g
     
  # Get the outlet flow rate

#%%Energy consumption 
class energycons:
    def energycalc(Qtot, QNaOH, Qin, QNaOH_1, QNaOH_2_add, QNaOH_2_st, dp, npump):
        Qtot=Qtot
        QNaOH=QNaOH
        Qin=Qin
        QNaOH_2_st=QNaOH_2_st
        QNaOH_2_add=QNaOH_2_add
        QNaOH_1=QNaOH_1        
        
        #Calculate the energy required for pumping in kW
        Epump_1=(Qin*dp+QNaOH_1*dp)*1e5/3600/(1000*npump) #(W)
        Epump_2=((Qin+QNaOH_1)*dp+(QNaOH_2_add+QNaOH_2_st)*dp)*1e5/3600/(1000*npump)       #(W)
        return Epump_1, Epump_2
                    

