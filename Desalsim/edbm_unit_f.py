import math 
import numpy as np
from . import scaleup
from .density_calc import density_calc 
from . import constants


#%%
#Molecular weight 
MW_Na=constants.MW_Na
MW_Cl=constants.MW_cl
MW_SO4=constants.MW_so4
MW_K=constants.MW_K
MW_Ca=constants.MW_Ca
MW_Mg=constants.MW_Mg
MW_HCO3=constants.MW_HCO3
MW_H=1.00784
MW_OH=17.008
MW_H2O=MW_H+MW_OH

#%%Calculations 
class EDBMCalc:
    """
    A class used to represent Mass Balance for Electrodialysis With Bipolar Membranes 
    ...

    Attributes
    ----------

    CNa_in, CCl_in, CK_in, CMg_in, CCa_in, CSO4_in : float
        Initial concentrations of various ions (g/l).
    CNa_out, CCl_out, CK_out, CMg_out, CCa_out, CSO4_out : float
        Outlet concentrations of various ions (g/l).
    Ci_s_in: float 
        Initial concentrations of various ions in salt channel (mol/l)
    Ci_a_in: float 
        Initial concentrations of various ions in acid channel (mol/l)
    Ci_b_in: float 
        Initial concentrations of various ions in base channel (mol/l)  
    EMF: float
        electromotive force (V) 
    KW_s_in: float
        Inlet ionic water product in salt channel 
    KW_a_in: float
        Inlet ionic water product in acid channel 
    KW_b_in: float
        Inlet ionic water product in base channel         
    M_h2o_a_in: float 
        Initial mass flow rate of water in acid channel (kg/h)
    M_h2o_b_in: float 
        Initial mass flow rate of water in base channel (kg/h)
    M_h2o_s_in: float 
        Initial mass flow rate of water in salt channel (kg/h)        
    N_trip: integer 
        Number of triplets of a channel 
    P: float
        Gross power needed (W)
    PM: float
        Molecular weight
    Q : float
        Flow rate (l/h)
    V_ext: float
        Voltage needed (V) 


    Methods
    -------
    flowrate(): 
        Calculates the flowrate in each channel 
    
    in_mass_flow_rates():
        Calculates the inlet mass flow rates of each ion, kg/h 
    
    acid_channel():
        Calculates flow rates and outlet concentration in Acid channel 
    
    base_channel():
        Calculates flow rates and outlet concentration in Base channel 
        
    Salt_channel():
        Calculates flow rates and outlet concentration in Salt channel 
    """
    def __init__(self, Qin, A, I_d, N, C_s_in, C_b_in, C_a_in, T ):
        #input conditions   
        F=96485.3 #Coulombs/mol
        self.A=A
        self.I_ext=self.A*I_d
        self.JA=3.6*self.I_ext/F
        self.T=T
        self.N_trip=N 
        #Initialization of parameters
        self.Qin=Qin
        
        Cc1, Cc2, Cc3, Cc4, Cc5, Cc6, Cc7, CH, COH= C_s_in
        Cb1, Cb2,Cb3,Cb4,Cb5,Cb6,Cb7,Cb8,Cb9=C_b_in
        Ca1,Ca2,Ca3,Ca4,Ca5,Ca6,Ca7,Ca8,Ca9=C_a_in
        self.CNa_s_in=Cc1/MW_Na #Na concentration (g/l) to (mol/l)
        self.CCl_s_in=Cc2/MW_Cl #Cl concentration (g/l) to (mol/l)
        self.CK_s_in=Cc3/MW_K #K concentration (g/l) to (mol/l)
        self.CMg_s_in=Cc4/MW_Mg #Mg concentration (g/l) to (mol/l)
        self.CCa_s_in=Cc5/MW_Ca #Ca concentration (g/l) to (mol/l)
        self.CSO4_s_in=Cc6/MW_SO4 #SO4 concentration (g/l) to (mol/l)
        self.CHCO3_s_in=Cc6/MW_HCO3 #HCO3 concentration (g/l) to (mol/l)      
        self.COH_s_in=COH/MW_OH
        self.CH_s_in=CH/MW_H
        self.PM_i=[MW_Na, MW_Cl, MW_K, MW_Mg, MW_Ca, MW_SO4, MW_HCO3, MW_H, MW_OH]
        
        self.d_a=density_calc(self.T,sum(C_a_in) )/1000
        self.d_b=density_calc(self.T,sum(C_b_in) )/1000
        self.d_s=density_calc(self.T,sum(C_s_in) )/1000
        
        #Create lists for initial concentration in each channel 
        self.Ci_s_in=[self.CNa_s_in,self.CCl_s_in,self.CK_s_in,self.CMg_s_in,self.CCa_s_in,self.CSO4_s_in,self.CHCO3_s_in, self.CH_s_in,self.COH_s_in]
        self.Ci_b_in=[Cb1/MW_Na, Cb2/ MW_Cl,Cb3/MW_K,Cb4/MW_Mg,Cb5/MW_Ca,Cb6/MW_SO4,Cb7/MW_HCO3,Cb8/MW_H,Cb9/MW_OH]
        self.Ci_a_in=[Ca1/MW_Na,Ca2/MW_Cl,Ca3/MW_K,Ca4/MW_Mg,Ca5/MW_Ca,Ca6/MW_SO4,Ca7/MW_HCO3,Ca8/MW_H,Ca9/MW_OH]
        
        #save initial concentration before recirculation: 
        self.Ci_s_in_0=self.Ci_s_in
        self.Ci_b_in_0=self.Ci_b_in
        self.Ci_a_in_0=self.Ci_a_in

        
    def flowrate(self):
        # Calculates the flowrate in each channel 
        self.Q1_s_in=self.Qin/self.N_trip #units: l/h
        self.Q1_b_in=self.Qin/self.N_trip #units: l/h
        self.Q1_a_in=self.Qin/self.N_trip #units: l/h

    def in_mass_flow_rates(self, ph_s):
        #Calculates the inlet mass flow rates of each ion, kg/h 
        #initialize lists 
        self.M_s_in=[[],[],[],[],[],[],[],[],[]]
        self.M_b_in=[[],[],[],[],[],[],[],[],[]]
        self.M_a_in=[[],[],[],[],[],[],[],[],[]]
                
        for i in range(0,9):
            self.M_s_in[i]=self.Q1_s_in*self.Ci_s_in[i]*self.PM_i[i]/1000 #units: kg/h
            self.M_b_in[i]=self.Q1_b_in*self.Ci_b_in[i]*self.PM_i[i]/1000 #units: kg/h
            self.M_a_in[i]=self.Q1_a_in*self.Ci_a_in[i]*self.PM_i[i]/1000 #units: kg/h
            
        #Calculate mass of water in the initial streams 
        self.M_h2o_a_in=self.Q1_a_in*self.d_a-sum(self.M_a_in)          
        self.M_h2o_b_in=self.Q1_b_in*self.d_b-sum(self.M_b_in) 
        self.M_h2o_s_in=self.Q1_s_in*self.d_s-sum(self.M_s_in) 
        
        #Calculate inlet ionic water product in each channel 
        self.KW_s_in=self.Ci_s_in[7]*self.Ci_s_in[8] 
        self.KW_a_in=self.Ci_a_in[7]*self.Ci_a_in[8]
        self.KW_b_in=self.Ci_b_in[7]*self.Ci_b_in[8]
    
        #Calculation of the inlet concentration of the Cl in salt channel (exception)
        if ph_s<7:
            self.M_s_in[1]=self.Q1_s_in*self.Ci_s_in[1]*self.PM_i[1]/1000+(self.M_s_in[7]/self.PM_i[7]-self.M_s_in[8]/self.PM_i[8])*(self.PM_i[1])
        elif ph_s>7:
            self.M_s_in[1]=self.Q1_s_in*self.Ci_s_in[1]*self.PM_i[1]/1000+(self.M_s_in[8]/self.PM_i[8]-self.M_s_in[7]/self.PM_i[7])*self.PM_i[1]
        

    def acid_channel(self):
        #Performs mass balance calculations for Acid channel 
        #Initialize ionic Mass and concentration in acid channel  
        self.M_a_out=[[],[],[],[],[],[],[],[],[]]
        self.Ci_a_out=[[],[],[],[],[],[],[],[],[]]
        
        #Calculation of outlet mass flow rate for H+
        self.M_a_out[7]=self.M_a_in[7]+self.JA*self.PM_i[7]
        
        #Calculation of outlet mass flow rate for Cl-
        self.M_a_out[1]=self.M_a_in[1]+self.JA*self.PM_i[1]
        #Calculation of outlet mass flow rate for OH-
        self.M_a_out[8]=self.M_a_in[8]
        #Calculation of outlet mass flow rate for Na+
        self.M_a_out[0]=self.M_a_in[0]
        #Calculation of outlet mass flow rate for water 
        self.M_h2o_a_out=self.M_h2o_a_in-0.5*self.JA*MW_H2O
        #Calculation total outlet mass flow rate 
        self.M_a_out_t=self.M_h2o_a_out+self.M_a_out[7]+self.M_a_out[1]+self.M_a_out[8]+self.M_a_out[0]

        #Calculation of outlet mass flow rate for other ionic species in channel 
        for i in range(2,7):
            self.M_a_out[i]=self.M_a_in[i]
            self.M_a_out_t=self.M_a_out_t+self.M_a_out[i]
        
        #Calculate of volumetric outlet flow rate
        self.Q1_a_out=self.M_a_out_t/1
        #Calculation of outlet concentration of single ions in channel 
        for i in range(0,9):
            self.Ci_a_out[i]=self.M_a_out[i]/(self.Q1_a_out*self.PM_i[i]/1000)

        
    def base_channel(self):
        #Performs mass balance calculations for Base channel 
        #Initialize ionic Mass and concentration in base channel  
        self.M_b_out=[[],[],[],[],[],[],[],[],[]]
        self.Ci_b_out=[[],[],[],[],[],[],[],[],[]]
       
        #Calculation of outlet mass flow rate for Na+
        self.M_b_out[0]=self.M_b_in[0]+self.JA*self.PM_i[0]
        #Calculation of outlet mass flow rate for OH-        
        self.M_b_out[8]=self.M_b_in[8]+self.JA*self.PM_i[8]
        #Calculation of outlet mass flow rate for Cl-
        self.M_b_out[1]=self.M_b_in[1]
        #Calculation of outlet mass flow rate for H+
        self.M_b_out[7]=self.M_b_in[7]
        #Calculation of outlet mass flow rate for water 
        self.M_h2o_b_out=self.M_h2o_b_in+0.5*self.JA*MW_H2O
        #Calculation total outlet mass flow rate 
        self.M_b_out_t=self.M_h2o_b_out+self.M_b_out[7]+self.M_b_out[1]+self.M_b_out[8]+self.M_b_out[0]
 
        #Calculation of outlet mass flow rate for other ionic species in channel 
        for i in range(2,7):
            self.M_b_out[i]=self.M_b_in[i]  
            self.M_b_out_t=self.M_b_out_t+self.M_b_out[i]
            
        d_s_new=density_calc(25, self.M_b_out_t)/1000
        #Calculate of volumetric outlet flow rate
        self.Q1_b_out=self.M_b_out_t/d_s_new #assume density is 1kg/l
        
        #Calculation of outlet concentration of single ions in channel 
        for i in range(0,9):
            self.Ci_b_out[i]=self.M_b_out[i]/(self.Q1_b_out*self.PM_i[i]/1000)
    
    def salt_channel(self, Cm_bp_H, Cm_bp_OH):
        #Performs mass balance calculations for Salt channel 
        #Initialize ionic Mass and concentration in salt channel  
        self.M_s_out=[[],[],[],[],[],[],[],[],[]]
        self.Ci_s_out=[[],[],[],[],[],[],[],[],[]]
        
        #Calculation of outlet mass flow rate for Na+
        self.M_s_out[0]=self.M_s_in[0]-(self.M_b_out[0]-self.M_b_in[0])
        #Calculation of outlet mass flow rate for Cl-
        self.M_s_out[1]=self.M_s_in[1]-(self.M_a_out[1]-self.M_a_in[1])
        #Calculation of outlet mass flow rate for OH-   
        self.M_s_out[8]=self.M_s_in[8]
        #Calculation of outlet mass flow rate for H+
        self.M_s_out[7]=self.M_s_in[7]
        #Calculation of outlet mass flow rate for water 
        self.M_h2o_s_out=self.M_h2o_s_in
        #Calculation total outlet mass flow rate 
        self.M_s_out_t=self.M_s_out[0]+self.M_s_out[1]+self.M_h2o_s_out#+self.M_s_out[8]+self.M_s_out[7]
        #Calculation of outlet mass flow rate for other ionic species in channel 
        for i in range(2,7):
            self.M_s_out[i]=self.M_s_in[i]
            self.M_s_out_t=self.M_s_out_t+self.M_s_out[i]

        #Calculate of volumetric outlet flow rate
        self.Q1_s_out=self.M_s_out_t/self.d_s
        
        #Calculation of outlet concentration of single ions in channel 
        for i in range(0,9):
            self.Ci_s_out[i]=self.M_s_out[i]/(self.Q1_s_out*self.PM_i[i]/1000)     

        
        #Calculation of voltage needed and electromotive force 
        self.c1_emf=math.log(((self.Ci_a_in[7]+self.Ci_a_out[7])/2)/Cm_bp_H)        
        self.c2_emf=-math.log(Cm_bp_OH/((self.Ci_b_in[8]+self.Ci_b_out[8])/2))        
        self.c3_emf=-math.log(((self.Ci_s_in[1]+self.Ci_s_out[1])/2)/((self.Ci_a_in[1]+self.Ci_a_out[1])/2))
        self.c4_emf=math.log(((self.Ci_b_in[0]+self.Ci_b_out[0])/2)/((self.Ci_s_in[0]+self.Ci_s_out[0])/2))  
        
            #EMF: electromotive force (V) 
        F=96485.3 #Coulombs/mol
        R_const=8.314462618 #kg⋅m2⋅s−2⋅K−1⋅mol−1
        # R_int=0.28 #ohm cm2
        R_int=45 #ohm cm2
        z=1
        self.EMF=((R_const*self.T/z/F)*(self.c1_emf+self.c2_emf+self.c3_emf+self.c4_emf))*self.N_trip   
    
            #V_ext: Voltage needed (V) 
        self.V_ext=self.EMF+((self.I_ext*R_int)/(self.A*10000))*self.N_trip
        
            #Calculate gross power needed (W)
        self.P=self.V_ext*self.I_ext         
