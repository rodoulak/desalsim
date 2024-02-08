import math 
import numpy as np
import scaleup 
import density_calc
import constants

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
    recycl_below1M(): 
        Recycling method for concentrations below 1M.
    recycl_1M():
        Recycling method for concentrations at 1M.
    """
    def __init__(self, Qin, Cc1, Cc2, Cc3, Cc4, Cc5, Cc6, Cc7, CH, COH, N,Cb1, Cb2,Cb3,Cb4,Cb5,Cb6,Cb7,Cb8,Cb9,Ca1,Ca2,Ca3,Ca4,Ca5,Ca6,Ca7,Ca8,Ca9):
        #input conditions   
        self.I_ext=A*I_d
        self.JA=3.6*self.I_ext/F
        #Initialization of parameters
        self.Qin=Qin
        self.CNa_s_in=Cc1/MW_Na #Na concentration (g/l) to (mol/l)
        self.CCl_s_in=Cc2/MW_Cl #Cl concentration (g/l) to (mol/l)
        self.CK_s_in=Cc3/MW_K #K concentration (g/l) to (mol/l)
        self.CMg_s_in=Cc4/MW_Mg #Mg concentration (g/l) to (mol/l)
        self.CCa_s_in=Cc5/MW_Ca #Ca concentration (g/l) to (mol/l)
        self.CSO4_s_in=Cc6/MW_SO4 #SO4 concentration (g/l) to (mol/l)
        self.CHCO3_s_in=Cc6/MW_HCO3 #HCO3 concentration (g/l) to (mol/l)      
        self.COH_s_in=COH/MW_OH
        self.PM_i=[MW_Na, MW_Cl, MW_K, MW_Mg, MW_Ca, MW_SO4, MW_HCO3, MW_H, MW_OH]
        self.N_trip=N 
        self.CH_s_in=CH/MW_H
        #Create lists for initial concentration in each channel 
        self.Ci_s_in=[self.CNa_s_in,self.CCl_s_in,self.CK_s_in,self.CMg_s_in,self.CCa_s_in,self.CSO4_s_in,self.CHCO3_s_in, self.CH_s_in,self.COH_s_in]
        self.Ci_b_in=[Cb1/MW_Na, Cb2/ MW_Cl,Cb3/MW_K,Cb4/MW_Mg,Cb5/MW_Ca,Cb6/MW_SO4,Cb7/MW_HCO3,1e-7,1e-7]
        self.Ci_a_in=[Ca1/MW_Na,Ca2/MW_Cl,Ca3/MW_K,Ca4/MW_Mg,Ca5/MW_Ca,Ca6/MW_SO4,Ca7/MW_HCO3,1e-7,1e-7]
        
        #save initial concentration before recirculation: 
        self.Ci_s_in_0=self.Ci_s_in
        self.Ci_b_in_0=self.Ci_b_in
        self.Ci_a_in_0=self.Ci_a_in

        
    def flowrate(self):
        # Calculates the flowrate in each channel 
        self.Q1_s_in=self.Qin/self.N_trip #units: l/h
        self.Q1_b_in=self.Qin/self.N_trip #units: l/h
        self.Q1_a_in=self.Qin/self.N_trip #units: l/h

    def in_mass_flow_rates(self):
        #Calculates the inlet mass flow rates of each ion, kg/h 
        #initialize lists 
        self.M_s_in=[[],[],[],[],[],[],[],[],[]]
        self.M_b_in=[[],[],[],[],[],[],[],[],[]]
        self.M_a_in=[[],[],[],[],[],[],[],[],[]]
                
        for i in range(0,9):
            self.M_s_in[i]=self.Q1_s_in*self.Ci_s_in[i]*self.PM_i[i]/1000 #units: kg/h
            self.M_b_in[i]=self.Q1_b_in*self.Ci_b_in[i]*self.PM_i[i]/1000 #units: kg/h
            self.M_a_in[i]=self.Q1_a_in*self.Ci_a_in[i]*self.PM_i[i]/1000 #units: kg/h
    
        self.M_h2o_a_in=self.Q1_a_in*d_a-sum(self.M_a_in)          
        self.M_h2o_b_in=self.Q1_b_in*d_b-sum(self.M_b_in) 
        self.M_h2o_s_in=self.Q1_s_in*d_s-sum(self.M_s_in) 
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
            
        #Calculate of volumetric outlet flow rate
        self.Q1_b_out=self.M_b_out_t/1 #assume density is 1kg/l
        
        #Calculation of outlet concentration of single ions in channel 
        for i in range(0,9):
            self.Ci_b_out[i]=self.M_b_out[i]/(self.Q1_b_out*self.PM_i[i]/1000)

    
    def salt_channel(self):
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
        self.M_s_out_t=self.M_s_out[0]+self.M_s_out[1]+self.M_s_out[8]+self.M_s_out[7]+self.M_h2o_s_out
        
        #Calculation of outlet mass flow rate for other ionic species in channel 
        for i in range(2,7):
            self.M_s_out[i]=self.M_s_in[i]
            self.M_s_out_t=self.M_s_out_t+self.M_s_out[i]
            
        #Calculate of volumetric outlet flow rate
        self.Q1_s_out=self.M_s_out_t/1.1
        
        #Calculation of outlet concentration of single ions in channel 
        for i in range(0,9):
            self.Ci_s_out[i]=self.M_s_out[i]/(self.Q1_s_out*self.PM_i[i]/1000)     

        
        #Calculation of voltage needed and electromotive force 
        self.c1_emf=math.log(((self.Ci_a_in[7]+self.Ci_a_out[7])/2)/Cm_bp_H)        
        self.c2_emf=math.log(Cm_bp_OH/((self.Ci_b_in[8]+self.Ci_b_out[8])/2))        
        self.c3_emf=math.log(((self.Ci_s_in[1]+self.Ci_s_out[1])/2)/((self.Ci_a_in[1]+self.Ci_a_out[1])/2))
        self.c4_emf=math.log(((self.Ci_b_in[0]+self.Ci_b_out[0])/2)/((self.Ci_s_in[0]+self.Ci_s_out[0])/2))  
        
            #EMF: electromotive force (V) 
        self.EMF=((R_const*T/z/F)*(self.c1_emf+self.c2_emf+self.c3_emf+self.c4_emf))*self.N_trip   
        
            #V_ext: Voltage needed (V) 
        self.V_ext=self.EMF+((self.I_ext*R_int)/(A*10000))*self.N_trip
        
            #Calculate gross power needed (W)
        self.P=self.V_ext*self.I_ext 
        
    def recycl_below1M(self):
        #Initialize the mass flow rate and concentration in each channel for recycling mode  
        self.M_s_in=self.M_s_out
        self.M_b_in=self.M_b_out
        self.M_a_in=self.M_a_out
        self.M_s_in_t=self.M_s_out_t
        self.M_a_in_t=self.M_a_out_t
        self.M_b_in_t=self.M_b_out_t
        self.Ci_s_in=self.Ci_s_out
        self.Ci_b_in=self.Ci_b_out
        self.Ci_a_in=self.Ci_a_out
        
        #Calculate density of outlet in each channel, kg/l
        self.d_a_out=density_calc.density_calc(25, sum(self.Ci_a_out))/1000
        self.d_b_out=density_calc.density_calc(25, sum(self.Ci_b_out))/1000
        self.d_s_out=density_calc.density_calc(25, sum(self.Ci_s_out))/1000
        

    def recycl_1M(self):
        #Calculate density of outlet in each channel, kg/l
        self.d_s_out=density_calc.density_calc(25, sum(self.Ci_s_out))/1000
        self.d_a_out=density_calc.density_calc(25, sum(self.Ci_a_out))/1000
        self.d_b_out=density_calc.density_calc(25, sum(self.Ci_b_out))/1000
        
        #Calculate recycling, outlet and new feed flow rate in each channel 
            #Salt channel 
        self.Q_s_r=rs*self.M_s_out_t/self.d_s_out
        self.Q_s_t_out=self.Q1_s_in/(1-rs)
        self.Q_s_out=self.Q1_s_in
        self.Q_s_in_t=self.Q1_s_in+self.Q_s_r
        self.M_s_r=self.Q_s_r/self.d_s_out
        self.M_s_t_out=self.Q_s_t_out/self.d_s_out
        self.M_s_in_t=self.Q1_s_in/d_in+self.Q_s_r/self.d_s_out

            #Acid channel 
        self.Q_a_r=r*self.M_a_out_t/self.d_a_out
        self.Q_a_t_out=self.Q_a_r/r
        self.Q_a_out=self.Q1_a_in
        self.Q_a_in_t=self.Q1_a_in+self.Q_a_r
        self.M_a_r=self.Q_a_r/self.d_a_out
        self.M_a_t_out=self.Q_a_t_out/self.d_a_out
        self.M_a_in_t=self.Q1_a_in/d_in+self.Q_a_r/self.d_a_out

            #Base channel 
        self.Q_b_r=r*self.M_b_out_t/self.d_b_out
        self.Q_b_t_out=self.Q1_b_in/(1-r)
        self.Q_b_out=self.Q1_b_in
        self.Q_b_in_t=self.Q1_b_in+self.Q_b_r
        self.M_b_r=self.Q_b_r/self.d_b_out
        self.M_b_t_out=self.Q_b_t_out/self.d_b_out
        self.M_b_in_t=self.Q1_b_in/d_in+self.Q_b_r/self.d_b_out       

        #Calculation of new concentration 
        for i in range(len(self.Ci_a_in)):
            self.Ci_s_in[i]=self.Ci_s_out[i]*self.M_s_r/self.M_s_in_t+self.Ci_s_in_0[i]*(self.M_s_in_t-self.M_s_r)/(self.M_s_in_t)      
            self.M_s_in[i]=self.Ci_s_in[i]*self.Q_s_in_t
            self.Ci_b_in[i]=self.Ci_b_out[i]*self.M_b_r/self.M_b_in_t+self.Ci_b_in_0[i]*(self.M_b_in_t-self.M_b_r)/(self.M_b_in_t) 
            self.M_b_in[i]=self.Ci_a_in[i]*self.Q_b_in_t
            self.Ci_a_in[i]=self.Ci_a_out[i]*self.Q_a_r/self.Q_a_in_t+self.Ci_a_in_0[i]*(self.Q1_a_in)/(self.Q_a_in_t) 
            self.M_a_in[i]=self.Ci_a_in[i]*self.Q_a_in_t

        #Calculation of final mass flow rates and ions mass flow rates
        self.M_b_out_f=self.Q_b_out*self.d_b_out
        self.M_a_out_f=self.Q_a_out*self.d_a_out
        self.M_s_out_f=self.Q_s_out*self.d_s_out
        
        self.M_s_in=[[],[],[],[],[],[],[],[],[]]
        self.M_b_in=[[],[],[],[],[],[],[],[],[]]
        self.M_a_in=[[],[],[],[],[],[],[],[],[]]
                
        for i in range(0,9):
            self.M_s_in[i]=self.Q_s_in_t*self.Ci_s_in[i]*self.PM_i[i]/1000 #units: kg/h
            self.M_b_in[i]=self.Q_b_in_t*self.Ci_b_in[i]*self.PM_i[i]/1000 #units: kg/h
            self.M_a_in[i]=self.Q_a_in_t*self.Ci_a_in[i]*self.PM_i[i]/1000 #units: kg/h

    
        self.M_h2o_a_in=self.Q_a_in_t*d_a-sum(self.M_a_in)          
        self.M_h2o_b_in=self.Q_b_in_t*d_b-sum(self.M_b_in) 
        self.M_h2o_s_in=self.Q_s_in_t*d_s-sum(self.M_s_in) 
        self.KW_s_in=self.Ci_s_in[7]*self.Ci_s_in[8]
        self.KW_a_in=self.Ci_a_in[7]*self.Ci_a_in[8]
        self.KW_b_in=self.Ci_b_in[7]*self.Ci_b_in[8]
    
        if ph_s<7:
            self.M_s_in[1]=self.Q_s_in_t*self.Ci_s_in[1]*self.PM_i[1]/1000+(self.M_s_in[7]/self.PM_i[7]-self.M_s_in[8]/self.PM_i[8])*(self.PM_i[1])
        elif ph_s>7:
            self.M_s_in[1]=self.Q_s_in_t*self.Ci_s_in[1]*self.PM_i[1]/1000+(self.M_s_in[8]/self.PM_i[8]-self.M_s_in[7]/self.PM_i[7])*self.PM_i[1]
            
#%% Example usage 
#constants
I_d=400 # The electricl current desnity Am2

#A=0.16 #active area of the membrane across which ion permeation occurs (A)

F=96485.3 #Coulombs/mol
R_const=8.314462618 #kg⋅m2⋅s−2⋅K−1⋅mol−1
# R_int=0.28 #ohm cm2
R_int=45 #ohm cm2
z=1
npump=0.8 #pump efficiency (units: -)
dp_f=0.1 #pressure drop for feed (units: bar)
dp_r=1.2 #pressure drop for recycling (units: bar)
dp_p=1 #pressure drop for extracting products (units: bar)

#Membrane characteristics

Cm_bp_H= 0.0000001 #mol/l 
Cm_bp_OH= 0.0000001 #mol/l 

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

T=20+273.15

#input data
r=1-0.4/(5-0.4) #recycling rate 
rs=1-1/(5-1) #recycling rate 
ph_s=4.71 #pH salt channel (units: -)
ph_b=7#pH base channel (units: -)
ph_a=7#pH acid channel (units: -)
d_a=1 #density acid channel (units: kg/l)
d_s=1 #density salt channel (units: kg/l)
d_b=1 #density base channel (units: kg/l)


#Feed concentration g/L
Cin_edbm=[13.44, 20.725, 1.146, 0, 0, 0.18]
d_in=density_calc.density_calc(25,sum(Cin_edbm))/1000

#Feed flow rate L/h
Q_in_edbm=1000

#Calculate water quantity in inflow 
Mw_in=Q_in_edbm/d_in 

#Set number of triplets 
N_trip=50

#Set membrane area based on the feed flow rate 
A=0.4#scaleup.scaleup(19.2, 300, Q_in_edbm)/3 

"or"
# N_trip=round(scaleup.scaleup(20, 300, Q_in_edbm),0)       
# A=0.4#19/40
# A_tot=19*Q_in_edbm/60 #scaleup.scaleup(0.19, 100, Q_in_edbm)#19*Q_in_edbm/60/3#
# N_trip=round(A_tot/A,0)/10

#Initialize concentration of Na in salt channel 
Cna_s=[]

#Create an instance of the EDBMCalc class with the defined parameters
edbm_dat=EDBMCalc(Q_in_edbm, Cin_edbm[0], Cin_edbm[1], Cin_edbm[2], Cin_edbm[3], Cin_edbm[4], Cin_edbm[5], 0,  10**(-ph_s), 3.01551E-11, round(N_trip,0), 0,0,0,0,0,0,0,10**(-ph_b), 10**(-(14-ph_b)), 0,0,0,0,0,0,0,10**(-ph_a), 10**(-(14-ph_a)))

# Call the necessary methods to calculate values
edbm_dat.flowrate()
edbm_dat.in_mass_flow_rates()
edbm_dat.acid_channel()
edbm_dat.base_channel()
edbm_dat.salt_channel()
Cna_s.append(edbm_dat.Ci_s_out[0])

#Loop for the recycling process 
while (edbm_dat.Ci_b_out[0]<1) and (edbm_dat.Ci_s_out[0]>=0.05):
    edbm_dat.recycl_below1M()
    edbm_dat.acid_channel()
    edbm_dat.base_channel()
    edbm_dat.salt_channel()
    Cna_s.append(edbm_dat.Ci_s_out[0])
    if edbm_dat.Ci_s_out[0]<0:

        break

# Loop for recycling at higher concentration
for i in range(1): 
    edbm_dat.recycl_1M()
    while (edbm_dat.Ci_b_out[0]<1) and (edbm_dat.Ci_s_out[0]>=0.05):
        edbm_dat.recycl_below1M()
        edbm_dat.acid_channel()
        edbm_dat.base_channel()
        edbm_dat.salt_channel()
        Cna_s.append(edbm_dat.Ci_s_out[0])
        
        if edbm_dat.Ci_s_out[0]<0:

              break

# Sum results

"Salt channel "
    #Concentration in salt channel 
Cbrine_out_t=sum(edbm_dat.Ci_s_out)
Cbrine_out=edbm_dat.Ci_s_out #mol/l
Cbrine_out_g=[Cbrine_out[0]*MW_Na, Cbrine_out[1]*MW_Cl, Cbrine_out[2]*MW_K, Cbrine_out[3]*MW_Mg, Cbrine_out[4]*MW_Ca, Cbrine_out[5]*MW_SO4] #g/l

    #Mass flow rate
M_s_out=edbm_dat.M_s_out_f*N_trip
print("M_s_out is "+str(M_s_out)+"kg/hr")
M_s_out_r=edbm_dat.M_s_r*N_trip
print("M_s_out recycling is "+str(M_s_out_r)+"kg/hr")
    #Volumetric flow rate 
Q_s_out=edbm_dat.Q_s_out*N_trip
print("Q_s_out is "+str(Q_s_out)+"l/hr")

"Base channel "
    #Concentration in base channel 
Cb_out=edbm_dat.Ci_b_out[0:6]
Cb_out_g=[Cb_out[0]*MW_Na, Cb_out[1]*MW_Cl, Cb_out[2]*MW_K, Cb_out[3]*MW_Mg, Cb_out[4]*MW_Ca, Cb_out[5]*MW_SO4]

    #Mass flow rate
M_b_out=edbm_dat.M_b_out_f*N_trip
print("M_b_out is "+str(M_b_out)+"kg/hr")
M_b_out_r=edbm_dat.M_b_r*N_trip

    #Volumetric flow rate 
Q_b_out=edbm_dat.Q_b_out*N_trip
print("Q_b_out is "+str(Q_b_out)+"l/hr")

    #Conversion to solid 
M_NaOH_out=Q_b_out*edbm_dat.Ci_b_out[0]*constants.MW_NaOH/1000 #kg/hr 

"Acid channel" 
    #Concentration in acid channel 
Ca_out=edbm_dat.Ci_a_out
Ca_out=edbm_dat.Ci_a_out[0:6]
Ca_out_g=[Ca_out[0]*MW_Na, Ca_out[1]*MW_Cl, Ca_out[2]*MW_K, Ca_out[3]*MW_Mg, Ca_out[4]*MW_Ca, Ca_out[5]*MW_SO4]

    #Mass flow rate 
M_a_out=edbm_dat.M_a_out_f*N_trip
print("M_a_out is "+str(M_a_out)+"kg/hr")
M_a_out_r=edbm_dat.M_a_r*N_trip
print("M_a_out recycling is "+str(M_a_out_r)+"kg/hr")

    #Volumetric flow rate 
Q_a_out=edbm_dat.Q_a_out*N_trip
print("Q_a_out is "+str(Q_a_out)+"l/hr")

    #Conversion to solid 
M_HCl_out=Q_a_out*constants.MW_HCl/1000 #kg/hr


#Calculate required amount of water for the operation mode 
Q_w_in=2*Q_in_edbm

#Calculate mass balance 
bal=(Q_in_edbm*d_in+2*Q_in_edbm)-M_a_out-M_s_out-M_b_out
print("balance edbm is "+str(bal))

#Energy consumption 

#Calculate energy consumption for pumping 
Ppump=(3*Q_in_edbm*dp_f+Q_s_out*dp_p+Q_b_out*dp_p+Q_a_out*dp_p+edbm_dat.Q_s_r*N_trip*dp_r+edbm_dat.Q_a_r*N_trip*dp_r+edbm_dat.Q_b_r*N_trip*dp_r)/1000/3600*1e5/npump #units: W -> l to m3 so /1000; bar to J 1e5N/m2*1J/m ; hr to 3660s  

#Total energy consumption 
E_el_Edbm=48*I_d*A/1000+Ppump/1000
SEC=(48*I_d*A)/(Q_b_out*(edbm_dat.Ci_b_out[0]-edbm_dat.Ci_b_in[0])*constants.MW_NaOH)

print("Total electrical consumption for EDBM is " + str(E_el_Edbm)+ " KW")
print("sec is "+str(SEC)+"kwh/kg naoh")

#%% sec based on calogero paper
#SEC= U*Id*A/(60*Qp*(Cp,out-Cp,in)*Mp) where -> U: electrical potential, Id: current density, Cp: product concentration, Mp: molar mass of desired product (g/mol)
SEC2=(48*I_d*A)/(Q_b_out*(edbm_dat.Ci_b_out[0]-edbm_dat.Ci_b_in[0])*constants.MW_NaOH)
print("sec2 is " +str(SEC2)) 
        