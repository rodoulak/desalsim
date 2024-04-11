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
        
        print("self.M_s_in is "+str((self.M_s_in)))
        print("self.M_a_in is "+str(self.M_a_in))
        print("self.M_b_in is "+str(self.M_b_in))
        print ("       ")
        self.M_h2o_a_in=self.Q1_a_in*d_a-sum(self.M_a_in)   
        print("M_h2o_a_in initial is "+str(self.M_h2o_a_in))
        self.M_h2o_b_in=self.Q1_b_in*d_b-sum(self.M_b_in) 
        print("M_h2o_b_in initial is "+str(self.M_h2o_b_in))
        self.M_h2o_s_in=self.Q1_s_in*d_s-sum(self.M_s_in) 
        print("M_h2o_s_in initial is "+str(self.M_h2o_s_in))
        print ("       ")
        #Calculate inlet ionic water product in each channel 
        self.KW_s_in=self.Ci_s_in[7]*self.Ci_s_in[8] 
        self.KW_a_in=self.Ci_a_in[7]*self.Ci_a_in[8]
        self.KW_b_in=self.Ci_b_in[7]*self.Ci_b_in[8]
    
        #Calculation of the inlet concentration of the Cl in salt channel (exception)
        if ph_s<7:
            self.M_s_in[1]=self.Q1_s_in*self.Ci_s_in[1]*self.PM_i[1]/1000+(self.M_s_in[7]/self.PM_i[7]-self.M_s_in[8]/self.PM_i[8])*(self.PM_i[1])
        elif ph_s>7:
            self.M_s_in[1]=self.Q1_s_in*self.Ci_s_in[1]*self.PM_i[1]/1000+(self.M_s_in[8]/self.PM_i[8]-self.M_s_in[7]/self.PM_i[7])*self.PM_i[1]
        
        self.M_s_in_0=self.Q1_s_in*d_in #units: kg/h
        self.M_b_in_0= self.Q1_a_in#units: kg/h
        self.M_a_in_0= self.Q1_b_in#units: kg/h

    def acid_channel(self):
        #Performs mass balance calculations for Acid channel 
        #Initialize ionic Mass and concentration in acid channel  
        self.M_a_out=[[],[],[],[],[],[],[],[],[]]
        self.Ci_a_out=[[],[],[],[],[],[],[],[],[]]
        
        #Calculation of outlet mass flow rate for H+
        self.M_a_out[7]=self.M_a_in[7]+self.JA*self.PM_i[7]
        
        #Calculation of outlet mass flow rate for Cl-
        self.M_a_out[1]=self.M_a_in[1]+self.JA*self.PM_i[1]
        print("self.JA*self.PM_i[1] is "+str(self.JA*self.PM_i[1]*N_trip))
        print("self.M_a_in[1] is "+str(self.M_a_in[1]*N_trip))
        print("self.M_a_out[1] is "+str(self.M_a_out[1]*N_trip))
        #Calculation of outlet mass flow rate for OH-
        self.M_a_out[8]=self.M_a_in[8]
        #Calculation of outlet mass flow rate for Na+
        self.M_a_out[0]=self.M_a_in[0]
        #Calculation of outlet mass flow rate for water 
        self.M_h2o_a_out=self.M_h2o_a_in-0.5*self.JA*MW_H2O
        #Calculation total outlet mass flow rate 
        self.M_a_out_t=self.M_h2o_a_out+self.M_a_out[7]+self.M_a_out[1]+self.M_a_out[8]+self.M_a_out[0]

        print("M_h2o_a_out is "+str(self.M_h2o_a_out*N_trip))
        print("self.M_a_out[0] is "+str(self.M_a_out[0]))
        print("self.M_a_out[8] is "+str(self.M_a_out[8]))
        print("self.M_a_in[1] is "+str(self.M_a_in[1]*N_trip))
        print("self.M_a_in[7] is "+str(self.M_a_in[7]*N_trip))
        print("self.M_a_out[1] is "+str(self.M_a_out[1]*N_trip))
        print("self.M_a_out[7] is "+str(self.M_a_out[7]*N_trip))
        #Calculation of outlet mass flow rate for other ionic species in channel 
        for i in range(2,7):
            self.M_a_out[i]=self.M_a_in[i]
            self.M_a_out_t=self.M_a_out_t+self.M_a_out[i]
        print("M_a_out_t is "+str(self.M_a_out_t*N_trip))
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
        print("M_h2o_b_out is "+str(self.M_h2o_b_out*N_trip))
        print("self.M_b_out[7] is "+str(self.M_b_out[7]))
        print("self.M_b_out[1] is "+str(self.M_b_out[1]))
        print("self.M_b_in[0] is "+str(self.M_b_in[0]*N_trip))
        print("self.M_b_out[0] is "+str(self.M_b_out[0]*N_trip))
        print("self.M_b_out[8] is "+str(self.M_b_out[8]*N_trip))
        print("sum of solids 1 in base channel is "+str((self.M_b_out[7]+self.M_b_out[1]+self.M_b_out[8]+self.M_b_out[0])*N_trip))
        #Calculation of outlet mass flow rate for other ionic species in channel 
        for i in range(2,7):
            self.M_b_out[i]=self.M_b_in[i]  
            self.M_b_out_t=self.M_b_out_t+self.M_b_out[i]
        
        print("sum of solids 2 in base channel is "+str(sum(self.M_b_out[2:7])*N_trip))
        print("M_b_out_t is "+str(self.M_b_out_t*N_trip))
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
        print("M_s_out[0] is "+str(round(self.M_s_out[0],2))
              +"M_s_in[0] is "+str(round(self.M_s_in[0],2))
              +"the difference is M_b_out[0]-M_b_in[0]"+str(round((self.M_b_out[0]-self.M_b_in[0]),2))
              )
        
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
        self.Q1_s_out=self.M_s_out_t/d_in
        
        #Calculation of outlet concentration of single ions in channel 
        for i in range(0,9):
            self.Ci_s_out[i]=self.M_s_out[i]/(self.Q1_s_out*self.PM_i[i]/1000)     

        if edbm_dat.Ci_s_out[0]>=0 and edbm_dat.Ci_s_out[1]>=0:
        #Calculation of voltage needed and electromotive force 
            self.c1_emf=math.log(((self.Ci_a_in[7]+self.Ci_a_out[7])/2)/Cm_bp_H)        
            self.c2_emf=-math.log(Cm_bp_OH/((self.Ci_b_in[8]+self.Ci_b_out[8])/2))        
            self.c3_emf=-math.log(((self.Ci_s_in[1]+self.Ci_s_out[1])/2)/((self.Ci_a_in[1]+self.Ci_a_out[1])/2))
            self.c4_emf=math.log(((self.Ci_b_in[0]+self.Ci_b_out[0])/2)/((self.Ci_s_in[0]+self.Ci_s_out[0])/2))  
            
                #EMF: electromotive force (V) 
            self.EMF=((R_const*T/z/F)*(self.c1_emf+self.c2_emf+self.c3_emf+self.c4_emf))*self.N_trip   
        
                #V_ext: Voltage needed (V) 
            self.V_ext=self.EMF+((self.I_ext*R_int)/(A*10000))*self.N_trip
            
                #Calculate gross power needed (W)
            self.P=self.V_ext*self.I_ext 
        
    # def backup_mode(self):
        
    def recycl_below1M(self):
        #Initialize the mass flow rate and concentration in each channel for recycling mode  
        self.M_s_in=self.M_s_out
        self.M_b_in=self.M_b_out
        self.M_a_in=self.M_a_out
        print("self.M_s_in recycling is "+str((self.M_s_in)))
        print("self.M_a_in recycling is "+str(self.M_a_in))
        print("self.M_b_in recycling is "+str(self.M_b_in))
        print ("       ")
        self.Ci_s_in=self.Ci_s_out
        self.Ci_b_in=self.Ci_b_out
        self.Ci_a_in=self.Ci_a_out
        
        self.Q1_a_in=self.Q1_a_out
        self.Q1_b_in=self.Q1_b_out
        self.Q1_s_in=self.Q1_s_out
        
        self.M_h2o_a_in=self.M_h2o_a_out#self.Q1_a_out*d_a-sum(self.M_a_in)   
        print("M_h2o_a_in recycling is "+str(self.M_h2o_a_in*N_trip))

        self.M_h2o_b_in=self.M_h2o_b_out#self.Q1_b_out*d_b-sum(self.M_b_in) 
        print("M_h2o_b_in recycling is "+str(self.M_h2o_b_in*N_trip))

        self.M_h2o_s_in=self.M_h2o_s_out
        print("M_h2o_s_in recycling is "+str(self.M_h2o_s_in*N_trip))
        print ("       ")
        
        #Calculate inlet ionic water product in each channel 
        self.KW_s_in=self.Ci_s_in[7]*self.Ci_s_in[8] 
        self.KW_a_in=self.Ci_a_in[7]*self.Ci_a_in[8]
        self.KW_b_in=self.Ci_b_in[7]*self.Ci_b_in[8]
    
        #Calculation of the inlet concentration of the Cl in salt channel (exception)
        if ph_s<7:
            self.M_s_in[1]=self.Q1_s_out*self.Ci_s_in[1]*self.PM_i[1]/1000+(self.M_s_in[7]/self.PM_i[7]-self.M_s_in[8]/self.PM_i[8])*(self.PM_i[1])
        elif ph_s>7:
            self.M_s_in[1]=self.Q1_s_out*self.Ci_s_in[1]*self.PM_i[1]/1000+(self.M_s_in[8]/self.PM_i[8]-self.M_s_in[7]/self.PM_i[7])*self.PM_i[1]
        
        

    def recycl_1M(self):
        #Calculate density of outlet in each channel, kg/l
        self.d_s_out=density_calc.density_calc(25, sum(self.Ci_s_out))/1000
        self.d_a_out=density_calc.density_calc(25, sum(self.Ci_a_out))/1000
        self.d_b_out=density_calc.density_calc(25, sum(self.Ci_b_out))/1000
        
        #Calculate recycling, outlet and new feed flow rate in each channel 
            #Salt channel 
        #Calculate recycling flow rate in salt stream 
        self.Q_s_r=self.M_s_out_t/self.d_s_out/(1/rs-1)
        #Calculate total outlet (including recycling stream) from salt stream 
        self.Q_s_t_out=self.M_s_out_t/self.d_s_out
        #self.Q_s_t_out=self.Q1_s_in*(1+rs)

        #Calculate total inflow rate (including recycling)
        self.Q_s_in_t=self.Q1_s_in+self.Q_s_r
        #Calculate final outlet flowrate in salt stream 
        self.Q_s_out=self.Q_s_in_t-self.Q_s_r
        #Calculate the mass flowrates 
        self.M_s_r=self.M_s_out_t/(1/rs-1)
        self.M_s_in_t=self.Q1_s_in*d_in+self.Q_s_r*self.d_s_out
        

            #Acid channel 
        #Calculate recycling flow rate in acid stream 
        self.Q_a_r=self.M_a_out_t/self.d_a_out/(1/r-1)
        #Calculate total outlet (including recycling stream) from acid stream
        self.Q_a_t_out=self.Q1_a_in+self.Q_a_r
        #self.Q_a_t_out=self.Q_a_r*(1+r)
        #Calculate total inflow rate (including recycling) to acid stream 
        self.Q_a_in_t=self.Q1_a_in+self.Q_a_r
        #Calculate final outlet flowrate in acid stream 
        self.Q_a_out=self.Q_a_in_t-self.Q_a_r
        #Calculate the mass flowrates 
        self.M_a_r=self.M_a_out_t/(1/r-1)
        self.M_a_in_t=self.Q1_a_in*d_in+self.Q_a_r*self.d_a_out


            #Base channel 
        #Calculate recycling flow rate in base stream 
        self.Q_b_r=self.M_b_out_t/self.d_b_out/(1/r-1)
        #Calculate total outlet (including recycling stream) from base stream
        self.Q_b_t_out=self.Q1_b_in+self.Q_b_r
        #self.Q_b_t_out=self.Q1_b_in*(1+r)
        #Calculate total inflow rate (including recycling) to base stream 
        self.Q_b_in_t=self.Q1_b_in+self.Q_b_r
        #Calculate final outlet flowrate in base stream 
        self.Q_b_out=self.Q_b_in_t-self.Q_b_r
        #Calculate the mass flowrates 
        self.M_b_r=self.M_b_out_t/(1/r-1)
        self.M_b_in_t=self.Q1_b_in*d_in+self.Q_b_r*self.d_b_out      

        #Calculation of new concentration 
        for i in range(len(self.Ci_a_in)):
            #self.Ci_s_in[i]=self.Ci_s_out[i]*self.M_s_r/self.M_s_in_t+self.Ci_s_in_0[i]*(self.M_s_in_t-self.M_s_r)/(self.M_s_in_t)      
            self.Ci_s_in[i]=self.Ci_s_out[i]*self.Q_s_r/self.Q_s_in_t+self.Ci_s_in_0[i]*(self.Q1_s_in)/(self.M_s_in_t)  
            self.M_s_in[i]=self.Ci_s_in[i]*self.Q_s_in_t
            self.Ci_b_in[i]=self.Ci_b_out[i]*self.Q_b_r/self.Q_b_in_t+self.Ci_b_in_0[i]*(self.Q1_b_in)/(self.M_b_in_t) 
            self.M_b_in[i]=self.Ci_a_in[i]*self.Q_b_in_t
            self.Ci_a_in[i]=self.Ci_a_out[i]*self.Q_a_r/self.Q_a_in_t+self.Ci_a_in_0[i]*(self.Q1_a_in)/(self.Q_a_in_t) 
            self.M_a_in[i]=self.Ci_a_in[i]*self.Q_a_in_t
        print("Ci_s_in new concentration is "+str((self.Ci_s_in)))
        print("Ci_a_in new concentration is "+str((self.Ci_a_in)))
        print("Ci_b_in new concentration is "+str((self.Ci_b_in)))
        print("              ")
        print("Mi_s_in new concentration is "+str((self.M_s_in)))
        print("Mi_a_in new concentration is "+str((self.M_a_in)))
        print("Mi_b_in new concentration is "+str((self.M_b_in)))
        print("              ")
        #Calculation of final mass flow rates and ions mass flow rates
        # self.M_b_out_f=self.Q_b_out*self.d_b_out
        # self.M_a_out_f=self.Q_a_out*self.d_a_out
        # self.M_s_out_f=self.Q_s_out*self.d_s_out
        
        self.M_s_in=[[],[],[],[],[],[],[],[],[]]
        self.M_b_in=[[],[],[],[],[],[],[],[],[]]
        self.M_a_in=[[],[],[],[],[],[],[],[],[]]
                
        for i in range(0,9):
            self.M_s_in[i]=self.Q_s_in_t*self.Ci_s_in[i]*self.PM_i[i]/1000 #units: kg/h
            self.M_b_in[i]=self.Q_b_in_t*self.Ci_b_in[i]*self.PM_i[i]/1000 #units: kg/h
            self.M_a_in[i]=self.Q_a_in_t*self.Ci_a_in[i]*self.PM_i[i]/1000 #units: kg/h
    
        self.M_h2o_a_in=self.M_a_in_t-sum(self.M_a_in)          
        self.M_h2o_b_in=self.M_b_in_t-sum(self.M_b_in) 
        self.M_h2o_s_in=self.M_s_in_t-sum(self.M_s_in) 
        
        self.KW_s_in=self.Ci_s_in[7]*self.Ci_s_in[8]
        self.KW_a_in=self.Ci_a_in[7]*self.Ci_a_in[8]
        self.KW_b_in=self.Ci_b_in[7]*self.Ci_b_in[8]
    
        if ph_s<7:
            self.M_s_in[1]=self.Q_s_in_t*self.Ci_s_in[1]*self.PM_i[1]/1000+(self.M_s_in[7]/self.PM_i[7]-self.M_s_in[8]/self.PM_i[8])*(self.PM_i[1])
        elif ph_s>7:
            self.M_s_in[1]=self.Q_s_in_t*self.Ci_s_in[1]*self.PM_i[1]/1000+(self.M_s_in[8]/self.PM_i[8]-self.M_s_in[7]/self.PM_i[7])*self.PM_i[1]
        
    def recycl_init(self):
        #Initialize the mass flow rate and concentration in each channel for recycling mode  
        # self.M_s_in=self.M_s_in
        # self.M_b_in=self.M_b_in
        # self.M_a_in=self.M_a_in_t
        self.Ci_s_in=self.Ci_s_in
        self.Ci_b_in=self.Ci_b_in
        self.Ci_a_in=self.Ci_a_in
    
    def recycl_final(self):
        self.M_s_out=self.M_s_out_t*(1-rs)
        self.M_a_out=self.M_a_out_t*(1-r)
        self.M_b_out=self.M_b_out_t*(1-r)
        
    def backup_mode(self):
        #acid channel
        self.previous_a=[self.Q1_a_out, self.M_a_out_t, self.Ci_a_out, self.M_a_out]
        self.previous_b=[self.Q1_b_out, self.M_b_out_t, self.Ci_b_out, self.M_b_out]
        self.previous_s=[self.Q1_s_out, self.M_s_out_t, self.Ci_s_out, self.M_s_out]  
        
    def restore_from_backup(self):
        # Restore previous iteration results
        self.Q1_a_out, self.M_a_out_t, self.Ci_a_out, self.M_a_out = self.previous_a
        self.Q1_b_out, self.M_b_out_t, self.Ci_b_out, self.M_b_out = self.previous_b
        self.Q1_s_out, self.M_s_out_t, self.Ci_s_out, self.M_s_out = self.previous_s
#%% Example usage 
#constants
I_d=400 # The electrical current desnity Am2

F=96485.3 #Coulombs/mol
R_const=8.314462618 #kg⋅m2⋅s−2⋅K−1⋅mol−1
# R_int=0.28 #ohm cm2
R_int=45 #ohm cm2
z=1
npump=0.8 #pump efficiency (units: -)
dp=1 #pressure drop (units: bar)

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
r=0.91 #recycling rate 
rs=0.75 #recycling rate 
ph_s=4.71 #pH salt channel (units: -)
ph_b=7#pH base channel (units: -)
ph_a=7#pH acid channel (units: -)
d_a=0.997 #density acid channel (units: kg/l)
#d_s=1 #density salt channel (units: kg/l)
d_b=0.997 #density base channel (units: kg/l)


#Feed concentration g/L
Cin_edbm=[17.0, 19.6, 0.303, 0, 0, 0.18]#[17.0, 19.6, 0.303, 0, 0, 0.18]#13.44
d_in=density_calc.density_calc(25,sum(Cin_edbm))/1000
d_s=d_in
#Feed flow rate L/h
Q_in_edbm=47000

#Calculate water quantity in inflow 
Mw_in=Q_in_edbm*d_in 

#Set number of triplets 
N_trip=50*47 #range: 20-200 triplets based on the inlet flow rate

#Set membrane area based on the feed flow rate, m2 
A=0.4 #range: 0.1-1

#Initialize concentration of Na in salt channel 
Cna_s=[]

#Create an instance of the EDBMCalc class with the defined parameters
edbm_dat=EDBMCalc(Q_in_edbm, Cin_edbm[0], Cin_edbm[1], Cin_edbm[2], Cin_edbm[3], Cin_edbm[4], Cin_edbm[5], 0,  10**(-ph_s), 3.01551E-11, N_trip, 0,0,0,0,0,0,0,10**(-ph_b), 10**(-(14-ph_b)), 0,0,0,0,0,0,0,10**(-ph_a), 10**(-(14-ph_a)))

# Call the necessary methods to calculate values
edbm_dat.flowrate()
edbm_dat.in_mass_flow_rates()
edbm_dat.acid_channel()
edbm_dat.base_channel()
edbm_dat.salt_channel()
Cna_s.append(edbm_dat.Ci_s_out[0])


Cbrine_out_t=sum(edbm_dat.Ci_s_out[0:6])
Cbrine_out=edbm_dat.Ci_s_out #mol/l
Cbrine_out_g=[Cbrine_out[0]*MW_Na, Cbrine_out[1]*MW_Cl, Cbrine_out[2]*MW_K, Cbrine_out[3]*MW_Mg, Cbrine_out[4]*MW_Ca, Cbrine_out[5]*MW_SO4] #g/l
print("salt channel Na concentration:"+str(round(Cbrine_out[0],2))+"M and "+str(round(Cbrine_out_g[0],2))+"g/l")
print("salt channel Cl concentration:"+str(round(Cbrine_out[1],2))+"M and "+str(round(Cbrine_out_g[1],2))+"g/l")
print ("M_s_out Na "+str(round(edbm_dat.M_s_out[0],2))+" and M_s_out Cl "+str(round(edbm_dat.M_s_out[1],2)))
print ("M_b_out Na "+str(round(edbm_dat.M_b_out[0],2)))
print ("M_a_out Cl "+str(round(edbm_dat.M_a_out[1],2)))
print("base channel Na concentration:"+str(round(edbm_dat.Ci_b_out[0],2))+"M")
print("acid channel Cl concentration:"+str(round(edbm_dat.Ci_a_out[1],2))+"M")
print("end of batch mode ")    
print("     ")

i=0
# #Loop for the recycling process (closed loop configuration)
print("edbm_dat.Ci_b_out[0] is "+str(edbm_dat.Ci_b_out[0]))
print("edbm_dat.Ci_s_out[0] is "+str(edbm_dat.Ci_s_out[0]))
print("edbm_dat.Ci_s_out[1] is "+str(edbm_dat.Ci_s_out[1]))
# previous_Cbrine_out_t=Cbrine_out_t
# previous_Cbrine_out=Cbrine_out
# previous_Cbrine_out_g=Cbrine_out_g
# previous_M_a_out=edbm_dat.M_a_out
# previous_M_a_out_t=edbm_dat.M_a_out_t
# previous_M_b_out=edbm_dat.M_b_out
# previous_M_b_out_t=edbm_dat.M_b_out_t
# previous_M_s_out=edbm_dat.M_s_out
# previous_M_s_out_t=edbm_dat.M_s_out_t
# previous_Q_s_out=edbm_dat.Q_s_out
edbm_dat.backup_mode()

while (edbm_dat.Ci_b_out[0]<1) or (edbm_dat.Ci_s_out[0]<=0) or (edbm_dat.Ci_s_out[1]<=0):
    i=i+1
    print("loop works "+str(i))    
    edbm_dat.recycl_below1M()
    edbm_dat.acid_channel()
    edbm_dat.base_channel()
    edbm_dat.salt_channel()
    Cna_s.append(edbm_dat.Ci_s_out[0])
    Cb_out=edbm_dat.Ci_b_out[0:6]
    Cb_out_g=[Cb_out[0]*MW_Na, Cb_out[1]*MW_Cl, Cb_out[2]*MW_K, Cb_out[3]*MW_Mg, Cb_out[4]*MW_Ca, Cb_out[5]*MW_SO4]
    print("base channel Na concentration "+str(round(Cb_out[0],2))+"M and "+str(round(Cb_out_g[0],2))+"g/l")
    if edbm_dat.Ci_s_out[0]<=0 or edbm_dat.Ci_s_out[1]<=0:
        print("negative values")
        edbm_dat.restore_from_backup()
        # Ci_b_out=edbm_dat.previous_b[2]
        # Ci_a_out=edbm_dat.previous_a[2]
        # Ci_s_out=edbm_dat.previous_s[2]
        # Q_a_out=edbm_dat.previous_a[0]
        # Q_b_out=edbm_dat.previous_b[0]
        # Q_s_out=edbm_dat.previous_s[0]
        # M_a_out_t=edbm_dat.previous_a[1]
        # M_b_out_t=edbm_dat.previous_b[1]
        # M_s_out_t=edbm_dat.previous_s[1]
        # M_a_out=edbm_dat.previous_a[3]
        # M_b_out=edbm_dat.previous_b[3]
        # M_s_out=edbm_dat.previous_s[3]
        print("M_s_out backup is "+str(M_s_out))
        break

Cbrine_out_t=sum(edbm_dat.Ci_s_out[0:6])
Cbrine_out=edbm_dat.Ci_s_out #mol/l
Cbrine_out_g=[Cbrine_out[0]*MW_Na, Cbrine_out[1]*MW_Cl, Cbrine_out[2]*MW_K, Cbrine_out[3]*MW_Mg, Cbrine_out[4]*MW_Ca, Cbrine_out[5]*MW_SO4] #g/l
print("salt channel Na concentration:"+str(round(Cbrine_out[0],2))+"M and "+str(round(Cbrine_out_g[0],2))+"g/l")
print("salt channel Cl concentration:"+str(round(Cbrine_out[1],2))+"M and "+str(round(Cbrine_out_g[1],2))+"g/l")
print ("M_s_out Na "+str(round(edbm_dat.M_s_out[0],2))+" and M_s_out Cl "+str(round(edbm_dat.M_s_out[1],2)))
    #Concentration in base channel 
Cb_out=edbm_dat.Ci_b_out[0:6]
Cb_out_g=[Cb_out[0]*MW_Na, Cb_out[1]*MW_Cl, Cb_out[2]*MW_K, Cb_out[3]*MW_Mg, Cb_out[4]*MW_Ca, Cb_out[5]*MW_SO4]
print("base channel Na concentration "+str(round(Cb_out[0],2))+"M and "+str(round(Cb_out_g[0],2))+"g/l")
print("     ")
Ca_out=edbm_dat.Ci_a_out[0:6]
Ca_out_g=[Ca_out[0]*MW_Na, Ca_out[1]*MW_Cl, Ca_out[2]*MW_K, Ca_out[3]*MW_Mg, Ca_out[4]*MW_Ca, Ca_out[5]*MW_SO4]
print("acid channel Cl concentration "+str(round(Ca_out[1],2))+"M and "+str(round(Ca_out_g[1],2))+"g/l")
M_s_out=edbm_dat.M_s_out_t*N_trip
print("M_s_out is "+str(M_s_out)+"kg/hr")
Q_s_out=edbm_dat.Q1_s_out*N_trip
print("Q_s_out is "+str(Q_s_out)+"l/hr")
M_b_out=edbm_dat.M_b_out_t*N_trip
print("M_b_out is "+str(M_b_out)+"kg/hr")
M_a_out=edbm_dat.M_a_out_t*N_trip
print("M_a_out is "+str(M_a_out)+"kg/hr")
print("     ")
bal=(Q_in_edbm*d_in+2*Q_in_edbm)-M_a_out-M_s_out-M_b_out
print("balance edbm is "+str(bal))
print("end of closed loop process ")    
print("     ")
# # Loop for recycling at higher concentration (Feed and Bleed Configuration)
# %%
edbm_dat.recycl_1M()
Cbrine_out_t=sum(edbm_dat.Ci_s_out[0:6])
Cbrine_out=edbm_dat.Ci_s_out #mol/l
Cbrine_out_g=[Cbrine_out[0]*MW_Na, Cbrine_out[1]*MW_Cl, Cbrine_out[2]*MW_K, Cbrine_out[3]*MW_Mg, Cbrine_out[4]*MW_Ca, Cbrine_out[5]*MW_SO4] #g/l
print("salt channel Na concentration:"+str(round(Cbrine_out[0],2))+"M and "+str(round(Cbrine_out_g[0],2))+"g/l")
print("salt channel Cl concentration:"+str(round(Cbrine_out[1],2))+"M and "+str(round(Cbrine_out_g[1],2))+"g/l")
print ("M_s_out Na "+str(round(edbm_dat.M_s_out[0],2))+" and M_s_out Cl "+str(round(edbm_dat.M_s_out[1],2)))
print("     ")
M_s_out=edbm_dat.M_s_out_t*N_trip
print("M_s_out is "+str(M_s_out)+"kg/hr")
Q_s_out=edbm_dat.Q_s_out*N_trip
print("Q_s_out is "+str(Q_s_out)+"l/hr")
M_b_out=edbm_dat.M_b_out_t*N_trip
print("M_b_out is "+str(M_b_out)+"kg/hr")
M_a_out=edbm_dat.M_a_out_t*N_trip
print("M_a_out is "+str(M_a_out)+"kg/hr")
print("     ")
bal=(Q_in_edbm*d_in+2*Q_in_edbm)-M_a_out-M_s_out-M_b_out
print("balance edbm is "+str(bal))
print("end of recycling 2 process ")    
print("     ")
#while (edbm_dat.Ci_b_out[0]<1) and (edbm_dat.Ci_s_out[0]>=0.2) and(edbm_dat.Ci_s_out[1]>=0.2) :
edbm_dat.recycl_init()
#edbm_dat.recycl_below1M()
edbm_dat.acid_channel()
edbm_dat.base_channel()
edbm_dat.salt_channel()
edbm_dat.recycl_final()
# Cna_s.append(edbm_dat.Ci_s_out[0])
        
        # #edbm_dat.recycl_1M()
        # if edbm_dat.Ci_s_out[0]<0:

        #       break

# Sum results

"Salt channel "
    #Concentration in salt channel 
Cbrine_out_t=sum(edbm_dat.Ci_s_out[0:6])
Cbrine_out=edbm_dat.Ci_s_out #mol/l
Cbrine_out_g=[Cbrine_out[0]*MW_Na, Cbrine_out[1]*MW_Cl, Cbrine_out[2]*MW_K, Cbrine_out[3]*MW_Mg, Cbrine_out[4]*MW_Ca, Cbrine_out[5]*MW_SO4] #g/l
print("salt channel Na concentration:"+str(round(Cbrine_out[0],2))+"M and "+str(round(Cbrine_out_g[0],2))+"g/l")
print("salt channel Cl concentration:"+str(round(Cbrine_out[1],2))+"M and "+str(round(Cbrine_out_g[1],2))+"g/l")

    #Mass flow rate
M_s_out=edbm_dat.M_s_out*N_trip
print("M_s_out is "+str(M_s_out)+"kg/hr")
# M_s_out_r=edbm_dat.M_s_r*N_trip
# print("M_s_out recycling is "+str(M_s_out_r)+"kg/hr")
    #Volumetric flow rate 
Q_s_out=edbm_dat.Q_s_out*N_trip
print("Q_s_out is "+str(Q_s_out)+"l/hr")

"Base channel "
    #Concentration in base channel 
Cb_out=edbm_dat.Ci_b_out[0:6]
Cb_out_g=[Cb_out[0]*MW_Na, Cb_out[1]*MW_Cl, Cb_out[2]*MW_K, Cb_out[3]*MW_Mg, Cb_out[4]*MW_Ca, Cb_out[5]*MW_SO4]
print("base channel Na concentration "+str(round(Cb_out[0],2))+"M and "+str(round(Cb_out_g[0],2))+"g/l")

#     #Mass flow rate
M_b_out=edbm_dat.M_b_out*N_trip
print("M_b_out is "+str(M_b_out)+"kg/hr")
# M_b_out_r=edbm_dat.M_b_r*N_trip
# print("M_b_out recycling is "+str(M_b_out_r)+"kg/hr")
#     #Volumetric flow rate 
Q_b_out=edbm_dat.Q_b_out*N_trip
print("Q_b_out is "+str(Q_b_out)+"l/hr")

#     #Conversion to solid 
M_NaOH_out=Q_b_out*edbm_dat.Ci_b_out[0]*constants.MW_NaOH/1000 #kg/hr 

"Acid channel" 
    #Concentration in acid channel 
Ca_out=edbm_dat.Ci_a_out
Ca_out=edbm_dat.Ci_a_out[0:6]
Ca_out_g=[Ca_out[0]*MW_Na, Ca_out[1]*MW_Cl, Ca_out[2]*MW_K, Ca_out[3]*MW_Mg, Ca_out[4]*MW_Ca, Ca_out[5]*MW_SO4]
print("acid channel Cl concentration "+str(round(Ca_out[1],2))+"M and "+str(round(Ca_out_g[1],2))+"g/l")

    #Mass flow rate 
M_a_out=edbm_dat.M_a_out*N_trip
print("M_a_out is "+str(M_a_out)+"kg/hr")
# M_a_out_r=edbm_dat.M_a_r*N_trip
# print("M_a_out recycling is "+str(M_a_out_r)+"kg/hr")

    #Volumetric flow rate 
Q_a_out=edbm_dat.Q_a_out*N_trip
#Q_a_out=edbm_dat.Q_a_out*N_trip
print("Q_a_out is "+str(Q_a_out)+"l/hr")

    #Conversion to solid 
M_HCl_out=Q_a_out*constants.MW_HCl/1000 #kg/hr

#Calculate required amount of water for the operation mode 
Q_w_in=2*Q_in_edbm

#Calculate mass balance 
bal=(Q_in_edbm*d_in+2*Q_in_edbm)-M_a_out-M_s_out-M_b_out
print("balance edbm is "+str(bal))
error_perc=abs(bal)/(Q_in_edbm*d_in+2*Q_in_edbm)*100
print("balance error percentage is "+str(round(error_perc,2))+"%")

#Energy consumption 
V_ext=edbm_dat.V_ext #xternal 
# #Calculate energy consumption for pumping 
# Ppump=(edbm_dat.Q_s_in_t*N_trip*dp+edbm_dat.Q_a_in_t*N_trip*dp+edbm_dat.Q_b_in_t*N_trip*dp)/1000/3600*1e5/npump #units: W -> l to m3 so /1000; bar to J 1e5N/m2*1J/m ; hr to 3660s  

#Calculate current efficiency 
Cb_in=[0]
CE=(Q_b_out)*(Cb_out[0]-Cb_in[0])*F/(3600*N_trip*I_d*A)*100 #%
print("current efficiency is "+str(CE))

CE2=(Q_b_out)*(edbm_dat.Ci_b_out[0]-edbm_dat.Ci_b_in[0])*F/(3600*N_trip*I_d*A)*100
print("current efficiency 2 is "+str(CE2))

#Total energy consumption 
# E_el_Edbm=V_ext*I_d*A/1000+Ppump/1000
SEC=(V_ext*I_d*A)/(Q_b_out*(edbm_dat.Ci_b_out[0]-edbm_dat.Ci_b_in[0])*constants.MW_NaOH)

# print("Total electrical consumption for EDBM is " + str(E_el_Edbm)+ " KW")
print("sec is "+str(SEC)+"kwh/kg naoh")
