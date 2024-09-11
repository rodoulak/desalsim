from desalsim.density_calc import density_calc 

#%%calculations 
class MEDCalculator:
    """
    A class used to represent Mass Balance for Multi-Effect Distillation 

    ...

    Attributes
    ----------
    Qf : float
        Flow rate (m^3/s).
    CNa_in, CCl_in, CK_in, CMg_in, CCa_in, CSO4_in : float
        Initial concentrations of various ions (g/l).

    Methods
    -------
    salinity_calc():
        Calculates the inflow salinity.
        
    mass_balance_med():
        Performs mass balance calculations.
        
    temperature_calc(T_s, T_N, T_cw_in, DT_loss):
        Calculates temperature-related parameters.
        
    performance_parameters():
        Calculates performance parameters.
        
    output_concentration():
        Calculates the concentrations in the output.
 
    """
    
    def __init__(self,Qf, Mf_med, CNa_in, CCl_in, CK_in, CMg_in, CCa_in, CSO4_in, T_in):
        #Initialize class attributes
        self.Qf = Qf
        self.Mf_med=Mf_med
        self.CNa_in = CNa_in
        self.CCl_in = CCl_in
        self.CK_in = CK_in
        self.CMg_in = CMg_in
        self.CCa_in = CCa_in
        self.CSO4_in = CSO4_in
        self.cons=[self.CNa_in,self.CCl_in,self.CK_in,self.CMg_in,self.CCa_in,self.CSO4_in]
        self.T_in=T_in
        self.d=density_calc(self.T_in, sum(self.cons))


    def salinity_calc(self):
        # Calculate inflow salinity
        self.salinity_in=sum(self.cons)
        self.xf=self.salinity_in/self.d*1000
        self.Mf=self.Mf_med/3600 #kg/s
        
        
    def mass_balance_med(self, Cb_out):
        # Perform mass balance calculations
        #Brine Concentration leaving effect n 
        self.xn=Cb_out
        
        #Calculate concentration factor 
        self.conc_f=self.xn/self.xf
        
        #Calculate brine flow rate of leaving effect n
        self.Bn=self.xf*self.Mf/self.xn #kg/s
        self.Qb=self.Bn*3600 #kg/hr
        
        #Calculate total distillate flow rate
        self.Mdist=self.Mf-self.Bn #kg/s
        self.Qdist=self.Mdist*3600 #kg/hr
        
        #Calculate distillate flow produces for each effect
        self.D1=self.Mdist/(1+self.lhv1/self.lhv2) #kg/s
        self.D2=self.D1*self.lhv1/self.lhv2 #kg/s

        
        
    def temperature_calc(self, DT_loss, T_N, T_s):
        # Calculate temperature-related parameters
        T3=69
        #Transfer coefficients of the effects,
        self.U1= 1.9695+0.012057*T3-0.000085989*T3**2+0.00000025651*T3**3
        self.U2=0.95*self.U1
        
        #Temperature drop in each effect 
        self.DT_t=T_s-T_N #oC
        self.DT1=self.DT_t/(self.U1*(1/self.U1+1/self.U2)) #oC
        self.DT2=self.DT1*self.U1/self.U2 #oC
        
        #Temperature profile of each effect 
        self.T1=T_s-self.DT1 #oC
        self.T2=self.T1-self.DT2 #oC
        
        #Latent heat of first effect 
        self.lhv1=2499.5698-0.9156*(self.T1-DT_loss)-0.048343*(self.T1-DT_loss)**2
        
        #Latent heat of second effect
        self.lhv2=2626.1
    
    def performance_parameters(self, lh_s,  T_cw_in):
        # Calculate performance parameters
        #Steam mass flow rate
        self.Ms=self.D1*self.lhv1/lh_s #kg/s
        
        #Gain Output Ratio
        self.GOR=self.Mdist/self.Ms
        
        #Condenser thermal load (Qc) kj/s ot kw
        self.Qc=self.D2*self.lhv2 
        
        #thermal load in the first effect (Q1) kw
        self.Q1= self.Ms*lh_s
        
        #Cooling water flow rate kg/s
        Cp_w=4182 # specific heat capacity of water (j/kgC)
        self.QCW= self.D2*self.lhv2*1000/(Cp_w*(self.T_in-T_cw_in))-self.Mf
        
        #Performance ratio (PR) 
        self.PR= self.Mdist/self.Mf
        
        #Condenser thermal load
        cp_sol=4184 # specific heat capacity of solution (j/kgC)
        self.Qsen=self.Mf*cp_sol*(self.T_in-25)
        
        #Total thermal load
        self.Q_Tot=self.Q1+self.Qsen/1000

    def output_concentration(self):
        # Calculate concentrations in the output, g/l
        self.CNa_out=self.CNa_in*self.Mf/self.Bn 
        self.CCl_out=self.CCl_in*self.Mf/self.Bn 
        self.CK_out=self.CK_in*self.Mf/self.Bn 
        self.CMg_out=self.CMg_in*self.Mf/self.Bn 
        self.CCa_out=self.CCa_in*self.Mf/self.Bn  
        self.CSO4_out=self.CSO4_in*self.Mf/self.Bn  
