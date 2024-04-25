from density_calc import density_calc
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
    
    def __init__(self,Qf, Mf_med, CNa_in, CCl_in, CK_in, CMg_in, CCa_in, CSO4_in):
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

    def salinity_calc(self):
        # Calculate inflow salinity
        self.salinity_in=sum(self.cons)
        self.xf=self.salinity_in/d*1000
        self.Mf=self.Mf_med/3600 #kg/s
        
        
    def mass_balance_med(self):
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

        
        
    def temperature_calc(self):
        # Calculate temperature-related parameters
        
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
    
    def performance_parameters(self):
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
        self.QCW= self.D2*self.lhv2*1000/(Cp_w*(T_in-T_cw_in))-self.Mf
        
        #Performance ratio (PR) 
        self.PR= self.Mdist/self.Mf
        
        #Condenser thermal load
        self.Qsen=self.Mf*cp_sol*(T_in-25)
        
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

#  #%%
# # #Example usage

# #Feed concentration
# components = ['Na', 'Cl', 'K', 'Mg', 'Ca', 'SO4']
# Cin_med = [10.36, 15.39, 0.36, 0.028, 0.02, 0.07]

# #Feed flow rate 
# Qf_med =1000 #l/hr

# #input conditions
# T=20
# #feed flow density 
# d=density_calc(T, sum(Cin_med)) 
# Mf_med=Qf_med*d/1000 #Mass flow rate (units: kg/hr)

# #assumptions:
# T_in=40 #(oC)
# T_N=45 #Temperature in the last effect (oC)
# N=2 #Number of effects (-)
# Cp_w=4182 # specific heat capacity of water (j/kgC)
# cp_sol=4184 # specific heat capacity of solution (j/kgC)
# T_cw_in=25 #intake cooling water temperature (oC)
# T_cw_out=35 #out cooling water temperature (oC)
# T_s=70 #steam temperature oC
# DT_loss=1 #temperature difference (oC)
# T3=69
# dp=0.1  # pressure drop (units: bar)
# dp_slurry=1 # pressure drop (units: bar)
# npump=0.8 #pump efficiency (units: -)
# Cb_out=200 #Brine Concentration leaving effect n (unit: g/l)
# Xr=5.5 # brine circulation flow rate (units: -)

# #latent heat of motive steam:
# if T_s<=55:
#     lh_s=2370
# elif T_s>55 and T_s<=60:
#     lh_s=2358
# elif (T_s>60) and (T_s<=65):
#     lh_s=2345
# elif (T_s>65) and (T_s<=70):
#     lh_s=2333
# elif (T_s>70) and (T_s<=75):
#     lh_s=2321


# # Create an instance of the MEDCalculator class
# med_dat = MEDCalculator(Qf_med, Mf_med, Cin_med[0], Cin_med[1], Cin_med[2], Cin_med[3], Cin_med[4], Cin_med[5])

# # Call methods to perform calculations
# med_dat.salinity_calc()
# med_dat.temperature_calc()
# med_dat.mass_balance_med()
# med_dat.performance_parameters()
# med_dat.output_concentration()

# # Sum results
# #Concentration of brine stream 
# Cconc_med = [med_dat.CNa_out, med_dat.CCl_out, med_dat.CK_out, med_dat.CMg_out, med_dat.CCa_out, med_dat.CSO4_out]

# #Brine flow rate 
# Qout_med=med_dat.Qb

# #Distillate water flow rate 
# Qprod_med=med_dat.Qdist

# #Calculate circulation flow rate 
# Qr=Xr*Qf_med

# print("Brine flow rate: "+ str(round(Qout_med,2))+"kg/hr")
# print("Total distillate flow rate: "+ str(round(Qprod_med,2))+"kg/hr") 
# # Calculate density for output concentration
# d_b = density_calc(45, sum(Cconc_med))

# print("Sum of output concentrations: " + str(round(sum(Cconc_med),2))+"g/l")
# print("-----------------------------------------")

# # Calculate mass balances
# bal_t=0
# for i in range(len(Cconc_med)):
#     bal_i = (Cin_med[i] * Qf_med /  1000) - (med_dat.Qb / (d_b / 1000) * Cconc_med[i]/1000)
#     bal_t=bal_t+bal_i
#     print("Mass balance for " + components[i] + ": " + str(round(bal_i,2)))
# error_bal=bal_t/(sum(Cin_med)* Qf_med /  1000)*100
# print("Ion balance error percentage is "+str(round(error_bal,2))+"%")
# print("-----------------------------------------")

# # Calculate energy consumption
# E_el_med = ((Qf_med * 3.5 + med_dat.QCW * 3600 * 2 + (Qr + med_dat.Qb) * 3.5 + med_dat.Qdist * 1) / (1000 * npump)) * 1e5 / 3600 / 1000  # kWh
# print("Electrical energy consumption: " + str(round(E_el_med,2)) + " kWh")

# SEC_el = E_el_med / (Qf_med / d)  # kWh/m3 feed
# print("Specific energy consumption (electrical) per m3 feed: " + str(round(SEC_el,2)) + " kWh/m3")

# SEC_el_prod = E_el_med / (med_dat.Qdist / 1000)  # kWh/m3 dist water
# print("Specific energy consumption (electrical) per m3 product (distilled water): " + str(round(SEC_el_prod,2)) + " kWh/m3")
# print("-----------------------------------------")

# E_th_med = med_dat.Q_Tot
# print("Total thermal energy consumption: " + str(round(E_th_med,2)) + " kW")

# SEC_th = E_th_med / (Qf_med / d)  # kWh_th/m3
# print("Specific energy consumption (thermal) per m3 feed: " + str(round(SEC_th,2)) + " kWh_th/m3")
# print("-----------------------------------------")

# #Calculate required cooling water 
# Qcw = med_dat.QCW * 3600 #units: kg/hr
# print("Cooling water flow rate: " + str(round(Qcw,2)) + " kg/hr")
# print("-----------------------------------------")

# # Chemical consumption
# Cchem = 0  # Placeholder for chemical consumption, update as needed

# # Print the results
# print("Total Chemical Consumption: " + str(round(Cchem,2))+"kg/hr")