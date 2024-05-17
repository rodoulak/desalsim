import numpy as np
import math
from desalsim.density_calc import density_calc 
from desalsim import constants 
from desalsim import scaleup

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

#%%
#user example 
# Constants
MWs = 58.44  # Molecular weight of NaCl (g/mol)
MWw=18.01528# Molecular weight of water (g/mol)

#Assumptions 
R = 0.002  # Resistance of rinse stream (ohm)
Rp = 0.015  # Resistance of polarization (ohm)
A = 1.1  # Active area of cell-pair (m^2)
F = 96485.3329  # Faraday constant (C/mol)
T=20+273.15 #Inlet/operating temperature (K)
dp=1 #pressure drop (units: bar)
npump=0.8 #pump efficiency (units: -)
rho_w=1000 # water desnity kg/m3

#model parameters 
D=1.61e-9 #Diffusion coefficient (m^2/s)
tcu=0.5
veloc=8.9e-7 #m2/s
h=0.5 #mm
Sh=18

#Input data 
Sc_i = 43.39  # Salinity at concentrate inlet (g/kg)
Sc_o = 200  # Salinity at concentrate outlet (g/kg)
Sd_o=20 # Salinity at doluate outlet (g/kg)
Sd_i=Sc_i # Salinity at diluate inlet (g/kg)
N = 50 # Number of computational cells per cell-pair
Vcp = 8  # Applied voltage (V)
Vel=2.1 #Voltage across the electrodes 
Ij=400 #Am2
Mem_eff=0.64 #Membrane efficiency
Ncp=1 #Number of identical parallel cell-pairs 

#Feed flow rate L/h
d_in_ed=(density_calc(T-273, Sc_i)/1000)
Qed_in=1000/d_in_ed
Qed_in_c=Qed_in/17/Ncp
Qed_in_d=Qed_in*16/17/Ncp

# Effective cell-pair area (m^2)
Acp =scaleup.scaleup(24, 1000, Qed_in) 
Acp_tot=Acp

# Initializations
Sc = np.zeros(N)
Sd = np.zeros(N)
Ns_c = np.zeros(N)
Ns_d = np.zeros(N)
Nw_c = np.zeros(N)
Nw_d = np.zeros(N)
Js = np.zeros(N)
Jw = np.zeros(N)
Mw_in_d_l=np.zeros(N)
Ms_d=np.zeros(N)
Mw_d=np.zeros(N)
M_d=np.zeros(N)
M_c=np.zeros(N)
Q_c=np.zeros(N)
Q_d=np.zeros(N)

# Set initial values salt stream 
Sc[0] = Sc_i
Ns_c[0]=Qed_in_c*d_in_ed*Sc[0]/MWs #mol/hr
Ms_in_c=Qed_in_c*d_in_ed*Sc_i/1000 #kg salr/hr
Mw_in_c=Qed_in_c*d_in_ed-Ms_in_c #kg water/hr
Nw_c[0]=Mw_in_c*1000/MWw #mol/hr
M_c[0]=Mw_in_c+Ms_in_c
Q_c[0]=Qed_in_c


# Set initial values diluate stream 
Csw= [17.17, 25.47, 0.57, 0.04, 0.03, 0.10]
Sd[0] = Sc_i #g/kg 
Ns_d[0]=Qed_in_d*d_in_ed*Sd[0]/(MWs) #mol/s
Ms_in_d=Qed_in_d*d_in_ed*Sd[0]/1000 #kg salr/hr
Mw_in_d=Qed_in_d*d_in_ed-Ms_in_d #kg water/hr
Nw_d[0]=Mw_in_d*1000/MWw #mol/hr
Mw_in_d_l[0]=Mw_in_d
Ms_d[0]=Ms_in_d
Mw_d[0]=Mw_in_d
M_d[0]=Mw_in_d+Ms_in_d
Q_d[0]=Qed_in_d

# Initialize the Acp_tot array
Acp_tot_j=Acp_tot/N
ed_em = ElectrodialysisCalc()

# Iterate over cells
for j in range(1, N):
    #Calculate salinity change 
    concentration_diff = Sc[j - 1] - Sd[j - 1]
    Sc[j] = Sc[j - 1] + (Sc_o - Sc_i) / (N - 1)
    Sd[j] = Sd[j - 1] + (Sd_o - Sd_i) / (N - 1)
    
    #Calculate net salt flux 
    Js[j] = (ElectrodialysisCalc.Ts_cp(Sd[j - 1]) * Ij / F - 
             (ElectrodialysisCalc.Ls_cp(Sc[j - 1], Sd[j - 1])) * concentration_diff)
    #calculate net water flux 
    Jw[j] = (ElectrodialysisCalc.Tw_cp(Sc[j - 1], Sd[j - 1]) * Ij / F +
             ElectrodialysisCalc.Lw_cp(Sc[j - 1]) * (ElectrodialysisCalc.p_osmo(Sc[j - 1], T, MWs) - 
                                                      ElectrodialysisCalc.p_osmo(Sd[j - 1], T, MWs)))
    
    #Calculate he total concentrate and dilute molar flow rates
    Ns_c[j] = Ns_c[j - 1] + Acp_tot_j * Js[j]
    Ns_d[j] = Ns_d[j - 1] - Acp_tot_j * Js[j]

    Nw_c[j] = Nw_c[j - 1] + Acp_tot_j * Jw[j]
    Nw_d[j] = Nw_d[j - 1] - Acp_tot_j * Jw[j]

    # Update the flow rates of the concentrate and dilute streams
    Q_c[j] = Nw_c[j] * MWw / (rho_w * (1 - Sc[j] / 1000))
    Q_d[j] = Nw_d[j] * MWw / (rho_w * (1 - Sd[j] / 1000))





Cc_na_f=Sc[N-1]/MWs*constants.MW_Na
Cc_cl_f=Sc[N-1]/MWs*constants.MW_cl
Sc_out=[Cc_na_f, Cc_cl_f]

#Calculate the concnetrate stream flow rate 
Mc=(Ns_c[N-1]*MWs/1000+Nw_c[N-1]*MWw/1000) #(kg/hr)
print("Mass flowrate concentrate stream is "+str(round(Mc,2))+ " kg/hr")

dc_out=density_calc(T-273, Sc[N-1])/1000 #(kg/l)
Qc=Mc/dc_out #concnetrate stream volume flow rate (l/hr)
i=2

for i in range(2,len(Csw)):
    Sc_out.append(Csw[i]*Qed_in_c/Qc)

print("Volume flowrate concentrate stream is "+str(round(Qc,2))+" l/hr")
print("The total effluent concentration concentrate stream  is " + str(round(Sc[N-1],2))+"g/kg")
print("-----------------------------------------")

#Calculations for diluate stream 
Md=(Ns_d[N-1]*MWs/1000+Nw_d[N-1]*MWw/1000) #mass flow rate (kg/hr)
print("Mass flowrate of diluate stream is "+str(round(Md,2))+" kg/hr")

Sd_f=Sd[N-1]
Cd_na_f=Sd_f/MWs*constants.MW_Na
Cd_cl_f=Sd_f/MWs*constants.MW_cl
dd_out=density_calc(T-273, Sd[N-1])/1000 #density of diluate stream
Qd=Md/dd_out #diluate stream volume flow rate (l/hr)

print("volume flowrate diluate stream is "+str(round(Qd,2))+" l/hr")
print("The total effluent concentration dilute is " + str(round(Sd[N-1],2))+"g/kg")

Sd_out=[Cd_na_f, Cd_cl_f]

for i in range(2,len(Csw)):
    Sd_out.append(Csw[i]*Qed_in/Qd)
print("-----------------------------------------")

#solid mass balance
bal=Qed_in-Md-Mc
bal=(Qed_in*sum(Csw) -Md*(sum(Sd_out))-Mc*Sc_o)/1000
print("Mass balance difference is "+str(round(bal,2)))
error_perc=abs(bal)/(Qed_in*sum(Csw))*100
print("Balance error percentage is "+str(round(error_perc,2))+"%")
print("-----------------------------------------")
#Energy consumption 
#power required 
Ws = 0
for j in range(N):
     Ws += Ij * Acp_tot_j * (Ncp * Vcp + Vel)
print("Power required is "+str(round(Ws/1000,2))+"KW")

#Calculate energy consumption for pumping 
Ppump_ed=(Qed_in_d*1+Qed_in_c*1+Qc*2+Qd*1)/1000/3600*1e5/npump
Eel_t_ed=Ws/1000+Ppump_ed/1000
print("Total energy consumption is "+str(round(Eel_t_ed,2))+"KW")
sec_ed=Eel_t_ed/(Qed_in/1000)
print("Specific energy consumption of Electrodialysis (ED) is "+str(round(sec_ed,2))+"KW/m3 feed")
