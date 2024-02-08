import numpy as np
import scaleup
import constants
import math
import density_calc

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
        """Osmotic pressure in concentrate compartment"""
        C1=S/MWs*constants.MW_Na/constants.MW_Na
        z1=1
        C2=S/MWs*constants.MW_cl/constants.MW_cl
        z2=1
        sum_Ci=sum([C1, C2])
        zi_2=[z1**2, z2**2]
        conc=[C1/constants.MW_Na, C2/constants.MW_cl] #mol/l
        mi=[C1*1000/(constants.MW_Na*1000*((1e+6-sum_Ci*1000)/1e+6)), C2*1000/(constants.MW_cl*1000*((1e+6-sum_Ci*1000)/1e+6))]
        mizi_2=[]
        for i in range(2) :
            mizi_2.append(mi[i]*zi_2[i])
            SI=sum(mizi_2)/2
            B=-348.662/(T)+6.72817-0.971307*math.log(T)
            C=40.5016/(T)-0.721404+0.103915*math.log(T)
            D=5321/(T)+233.76-0.9297*(T)+0.001417*(T)**2-0.0000008292*(T)**3
            S=1.17202*(sum(mizi_2)/sum(mi))*0.9982**0.5*(23375.556/(D*(T)))**1.5
            fi=1-S/(3.375*SI)*((1+1.5*SI**0.5)-2*math.log(1+1.5*SI**0.5)-1/(1+1.5*SI**0.5))+B*sum(mi)/2+C*(sum(mi)/2)**2
            sum_conc=sum(conc)
            p_osmo_psig=1.205*fi*(T)*sum(mi)
            p_osmo=p_osmo_psig/14.3 #bar
        return p_osmo
    
    def p_osmo2(S):
        """alternative calculation for Osmotic pressure in concentrate compartment"""
        C1=S/MWs*constants.MW_Na/constants.MW_Na
        C2=S/MWs*constants.MW_cl/constants.MW_cl
        sum_Ci=sum([C1, C2])
        p_osmo=0.0831446261815324*sum_Ci*T
        return p_osmo
    
    def dC(Ts_cp):
        """Change in concentration"""
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
d_in_ed=(density_calc.density_calc(T-273, Sc_i)/1000)
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
    print(f"Ts_cp={ElectrodialysisCalc.Ts_cp(Sd[j - 1])}, Ij={Ij}, F={F}, "
          f"Ls_cp={ElectrodialysisCalc.Ls_cp(Sc[j - 1], Sd[j - 1])}, Sc[j-1]={Sc[j - 1]}, Sd[j-1]={Sd[j - 1]}")
    print("Js is " + str(Js[j]))
    #calculate net water flux 
    Jw[j] = (ElectrodialysisCalc.Tw_cp(Sc[j - 1], Sd[j - 1]) * Ij / F +
             ElectrodialysisCalc.Lw_cp(Sc[j - 1]) * (ElectrodialysisCalc.p_osmo2(Sc[j - 1]) - 
                                                      ElectrodialysisCalc.p_osmo2(Sd[j - 1])))
    
    #Calculate he total concentrate and dilute molar flow rates
    Ns_c[j] = Ns_c[j - 1] + Acp_tot_j * Js[j]
    Ns_d[j] = Ns_d[j - 1] - Acp_tot_j * Js[j]

    Nw_c[j] = Nw_c[j - 1] + Acp_tot_j * Jw[j]
    Nw_d[j] = Nw_d[j - 1] - Acp_tot_j * Jw[j]

    # Update the flow rates of the concentrate and dilute streams
    Q_c[j] = Nw_c[j] * MWw / (rho_w * (1 - Sc[j] / 1000))
    Q_d[j] = Nw_d[j] * MWw / (rho_w * (1 - Sd[j] / 1000))
    print("sd for j " + str(j) + " is " + str(Sd[j]))

print("effluent concentration dilute is " + str(Sd))
print("effluent concentration brine is " + str(Sc))


Cc_na_f=Sc[N-1]/MWs*constants.MW_Na
Cc_cl_f=Sc[N-1]/MWs*constants.MW_cl
Sc_out=[Cc_na_f, Cc_cl_f]

#Calculate the concnetrate stream flow rate 
Mc=(Ns_c[N-1]*MWs/1000+Nw_c[N-1]*MWw/1000) #(kg/hr)
print("Mass flowrate concentrate stream is "+str(Mc))
dc_out=density_calc.density_calc(T-273, Sc[N-1])/1000 #(kg/l)
Qc=Mc/dc_out #concnetrate stream volume flow rate (l/hr)
i=2

for i in range(2,len(Csw)):
    Sc_out.append(Csw[i]*Qed_in_c/Qc)

print("volume flowrate concentrate stream is "+str(Qc))

#Calculations for diluate stream 
Md=(Ns_d[N-1]*MWs/1000+Nw_d[N-1]*MWw/1000) #mass flow rate (kg/hr)
print("Mass flowrate diluate stream is "+str(Md))
Sd_f=Sd[N-1]
Cd_na_f=Sd_f/MWs*constants.MW_Na
Cd_cl_f=Sd_f/MWs*constants.MW_cl
dd_out=density_calc.density_calc(T-273, Sd[N-1])/1000 #density of diluate stream
Qd=Md/dd_out #diluate stream volume flow rate (l/hr)
print("volume flowrate diluate stream is "+str(Qd))

Sd_out=[Cd_na_f, Cd_cl_f]

for i in range(2,len(Csw)):
    Sd_out.append(Csw[i]*Qed_in/Qd)
    

#solid mass balance
bal=Qed_in-Md-Mc
bal=(Qed_in*sum(Csw) -Md*(sum(Sd_out))-Mc*Sc_o)/1000
print("balance is "+str(bal))

#Energy consumption 
#power required 
Ws = 0
for j in range(N):
     Ws += Ij * Acp_tot_j * (Ncp * Vcp + Vel)
print("Ws is "+str(Ws))

#Calculate energy consumption for pumping 
Ppump_ed=(Qed_in_d*1+Qed_in_c*1+Qc*2+Qd*1)/1000/3600*1e5/npump
Eel_t_ed=Ws/1000+Ppump_ed/1000
print("total energy consumption is "+str(Eel_t_ed))
sec_ed=Eel_t_ed/(Qed_in/1000)
print("sec ed is "+str(sec_ed))


