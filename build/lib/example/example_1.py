#Import functions 
from nanofiltration_unit_f import OsmoticPressure
from nanofiltration_unit_f import molarity
from nanofiltration_unit_f import NFMass
from nanofiltration_unit_f import NfEnergy

from mfpfr_unit_f import MFPFRCALC
from mfpfr_unit_f import HClAddition
from mfpfr_unit_f import energycons

from ed_unit_f import ElectrodialysisCalc

from edbm_unit_f import EDBMCalc

from economic_f import revenue
from economic_f import econom

from density_calc import density_calc
import constants
import scaleup
import numpy as np
import pandas as pd
import math


#%%Constants

#PRICES
ele_price=0.18 # euro/kwh
steam_price=0.01 # euro/kwh
HCl_price=0.000632058 #euro/l
NaOH_price=0.33 #euro/l
water_price=0.001 # euro/kg
Na2SO4_price=0.05 #euro/kg
NaCl_price=0.066 #euro/kg
MGOH_price=1 #euro/kg
CaOH_price=0.125 #euro/kg


#Molecular weight 
    #molecular weight
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
MW_values = [MW_Na, MW_Cl, MW_K, MW_Mg, MW_Ca, MW_SO4]


MW_MgOH=58.3197
MW_CaOH=74.09
MW_Na2SO4=142.039458
MW_HCl=36.458
MW_NaOH=39.997
MW_NaCl=58.442729
MW_Na2SO4=142.04  
MW_Na2SO4_10h20=322.192258                                                                                                   #Molecular weight of sodium sulfate : g/mol
MW_KCl=74.5513

#constants 
d_Na2SO4=2.66*1000                                                                                                  #Density of sodium sulfate (25oC) : kg/m^3 (source:???) 
Cp_Na2SO4=128.2/MW_Na2SO4*1000   
d_NaCl=2.16*1000                                                                                                    #Density of NaCl (25oC) : kg/m^3 (source:https://pubchem.ncbi.nlm.nih.gov/compound/sodium_chloride#section=Solubility) 
Cp_NaCl=(50.5*MW_NaCl)*1000    
d_w=1*1000                                                                                                      #Density of water (25oC) : kg/m^3 (source:Wikipedia) 
Cp_w=4.1813*1000   
d_KCl=1.98*1000

#input data
#emissions
CO2_el=0.275 # kg/kwh
CO2_st=0#.255 # kg/kwh

#input data
T=20+273 #Operating temperature (units: K)
#Feed concentration
components = ['Na', 'Cl', 'K', 'Mg', 'Ca', 'SO4']
Ci_in = [12.33, 21.67, 0.45, 1.39, 0.45, 3.28]
z_values = [1, -1, 1, 2, 2, -2]
c_values = [Ci / 1000 for Ci in Ci_in]
mg_in = sum(c_values)
#Feed flow density 
d_in = density_calc(T-273, mg_in)  # kg/m3

#Feed flowrate
Qsw = 3000 / 24 * d_in #m3/d


#%%
#Calculations
"""--------Nanofiltration--------"""

Qf_nf = Qsw  # kg/hr
#Constants
R=8.314 #gas constant (units: J / mol·K)

#Asuumptions  
rjr_values = [0.16, 0.29, 0.21, 0.98, 0.95, 0.98] #Ions rejection rates based on membrane characteristics (units: -)
Wrec = 0.7 # Water recovery based on membrane characteristics (units: -)
n=0.8 #pump efficiency (units: -)
dp=2 # pressure drop (units: bar)

# Create molarity objects and calculate meq for each component
molarity_objects = [molarity(MW, z, Ci) for MW, z, Ci in zip(MW_values, z_values, Ci_in)]
meq_values = [m.calculate_meq() for m in molarity_objects]

# Function to create NFMass objects for different components
def create_nfmass_objects(components, C_in, rjr_values, Wrec, Qf):
    return [NFMass(comp, Ci, rjr, Wrec, Qf) for comp, Ci, rjr in zip(components, C_in, rjr_values)]

# Create NFMass objects for different components
nfmass_objects = create_nfmass_objects(components, Ci_in, rjr_values, Wrec, Qf_nf)

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
P_osmo_f = OsmoticPressure(Ci_in, z_values, T).osmotic_pressure_calculation()
P_osmo_p = OsmoticPressure(Cperm, z_values, T).osmotic_pressure_calculation()
P_osmo_c = OsmoticPressure(Cconc, z_values, T).osmotic_pressure_calculation()

d_p=density_calc(T-273, sum(Cperm))

#Calculate Energy consumption 
nf_energy=NfEnergy(P_osmo_c, P_osmo_f, P_osmo_p, dp, d_p, Qperm, Qf_nf, d_in,n)
result=nf_energy.calculate_energy_consumption()
E_el_nf = nf_energy.E_el_nf
for key, value in result.items():
        print(f"{key}: {value}")
        
QHCl_nf=0
Qantsc_nf=0
    
#%%
"""--------MF-PFR--------"""

#constants
kps_MgOH=5.61*0.000000000001 #product solublity of Mg(OH)2
kps_CaOH=5.5*0.000001 #product solublity of Ca(OH)2
d_mgoh_2=2.34 #Mg(OH)2 density (units: kg/l)
d_caoh_2=2.211 #Ca(OH)2 density (units: kg/l)

#assumptions
dp=0.5 # pressure drop (units: bar)
dp_HCl=0.3# pressure drop HCl solution (units: bar)
npump=0.8 #pump efficiency (units: -)

# Define the input parameters
T=20+273.15 #Temperature (units: K)
C_NaOH = [1, 1]  # Concentration of NaOH solution for step 1 and step 2 in MOL/L
conv = [95, 93]  # Conversion rate for step 1 and step 2
Cin_mfpfr = Cconc  # Concentrations of [Na, Cl, K, Mg, Ca, SO4]
Qin_mfpfr = Qconc # Flow rate in l/hr

# Calculate the density of the input
d_in = density_calc(25, sum(Cin_mfpfr)) / 1000

# Create an instance of the inputpar class with the defined parameters
mfpfr_dat = MFPFRCALC(Qin_mfpfr, Cin_mfpfr, *C_NaOH, *conv)

# Call the calc_step1 and calc_step2 methods to calculate the necessary values
mfpfr_dat.calc_step1(kps_MgOH, d_mgoh_2)
mfpfr_dat.calc_step2(d_mgoh_2, d_caoh_2 )
ph_2=mfpfr_dat.ph_2

# Calculate the total outlet concentration and the concentration of sulfate ions
Cour_mfpfr = sum([mfpfr_dat.CNa_out_2, mfpfr_dat.CCa_out_2, mfpfr_dat.CCl_out_2, mfpfr_dat.CK_out_2, mfpfr_dat.CMg_out_2, mfpfr_dat.CSO4_out_2]) # mol/l
CSO4_out_2 = mfpfr_dat.CSO4_out_2 # mol/l

# Calculate the concentration of sodium chloride ions
CNa_out_2 = mfpfr_dat.CNa_out_2
Cnacl_out = CNa_out_2 - 2 * CSO4_out_2

# Get the molar masses of the compounds
M_MgOH2_1 = mfpfr_dat.M_MgOH2_1
M_CaOH2 = mfpfr_dat.M_CaOH2_2
M_MgOH2 = mfpfr_dat.M_MgOH2_2

# Create a list of the outlet concentrations in mol/l
Cout_all_m = [mfpfr_dat.CNa_out_2, mfpfr_dat.CCl_out_2, mfpfr_dat.CK_out_2, mfpfr_dat.CMg_out_2, mfpfr_dat.CCa_out_2, mfpfr_dat.CSO4_out_2]

# Calculate the outlet concentrations in g/l
MW = [MW_Na, MW_Ca, MW_Cl, MW_K, MW_Mg, MW_SO4]
Cout_mfpfr_g = [Cout_all_m[i] * MW[i] for i in range(len(Cout_all_m))] # g/l

#Outlet flow rate 
Qout_2 = mfpfr_dat.Qout_2

# Calculate the chemical consumption of NaOH
QNAOH = mfpfr_dat.QNaOH_1 + mfpfr_dat.QNaOH_2_add + mfpfr_dat.QNaOH_2_st # convert to kg

#hcl to decrease ph to 7
HCl_conc=1 #l of HCl 1M
unit = HClAddition(Qout_2, Cout_all_m, MW_Cl, ph_2, HCl_conc)

# Call the calculate_HCl_addition method
QHCl, Cout_mfpfr_g = unit.calculate_HCl_addition(Cout_mfpfr_g)

# Print the volume of HCl added and the outlet concentration of chloride ions
print(f"HCl flow rate is {round(QHCl,2)} l/hr")
print(f"NaOH flow rate is {round(QNAOH,2)} l/hr")
print("-----------------------------------------")

#Calculate final outlet flow rate
Qout_f=Qout_2+QHCl #l/h
d_out_s=density_calc(25,round(sum(Cout_mfpfr_g),2))/1000 #kg/m3
M_mfpfr_out=Qout_f*d_out_s #kg/h

print("Mg(OH)2 mass flow rate is "+str(round(M_MgOH2_1,2))+"kg/hr")
print("Ca(OH)2 mass flow rate is "+str(round(M_CaOH2,2))+"kg/hr")
print("Total effluent flow rate is "+str(round(M_mfpfr_out,2))+"kg/hr")
print("Total effluent flow rate is "+str(round(Qout_f,2))+"kg/hr")
print("-----------------------------------------")

#electicity consumption
    # Create an instance of the inputpar class with the defined parameters
Epump_1, Epump_2=energycons.energycalc(mfpfr_dat.Qout_2, QNAOH, Qin_mfpfr, mfpfr_dat.QNaOH_1, mfpfr_dat.QNaOH_2_add, mfpfr_dat.QNaOH_2_st, dp, npump)

    #Electricity consumption for pumping , KWh
E_el_mfpf=(Epump_1+Epump_2+(QHCl*dp_HCl)*1e5/3600/(1000*npump))/1000
print("Total electricity energy consumption is "+str(round(E_el_mfpf,2))+ " KW")

    #Electricity consumption for filtration unit 
E_fil=scaleup.scaleup(0.5, 0.3*1000, M_mfpfr_out)

    #Total electricity consumption, KWh
E_el_mfpf=E_el_mfpf+E_fil

    #Specific energy consumption per kg of Mg(OH)2, KWh/kg of Mg(OH)2
SEC_el_prod=(E_el_mfpf)/(M_MgOH2)
print("Specific energy consumption per product is "+str(round(SEC_el_prod,2))+" KWh/kg product ")

    #Specific energy consumption per feed, KWh/m3 of feed
SEC_el_feed=(E_el_mfpf)/(Qin_mfpfr/1000)
print("Specific energy consumption per brine intake is "+str(round(SEC_el_feed,2))+" KWh/m3 of feed ")

#%%
"""--------Electrodialysis--------"""
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
Sc_i = sum(Cperm)  # Salinity at concentrate inlet (g/kg)
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
Qed_in=Qperm
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
    print("sd for j " + str(j) + " is " + str(Sd[j]))


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

#%%
"""--------Electrodialysis with Bipolar Membranes--------"""
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
T_in=25

#Membrane characteristics

Cm_bp_H= 0.0000001 #mol/l 
Cm_bp_OH= 0.0000001 #mol/l 


#input data
r=1-0.4/(5-0.4) #recycling rate 
rs=1-1/(5-1) #recycling rate 
ph_s=4.71 #pH salt channel (units: -)
ph_b=7#pH base channel (units: -)
ph_a=7#pH acid channel (units: -)
d_a=1 #density acid channel (units: kg/l)
d_s=1 #density salt channel (units: kg/l)
d_b=1 #density base channel (units: kg/l)

#Feed flow rate L/h
Q_in_edbm=M_mfpfr_out+Mc

#Feed concentration g/L
Cin_edbm=sum(Cout_mfpfr_g)*M_mfpfr_out/Q_in_edbm+Sc_o*Mc/Q_in_edbm
C_in_mix=[]
for i in range(len(Cconc)):
    C_in_mix.append(Cout_mfpfr_g[i]*M_mfpfr_out/Q_in_edbm+Sc_out[i]*Mc/Q_in_edbm)
    
Cin_edbm=C_in_mix #[Na, Cl, K, Mg, Ca, SO4 ]
Cin_edbm.extend([0, 10**(-ph_s), 3.01551E-11])
d_in=density_calc(T_in,sum(Cin_edbm))/1000

C_b_in=[0,0,0,0,0,0,0,10**(-ph_b), 10**(-(14-ph_b))]
C_a_in=[0,0,0,0,0,0,0,10**(-ph_a), 10**(-(14-ph_a))]
#Calculate water quantity in inflow 
Mw_in=Q_in_edbm/d_in 

#Set number of triplets 
N_trip=50*47

#Set membrane area based on the feed flow rate 
A=0.4

#Initialize concentration of Na in salt channel 
Cna_s=[]

#Create an instance of the EDBMCalc class with the defined parameters
edbm_dat=EDBMCalc(Q_in_edbm, A, I_d, N_trip, Cin_edbm, C_b_in, C_a_in, T_in )

# Call the necessary methods to calculate values
edbm_dat.flowrate()
edbm_dat.in_mass_flow_rates(ph_s)
edbm_dat.acid_channel()
edbm_dat.base_channel()
edbm_dat.salt_channel(Cm_bp_H, Cm_bp_OH)
Cna_s.append(edbm_dat.Ci_s_out[0])

# # Sum results

"Salt channel "
    #Concentration in salt channel 
Cbrine_out_t=sum(edbm_dat.Ci_s_out)
Cbrine_out=edbm_dat.Ci_s_out #mol/l
Cbrine_out_g=[Cbrine_out[0]*MW_Na, Cbrine_out[1]*MW_Cl, Cbrine_out[2]*MW_K, Cbrine_out[3]*MW_Mg, Cbrine_out[4]*MW_Ca, Cbrine_out[5]*MW_SO4] #g/l

     #Mass flow rate
M_s_out=edbm_dat.M_s_out_t*N_trip
print("Salt channel: Mass flow rate out is "+str(round(M_s_out,2))+"kg/hr")
    #Volumetric flow rate 
Q_s_out=edbm_dat.Q1_s_out*N_trip
print("Salt channel: Volumetric flow rate out is "+str(round(Q_s_out,2))+"l/hr")
print("Na concentration:"+str(round(Cbrine_out[0],2))+"M and "+str(round(Cbrine_out_g[0],2))+"g/l")
print("Cl concentration:"+str(round(Cbrine_out[1],2))+"M and "+str(round(Cbrine_out_g[1],2))+"g/l")
print("-----------------------------------------")

"Base channel "
    #Concentration in base channel 
Cb_out=edbm_dat.Ci_b_out[0:6]
Cb_out_g=[Cb_out[0]*MW_Na, Cb_out[1]*MW_Cl, Cb_out[2]*MW_K, Cb_out[3]*MW_Mg, Cb_out[4]*MW_Ca, Cb_out[5]*MW_SO4]

     #Mass flow rate
M_b_out=edbm_dat.M_b_out_t*N_trip
print("Base channel: Mass flow rate out is "+str(round(M_b_out,2))+"kg/hr")

     #Volumetric flow rate 
Q_b_out=edbm_dat.Q1_b_out*N_trip
print("Base channel: Volumetric flow rate out is "+str(round(Q_b_out,2))+"l/hr")

print("Na concentration "+str(round(Cb_out[0],2))+"M and "+str(round(Cb_out_g[0],2))+"g/l")
    #Conversion to solid 
M_NaOH_out=Q_b_out*edbm_dat.Ci_b_out[0]*constants.MW_NaOH/1000 #kg/hr 
print("-----------------------------------------")

"Acid channel" 
    #Concentration in acid channel 
Ca_out=edbm_dat.Ci_a_out
Ca_out=edbm_dat.Ci_a_out[0:6]
Ca_out_g=[Ca_out[0]*MW_Na, Ca_out[1]*MW_Cl, Ca_out[2]*MW_K, Ca_out[3]*MW_Mg, Ca_out[4]*MW_Ca, Ca_out[5]*MW_SO4]

    #Mass flow rate 
M_a_out=edbm_dat.M_a_out_t*N_trip
print("Acid channel: Mass flow rate out is "+str(round(M_a_out,2))+"kg/hr")

    #Volumetric flow rate 
Q_a_out=edbm_dat.Q1_a_out*N_trip
print("Acid channel: Volumetric flow rate out is "+str(round(Q_a_out,2))+"l/hr")
print("Cl concentration "+str(round(Ca_out[1],2))+"M and "+str(round(Ca_out_g[1],2))+"g/l")
print("-----------------------------------------")

    #Conversion to solid 
M_HCl_out=Q_a_out*constants.MW_HCl/1000 #kg/hr

#Calculate required amount of water for the operation mode 
Q_w_in=2*Q_in_edbm

#Calculate mass balance 
bal=(Q_in_edbm*d_in+2*Q_in_edbm)-M_a_out-M_s_out-M_b_out
print("Mass balance difference is "+str(round(bal,2)))
error_perc=abs(bal)/(Q_in_edbm*d_in+2*Q_in_edbm)*100
print("Balance error percentage is "+str(round(error_perc,2))+"%")

#Energy consumption 
V_ext=edbm_dat.V_ext #xternal 
#Calculate energy consumption for pumping 
Ppump=(edbm_dat.Q1_s_in*N_trip*dp+edbm_dat.Q1_a_in*N_trip*dp+edbm_dat.Q1_b_in*N_trip*dp)/1000/3600*1e5/npump #units: W -> l to m3 so /1000; bar to J 1e5N/m2*1J/m ; hr to 3660s  

#Calculate current efficiency 
Cb_in=[0]
CE=(Q_b_out)*(Cb_out[0]-Cb_in[0])*F/(3600*N_trip*I_d*A)*100 #%
print("Current efficiency is "+str(round(CE,2))+"%")
print("-----------------------------------------")
#Total energy consumption 
E_el_Edbm=V_ext*I_d*A/1000+Ppump/1000
SEC=(V_ext*I_d*A)/(Q_b_out*(edbm_dat.Ci_b_out[0]-edbm_dat.Ci_b_in[0])*constants.MW_NaOH)

print("Total electrical consumption for EDBM is " + str(round(E_el_Edbm,2))+ " KW")
print("Specific energy consumption is "+str(round(SEC,2))+"kwh/kg NaOH")

#%%
#Concentration of saline solution 
Q_out_b=Qd+Q_s_out
Cout_b=[]
for i in range(len(Cbrine_out_g)):
    Cout_b.append(Cbrine_out_g[i]*Q_s_out/Q_out_b+Sd_out[i]*Qd/Q_out_b)
Cout_b_t=sum(Cout_b)

#%%Indicators:
class indic:
    
    def __init__(self, tech, Qout, Qin, Qprod_1, prod_1_name, Qprod_2, prod_2_name, E_el, E_th, Cin, Cout, chem1,chem2):
        self.tech=tech
        self.Qout=Qout
        self.Qin=Qin
        self.E_el=E_el
        self.E_th=E_th
        self.Cin=Cin
        self.Cout=Cout
        self.Qprod_1=Qprod_1
        self.Qprod_2=Qprod_2
        self.prod_1_name=prod_1_name
        self.prod_2_name=prod_2_name
        self.chem=chem1+chem2
        
        
    def techn_indi(self):        
        if self.prod_1_name=="water": 
            self.Water_rr=self.Qprod_1/self.Qin*100 # %
            self.Qwater=self.Qprod_1
        elif self.prod_2_name=="water":
            self.Water_rr=self.Qprod_2/self.Qin*100 # %
            self.Qwater=self.Qprod_2
        else:
            self.Qwater=0
            self.Water_rr=0
            print("no water production by " + self.tech)
        
#%%Sumsummarize results  
tec1=indic("NF", Qconc, Qf_nf, Qperm, "none", 0, "none", E_el_nf, 0, Ci_in, Cconc, QHCl_nf, Qantsc_nf)   
tec1.techn_indi()

tec2=indic("MF-PFR", Qout_2, Qconc, M_MgOH2_1,"Mg(OH)2", M_CaOH2, "Ca(OH)2", E_el_mfpf, 0, Cconc, Cout_mfpfr_g, QNAOH, QHCl)
tec2.techn_indi()

tec3=indic("ED", Qc, Qed_in, 0, "none", 0, "none",  Eel_t_ed, 0, Sc_i, Sd_out, 0,0)
tec3.techn_indi()

tec4=indic("EDBM", Q_s_out, Q_in_edbm, Q_b_out,"NaOH", Q_a_out, "HCl", E_el_Edbm, 0, Cin_edbm, Cbrine_out_g,0,0 )
tec4.techn_indi()

#
##list results 
tec_names=[tec1.tech,tec2.tech, tec3.tech,  tec4.tech]
E_el_all=[tec1.E_el,tec2.E_el, tec3.E_el, tec4.E_el]
E_th_all=[tec1.E_th,tec2.E_th, tec3.E_th,  tec4.E_th]
Cout_all=[tec1.Cout,tec2.Cout, tec3.Cout,  tec4.Cout]
Qout_all=[tec1.Qout, tec2.Qout, tec3.Qout,  tec4.Qout]
Qchem_all=[tec1.chem,tec2.chem,tec3.chem, tec4.chem]

#REQUIRED AMOUNT OF CHEMICALS 
Qnaoh_ex=Q_b_out-QNAOH
if Qnaoh_ex<0:
    Qnaoh_s=0
    Qnaoh_need=-Qnaoh_ex
else:
    Qnaoh_s=Qnaoh_ex
    Qnaoh_need=0
print("remaining amount of naoh is "+str(Qnaoh_ex) + "kg/hr")
Qhcl_ex=Q_a_out-QHCl
print("remaining amount of hcl is "+str(Qhcl_ex) + "kg/hr")
if Qhcl_ex<0:
    Qhcl_s=0
    Qhcl_need=-Qhcl_ex
else:
    Qhcl_s=Qhcl_ex
    Qhcl_need=0
#naoh need solid
Mnaoh_need=Qnaoh_need*constants.MW_NaOH/1000
Mhcl_need=Qhcl_need*constants.MW_HCl/1000

#
sum_el_en=sum(E_el_all)
sum_th_en=sum(E_th_all)
#%%
#Technical 
#T-1) quantity of water production 
Qw_tot=tec1.Qwater+tec2.Qwater+ tec3.Qwater+ tec4.Qwater
abs_Qw=Qw_tot-Q_w_in
print("absolute water production is "+str(abs_Qw))
#T-1b) Water recovery  
Qsw=Qsw
rec=Qw_tot/Qsw*100 #%
#T-2)energy performance: specific electrical energy consumption #kwh/kg of desalinated water 
SEC=sum(E_el_all)/Qw_tot
print("specific electrical consumption is " +str(SEC)+ " Kwh/kg of desalinated water")
Tot_el=sum(E_el_all) #kWh 
Tot_th=sum(E_th_all) #kWh

#T-3)Brine production -> kg of brine/ kg of seawater  
Qbrine_f=0
for i in range(len(tec_names)):
    if tec_names[i]=="EDBM" or tec_names[i]=="ED":
        if (sum(Cout_all[i])>sum(Csw)).any():
            Qbrine_f=Qbrine_f+Qout_all[i]
            print("for "+str(tec_names[i])+" brine concentration is "+str(sum(Cout_all[i])))
            print("for "+str(tec_names[i])+" brine production is "+str(Qout_all[i]))
        else:
            Qbrine_f=Qbrine_f+0
Q_br_prod=Qbrine_f/Qsw
print("brine production is "+str(Q_br_prod)+" kg of brine/kg of seawater")
#------------------------------------------------------------------------------
#%%
#Environmental 
#Env-1)Climate change impact ->Carbon dioxide emission kg co2 
emis=(sum(E_el_all)*CO2_el)+(sum(E_th_all)*CO2_st)
emis_t=((sum(E_el_all)*CO2_el)+(sum(E_th_all)*CO2_st))*constants.hr
print("Total specific Carbon dioxide emission "+str(emis_t)+" kg co2")
#Env-2) Resource utilization -> chemical consumption in kg of chemical/kg of seawater 
Chem_cons=(sum(Qchem_all))/Qsw #check units and convert ml/l to kg with densities 
print("Total chemical consumption is "+str(Chem_cons)+" kg of chemical/kg of seawater ")
#Env-3) Water footprint/ water use 
Qw_use=Q_w_in+Qnaoh_need+Qhcl_need

#%%
#Economic 
#Input data 
hr=24*300 #hours/year
lf= 20 #years
r=0.06 # interest rate
inf=0.02 # inflation rate 

#prices 
el_pr=0.253 #euro/kwh
s_pr=0 #euro/kwh
co2_em_c=23 #euro/ton co2 
nacl_pr=66/1000 # euro/kg 
hcl_pr=1040.0/180 #euro/l 1M solution 
w_pr=0.001 #euro/kg
na2so4_pr=140*0.83/1000 #euro/kg
mgoh2_pr=1000/1000 #euro/kg
caoh2_pr=125/1000#euro/kg
naoh_pr=6840/950#euro/L 1M NaOH solution  
antisc_pr=0.002 #euro/ml
cw_pr=0.118/1000#euro/kg
m_salt_pr=5/1000#euro/kg

#density
d_naoh=1.04 #kg/l for 1M solution 
d_hcl=1.0585#kg/l for 1M solution 

#Assumptions 
main_c_percent=0.03  # % of the fixed-capital investment
oper_sup_c_percent=0.05  # % of maintenance 
oper_lab_c_percent=0.15    # % of annual OPEX
super_c_percent=0.15 # % of operating labor 
lab_c_percent=0.15 # % of operating labor 
pat_c_percent=0.03 # % of annual OPEX
fix_char_percent=0.05 # % of annual OPEX
over_c_percent=0.05 # %of annual OPEX
norm_factor =(1 -oper_lab_c_percent-oper_lab_c_percent*super_c_percent- oper_lab_c_percent*lab_c_percent-pat_c_percent-fix_char_percent-over_c_percent)
economic_assumptions=[main_c_percent,oper_sup_c_percent, oper_lab_c_percent, super_c_percent, lab_c_percent, pat_c_percent, 
                      fix_char_percent, over_c_percent, norm_factor ]

CAPEX=0
OPEX=0
CO2_c=0
OPEX_with_ext=0
Mf_basic_sc=[]
capex_list=[]
opex_list=[]
eq_c=[constants.eq_c_nf, constants.eq_c_ed, constants.eq_c_mfpfr, constants.eq_c_edbm]   
el_conc=E_el_all
s_conc=E_th_all
chem1_conc=[Qantsc_nf,0,Mnaoh_need,0]
chem1_pr=[0.002,0,constants.naoh_pr_s,0]
chem2_conc=[0,0,Mhcl_need,0]
chem2_pr=[0.002,0,constants.hcl_pr_s,0]
wat_conc=[0,0, Qnaoh_need, Q_w_in]
cw_conc=[0,0, 0, 0]
Mf_basic_sc=[constants.Mf_basic_sc[0],constants.Mf_basic_sc[8], constants.Mf_basic_sc[3], constants.Mf_basic_sc[6]]

Mf_sce=[Qsw,   Qed_in, Qin_mfpfr, Q_in_edbm]
for i in range(len(eq_c)):
    if Mf_basic_sc[i]!=Mf_sce[i]:
        eq_c[i]=scaleup.scaleup_eq(eq_c[i],Mf_basic_sc[i],Mf_sce[i],tec_names[i])
        
#Cost calculation 
for i in range(len(eq_c)):
    total_econom=econom(eq_c[i], el_conc[i], s_conc[i], chem1_conc[i], chem1_pr[i],chem2_conc[i], chem2_pr[i], cw_conc[i], wat_conc[i])
    total_econom.capex_calc()
    total_econom.opex_calc(hr, el_pr, s_pr, cw_pr, w_pr, economic_assumptions)
    CAPEX=total_econom.t_capital_inv
    OPEX=total_econom.opex
    opex_list.append(total_econom.opex)
    capex_list.append(total_econom.t_capital_inv)    
print("Tottal operating cost (OPEX) is "+str(OPEX)+ " Euro/year") 
print("Total investment cost (CAPEX) of system is " + str(round(CAPEX))+ " Euro")    
print("-----------------------------------------")

#amortisation factor 
a=(constants.r*(1+constants.r)**constants.lf)/((1+constants.r)**constants.lf-1)

"""Revenue calculation"""
#Input data 
reve_t=0
reve_list=[]
prd=[10,2,3,0,0,0]    
prd_name= ["Water", "NaCl", "Mg(OH)2", "Na2SO4", "NaOH", "HCl"]   

#Revenue calculation
for i in range(len(prd)):
    rev_calc=revenue(prd[i], prd_name[i])    
    rev_calc.rev(hr, w_pr, nacl_pr, mgoh2_pr,na2so4_pr, naoh_pr, hcl_pr)
    print("Revenues from selling product " + prd_name[i]+" are " + str(round(rev_calc.rev_prd,2))+" Euro/year")
    reve_t = reve_t+rev_calc.rev_prd
    reve_list.append(rev_calc.rev_prd)
    
#%%Present results
#Summary tables 
    #Table C: Ion concentration 
sum_table_C={'F1: Seawater': Ci_in,
           'F2: NF permeate':Cperm, 
           'F3: NF concentrate': Cconc,
           'F4: ED saline solution': Sd_out, 
           'F5: ED to EDBM': Sc_out,
           'F6: MF-PFR out': Cout_mfpfr_g, 
           'F7: MgOH2': [0,0,0,0,0,0],
           'F8: CaOH2': [0,0,0,0,0,0],
           'F9: EDBM: ACID': Ca_out_g[0:6],
           'F10: EDBM: BASE ':Cb_out_g[0:6],
           'F11: EDBM: SALT':Cbrine_out_g[0:6],
           }
sum_table_C=pd.DataFrame(sum_table_C)
print(sum_table_C)

    #Table D: Density
sum_table_d={'F1: Seawater': density_calc(25,sum(Ci_in))/1000,
           'F2: NF permeate':density_calc(25,sum(Cperm))/1000, 
           'F3: NF concentrate': density_calc(25,sum(Cconc))/1000,
           'F4: ED saline solution': density_calc(25,sum(Sd_out))/1000, 
           'F5: ED to EDBM': density_calc(25,sum(Sc_out))/1000,
           'F6: MF-PFR out':density_calc(25,sum( Cout_mfpfr_g))/1000, 
           'F7: MgOH2': 1,
           'F8: CaOH2': 1,
           'F9: EDBM: ACID': density_calc(25,sum(Ca_out_g[0:6]))/1000,
           'F10: EDBM: BASE ': density_calc(25,sum(Cb_out_g[0:6]))/1000,
           'F11: EDBM: SALT': density_calc(25,sum(Cbrine_out_g[0:6]))/1000,
           }
sum_table_d=pd.DataFrame(sum_table_d, index=['density'])
print(sum_table_d)

        #Table F: Flow rates  
sum_table_f={'F1: Seawater': round(Qf_nf,2),
           'F2: NF permeate':round(Qperm,2), 
           'F3: NF concentrate': round(Qconc,2),
           'F4: ED diluate': round(  Md,2), 
           'F5: ED to EDBM': round(Mc,2),
           'F6: MF-PFR out': round(M_mfpfr_out,2), 
           'F7: MgOH2': round(M_MgOH2_1,2),
           'F8: CaOH2': round(M_CaOH2,2),
           'F9: EDBM: ACID':round(Q_a_out,2),
           'F10: EDBM: BASE ':round(Q_b_out,2),
           'F11: EDBM: SALT':round(Q_s_out,2),
           'F12: EDBM: water in':round(Q_w_in,2),
           'F13: MF-PFR NaOH': round(QNAOH,2),
           'F14: MF-PFR HCl': round(QHCl,2),
           'F15: NF antiscalant':round(Qantsc_nf,2)
           }
sum_table_f=pd.DataFrame(sum_table_f, index=['flow rate'])
print(sum_table_f)

#sankey diagram 
import plotly.graph_objects as go
from plotly.offline import plot

#Figure 1: Sankey diagram for Example: Mass flow rates
fig = go.Figure(data=[go.Sankey(
    node = dict(
      pad = 15,
      thickness = 20,
      line = dict(color = "black", width = 0.5),
      label = ["Seawater", "Water", "NaCl",  "NF", "ED", "MF-PFR", "EDBM", "HCl", "NaOH", "Mg(OH)\u2082", "Ca(OH)\u2082", "Utilities", "Saline solution"],
      color = "blue"
    ),
    link = dict(
      source = [0,3,3,4,4,5,5,5,6,6,6,11,11], 
      target = [3,4,5,6,12,6,9,10,7,8,12,5,6],
      value = [Qsw, Qperm, Qconc,Mc, Md,M_mfpfr_out, M_MgOH2_1,M_CaOH2,Q_a_out, Q_b_out,Q_s_out, QNAOH+QHCl,Q_w_in ]
  ))])
color_for_nodes = ["lightsteelblue","darkcyan","maroon", "midnightblue", "midnightblue", "midnightblue", "maroon"]
fig.update_traces(node_color = color_for_nodes)
fig.update_layout(title_text="Sankey diagram for Example: Mass flow rates ", font_size=15)
fig.show()
plot(fig)

#%% Results to excel 
dfel=pd.DataFrame(E_el_all ,tec_names)
dfeth=pd.DataFrame(E_th_all, tec_names)
dfind=(Qw_tot, abs_Qw, rec, Tot_el, Tot_th, emis_t, Q_w_in,  OPEX, CAPEX)
dfprod=(Qw_tot, M_MgOH2_1, M_CaOH2,  M_b_out, M_a_out)

dfprodn=("Water (kg/hr)",  "Mg(OH)2(kg/hr)", "Ca(OH)2(kg/hr)",  "1M NaOH (kg/hr)","1M HCl(kg/hr)")
ind=np.array(["Water production", "absolute water production", "water recovery","Total electrical consumption", "Total thermal energy consumption", "Carbon dioxide emission kg co2/year",
              "Water footprint","OPEX", "CAPEX"])
units=pd.DataFrame(["kg/hr", "kg/hr", "%","KWh", "KWh"," kg co2/year","kg/hr","€/year", "€"])
dfind_t=pd.DataFrame(data=dfind, index=ind)
dfprodn=pd.DataFrame(data=dfprod, index=dfprodn)
with pd.ExcelWriter('results_example.xlsx') as writer:
    sum_table_f.to_excel(writer,sheet_name="example")
    sum_table_C.to_excel(writer,sheet_name="example", startcol=0, startrow=6, header=True)
    sum_table_d.to_excel(writer,sheet_name="example", startcol=0, startrow=2, header=True)
    dfel.to_excel(writer,sheet_name="example", startcol=0, startrow=14, header=True)
    dfeth.to_excel(writer,sheet_name="example", startcol=2, startrow=14, header=True)
    dfprodn.to_excel(writer,sheet_name="example", startcol=0, startrow=23, header=True)
    #dfprod.to_excel(writer,sheet_name="water_mining_scenario", startcol=1, startrow=24, header=True)
    dfind_t.to_excel(writer,sheet_name="indicators")
    units.to_excel(writer,sheet_name="indicators",startcol=2, startrow=0, header=True)