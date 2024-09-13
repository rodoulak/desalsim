#Scenario 2
##treatment train: NF-> MED-> THERMAL CRYST                 
from desalsim.nanofiltration_unit_f import OsmoticPressure
from desalsim.nanofiltration_unit_f import NFMass
from desalsim.nanofiltration_unit_f import NfEnergy

from desalsim.med_unit_f import MEDCalculator

from desalsim.thermal_cryst_f import thermal_calc
from desalsim.thermal_cryst_f import conc_cal
from desalsim.thermal_cryst_f import calculate_energy


from desalsim import constants

from desalsim.economic_f import revenue
from desalsim.economic_f import econom

from desalsim import scaleup
from desalsim.density_calc import density_calc
import numpy as np
import pandas as pd
#%%
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
MW_CaSO4=constants.MW_Ca+constants.MW_so4

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

# Function to create NFMass objects for different components
def create_nfmass_objects(components, C_in, rjr_values, Wrec, Qf):
    return [NFMass(comp, Ci, rjr, Wrec, Qf) for comp, Ci, rjr in zip(components, C_in, rjr_values)]

# Create NFMass objects for different components
nfmass_objects = create_nfmass_objects(components, Ci_in, rjr_values, Wrec, Qf_nf)

Cconc = [nf_mass.Cconci for nf_mass in nfmass_objects]
Cperm = [nf_mass.Cpermi for nf_mass in nfmass_objects]
Qperm = nfmass_objects[0].Qperm  # kg/hr
Qconc= nfmass_objects[0].Qconc # kg/hr

# Calculate Osmotic Pressure
P_osmo_f = OsmoticPressure(c_values, z_values, T).osmotic_pressure_calculation()
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
#%%calculation
"""--------MED--------"""

#Feed concentration
components = ['Na', 'Cl', 'K', 'Mg', 'Ca', 'SO4']
Cin_med = Cperm

#Feed flow rate 
Qf_med =Qperm

#input conditions
T=20
#feed flow density 
d=density_calc(T, sum(Cin_med)) 
Mf_med=Qf_med*d/1000 #Mass flow rate (units: kg/hr)

#assumptions:
T_in=40 #(oC)
T_N=45 #Temperature in the last effect (oC)
N=2 #Number of effects (-)
Cp_w=4182 # specific heat capacity of water (j/kgC)
cp_sol=4184 # specific heat capacity of solution (j/kgC)
T_cw_in=25 #intake cooling water temperature (oC)
T_cw_out=35 #out cooling water temperature (oC)
T_s=70 #steam temperature oC
DT_loss=1 #temperature difference (oC)
dp=0.1  # pressure drop (units: bar)
dp_slurry=1 # pressure drop (units: bar)
npump=0.8 #pump efficiency (units: -)
Cb_out=200 #Brine Concentration leaving effect n (unit: g/l)
Xr=5.5 # brine circulation flow rate (units: -)

#latent heat of motive steam:
if T_s<=55:
    lh_s=2370
elif T_s>55 and T_s<=60:
    lh_s=2358
elif (T_s>60) and (T_s<=65):
    lh_s=2345
elif (T_s>65) and (T_s<=70):
    lh_s=2333
elif (T_s>70) and (T_s<=75):
    lh_s=2321

# Create an instance of the MEDCalculator class
med_dat = MEDCalculator(Qf_med, Mf_med, Cin_med[0], Cin_med[1], Cin_med[2], Cin_med[3], Cin_med[4], Cin_med[5], T_in)

# Call methods to perform calculations
med_dat.salinity_calc()
med_dat.temperature_calc(DT_loss, T_N, T_s)
med_dat.mass_balance_med(Cb_out)
med_dat.performance_parameters(lh_s,  T_cw_in)
med_dat.output_concentration()

# Sum results
#Concentration of brine stream 
Cconc_med = [med_dat.CNa_out, med_dat.CCl_out, med_dat.CK_out, med_dat.CMg_out, med_dat.CCa_out, med_dat.CSO4_out]

#Brine flow rate 
Qout_med=med_dat.Qb

#Distillate water flow rate 
Qprod_med=med_dat.Qdist

#Calculate circulation flow rate 
Qr=Xr*Qf_med

# Calculate density for output concentration
d_b = density_calc(45, sum(Cconc_med))

print("Sum of output concentrations: " + str(round(sum(Cconc_med),2))+"g/l")
print("-----------------------------------------")

# Calculate energy consumption
E_el_med = ((Qf_med * 3.5 + med_dat.QCW * 3600 * 2 + (Qr + med_dat.Qb) * 3.5 + med_dat.Qdist * 1) / (1000 * npump)) * 1e5 / 3600 / 1000  # kWh
print("Electrical energy consumption: " + str(round(E_el_med,2)) + " kWh")

SEC_el = E_el_med / (Qf_med / d)  # kWh/m3 feed
print("Specific energy consumption (electrical) per m3 feed: " + str(round(SEC_el,2)) + " kWh/m3")

SEC_el_prod = E_el_med / (med_dat.Qdist / 1000)  # kWh/m3 dist water
print("Specific energy consumption (electrical) per m3 product (distilled water): " + str(round(SEC_el_prod,2)) + " kWh/m3")
print("-----------------------------------------")

E_th_med = med_dat.Q_Tot
print("Total thermal energy consumption: " + str(round(E_th_med,2)) + " kW")

SEC_th = E_th_med / (Qf_med / d)  # kWh_th/m3
print("Specific energy consumption (thermal) per m3 feed: " + str(round(SEC_th,2)) + " kWh_th/m3")
print("-----------------------------------------")

#Calculate required cooling water 
Qcw_med = med_dat.QCW * 3600 #units: kg/hr
print("Cooling water flow rate: " + str(round(Qcw_med,2)) + " kg/hr")
print("-----------------------------------------")

# Chemical consumption
Cchem_med = 0  # Placeholder for chemical consumption, update as needed

# Print the results
print("Total Chemical Consumption: " + str(round(Cchem_med,2))+"kg/hr")

#%%
#calculation thermal crystallizer 
"""--------TCr--------"""
C_in_mix=[]

#constants
R=8.31446261815324 #gas constant
d_salt=1.974463 #density_salt (NaCl)=1974.5 kg/m3

    # Thermodynamics 
Cp_f=3.14 # Feed specific heat capacity (units: KJ* Kg*oC)
CP_cw=4.187 # Water specific heat capacity (units: KJ* Kg*oC)
UA=45990
LHV_v=2107.92 # kj/kg (gathered from steam table)
LHV_s=2357.69 # kj/kg (gathered from steam table)

#Assumptions
dp=0.1 # pressure drop (units: bar)
dp_slurry=3.5 #pressure drop for slurry streams (units: bar)
dp_f=3.5 #pressure drop for feed (units: bar)
dp_w=1 #pressure drop for qater streams (units: bar)
dp_cw=2 #pressure drop for cooling water (units: bar) 
salt_mois=20 # % moinsture in salt stream (units: %)
npump=0.8 #pump efficiency (units: -)


#Input parameters
T_in=40 #oC
T_top=60 #oC
T_cw_f=25 #oC
T_cw_o=40 #oC
T_op=60 #oC

#Feed flowrate
Qf_tcr=Qconc+ Qout_med#kg/h

#Feed concentration 

Cf_tcr=sum(Cconc_med)*Qout_med/Qf_tcr+sum(Cconc)*Qconc/Qf_tcr
for i in range(len(Cconc)):
    C_in_mix.append(Cconc[i]*Qout_med/Qf_tcr+Cconc[i]*Qconc/Qf_tcr)
    
Cf_tcr_in=C_in_mix #[Na, Cl, K, Mg, Ca, SO4 ]
Cf_s=sum(Cf_tcr_in) #total salt concentration (g/L)
d_sol=density_calc(T_in, Cf_tcr)/1000 #density of feed
Cf_caso4=Cf_tcr_in[4]*MW_CaSO4/MW_Ca #CaSO4 concentration in feed solution (g/L)

#Calculation 
th_cryst_dat=thermal_calc(T_op, Qf_tcr, Cf_tcr, Cf_caso4, T_in, Cf_tcr_in, salt_mois, LHV_v, LHV_s, T_cw_o, T_cw_f)
th_cryst_dat.mass_bal_cryst()
th_cryst_dat.heat_bal_cryst()

#Recovered salt flow rate (kg/h)
M_Nacl=th_cryst_dat.solid_mass
print("Mass flowrate of recovered salt is "+str(round(M_Nacl,2))+"kg/hr")

# Calculation of the concentration of different ions in the solution
th_cryst_dat_2=conc_cal(Qf_tcr, M_Nacl , 'Na',Cf_tcr_in[0], 'cl',Cf_tcr_in[1],'k', Cf_tcr_in[2], 'mg', Cf_tcr_in[3], 'ca', Cf_tcr_in[4], 'so4', Cf_tcr_in[5], T_in)
th_cryst_dat_2.molarity()
th_cryst_dat_2.conc_salt_comp()
th_cryst_dat_2.salt_conc()

#Evaporation mass flow rate (water recovery), (kg/h)
Q_evap_mass_tcr=th_cryst_dat.ev_mass
print("Mass flowrate of recovered water is "+str(round(Q_evap_mass_tcr,2))+"kg/hr")
print("-----------------------------------------")
#Salt stream concentration (g/L)
Csalt_out=[th_cryst_dat_2.CNa,th_cryst_dat_2.CCl, th_cryst_dat_2.CK, th_cryst_dat_2.CMg, th_cryst_dat_2.CCa, th_cryst_dat_2.CSO4]

#Cooling water requirements (kg/h)
Qcw_tcr=th_cryst_dat.cw_mass
print("Mass flowrate of required cooling water is "+str(round(Qcw_tcr,2))+"kg/hr")
print("-----------------------------------------")
#Calculate energy consumption 
heat_req=th_cryst_dat.heat_req
E_el_th_Cr, E_th_th_Cr, SEC_el_f, SEC_el_NaCl, SEC_el_w = calculate_energy(Qf_tcr, Q_evap_mass_tcr, Qcw_tcr, M_Nacl, heat_req, d_sol, dp_f, dp_w, dp_slurry, dp_cw, npump)
print(f"SEC_el_prod is {SEC_el_w} KWh/m3 product")
print(f"SEC_el_prod2 is {SEC_el_NaCl} KWh/kg of NaCl product")
print("SEC_th_prod is "+str( round((E_th_th_Cr/(Qf_tcr/1000)),2))+" KWh/m3 intake")

# Chemical consumption
Cchem_tcr = 0  # Placeholder for chemical consumption, update as needed


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
        
#%% Sumsummarize results   

tec1=indic("NF", Qconc, Qf_nf, Qperm, "none", 0, "none", E_el_nf, 0, Ci_in, Cconc, QHCl_nf, Qantsc_nf)   
tec1.techn_indi()

tec2=indic("MED", Qout_med, Qf_med, Qprod_med, "water", 0, "none",  E_el_med, E_th_med, Cin_med, Cconc_med, Cchem_med,0)
tec2.techn_indi()

tec3=indic("Thermal Cryst", 0, Qf_tcr, Q_evap_mass_tcr,"water", M_Nacl, "NaCl", E_el_th_Cr, E_th_th_Cr, Cf_tcr, 0,0,0 )
tec3.techn_indi()


##list results 
tec_names=[tec1.tech, tec2.tech, tec3.tech]
E_el_all=[tec1.E_el, tec2.E_el, tec3.E_el]
E_th_all=[tec1.E_th, tec2.E_th, tec3.E_th]
Cout_all=[tec1.Cout, tec2.Cout, tec3.Cout]
Qout_all=[tec1.Qout, tec2.Qout, tec3.Qout]
Qchem_all=[tec1.chem,tec2.chem, tec3.chem]

sum_el_en=sum(E_el_all)
sum_th_en=sum(E_th_all)
#%%
Qsw=Qsw
#Technical indicators 
#T-1) quantity of water production 
Qw_tot=tec1.Qwater+tec2.Qwater+ tec3.Qwater
#print("total water production is " + str(Qw_tot))
abs_Qw=Qw_tot
#T-1b) Water recovery  
rec=Qw_tot/Qsw*100 #%
#T-2) Specific electrical energy consumption #kwh/kg of desalinated water 
SEC=sum(E_el_all)/Qsw
print("specific electrical consumption is " +str(SEC)+ " Kwh/kg of desalinated water")

#T-3)Brine production -> kg of brine/ kg of seawater 
Qbrine_f=0
for i in range(len(tec_names)):
    if tec_names[i]=="Thermal Cryst":
        if Cout_all[i]>(constants.Csw):
            Qbrine_f=Qbrine_f+Qout_all[i]
        else:
            Qbrine_f=Qbrine_f+0
Q_br_prod=Qbrine_f/Qsw
print("brine production is "+str(Q_br_prod)+" kg of brine/kg of seawater")
#------------------------------------------------------------------------------
#%%
#Environmental 
#Env-1)Climate change impact ->Carbon dioxide emission kg co2
emis=(sum(E_el_all)*constants.CO2_el)
emis_t=((sum(E_el_all)*constants.CO2_el)+(sum(E_th_all)*constants.CO2_st))*constants.hr
print("Total specific Carbon dioxide emission "+str(emis_t)+" kg co2")
#Env-2) Resource utilization -> chemical consumption in kg of chemical/kg of seawater 
Chem_cons=(sum(Qchem_all))/Qsw #check units and convert ml/l to kg with densities 
print("Total chemical consumption is "+str(Chem_cons)+" kg of chemical/kg of seawater ")

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

#Assumptions for CAPEX 
inst_percent=0.25
buildings_percent=0.2
land_percent=0.06
indirect_c_percent=0.15
workinf_c_percent=0.2
capex_assumptions=[inst_percent,buildings_percent, land_percent, indirect_c_percent, workinf_c_percent]


#Assumptions for OPEX
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

# Initialize values
CAPEX=0
OPEX=0
CO2_c=0

 # Create list with capital cost for each unit.
capex_list=[]
 # Create list with operating for each unit.
opex_list=[]

#summarise results for utilities 
el_conc=[E_el_nf,E_el_med,E_el_th_Cr] 
s_conc=[0,E_th_med,E_th_th_Cr]
chem1_conc=[Qantsc_nf,Cchem_med,Cchem_tcr]
chem1_pr=[0.002,0,0]
chem2_conc=[0,0,0]
chem2_pr=[0,0,0]
cw_conc=[0,Qcw_med, Qcw_tcr]
wat_conc=[0,0, 0]

if sum(wat_conc)==0:
    Q_w_in=0

# Calculation of scale-up equipment cost 
# Import equipment cost for reference scenario from constants function
eq_c=[constants.eq_c_nf, constants.eq_c_med, constants.eq_c_th_cr]   
    # Capacity of reference scenario
Mf_basic_sc=[constants.eq_c_nf,constants.eq_c_med,constants.eq_c_th_cr]
Mf_sce=[Qsw, Qperm,Qf_tcr]
    
    # Calculation of the new equipment cost
for i in range(len(eq_c)):
    if Mf_basic_sc[i]!=Mf_sce[i]:
        eq_c[i]=scaleup.scaleup(eq_c[i],Mf_basic_sc[i],Mf_sce[i])
        
#Cost calculation 
for i in range(len(eq_c)):
    total_econom=econom(eq_c[i], el_conc[i], s_conc[i], chem1_conc[i], chem1_pr[i],chem2_conc[i], chem2_pr[i], cw_conc[i], wat_conc[i])
    total_econom.capex_calc(capex_assumptions)
    total_econom.opex_calc(hr, el_pr, s_pr, cw_pr, w_pr, economic_assumptions)
    CAPEX=total_econom.t_capital_inv
    OPEX=total_econom.opex
    opex_list.append(total_econom.opex)
    capex_list.append(total_econom.t_capital_inv)    
print("Total operating cost (OPEX) is "+str(OPEX)+ " Euro/year") 
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

#%%present results
#Summary tables 
    #Table C: Ion concentration 
sum_table_C={'F1: Seawater': Ci_in,
           'F2: NF permeate': Cperm, 
           'F3: NF concentrate': Cconc,
           'F4: MED brine': Cconc_med, 
           'F5: MED water': [0,0,0,0,0,0],
           'F6: Thermal Cryst water': [0,0,0,0,0,0],
           'F7: Thermal Cryst salt': Csalt_out,
           }
sum_table_C=pd.DataFrame(sum_table_C)
print(sum_table_C)

    #Table D: Density
sum_table_d={'F1: Seawater': density_calc(25,sum(Ci_in))/1000,
           'F2: NF permeate': density_calc(25,sum(Cperm))/1000, 
           'F3: NF concentrate': density_calc(25,sum(Cconc))/1000,
           'F4: MED brine': density_calc(25,sum(Cconc_med))/1000, 
           'F5: MED water': 1,
           'F6: Thermal Cryst water': 1,
           'F7: Thermal Cryst salt': density_calc(25,sum(Csalt_out))/1000,
           }
sum_table_d=pd.DataFrame(sum_table_d, index=['density'])
print(sum_table_d)

        #Table F: Flow rates  
sum_table_f={'F1: Seawater': round(Qsw,2),
           'F2: NF permeate': round(Qperm,2), 
           'F3: NF concentrate': round(Qconc,2),
           'F4: MED brine': round(Qout_med,2), 
           'F5: MED water': round(Qprod_med,2),
           'F6: Thermal Cryst water': round(Q_evap_mass_tcr,2),
           'F7: Thermal Cryst salt':round( M_Nacl,2),
           'F8: NF antiscalant': round(Qantsc_nf,2), 
           'F9: MED Cooling water': round(Qcw_med,2), 
           'F10: Thermal Cryst Cooling water': round(Qcw_tcr,2), 
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
      label = ["Seawater", "Water", "NaCL",  "NF", "MED", "Cryst", "Mixed Salt"],
      color = "blue"
    ),
    link = dict(
      source = [0, 3,3,4,4, 5,5], 
      target = [3,4,5,5,1,1,6],
      value = [Qsw, Qperm, Qconc, Qout_med, Qprod_med, Q_evap_mass_tcr, M_Nacl]
  ))])
color_for_nodes = ["lightsteelblue","darkcyan","maroon", "midnightblue", "midnightblue", "midnightblue", "maroon"]
fig.update_traces(node_color = color_for_nodes)
fig.update_layout(title_text="Sankey diagram for Example: Mass flow rates ", font_size=15)
fig.show()
plot(fig)

#%% Results to excel 
dfel=pd.DataFrame(E_el_all ,tec_names)
dfeth=pd.DataFrame(E_th_all, tec_names)
dfind=(Qw_tot, abs_Qw, rec, sum_el_en, sum_th_en, emis_t, Q_w_in,  OPEX, CAPEX)
dfprod=(Qw_tot,M_Nacl, 0, 0, 0, 0, 0)

dfprodn=("Water (kg/hr)", "NaCl (kg/hr)", "Mg(OH)\u2082 (kg/hr)", "Ca(OH)\u2082 (kg/hr)", "Na\u2082SO\u2084(kg/hr)", "1M NaOH (kg/hr)","1M HCl(kg/hr)")
ind=np.array(["Water production", "absolute water production", "water recovery","Total electrical consumption", "Total thermal energy consumption", "Carbon dioxide emission kg co2/year",
              "Water footprint","OPEX", "CAPEX"])
units=pd.DataFrame(["kg/hr", "kg/hr", "%","KWh", "KWh"," kg CO\u2082/year","kg/hr","€/year", "€"])
dfind_t=pd.DataFrame(data=dfind, index=ind)
dfprodn=pd.DataFrame(data=dfprod, index=dfprodn)
with pd.ExcelWriter('results_example .xlsx') as writer:
    sum_table_f.to_excel(writer,sheet_name="example ")
    sum_table_C.to_excel(writer,sheet_name="example ", startcol=0, startrow=4, header=True)
    dfel.to_excel(writer,sheet_name="example ", startcol=0, startrow=14, header=True)
    dfeth.to_excel(writer,sheet_name="example ", startcol=2, startrow=14, header=True)
    dfprodn.to_excel(writer,sheet_name="example ", startcol=0, startrow=23, header=True)
    dfind_t.to_excel(writer,sheet_name="indicators")
    units.to_excel(writer,sheet_name="indicators",startcol=2, startrow=0, header=True)