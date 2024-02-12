#Scenario 2
##treatment train: NF-> MED-> THERMAL CRYST                 
from nanofiltration_unit_f import OsmoticPressure
from nanofiltration_unit_f import molarity
from nanofiltration_unit_f import NFMass
from nanofiltration_unit_f import NfEnergy

from med_unit_f import MEDCalculator

from thermal_cryst_f import thermal_calc
from thermal_cryst_f import conc_cal
from thermal_cryst_f import calculate_energy

from selected_indicators import indic

import constants

from economic_f import revenue
from economic_f import econom

import scaleup
import density_calc
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
Qconc= nfmass_objects[0].Qconc # kg/hr

# Calculate Osmotic Pressure
P_osmo_f = OsmoticPressure(c_values, z_values, T).osmotic_pressure_calculation()
P_osmo_p = OsmoticPressure(Cperm, z_values, T).osmotic_pressure_calculation()
P_osmo_c = OsmoticPressure(Cconc, z_values, T).osmotic_pressure_calculation()

d_p=density_calc(T-273, sum(Cperm))

#Calculate Energy consumption 
E_el_nf=NfEnergy(P_osmo_c, P_osmo_f, P_osmo_p, dp, d_p, Qperm, Qf_nf, d_in,n)
result=E_el_nf.calculate_energy_consumption()
for key, value in result.items():
        print(f"{key}: {value}")
        
QHCl_nf=0
Qantsc_nf=0
#%%calculation MED 
#constants 
#conditions
T=20+273.15
d=1018.494

#assumptions:
conc_f=8.5
T_in=40 #(oC)
T_N=45 #Temperature in the last effect (oC)
N=2 #Number of effects (-)
Cp_w=4182 # specific heat capacity of water (j/kgC)
cp_sol=4184
T_cw_in=25 #intake cooling water temperature (oC)
T_cw_out=35 #out cooling water temperature (oC)
T_s=70 #steam temperature oC
DT_loss=1 #oC
T3=69
#Ms=819
dp=0.1
dp_slurry=1
npump=0.8
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

#calculations 
Qf_med=Qperm
#Cin_med=sum([Comp1.Cpermi, Comp2.Cpermi,Comp3.Cpermi, Comp4.Cpermi, Comp5.Cpermi, Comp6.Cpermi])
Cin_med_t=sum(Cperm)
Cin_med=Cperm
print("Cin_ med is " + str(Cin_med))
med_dat=MEDCalculator(Qf_med,Cin_med[0],Cin_med[1],Cin_med[2],Cin_med[3],Cin_med[4],Cin_med[5])
#med_dat=med_calc(Qf_med, Comp1.Cpermi, Comp2.Cpermi,Comp3.Cpermi, Comp4.Cpermi, Comp5.Cpermi, Comp6.Cpermi)
med_dat.sal_calc()
med_dat.temp_calc()
med_dat.mass_bal_med()
med_dat.performance_par()
med_dat.out_conc()
#sum results 
Cconc=[med_dat.CNa_out, med_dat.CCl_out, med_dat.CK_out,med_dat.CMg_out,med_dat.CCa_out,med_dat.CSO4_out]
print(sum(Cconc))
Cconc_med_all=Cconc
Cconc_med=sum(Cconc)
Qout_med=med_dat.Qb
Qprod_med=med_dat.Qdist

Qr=5.5*Qf_med
#assumptions for dp: dp=2 bar for cw 
#dp for feed is the same as the recirculation: 3.5bar 
#dp for distillate water 1bar 
E_el_med=((Qf_med*3.5+med_dat.QCW*3600*2+(Qr+med_dat.Qb)*3.5+med_dat.Qdist*2)//(1000*npump))*1e5/3600/1000 #kwh
print("E_el_med is " +str(E_el_med)+ " kWh")
Qcw_med=med_dat.QCW*3600
print("Cooling water flow rate is "+str(Qcw_med)+" kg/hr")
E_th_med=med_dat.Q_Tot
print("Total thermal energy consumption is "+str(E_th_med)+" KW")
Cchem=0
#%%
#calculation thermal crystallizer 
C_in_mix=[]
dp=0.1
dp_slurry=1
#Input parameters
T_in=25 #oC
T_top=60 #oC
T_cw_f=25 #oC
T_cw_o=40 #oC
T_op=60
Qth_cr_in=Qout_med+nanofiltration_unit.Qconc #kg/h
Cf_s=Cconc_med*Qout_med/Qth_cr_in+nanofiltration_unit.Cout_conc*nanofiltration_unit.Qconc/Qth_cr_in
for i in range(len(Cconc)):
    C_in_mix.append(Cconc[i]*Qout_med/Qth_cr_in+nanofiltration_unit.Ci_out_c[i]*nanofiltration_unit.Qconc/Qth_cr_in)
Cf_caso4=0.075132776
th_cryst_dat=thermal_calc( Qth_cr_in, Cf_s, Cf_caso4, T_in, C_in_mix )
th_cryst_dat.mass_bal_cryst()
th_cryst_dat.heat_bal_cryst()
print(th_cryst_dat.ev_mass)
print(th_cryst_dat.solid_mass)
print(th_cryst_dat.nacl_sat_wt)
print("bpe is " + str(th_cryst_dat.BPE))
print ("hear required is "+str(th_cryst_dat.heat_req))
print ("steam temp is "+str(th_cryst_dat.T_s))
print ("steam mass is "+str(th_cryst_dat.steam_mass))
print ("cw mass is "+str(th_cryst_dat.cw_mass))
x=th_cryst_dat.solid_mass

th_cryst_dat_2=conc_cal(Qth_cr_in, x , 'Na',C_in_mix[0], 'cl',C_in_mix[1],'k', C_in_mix[2], 'mg', C_in_mix[3], 'ca', C_in_mix[4], 'so4', C_in_mix[5])
th_cryst_dat_2.molarity()
th_cryst_dat_2.conc_salt_comp()
th_cryst_dat_2.salt_conc()
M_Nacl=th_cryst_dat.solid_mass
Qcw_th=th_cryst_dat.cw_mass
Q_evap_mass=th_cryst_dat.ev_mass
Csalt_out_th=[th_cryst_dat_2.CNa,th_cryst_dat_2.CCl, th_cryst_dat_2.CK, th_cryst_dat_2.CMg, th_cryst_dat_2.CCa, th_cryst_dat_2.CSO4]
E_el_th_Cr=((Qth_cr_in*3.5+th_cryst_dat.ev_mass*1+5.5*th_cryst_dat.nacl_sat_wt*3.5+2*Qcw_th*2)/(1000*npump))*1e5/3600/1000#kwh

Qchem=0
Qchem=0
E_th_th_Cr=th_cryst_dat.heat_req/3600 #kwh
#%%indicators
#%%main 

#Call units' functions 
#def __init__(self, tech, Qout, Qin, Qprod_1, prod_1_name, Qprod_2, prod_2_name, E_el, E_th, Cin, Cout, chem1, chem2)       
#print(E_el_nf)
tec1=indic("NF", nanofiltration_unit.Qconc, nanofiltration_unit.Qsw, nanofiltration_unit.Qperm, "none", 0, "none", nanofiltration_unit.E_el_nf, 0, nanofiltration_unit.Cin_conc, nanofiltration_unit.Cout_conc, nanofiltration_unit.Qhcl,nanofiltration_unit.Qantsc)     
tec1.techn_indi()
tec1.env_ind()
print("The CO2 emissions based on the energy consumption of "+ tec1.tech +" are "+str(round(tec1.emis))+" kgCO2/m3 of product")
#
tec2=indic("MED", Qout_med, Qf_med, Qprod_med, "water", 0, "none",  E_el_med, E_th_med, Cin_med, Cconc_med, 0,0)
tec2.techn_indi()
tec2.env_ind()
print("The CO2 emissions based on the energy consumption of "+ tec2.tech +" are "+str(round(tec2.emis))+" kgCO2/m3 of product")

tec3=indic("Thermal Cryst", 0, Qth_cr_in, th_cryst_dat.ev_mass,"water", th_cryst_dat.nacl_sat_wt, "NaCl", E_el_th_Cr, E_th_th_Cr, Cf_s, 0,0,0 )
tec3.techn_indi()
tec3.env_ind()

##list results 
tec_names=[tec1.tech, tec2.tech, tec3.tech]
E_el_all=[tec1.E_el, tec2.E_el, tec3.E_el]
E_th_all=[tec1.E_th, tec2.E_th, tec3.E_th]
Cout_all=[tec1.Cout, tec2.Cout, tec3.Cout]
Qout_all=[tec1.Qout, tec2.Qout, tec3.Qout]
Qchem_all=[tec1.chem,tec2.chem, tec3.chem]
#%%
Qsw=nanofiltration_unit.Qsw
#Technical 
#T-1) quantity of water production 
Qw_tot=tec1.Qwater+tec2.Qwater+ tec3.Qwater
#print("total water production is " + str(Qw_tot))
abs_Qw=Qw_tot
#T-1b) Water recovery  
rec=Qw_tot/Qsw*100 #%
#T-2)energy performance
#T-2a) specific electrical energy consumption #kwh/kg of desalinated water 
SEC=sum(E_el_all)/Qsw
print("specific electrical consumption is " +str(SEC)+ " Kwh/kg of desalinated water")
#2b) share of renewable energy %
#2cb) Energy consumption 
Tot_el=sum(E_el_all) #kWh 
Tot_th=sum(E_th_all) #kWh
#T-3)efficiency 
#T-3a) Resource efficiency
Rr=(Qw_tot)/(Qsw)
Rr2=(Qw_tot)/(Qsw)
print("resource efficiency of the system is "+str(Rr))
#T-3b) Salt recovery 
Salt_r=0
#T-3b)brine production -> kg of brine/ kg of seawater 
Qbrine_f=0
for i in range(len(tec_names)):
    if tec_names[i]=="Thermal Cryst":
        if Cout_all[i]>(constants.Csw):
            Qbrine_f=Qbrine_f+Qout_all[i]
        else:
            Qbrine_f=Qbrine_f+0
Q_br_prod=Qbrine_f/nanofiltration_unit.Qsw
print("brine production is "+str(Q_br_prod)+" kg of brine/kg of seawater")
#------------------------------------------------------------------------------
#%%
#Environmental 
#Env-1)Climate change impact ->Carbon dioxide emission kg co2
emis=(sum(E_el_all)*constants.CO2_el)
emis_t=((sum(E_el_all)*constants.CO2_el)+(sum(E_th_all)*constants.CO2_st))*constants.hr
print("Total specific Carbon dioxide emission "+str(emis_t)+" kg co2")
#Env-2) Resource utilization -> chemical consumption in kg of chemical/kg of seawater 
Chem_cons=(sum(Qchem_all))/nanofiltration_unit.Qsw #check units and convert ml/l to kg with densities 
print("Total chemical consumption is "+str(Chem_cons)+" kg of chemical/kg of seawater ")
#Env-3) The aquatic eco-toxic impact of brine disposal -> ECO-TOXICITY 
#C_brine=edbm_unit_recycling.Cbrine_out #[ CNa, CCl, CK, CMg, CCa, CSO4, Chco3, CH, COH]
#ETP=0
#CF=4.62e-1
#for i in range(len(C_brine)-2):
#    ETP=ETP+C_brine[i]*CF
#print("The aquatic eco-toxic impact of brine disposal is "+str(ETP)+"")
#Env-4) Environmental prodection -> minimization of waste in kg of brine/CO2 emissions 
Min_waste=(Qout_med+nanofiltration_unit.Qconc-Qbrine_f)/((sum(E_el_all)*constants.CO2_el)+(sum(E_th_all)*constants.CO2_st))
print("Emissions from minimization of waste is "+str(Min_waste)+" kg of brine/CO2 emissions ")
#Env-5) Resource recovery efficiency -> Resource efficiency (%)

print("Resource recovery efficiency is "+str()+" ")
print(""+str()+"")


#%%
#Economic 
CAPEX=0
OPEX=0
CO2_c=0
OPEX_with_ext=0
eq_c=[constants.eq_c_nf, constants.eq_c_med, constants.eq_c_th_cr]   
opex_list=[]
capex_list=[]
Mf_basic_sc=[2290.858979,1695.830673,199.5094909]
Mf_sce=[nanofiltration_unit.Qsw, nanofiltration_unit.Qperm,Qth_cr_in]
el_conc=[nanofiltration_unit.E_el_nf,E_el_med,E_el_th_Cr] 
s_conc=[0,E_th_med,E_th_th_Cr]
chem1_conc=[nanofiltration_unit.Qantsc,0,0]
chem1_pr=[0.002,0,0]
chem2_conc=[0,0,0]
chem2_pr=[0,0,0]
cw_conc=[0,Qcw_med, Qcw_th]
wat_conc=[0,0, 0]
if sum(wat_conc)==0:
    Q_w_in=0
    
Mf_sce=[nanofiltration_unit.Qsw, nanofiltration_unit.Qperm,Qth_cr_in]
for i in range(len(eq_c)):
    if Mf_basic_sc[i]!=Mf_sce[i]:
        eq_c[i]=scaleup.scaleup_eq(eq_c[i],Mf_basic_sc[i],Mf_sce[i],tec_names[i])
        
   
for i in range(len(eq_c)):
    total_econom=econom(eq_c[i], el_conc[i], s_conc[i], chem1_conc[i], chem1_pr[i],chem2_conc[i], chem2_pr[i], cw_conc[i],wat_conc[i])
    total_econom.capex_calc()
    total_econom.opex_calc()
    CAPEX=CAPEX+total_econom.t_capital_inv
    OPEX=OPEX+total_econom.opex
    opex_list.append(total_econom.opex)
    print("opex for i is "+str(OPEX))
    total_econom.carbon_calc() 
    CO2_c=CO2_c+total_econom.carbon_c
    OPEX_with_ext=OPEX_with_ext+total_econom.opex
    
    capex_list.append(total_econom.t_capital_inv)
    
print("total capex of system: " + str(round(CAPEX))+ " euro")
print("total OPEX without externalities of system: " + str(round(OPEX))+ " euro/year") 
print("total OPEX with externalities of system: " + str(round(OPEX_with_ext))+ " euro/year") 
print("total carbon tax cost: " + str(round(CO2_c))+ " euro/year")     

#amortisation factor 
a=(constants.r*(1+constants.r)**constants.lf)/((1+constants.r)**constants.lf-1)
reve_t=0
prd=[Qprod_med+Q_evap_mass, M_Nacl]    
prd_name= ["water", "nacl"]   
reve_list=[]
for i in range(len(prd)):
    rev_calc=revenue(prd[i], prd_name[i])    
    rev_calc.rev()
    print("rev_calc.rev_prd for " + prd_name[i]+" is " + str(rev_calc.rev_prd))
    reve_t = reve_t+rev_calc.rev_prd
    reve_list.append(rev_calc.rev_prd)
#EC-1) Product value -> levelized cost 
print(""+str()+"")
#EC-2)impact of externalities on the production cost 
#EC-2a)production cost with externalities in euro/m3 of seawater 
Pr_c_ext=OPEX_with_ext/Qsw
print("production cost with externalities is "+str()+" euro/m3 of seawater")
#EC-2b)production cost without externalities in euro/m3 of seawater
Pr_c=OPEX/Qsw
print("production cost without externalities is "+str(Pr_c)+" euro/m3 of seawater")
#EC-3) profitability -> production efficiency in euro/euro 
prd_eff=(reve_t)/(CAPEX*a+OPEX)
print("production efficiency is "+str(prd_eff)+"")
#EC-4) profitability -> economic margin 
Rev=reve_t
Ec_m= (CAPEX*a+OPEX-reve_t)/(Qw_tot*constants.hr/1000)
print("economic margin is "+str(Ec_m)+"euro/m3 of water")
OPEX=OPEX
CAPEX=CAPEX
oper_c_t=CAPEX*a+OPEX
print("CAPEX_sc3 is "+str(CAPEX)+"")
print(""+str()+"")
#%% Levelized cost (joint cost)
#product 1: water -> whole system 
LPCw=(oper_c_t-sum(reve_list[1:]))/(Qw_tot/1000*constants.hr)
print("LPCw is "+str(LPCw) +" euro/m3")
#product 2: NaCl -> water chain (NF->MED->TCr)
C_w_chain=capex_list[0]+capex_list[1]+capex_list[2]
O_w_chain=opex_list[0]+opex_list[1]+opex_list[2]
LPCnacl_1=(C_w_chain*a+O_w_chain-reve_list[0])/(M_Nacl*constants.hr)
print("LPCnacl is "+str(LPCnacl_1)+" euro/kg")
LPCnacl_2=(oper_c_t-reve_list[0])/(M_Nacl*constants.hr)
print("LPCnacl is "+str(LPCnacl_2)+" euro/kg")
#product 3: Mg(OH)2 -> whole process chain 
LPCmg_1=0
LPCmg_2=0
#product 3: Na2so4 -> whole process chain 
LPCna2so4_1=0
LPCna2so4_2=0
#product 5: HCl -> whole process chain 
LPChcl_1=0
LPChcl_2=0
#product 6: NaOH -> whole process chain 

LPCnaoh_1=0
LPCnaoh_2=0
#%%present results

#%%
sum_el_en=sum(E_el_all)
sum_th_en=sum(E_th_all)

sum_table_C={'F1: Seawater': nanofiltration_unit.Ci_in,
           'F2: NF permeate': nanofiltration_unit.Ci_out_p, 
           'F3: NF concentrate': nanofiltration_unit.Ci_out_c,
           'F4: MED brine': Cconc_med_all, 
           'F5: MED water': [0,0,0,0,0,0],
           'F6: Thermal Cryst water': [0,0,0,0,0,0],
           'F7: Thermal Cryst salt': Csalt_out_th,
           }
sum_table_C=pd.DataFrame(sum_table_C)
print(sum_table_C)
sum_table_d={'F1: Seawater': density_calc.density_calc(25,sum(nanofiltration_unit.Ci_in))/1000,
           'F2: NF permeate': density_calc.density_calc(25,sum(nanofiltration_unit.Ci_out_p))/1000, 
           'F3: NF concentrate': density_calc.density_calc(25,sum(nanofiltration_unit.Ci_out_c))/1000,
           'F4: MED brine': density_calc.density_calc(25,sum(Cconc_med_all))/1000, 
           'F5: MED water': 1,
           'F6: Thermal Cryst water': 1,
           'F7: Thermal Cryst salt': density_calc.density_calc(25,sum(Csalt_out_th))/1000,
           }
sum_table_d=pd.DataFrame(sum_table_d, index=['density'])
print(sum_table_d)


sum_table_f={'F1: Seawater': round(nanofiltration_unit.Qsw,2),
           'F2: NF permeate': round(nanofiltration_unit.Qperm,2), 
           'F3: NF concentrate': round(nanofiltration_unit.Qconc,2),
           'F4: MED brine': round(Qout_med,2), 
           'F5: MED water': round(Qprod_med,2),
           'F6: Thermal Cryst water': round(Q_evap_mass,2),
           'F7: Thermal Cryst salt':round( M_Nacl,2),
           'F8: NF antiscalant': round(nanofiltration_unit.Qantsc,2), 
           'F9: MED Cooling water': round(Qcw_med,2), 
           'F10: Thermal Cryst Cooling water': round(Qcw_th,2), 
           }
sum_table_f=pd.DataFrame(sum_table_f, index=['flow rate'])
print(sum_table_f)

#sankey diagram 
import plotly.graph_objects as go
from plotly.offline import plot

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
      value = [nanofiltration_unit.Qsw, nanofiltration_unit.Qperm, nanofiltration_unit.Qconc, Qout_med, Qprod_med, Q_evap_mass, M_Nacl]
  ))])
color_for_nodes = ["lightsteelblue","darkcyan","maroon", "midnightblue", "midnightblue", "midnightblue", "maroon"]
fig.update_traces(node_color = color_for_nodes)
fig.update_layout(title_text="Sankey diagram for Scenario 3: Water recovery ", font_size=15)
fig.show()

plot(fig)

fig = go.Figure(data=[go.Sankey(
    node = dict(
      pad = 15,
      thickness = 20,
      line = dict(color = "black", width = 0.5),
      label = [ "NF", "MED", "Cryst", "total"],
      color = "blue"
    ),
    link = dict(
      source = [3, 3,3], 
      target = [0,1,2],
      value = [E_el_all[0], E_el_all[1], E_el_all[2]]
  ))])
color_for_nodes = ["lightsteelblue","darkcyan","maroon", "midnightblue", "midnightblue", "midnightblue", "maroon"]
fig.update_traces(node_color = color_for_nodes)
fig.update_layout(title_text="Sankey diagram for Scenario 3: Water recovery, E_el ", font_size=15)
fig.show()

plot(fig)

# fig = go.Figure(data=[go.Sankey(
#     node = dict(
#       pad = 15,
#       thickness = 20,
#       line = dict(color = "black", width = 0.5),
#       label = [ "NF", "MED", "Cryst", "total"],
#       color = "blue"
#     ),
#     link = dict(
#       source = [3, 3,3], 
#       target = [0,1,2],
#       value = [E_th_all[0], E_th_all[1], E_th_all[2]]
#   ))])
# color_for_nodes = ["lightsteelblue","darkcyan","maroon", "midnightblue", "midnightblue", "midnightblue", "maroon"]
# fig.update_traces(node_color = color_for_nodes)
# fig.update_layout(title_text="Sankey diagram for Scenario 3: Water recovery, E_th ", font_size=15)
# fig.show()

# plot(fig)

#%%export to excel
dfel=pd.DataFrame(E_el_all ,tec_names)
dfeth=pd.DataFrame(E_th_all, tec_names)
dfind=(Qw_tot, abs_Qw, rec, Tot_el, Tot_th, emis_t, Q_w_in,  OPEX, CAPEX, Ec_m, Pr_c)
dfprod=(Qw_tot,M_Nacl, 0, 0, 0, 0, 0)

dfprodn=("Water (kg/hr)", "NaCl (kg/hr)", "Mg(OH)2(kg/hr)", "Ca(OH)2(kg/hr)", "Na2SO4(kg/hr)", "1M NaOH (kg/hr)","1M HCl(kg/hr)")
ind=np.array(["Water production", "absolute water production", "water recovery","Total electrical consumption", "Total thermal energy consumption", "Carbon dioxide emission kg co2/year",
              "Water footprint","OPEX", "CAPEX","economic margin in €/m\u00B3 of water", "Production cost"])
units=pd.DataFrame(["kg/hr", "kg/hr", "%","KWh", "KWh"," kg co2/year","kg/hr","€/year", "€","€/m\u00B3 of water", "euro/m\u00B3 of seawater"])
dfind_t=pd.DataFrame(data=dfind, index=ind)
dfprodn=pd.DataFrame(data=dfprod, index=dfprodn)
with pd.ExcelWriter('results_water_scenario.xlsx') as writer:
    sum_table_f.to_excel(writer,sheet_name="water_scenario")
    sum_table_C.to_excel(writer,sheet_name="water_scenario", startcol=0, startrow=4, header=True)
    dfel.to_excel(writer,sheet_name="water_scenario", startcol=0, startrow=14, header=True)
    dfeth.to_excel(writer,sheet_name="water_scenario", startcol=2, startrow=14, header=True)
    dfprodn.to_excel(writer,sheet_name="water_scenario", startcol=0, startrow=23, header=True)
    #dfprod.to_excel(writer,sheet_name="water_mining_scenario", startcol=1, startrow=24, header=True)
    dfind_t.to_excel(writer,sheet_name="indicators")
    units.to_excel(writer,sheet_name="indicators",startcol=2, startrow=0, header=True)