import matplotlib.pyplot as plt
import mfpfr_unit2 
import nanofiltration_unit
import med_unit 
import thermal_cryst
import economic 
import edbm_unit_recycling
#from tabulate import tabulate

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

hours=8

#Molecular weight 
MW_Na=22.99
MW_cl=35.453
MW_so4=32.064+4*15.999
MW_K=39.102
MW_Ca=40.08
MW_Mg=24.31
MW_HCO3=1.008+12.011+3*15.999
MW_MgOH=58.3197
MW_CaOH=74.09
MW_Na2SO4=142.039458
MW_HCl=36.458
MW_NaOH=39.997
MW_NaCl=58.442729
 
#emissions
CO2_el=0.275 # kg/kwh
CO2_st=0.255 # kg/kwh

#input data
Csw=39.1 #g/l 


#%%Indicators:
#----------------------------------------------------------------------------------------------------------------------------------
#Technological indicators  

#Quantity of water produced: Kg/hr or t/y
#Quality of water produced (purity of the products): 	%
#energy consumption electrical and thermal(Quantifies the efficiency of transforming energy ):	kWh
#Product water recovery ratio (Water_rr):		%
#Total brine reduction: 		t/yr or %
#Minimization of waste (Reduction of waste (brine) by product or by emissions): 	Kg/m3 of water and/or kg/CO2 emissions 
#----------------------------------------------------------------------------------------------------------------------------------

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
            
       
#----------------------------------------------------------------------------------------------------------------------------------  
#Carbon dioxide emission	Based on the energy consumption 	kg CO2-Equ/m3 product water
#Resource utilisation (resource usage): 	Water usage	l/use
#                                       	Water loss	l/use
#	                                      Energy use	KWh/use
#	                                      Chemical use 	l/use
#Brine discharges		m3/m3 of seawater or m3/m3 of fresh water
#Energy consumption 		KWh/m3 of seawater
#Recycling/reuse rate 	Number (or volume, or mass) fraction of resource recovery involving collection and treatment of a waste product for use as a raw material	%
#Resource efficiency 	Number (or volume, or mass) ratio of useful material output and material input	%
#Energy intensity 	Total energy used per units produced	KWh/ m3 of fresh water or KWh/kg of salt (for each salt)
            
#----------------------------------------------------------------------------------------------------------------------------------            

        
#class energy_indicator:
#    def __init__(self, El_el_nf, El_el_med, El_el_therm, El_el_mfpf, El_el_efc, El_el_edbm, Qwater_tot):
#        self.El_el_nf=El_el_nf
#        self.El_el_med=El_el_med
#        self.El_el_therm=El_el_therm
#        self.El_el_mfpf=El_el_mfpf
#        self.El_el_efc=El_el_efc
#        self.El_el_edbm=El_el_edbm
#    
#    def ene_ind(self):
#        self.E_el_tot=self.El_el_nf+self.El_el_med+ self.El_el_therm+ self.El_el_mfpf+ self.El_el_efc+ self.El_el_edbm
#        #self.E_th_tot=E_th
#        self.E_int=self.E_el_tot/self.Qwater_tot #energy intensity for water 
#        
    def env_ind(self):
        self.emis=(self.E_el*CO2_el+self.E_th*CO2_st)/self.Qprod_1*1000 #Carbon dioxide emission: kg CO2-Equ/m3 product water
        #self.emis_t=(self.E_el_tot*CO2_el+self.E_th_tot*CO2_st)/self.Qprod_1*1000
        self.brine_disc=self.Qout/self.Qin
#%%main 
#Scenario 1 (WATER MINING - CASE study)
##treatment train: NF-> MED-> THERMAL CRYST
                  #->MFPFR->EFC->EDBM 
#Call units' functions 
#def __init__(self, tech, Qout, Qin, Qprod_1, prod_1_name, Qprod_2, prod_2_name, E_el, E_th, Cin, Cout)       
#print(E_el_nf)
tec1=indic("NF", nanofiltration_unit.Qconc, nanofiltration_unit.Qsw, nanofiltration_unit.Qperm, "none", 0, "none", nanofiltration_unit.E_el_nf, 0, nanofiltration_unit.Cin_conc, nanofiltration_unit.Cout_conc, nanofiltration_unit.Qhcl,nanofiltration_unit.Qantsc)     
tec1.techn_indi()
tec1.env_ind()
print("The CO2 emissions based on the energy consumption of "+ tec1.tech +" are "+str(round(tec1.emis))+" kgCO2/m3 of product")
#
tec2=indic("MED", med_unit.Qout_med, med_unit.Qf_med, med_unit.Qprod_med, "water", 0, "none",  med_unit.E_el_med, 300, med_unit.Cin_med, med_unit.Cconc_med, 0,0)
tec2.techn_indi()
tec2.env_ind()
#print("The specific energy consumption of "+ tec2.tech +"is "+str(tec2.SEC)+" kwh/ton")
print("The CO2 emissions based on the energy consumption of "+ tec2.tech +" are "+str(round(tec2.emis))+" kgCO2/m3 of product")

tec3=indic("Thermal Cryst", 226.6711827, 377.7853045, 17.9520953,"water", 133.1620265, "NaCl", 9.229688143, 0, 1, 2,0,0 )
tec3.techn_indi()
tec3.env_ind()
##
tec4=indic("MF-PFR", 226.6711827, 377.7853045, mfpfr_unit2.M_MgOH2_1,"Mg(OH)2", 133.1620265, "Ca(OH)2", mfpfr_unit2.E_el_mfpf, 0, 1, 2, mfpfr_unit2.QNAOH, mfpfr_unit2.QHCl)
tec4.techn_indi()
tec4.env_ind()
#print("The specific energy consumption of "+ tec4.tech +"is "+str(tec4.SEC)+" kwh/ton")
##
tec5=indic("EFC", 226.6711827, 377.7853045, 17.9520953,"Na2SO4", 133.1620265, "water", 9.229688143, 0, 1, 2,0,0 )
tec5.techn_indi()
tec5.env_ind()
##
tec6=indic("EDBM", 226.6711827, 377.7853045, 17.9520953,"NaOH", 133.1620265, "HCl", 9.229688143, 0, 1, 2,0,0 )
tec6.techn_indi()
tec6.env_ind()
#
##list results 
tec_names=[tec1.tech, tec2.tech, tec3.tech, tec4.tech, tec5.tech, tec6.tech]
E_el_all=[tec1.E_el, tec2.E_el, tec3.E_el, tec4.E_el, tec5.E_el, tec6.E_el]
E_th_all=[tec1.E_th, tec2.E_th, tec3.E_th, tec4.E_th, tec5.E_th, tec6.E_th]
Cout_all=[tec1.Cout, tec2.Cout, tec3.Cout, tec4.Cout, tec5.Cout, tec6.Cout]
Qout_all=[tec1.Qout, tec2.Qout, tec3.Qout, tec4.Qout, tec5.Qout, tec6.Qout]
Qchem_all=[tec1.chem,tec2.chem, tec3.chem,tec4.chem, tec5.chem,tec6.chem]
#%%
#Technical 
#T-1) quantity of water production 
Qw_tot=tec1.Qwater+tec2.Qwater+ tec3.Qwater+ tec4.Qwater+ tec5.Qwater+ tec6.Qwater
#print("total water production is " + str(Qw_tot))

#T-2)energy performance
#T-2a) specific electrical energy consumption #kwh/kg of desalinated water 
SEC=sum(E_el_all)/Qw_tot
print("specific electrical consumption is " +str(SEC)+ " Kwh/kg of desalinated water")
#2b) share of renewable energy %

#T-3)efficiency 
#brine production -> kg of brine/ kg of seawater 
Qbrine_f=0
for i in range(len(tec_names)):
    if tec_names[i]=="EDBM":
        if Cout_all[i]>Csw:
            Qbrine_f=Qbrine_f+Qout_all[i]
        else:
            Qbrine_f=Qbrine_f+0
Q_br_prod=Qbrine_f/nanofiltration_unit.Qsw
print("brine production is "+str(Q_br_prod)+" kg of brine/kg of seawater")
#------------------------------------------------------------------------------
#%%
#Environmental 
#Env-1)Climate change impact ->Carbon dioxide emission kg co2/kg desalinated water 
emis_t=((sum(E_el_all)*CO2_el)+(sum(E_th_all)*CO2_st))/Qw_tot
print("Total specific Carbon dioxide emission "+str(emis_t)+" kg co2/kg desalinated water")
#Env-2) Resource utilization -> chemical consumption in kg of chemical/kg of seawater 
Chem_cons=(sum(Qchem_all))/nanofiltration_unit.Qsw #check units and convert ml/l to kg with densities 
print("Total chemical consumption is "+str(Chem_cons)+" kg of chemical/kg of seawater ")
#Env-3) The aquatic eco-toxic impact of brine disposal -> ECO-TOXICITY 
C_brine=edbm_unit_recycling.Cbrine_out #[ CNa, CCl, CK, CMg, CCa, CSO4, Chco3, CH, COH]
ETP=0
CF=4.62e-1
for i in range(len(C_brine)-2):
    ETP=ETP+C_brine[i]*CF
print("The aquatic eco-toxic impact of brine disposal is "+str(ETP)+"")
#Env-4) Environmental prodection -> minimization of waste in kg of brine/CO2 emissions 
Min_waste=(med_unit.Qout_med+nanofiltration_unit.Qconc-Qbrine_f)/((sum(E_el_all)*CO2_el)+(sum(E_th_all)*CO2_st))
print("Emissions from minimization of waste is "+str(Min_waste)+" kg of brine/CO2 emissions ")
#Env-5) Resource recovery efficiency -> Resource efficiency (%)

print("Resource recovery efficiency is "+str()+" ")
print(""+str()+"")


#%%
#Economic 
#EC-1) Product value -> levelized cost 
print(""+str()+"")
#EC-2)impact of externalities on the production cost 
#EC-2a)production cost with externalities in euro/m3 of seawater 
Pr_c_ext=economic.OPEX_with_ext/Qw_tot
print("production cost with externalities is "+str()+" euro/m3 of seawater")
#EC-2b)production cost without externalities in euro/m3 of seawater
Pr_c=economic.OPEX/Qw_tot
print("production cost without externalities is "+str(Pr_c)+" euro/m3 of seawater")
#EC-3) profitability -> production efficiency in euro/euro 
prd_eff=(economic.reve_t)/(economic.CAPEX*economic.a+economic.OPEX)
print("production efficiency is "+str(prd_eff)+"")
print(""+str()+"")
print(""+str()+"")
#%%present results
#print(tabulate(E_el_all, headers=tec_names, tablefmt="grid", showindex="always"))
# Plot bar chart with data points
#plt.bar(tec_names, SEC_all)
## Display the plot
#plt.ylabel('Specific energy consumption (kwh/m3)')
#plt.xlabel('Technology')
#plt.show()

#%%scenario 2 : RO BRINE AS FEED 
##treatment train: NF-> MED-> THERMAL CRYST
                  #->MFPFR->EFC->EDBM 

#%%Scenario 3: water recovery
##treatment train: NF-> MED-> THERMAL CRYST     
                  #nf concentrate + med brine to thermal cryst 
