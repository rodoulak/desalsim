# -*- coding: utf-8 -*-
import nanofiltration_unit
import med_unit
import thermal_cryst
import mfpfr_unit2
import edbm_unit_recycling
import constants 
import scaleup
import efc_conc_step
import efc_unit

#%%constants
hr=24*300 #hours/year
lf= 20 #years
r=0.06 # interest rate
inf=0.02 # inflation rate 
co2_el_f=0.275 # emission factor electricity 
co2_s_f=0.633 # emission factor steam

#%%prices 
el_pr=0.253 #euro/kwh
s_pr=0#0.0083 #euro/kwh
co2_em_c=23 #euro/ton co2 
nacl_pr=66/1000 # euro/kg (66euro/ton)

hcl_pr=1040.0/180 #euro/l 1M solution (https://www.sigmaaldrich.com/NL/en/product/mm/100313?gclid=CjwKCAjw9pGjBhB-EiwAa5jl3EARJcEZTu2dTg76NbdyFh3fTORweHxqZex8qeIcxNPsdIFTOlF0sRoCnh4QAvD_BwE&gclsrc=aw.ds)
hcl_0_5_pr=  0#125/1000*0.5
w_pr=0.001 #euro/kg
na2so4_pr=140*0.83/1000 #euro/kg
mgoh2_pr=1000/1000 #euro/kg
caoh2_pr=125/1000#euro/kg
naoh_pr=6840/950#euro/L 1M NaOH solution  (https://www.sigmaaldrich.com/NL/en/search/1m-naoh-solution?focus=products&page=1&perpage=30&sort=relevance&term=1m%20naoh%20solution&type=product_name)
naoh_pr_s=330/1000 #euro/kg
hcl_pr_s=125/1000 #euro/kg
naoh_0_5_pr=66/1000*0.5
antisc_pr=0.002 #euro/ml
cw_pr=0.118/1000#euro/kg
m_salt_pr=5/1000#euro/kg


#density
d_naoh=1.04 #kg/l for 1M solution 
d_hcl=1.0585#kg/l for 1M solution 

#%%
#input data
#REQUIRED AMOUNT OF CHEMICALS 
Qnaoh_ex=edbm_unit_recycling.Q_b_out-mfpfr_unit2.QNAOH
print("edbm_unit_recycling.Q_b_out is "+str(edbm_unit_recycling.Q_b_out))
print("mfpfr_unit2.QNAOH is "+str(mfpfr_unit2.QNAOH))
if Qnaoh_ex<0:
    Qnaoh_s=0
    Qnaoh_need=-Qnaoh_ex
else:
    Qnaoh_s=Qnaoh_ex
    Qnaoh_need=0
print("remaining amount of naoh is "+str(Qnaoh_ex) + "kg/hr")
Qhcl_ex=edbm_unit_recycling.Q_a_out-mfpfr_unit2.QHCl
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
#%%
#symbols:
#hr-> hours
#lf -> lifetime 
#el -> electricity
#s-> steam
#pr -> price 
#c-> cost
#eq-> equipment 
#conc -> consumption 
#wat -> water 
#cw -> cooling water 
#E -> Energy 
#prd -> product mass rate 

#%%constants
class econom:
    def __init__(self, eq_c, el_conc, s_conc, chem1_conc, chem1_pr,chem2_conc, chem2_pr, cw_conc, wat_conc):
        self.eq_c= eq_c
        self.el_conc=el_conc
        self.s_conc=s_conc
        self.chem1_conc=chem1_conc
        self.chem1_pr=chem1_pr
        self.chem2_conc=chem2_conc
        self.chem2_pr=chem2_pr
        self.cw_conc=cw_conc
        self.wat_conc=wat_conc
        
    def capex_calc(self):
        self.inst_c=0.25*self.eq_c  #installation cost 
        self.buil_c=0.2*self.eq_c   #building, process and auxillary cost 
        self.land_c=0.06*self.eq_c   # land cost 
        self.hard_c=self.eq_c+self.inst_c # hardware costs
        self.dir_c= self.hard_c+self.buil_c+ self.land_c #direct costs 
        self.ind_c=0.15*self.dir_c #indirect costs 
        self.fix_c=self.dir_c+self.ind_c #fixed-capital investment 
        self.work_c=0.2*self.fix_c # working capital 
        self.t_capital_inv=self.fix_c+self.work_c # total capital investment 
        #print("total capital investment: " + str(self.t_capital_inv)+ " euro ")
        
    def opex_calc(self):
        self.E_el=self.el_conc*hr*el_pr
        self.E_th=self.s_conc*hr*s_pr
        self.t_E_c=self.E_el+self.E_th
        print("t_E_c "+str(self.t_E_c))
        self.chem_c=self.chem1_conc*hr*self.chem1_pr+self.chem2_conc*hr*self.chem2_pr
        print("chem_c"+str(self.chem_c))
        self.cw_c=self.cw_conc*hr*cw_pr
        self.wat_c=self.wat_conc*hr*w_pr
        self.main_c=0.03*self.fix_c #maintenance cost 
        print("fix cost"+str(self.fix_c))
        print("main_c "+str(self.main_c))
        self.oper_sup_c=0.05*self.main_c # operating suppliers cost 
        print("oper_sup_c "+str(self.oper_sup_c))
        self.oper_lab_c=(self.t_E_c+self.chem_c+self.wat_c+self.cw_c+self.main_c+self.oper_sup_c)/0.675*0.15 #operating labor 
        print("oper_lab_c "+str(self.oper_lab_c))
        self.super_c=0.15*self.oper_lab_c # direct supervisory and clerical labor 
        self.lab_c=0.15*self.oper_lab_c#Laboratory charges
        self.pat_c=(self.t_E_c+self.chem_c+self.wat_c+self.cw_c+self.main_c+self.oper_sup_c)/0.675*0.03     #Patents and royalties
        self.fix_char=(self.t_E_c+self.chem_c+self.wat_c+self.cw_c+self.main_c+self.oper_sup_c)/0.675*0.05  #Fixed charges
        self.over_c=(self.t_E_c+self.chem_c+self.wat_c+self.cw_c+self.main_c+self.oper_sup_c)/0.675*0.05 #plant overhead costs 
        self.opex=self.t_E_c+self.chem_c+self.wat_c+self.cw_c+self.main_c+self.oper_sup_c+self.oper_lab_c+self.super_c+self.lab_c+self.pat_c+self.fix_char+self.over_c
        self.opex2=(self.t_E_c+self.chem_c+self.wat_c+self.cw_c+self.main_c+self.oper_sup_c)/0.675
        #print("total operating cost (opex): " + str(self.opex)+ " euro/year")
    
    def carbon_calc(self):
        self.carbon_c=(self.el_conc*co2_el_f+self.s_conc*co2_s_f)*hr*co2_em_c/1000
        self.E_el=self.el_conc*hr*el_pr
        self.E_th=self.s_conc*hr*s_pr
        self.t_E_c=self.E_el+self.E_th
        self.chem_c=self.chem1_conc*hr*self.chem1_pr+self.chem2_conc*hr*self.chem2_pr
        self.cw_c=self.cw_conc*hr*cw_pr
        self.wat_c=self.wat_conc*hr*w_pr
        self.main_c=0.07*self.fix_c #maintenance cost 
        self.oper_sup_c=0.15*self.main_c # operating suppliers cost 
        self.oper_lab_c=(self.t_E_c+self.chem_c+self.wat_c+self.cw_c+self.main_c+self.oper_sup_c+self.carbon_c)/0.255*0.15 #operating labor 
        self.super_c=0.15*self.oper_lab_c # direct supervisory and clerical labor 
        self.lab_c=0.15*self.oper_lab_c#Laboratory charges
        self.pat_c=(self.t_E_c+self.chem_c+self.wat_c+self.cw_c+self.main_c+self.oper_sup_c)/0.255*0.03     #Patents and royalties
        self.fix_char=self.oper_lab_c  #Fixed charges
        self.over_c=(self.t_E_c+self.chem_c+self.wat_c+self.cw_c+self.main_c+self.oper_sup_c+self.carbon_c)/0.255*0.1 #plant overhead costs 
        self.opex=self.t_E_c+self.chem_c+self.wat_c+self.cw_c+self.main_c+self.oper_sup_c+self.oper_lab_c+self.super_c+self.lab_c+self.pat_c+self.fix_char+self.over_c+self.carbon_c

    
#%%
class revenue:
    def __init__(self, prd, prd_name):
        self.prd=prd
        self.prd_name=prd_name
        
        
    def rev(self):
        if self.prd_name=="water":
            self.rev_prd=self.prd*w_pr*hr #euro/year
        elif self.prd_name=="nacl":
            self.rev_prd=self.prd*nacl_pr*hr #euro/year
        elif self.prd_name=="mgoh2":
            self.rev_prd=self.prd*mgoh2_pr*hr #euro/year
        elif self.prd_name=="na2so4":
            self.rev_prd=self.prd*na2so4_pr*hr #euro/year
        elif self.prd_name=="naoh":
            self.prd=self.prd*d_naoh
            self.rev_prd=self.prd*naoh_pr*hr #euro/year
        elif self.prd_name=="naoh_0.5":
            self.prd=self.prd*(constants.MW_NaOH/2/1000)
            self.rev_prd=self.prd*naoh_pr*hr #euro/year
        elif self.prd_name=="hcl":
            self.prd=self.prd*d_hcl
            self.rev_prd=self.prd*hcl_pr*hr #euro/year
        elif self.prd_name=="hcl_0.5":
            self.prd=self.prd*(constants.MW_HCl/2/1000)
            self.rev_prd=self.prd*hcl_pr*hr #euro/year
#%%Cost calculation 
        #eq_c, el_conc, s_conc, chem_conc, chem_pr, wat_conc
#total_capex_nf=econom(160000, 23, 0, 5, 0.002, 0)
#total_capex_nf.capex_calc()
#total_capex_nf.opex_calc()
#total_capex_nf.carbon_calc()
#        
#total_capex_med=econom(220000, 12, 530.7575, 0, 0, 42088.1547154885)
#total_capex_med.capex_calc()
#total_capex_med.opex_calc()
#total_capex_med.carbon_calc()        
#
#total_capex_cryst=econom(94800, 9.95, 0, 0, 0, 3000.708)
#total_capex_cryst.capex_calc()
#total_capex_cryst.opex_calc()
#total_capex_cryst.carbon_calc()   
#
#total_capex_mfpfr=econom(80000, 9.182285752, 0, 350.5882985, 0.0007029, 0)
#total_capex_mfpfr.capex_calc()
#total_capex_mfpfr.opex_calc()
#total_capex_mfpfr.carbon_calc()  
#
#total_capex_efc=econom(72579.59, 9.229688143, 0, 0, 0, 0)
#total_capex_efc.capex_calc()
#total_capex_efc.opex_calc()
#total_capex_efc.carbon_calc()  
#
#total_capex_edbm=econom(90000, 3.90838595, 0, 0, 0, 453.3423653) #change chemical consumption and cost 
#total_capex_edbm.capex_calc()
#total_capex_edbm.opex_calc()
#total_capex_edbm.carbon_calc()    
#        
#CAPEX=total_capex_nf.t_capital_inv+total_capex_med.t_capital_inv+total_capex_cryst.t_capital_inv+total_capex_mfpfr.t_capital_inv+total_capex_efc.t_capital_inv+total_capex_edbm.t_capital_inv
#print("total capex of system: " + str(round(CAPEX))+ " euro")        
#OPEX=total_capex_nf.opex+total_capex_med.opex+total_capex_cryst.opex+total_capex_mfpfr.opex+total_capex_efc.opex+total_capex_edbm.opex
#print("total OPEX of system: " + str(round(OPEX))+ " euro/year")     
#CO2_c=total_capex_nf.carbon_c+total_capex_med.carbon_c+total_capex_cryst.carbon_c+total_capex_mfpfr.carbon_c+total_capex_efc.carbon_c+total_capex_edbm.carbon_c
#print("total carbon tax cost: " + str(round(CO2_c))+ " euro/year")     

#%%cost calculation vs2
CAPEX=0
OPEX=0
OPEX2=0
CO2_c=0
OPEX_with_ext=0
Mf_basic_sc=[]
opex_list=[]
capex_list=[]
tec_names=["NF", "MED", "CRY", "MFPFR", "NF-2", "EFC", "EDBM"]
eq_c=[constants.eq_c_nf, constants.eq_c_med, constants.eq_c_th_cr, constants.eq_c_mfpfr, constants.eq_c_nf_conc,constants.eq_c_efc, constants.eq_c_edbm]
el_conc=[nanofiltration_unit.E_el_nf,med_unit.E_el_med,thermal_cryst.E_el_th_Cr,mfpfr_unit2.E_el_mfpf,efc_conc_step.E_el_nf_efc, efc_unit.E_t,edbm_unit_recycling.P_t] 
s_conc=[0,med_unit.E_th_med,thermal_cryst.E_th_th_Cr,0,0,0,0 ]
chem1_conc=[nanofiltration_unit.Qantsc,0,0,Mnaoh_need,0,0,0]
chem1_pr=[0.002,0,0,naoh_pr_s,0,0,0]
chem2_conc=[0,0,0,Mhcl_need,0,0,0]
chem2_pr=[0.002,0,0,hcl_pr_s,0,0,0]
wat_conc=[0,0, 0, Qnaoh_need,0,0, edbm_unit_recycling.Q_w_in ]
cw_conc=[0,med_unit.Qcw, thermal_cryst.Qcw, 0,0,0, 0 ]
Mf_basic_sc.append(constants.Mf_basic_sc[0])
Mf_basic_sc.extend(constants.Mf_basic_sc)
Mf_sce=[nanofiltration_unit.Qsw, med_unit.Qf_med, thermal_cryst.Qf, mfpfr_unit2.Qin_mfpfr,mfpfr_unit2.Qout_2, efc_conc_step.Qconc, efc_unit.Mtot]
for i in range(len(eq_c)):
    if Mf_basic_sc[i]!=Mf_sce[i]:
        eq_c[i]=scaleup.scaleup_eq(eq_c[i],Mf_basic_sc[i],Mf_sce[i],tec_names[i])
for i in range(len(eq_c)):
    total_econom=econom(eq_c[i], el_conc[i], s_conc[i], chem1_conc[i], chem1_pr[i],chem2_conc[i], chem2_pr[i], cw_conc[i], wat_conc[i])
    total_econom.capex_calc()
    total_econom.opex_calc()
    CAPEX=CAPEX+total_econom.t_capital_inv
    OPEX=OPEX+total_econom.opex
    OPEX2=OPEX2+total_econom.opex2
    print("opex is "+str(total_econom.opex))
    print("opex 2 is "+str(total_econom.opex2))
    opex_list.append(total_econom.opex)
    total_econom.carbon_calc() 
    CO2_c=CO2_c+total_econom.carbon_c
    OPEX_with_ext=OPEX_with_ext+total_econom.opex
    capex_list.append(total_econom.t_capital_inv)
    
CAPEX=CAPEX.item()
OPEX=OPEX[0][0] 
# print("opex is "+str(OPEX)) 
# print("total capex of system: " + str(round(CAPEX))+ " euro")
# print("total OPEX without externalities of system: " + str(round(OPEX))+ " euro/year") 
# print("total OPEX with externalities of system: " + str(round(OPEX_with_ext))+ " euro/year") 
# print("total carbon tax cost: " + str(round(CO2_c))+ " euro/year")     

#amortisation factor 
a=(r*(1+r)**lf)/((1+r)**lf-1)
print("a is " +str(a))
#%%Revenue 
#REQUIRED AMOUNT OF CHEMICALS 
Qnaoh_ex=edbm_unit_recycling.Q_b_out-mfpfr_unit2.QNAOH
if Qnaoh_ex<0:
    Qnaoh_s=0
    Qnaoh_need=-Qnaoh_ex
else:
    Qnaoh_s=Qnaoh_ex
    Qnaoh_need=0
print("remaining amount of naoh is "+str(Qnaoh_ex) + "kg/hr")
Qhcl_ex=edbm_unit_recycling.Q_a_out-mfpfr_unit2.QHCl
print("remaining amount of hcl is "+str(Qhcl_ex) + "kg/hr")
if Qhcl_ex<0:
    Qhcl_s=0
    Qhcl_need=-Qhcl_ex
else:
    Qhcl_s=Qhcl_ex
    Qhcl_need=0
reve_t=0
reve_list=[]
prd=[med_unit.Qprod_med+thermal_cryst.Q_evap_mass+efc_unit.M_ice, thermal_cryst.M_Nacl, mfpfr_unit2.M_MgOH2_1, efc_unit.M_compound1_cr, Qnaoh_s*constants.MW_NaOH/1000, Qhcl_s*constants.MW_HCl/1000]    
prd_name= ["water", "nacl", "mgoh2", "na2so4", "naoh", "hcl"]   
for i in range(len(prd)):
    rev_calc=revenue(prd[i], prd_name[i])    
    rev_calc.rev()
    print("rev_calc.rev_prd for " + prd_name[i]+" is " + str(rev_calc.rev_prd))
    reve_t = reve_t+rev_calc.rev_prd
    reve_list.append(rev_calc.rev_prd)