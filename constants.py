# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 16:10:12 2022

@author: rodoulaktori
"""
MW_Na=22.99
MW_cl=35.453
MW_so4=32.064+4*15.999
MW_K=39.102
MW_Ca=40.08
MW_Mg=24.31
MW_HCO3=1.008+12.011+3*15.999

#PRICES
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


#%%

#Molecular weight 
MW_Na=22.99
MW_cl=35.453
MW_so4=32.064+4*15.999
MW_K=39.0983
MW_Ca=40.078
MW_Mg=24.305
MW_HCO3=1.008+12.011+3*15.999
MW_MgOH=58.3197
MW_CaOH=74.09
MW_Na2SO4=142.039458
MW_HCl=36.458
MW_NaOH=39.997
MW_NaCl=58.442729
MW_na2so410h20=322.192258

#emissions
CO2_el=0.275 # kg/kwh
CO2_st=0# kg/kwh

#input data
Csw=39.1 #g/l 

#economic
hr=24*300 #hours/year
lf= 20 #years
r=0.06 # interest rate
inf=0.02 # inflation rate 
co2_el_f=0.275 # emission factor electricity 
co2_s_f=0.633 # emission factor steam


#cost from pilot data
eq_c_nf=160000
eq_c_med=220000
eq_c_th_cr=94800
eq_c_mfpfr=80000
eq_c_efc=72579.59
eq_c_edbm= 90000   
eq_c_ev_pd=10000
eq_c_nf_conc=7000
eq_c_ed=123437.5

#capacity of basic scenario
#Mf_basic_sc=nf, med, thermal cryst, mfpfr, nf conc step, efc, edbm, ev_ponds ]
Mf_basic_sc=[2290.858979,1695.830673,199.5094909, 595.0283063,300,100, 229.0522085, 253, 55130.92]