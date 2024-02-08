# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 15:11:14 2022

@author: rodoulaktori
"""
#from nanofiltration_unit import molarity
#from nanofiltration_unit import nfmass
from nanofiltration_unit_f import Osmotic_pressure
import mfpfr_unit2
import constants
import density_calc
Cpermi=[]
Cconci=[]
Qf=mfpfr_unit2.Mout_2  
print("qf is " + str(Qf))      
Ci_in=mfpfr_unit2.Cout_mfpfr_g
#Na=molarity(constants.MW_Na, 1, Ci_in[0])


mg_in=sum(Ci_in)
print("mg_in is " + str(mg_in))  
d_in=density_calc.density_calc(25, mg_in)
print("d_in is " + str(d_in)) #kg/m3 
     
##(self, comp, Cfeedi, rjr, Wrec, Qf):        
rjr=[0.081847925, 0.05, 0.149762405, 0.88, 0.776183278, 0.96]
Wrec=0.6
Qperm= Wrec*Qf
Qconc=Qf-Qperm
for i in range(len(Ci_in)):
    Cpermi.append((1-rjr[i])*Ci_in[i])
    Cconci.append((Qf*Ci_in[i]-Qperm*Cpermi[i])/Qconc)

#Comp1=nfmass('Na', Ci_in[0], rjr_1[0],Wrec_1, Qf)    
#Comp2=nfmass('Cl', Ci_in[1], rjr_1[1], Wrec_1, Qf) 
#Comp3=nfmass('K', Ci_in[2], rjr_1[2], Wrec_1, Qf) 
#Comp4=nfmass('Mg', Ci_in[3],rjr_1[3], Wrec_1, Qf) 
#Comp5=nfmass('Ca', Ci_in[4], rjr_1[4], Wrec_1, Qf) 
#Comp6=nfmass('SO4', Ci_in[5], rjr_1[5], Wrec_1, Qf) 
#
#
#Comp1.permcal()
#Comp2.permcal()
#Comp3.permcal()
#Comp4.permcal()
#Comp5.permcal()
#Comp6.permcal()
#
##results 
Cconc=Cconci
print("concentration concentrate stream "+str(Cconc))

mg_out=sum(Cconc)
print("mg_in is " + str(mg_out))  
d_out=density_calc.density_calc(25, mg_out)
print("d_in is " + str(d_out)) #kg/m3 
d_p=density_calc.density_calc(25, sum(Cpermi))

Cperm=Cpermi
#print("concentration permeate stream #2st pass  #0 step "+str(Cperm))
##print("                                              ")
#Qperm=Comp1.Qperm
print("Qperm"+str(Qperm))

#Qconc=Comp1.Qconc
print("Qconc is "+str(Qconc))
P_osmo_f=Osmotic_pressure(Ci_in[0], 1,Ci_in[1], -1, Ci_in[2],1, Ci_in[3],2, Ci_in[4], 2, Ci_in[5], -2) #print("osmotic pressure of feed"+str(P_osmo_f.p_osmo))
P_osmo_p=Osmotic_pressure(Cpermi[0],1, Cpermi[1], -1,Cpermi[2],1, Cpermi[3],2, Cpermi[4], 2,Cpermi[5], -2) #print(P_osmo_p.p_osmo)
P_osmo_c=Osmotic_pressure(Cconci[0],1, Cconci[1],-1, Cconci[2],1, Cconci[3],2, Cconci[4], 2,Cconci[5], -2)
#
#mass balances 
for i in range(len(Cconc)):
    bal_i=Ci_in[i]*Qf-(Qperm*Cperm[i]+Qconc*Cconc[i])
    print(" for "+ str(i)+ " mass balance equals: " + str(bal_i))

#electircal consumption 
dp=2
Papplied=(P_osmo_c.p_osmo+P_osmo_f.p_osmo)/2-P_osmo_p.p_osmo+dp #print("applied pump pressure 2 pass " +str(Papplied_2))

Ppump=Papplied*Qperm/1000*1e5/3600 #print("pump pressure 2 pass " +str(Ppump_2)) #W

E_el_nf_efc=Ppump/1000/0.8 #kw

Qconc_m3=Qconc/d_out #m3/hr