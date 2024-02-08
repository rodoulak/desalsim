import math
import numpy as np
from density_calc import density_calc
import scaleup


#%%constants
R=8.31446261815324 #gas constant
T=20+273.15 #K
MW_Na=22.99
MW_cl=35.453
MW_so4=32.064+4*15.999
MW_K=39.102
MW_Ca=40.08
MW_Mg=24.31
MW_HCO3=1.008+12.011+3*15.999
#assumptions 
dp=2

#molarity 
class molarity:
    def __init__(self, MW, zi, Ci):
        self.MW=MW
        self.zi=zi
        self.Ci=Ci
        self.meq=self.Ci*1000*self.zi/self.MW
        #self.mi=self.Ci/(self.MW*(1000-self.sum_conc)/1000                    #CHECK AGAIN THE UNITS 
        
#%%
class nfmass:
    def __init__(self, comp, Cfeedi, rjr, Wrec, Qf):
        self.comp=comp
        self.Cfeedi=Cfeedi #ion concentration g/l
        self.rjr=rjr #rejection of each ion of the membrane 
        self.Wrec = Wrec # water recovery first pass 
        self.Qf=Qf
        
                       
    def permcal(self):
        self.Qperm= self.Wrec*self.Qf
        self.Qconc=Qf-self.Qperm
        self.Cpermi = (1-self.rjr)*self.Cfeedi
        self.Cconci = (Qf*self.Cfeedi-self.Qperm*self.Cpermi)/self.Qconc
        print("Qconc is "+str(self.Qconc))
        print("self.Qperm is "+str(self.Qperm))
#        print(" ion concenctration of " + self.comp+ " in permeate stream is"+" "+ str(self.Cpermi))
#        print(" ion concenctration of " + self.comp+ " in concentrate stream is"+" "+str(self.Cconci))

#%%Energy consumption 
class NfEnergy:   
        pass 
#%% osmotic pressure 
class Osmotic_pressure:
    
        def __init__(self, C1, z1,  C2, z2,  C3, z3,  C4, z4, C5, z5, C6, z6):
            self.C1=C1
            self.z1=z1
            self.C2=C2
            self.z2=z2
            self.C3=C3
            self.z3=z3
            self.C4=C4
            self.z4=z4
            self.C5=C5
            self.z5=z5
            self.C6=C6
            self.z6=z6
            self.sum_Ci=sum([C1, C2, C3, C4, C5, C6])
            self.zi_2=[z1**2, z2**2, z3**2, z4**2, z5**2, z6**2]
            self.conc=[C1/MW_Na, C2/MW_cl, C3/MW_K, C4/MW_Mg, C5/MW_Ca, C6/MW_so4] #mol/l
            self.mi=[C1*1000/(MW_Na*1000*((1e+6-self.sum_Ci*1000)/1e+6)), C2*1000/(MW_cl*1000*((1e+6-self.sum_Ci*1000)/1e+6)), C3*1000/(MW_K*1000*((1e+6-self.sum_Ci*1000)/1e+6)), C4*1000/(MW_Mg*1000*((1e+6-self.sum_Ci*1000)/1e+6)), C5*1000/(MW_Ca*1000*((1e+6-self.sum_Ci*1000)/1e+6)), C6*1000/(MW_so4*1000*((1e+6-self.sum_Ci*1000)/1e+6))]
            self.mizi_2=[]
            for i in range(6) :
                self.mizi_2.append(self.mi[i]*self.zi_2[i])
            self.SI=sum(self.mizi_2)/2
            self.B=-348.662/(T)+6.72817-0.971307*math.log(T)
            self.C=40.5016/(T)-0.721404+0.103915*math.log(T)
            self.D=5321/(T)+233.76-0.9297*(T)+0.001417*(T)**2-0.0000008292*(T)**3
            self.S=1.17202*(sum(self.mizi_2)/sum(self.mi))*0.9982**0.5*(23375.556/(self.D*(T)))**1.5
            self.fi=1-self.S/(3.375*self.SI)*((1+1.5*self.SI**0.5)-2*math.log(1+1.5*self.SI**0.5)-1/(1+1.5*self.SI**0.5))+self.B*sum(self.mi)/2+self.C*(sum(self.mi)/2)**2
            self.sum_conc=sum(self.conc)
            self.p_osmo_psig=1.205*self.fi*(T)*sum(self.mi)
            self.p_osmo=self.p_osmo_psig/14.3 #bar
           
#check units 

#%% 
#Input data    
    
Na=molarity(MW_Na, 1, 11.9)

Ci_in=[12.33, 21.67, 0.45, 1.39, 0.45, 3.28]   
mg_in=sum(Ci_in)
print("mg_in is " + str(mg_in))  
d_in=density_calc(25, mg_in)  #kg/m3 
print("d_in is " + str(d_in)) #kg/m3 
Qsw = 3000/24*d_in   #3000m3/d #149612.002308333  kg/hr  #2290.858979 #m3/hr
Qf=Qsw   #kg/hr
rjr_1=[0.16, 0.29, 0.21, 0.98, 0.95,0.98]
rjr_2=[0.12, 0.16, 0.09, 0.94, 0.85,0.95]
Wrec_1=0.7
Wrec_2=0.9

#1st pass  #0 step           
Comp1=nfmass('Na', Ci_in[0], rjr_1[0],Wrec_1, Qf)    
Comp2=nfmass('Cl', Ci_in[1], rjr_1[1], Wrec_1, Qf) 
Comp3=nfmass('K', Ci_in[2], rjr_1[2], Wrec_1, Qf) 
Comp4=nfmass('Mg', Ci_in[3],rjr_1[3], Wrec_1, Qf) 
Comp5=nfmass('Ca', Ci_in[4], rjr_1[4], Wrec_1, Qf) 
Comp6=nfmass('SO4', Ci_in[5], rjr_1[5], Wrec_1, Qf) 


Comp1.permcal()
Comp2.permcal()
Comp3.permcal()
Comp4.permcal()
Comp5.permcal()
Comp6.permcal()
Cconc=[Comp1.Cconci, Comp2.Cconci, Comp3.Cconci, Comp4.Cconci, Comp5.Cconci, Comp6.Cconci]
#print("concentration concentrate stream #2st pass  #0 step "+str(Cconc))
Cin_conc=(Comp1.Cfeedi+ Comp2.Cfeedi+Comp3.Cfeedi+Comp4.Cfeedi+Comp5.Cfeedi+Comp6.Cfeedi)
print("Cin is " + str(Cin_conc))
Cperm=[Comp1.Cpermi, Comp2.Cpermi, Comp3.Cpermi, Comp4.Cpermi, Comp5.Cpermi, Comp6.Cpermi]
#print("concentration permeate stream #2st pass  #0 step "+str(Cperm))
#print("                                              ")

P_osmo_f=Osmotic_pressure(Comp1.Cfeedi, 1,Comp2.Cfeedi, -1, Comp3.Cfeedi,1, Comp4.Cfeedi,2, Comp5.Cfeedi, 2, Comp6.Cfeedi, -2) #print("osmotic pressure of feed"+str(P_osmo_f.p_osmo))
P_osmo_p=Osmotic_pressure(Comp1.Cpermi,1, Comp2.Cpermi, -1,Comp3.Cpermi,1, Comp4.Cpermi,2, Comp5.Cpermi, 2,Comp6.Cpermi, -2) #print(P_osmo_p.p_osmo)
P_osmo_c=Osmotic_pressure(Comp1.Cconci,1, Comp2.Cconci,-1, Comp3.Cconci,1, Comp4.Cconci,2, Comp5.Cconci, 2,Comp6.Cconci, -2)

#print("Permeate flow rate: "+ str(round(Comp1.Qperm,1))+" m3/h")   
#print("Concentrate flow rate: "+ str(round(Comp1.Qconc,1))+" m3/h")   
#
#print("                                              ")
#second pass #0 step

Qf=Comp1.Qperm

Comp1=nfmass('Na', Comp1.Cpermi, rjr_2[0], Wrec_2, Qf)    
Comp2=nfmass('Cl', Comp2.Cpermi,rjr_2[1], Wrec_2, Qf) 
Comp3=nfmass('K', Comp3.Cpermi, rjr_2[2], Wrec_2, Qf) 
Comp4=nfmass('Mg', Comp4.Cpermi, rjr_2[3], Wrec_2, Qf) 
Comp5=nfmass('Ca', Comp5.Cpermi, rjr_2[4], Wrec_2, Qf) 
Comp6=nfmass('SO4', Comp6.Cpermi, rjr_2[5], Wrec_2, Qf) 

Comp1.permcal()
Comp2.permcal()
Comp3.permcal()
Comp4.permcal()
Comp5.permcal()
Comp6.permcal()
Cconc=[Comp1.Cconci, Comp2.Cconci, Comp3.Cconci, Comp4.Cconci, Comp5.Cconci, Comp6.Cconci]
#print("concentration concentrate stream #2st pass  #0 step "+str(Cconc))

Cperm=[Comp1.Cpermi, Comp2.Cpermi, Comp3.Cpermi, Comp4.Cpermi, Comp5.Cpermi, Comp6.Cpermi]
#print("concentration permeate stream #2st pass  #0 step "+str(Cperm))
#print("                                              ")

P_osmo_f=Osmotic_pressure(Comp1.Cfeedi, 1,Comp2.Cfeedi, -1, Comp3.Cfeedi,1, Comp4.Cfeedi,2, Comp5.Cfeedi, 2, Comp6.Cfeedi, -2) #print(P_osmo_f.p_osmo)
P_osmo_p=Osmotic_pressure(Comp1.Cpermi,1, Comp2.Cpermi, -1,Comp3.Cpermi,1, Comp4.Cpermi,2, Comp5.Cpermi, 2,Comp6.Cpermi, -2) #print(P_osmo_p.p_osmo)
P_osmo_c=Osmotic_pressure(Comp1.Cconci,1, Comp2.Cconci,-1, Comp3.Cconci,1, Comp4.Cconci,2, Comp5.Cconci, 2,Comp6.Cconci, -2)

#print("Permeate flow rate: "+ str(round(Comp1.Qperm,1))+" m3/h")   
#print("Concentrate flow rate: "+ str(round(Comp1.Qconc,1))+" m3/h")   
#print("end of 0 step")
#print("                                              ")

Qconc_to_st2=Comp1.Qconc

##---------------------------------------------------------------
#1st pass  #1 step 
Qf=Qsw+Qconc_to_st2
#print("qconc: "+str(Comp1.Qconc))
#print(Qf)
Comp1=nfmass('Na', 11.9*Qsw/Qf+Comp1.Cconci*Qconc_to_st2/Qf, rjr_1[0],Wrec_1, Qf)    
Comp2=nfmass('Cl', 21.8*Qsw/Qf+Comp2.Cconci*Qconc_to_st2/Qf, rjr_1[1],Wrec_1, Qf) 
Comp3=nfmass('K', 0.400*Qsw/Qf+Comp3.Cconci*Qconc_to_st2/Qf, rjr_1[2],Wrec_1, Qf) 
Comp4=nfmass('Mg', 1.400*Qsw/Qf+Comp4.Cconci*Qconc_to_st2/Qf, rjr_1[3],Wrec_1, Qf) 
Comp5=nfmass('Ca', 0.400*Qsw/Qf+Comp5.Cconci*Qconc_to_st2/Qf, rjr_1[4],Wrec_1, Qf) 
Comp6=nfmass('SO4', 3.200*Qsw/Qf+Comp6.Cconci*Qconc_to_st2/Qf, rjr_1[5],Wrec_1, Qf) 
Comp1.permcal()
Comp2.permcal()
Comp3.permcal()
Comp4.permcal()
Comp5.permcal()
Comp6.permcal()

#second pass #1 step
qc_out=Comp1.Qconc
Cconc_out=[Comp1.Cconci, Comp2.Cconci, Comp3.Cconci, Comp4.Cconci, Comp5.Cconci, Comp6.Cconci]
Qf=Comp1.Qperm

Comp1=nfmass('Na', Comp1.Cpermi, rjr_2[0], Wrec_2, Qf)    
Comp2=nfmass('Cl', Comp2.Cpermi, rjr_2[1], Wrec_2, Qf) 
Comp3=nfmass('K', Comp3.Cpermi, rjr_2[2], Wrec_2, Qf) 
Comp4=nfmass('Mg', Comp4.Cpermi, rjr_2[3], Wrec_2, Qf) 
Comp5=nfmass('Ca', Comp5.Cpermi, rjr_2[4], Wrec_2, Qf) 
Comp6=nfmass('SO4', Comp6.Cpermi, rjr_2[5], Wrec_2, Qf) 

Comp1.permcal()
Comp2.permcal()
Comp3.permcal()
Comp4.permcal()
Comp5.permcal()
Comp6.permcal()
#print("Permeate flow rate: "+ str(Comp1.Qperm)+" m3/h")  
#print("Permeate flow rate: "+ str(round(Comp1.Qperm,1))+" m3/h")   
#print("Concentrate flow rate: "+ str(round(Comp1.Qconc,1))+" m3/h") 

##total energy
#nf_E_el=(Ppump_1+Ppump_2)/Qsw
Qconc_to_st2=Comp1.Qconc
#
##-------------------------------------------------------------------------------------------------------
##1st pass  # step 
Qf=Qsw+Qconc_to_st2
#print("qconc: "+str(Comp1.Qconc))
#print(Qf)
Comp1=nfmass('Na', 11.9*Qsw/Qf+Comp1.Cconci*Qconc_to_st2/Qf, rjr_1[0],Wrec_1, Qf)    
Comp2=nfmass('Cl', 21.8*Qsw/Qf+Comp2.Cconci*Qconc_to_st2/Qf, rjr_1[1],Wrec_1, Qf) 
Comp3=nfmass('K', 0.400*Qsw/Qf+Comp3.Cconci*Qconc_to_st2/Qf, rjr_1[2],Wrec_1, Qf) 
Comp4=nfmass('Mg', 1.400*Qsw/Qf+Comp4.Cconci*Qconc_to_st2/Qf, rjr_1[3],Wrec_1, Qf) 
Comp5=nfmass('Ca', 0.400*Qsw/Qf+Comp5.Cconci*Qconc_to_st2/Qf, rjr_1[4],Wrec_1, Qf) 
Comp6=nfmass('SO4', 3.200*Qsw/Qf+Comp6.Cconci*Qconc_to_st2/Qf, rjr_1[5],Wrec_1, Qf) 
Comp1.permcal()
Comp2.permcal()
Comp3.permcal()
Comp4.permcal()
Comp5.permcal()
Comp6.permcal()

#print("Permeate flow rate: "+ str(Comp1.Qperm)+" m3/h")  
#print("Permeate flow rate: "+ str(round(Comp1.Qperm,1))+" m3/h")   
#print("Concentrate flow rate: "+ str(round(Comp1.Qconc,1))+" m3/h")   
#print("                                                              ")
#
#applied pressure calculations 
P_osmo_f=Osmotic_pressure(Comp1.Cfeedi, 1,Comp2.Cfeedi, -1, Comp3.Cfeedi,1, Comp4.Cfeedi,2, Comp5.Cfeedi, 2, Comp6.Cfeedi, -2) #print(P_osmo_f.p_osmo)
P_osmo_p=Osmotic_pressure(Comp1.Cpermi,1, Comp2.Cpermi, -1,Comp3.Cpermi,1, Comp4.Cpermi,2, Comp5.Cpermi, 2,Comp6.Cpermi, -2) #print(P_osmo_p.p_osmo)
P_osmo_c=Osmotic_pressure(Comp1.Cconci,1, Comp2.Cconci,-1, Comp3.Cconci,1, Comp4.Cconci,2, Comp5.Cconci, 2,Comp6.Cconci, -2) #print(P_osmo_c.p_osmo)

Ci_out_p=[Comp1.Cpermi, Comp2.Cpermi,Comp3.Cpermi, Comp4.Cpermi, Comp5.Cpermi, Comp6.Cpermi]
d_p=density_calc(25, sum(Ci_out_p))
Papplied_1=(P_osmo_c.p_osmo+P_osmo_f.p_osmo)/2-P_osmo_p.p_osmo+dp
Ppump_1=Papplied_1*Comp1.Qperm/d_p* 1e5/3600#*1e5

#sum results 
Cconc=[Comp1.Cconci, Comp2.Cconci, Comp3.Cconci, Comp4.Cconci, Comp5.Cconci, Comp6.Cconci]
Qconc_1=Comp1.Qconc #print("concentration concentrate stream #1st pass  #2 step "+str(Cconc))
print("Qconc_1 is "+ str(Qconc_1))

Cperm=[Comp1.Cpermi, Comp2.Cpermi, Comp3.Cpermi, Comp4.Cpermi, Comp5.Cpermi, Comp6.Cpermi] #print("concentration permeate stream #1st pass  #2 step "+str(Cperm))

#second pass #2 step

Qf=Comp1.Qperm

Comp1=nfmass('Na', Comp1.Cpermi, rjr_2[0], Wrec_2, Qf)    
Comp2=nfmass('Cl', Comp2.Cpermi, rjr_2[1], Wrec_2, Qf) 
Comp3=nfmass('K', Comp3.Cpermi, rjr_2[2], Wrec_2, Qf) 
Comp4=nfmass('Mg', Comp4.Cpermi, rjr_2[3], Wrec_2, Qf) 
Comp5=nfmass('Ca', Comp5.Cpermi, rjr_2[4], Wrec_2, Qf) 
Comp6=nfmass('SO4', Comp6.Cpermi, rjr_2[5], Wrec_2, Qf) 

Comp1.permcal()
Comp2.permcal()
Comp3.permcal()
Comp4.permcal()
Comp5.permcal()
Comp6.permcal()
#print("Permeate flow rate: "+ str(Comp1.Qperm)+" m3/h")  
#print("Permeate flow rate: "+ str(round(Comp1.Qperm,1))+" m3/h")   
#print("Concentrate flow rate: "+ str(round(Comp1.Qconc,1))+" m3/h") 
#

##applied pressure calculations 
P_osmo_f=Osmotic_pressure(Comp1.Cfeedi, 1,Comp2.Cfeedi, -1, Comp3.Cfeedi,1, Comp4.Cfeedi,2, Comp5.Cfeedi, 2, Comp6.Cfeedi, -2) #print(P_osmo_f.p_osmo)
P_osmo_p=Osmotic_pressure(Comp1.Cpermi,1, Comp2.Cpermi, -1,Comp3.Cpermi,1, Comp4.Cpermi,2, Comp5.Cpermi, 2,Comp6.Cpermi, -2) #print(P_osmo_p.p_osmo)
P_osmo_c=Osmotic_pressure(Comp1.Cconci,1, Comp2.Cconci,-1, Comp3.Cconci,1, Comp4.Cconci,2, Comp5.Cconci, 2,Comp6.Cconci, -2) #print(P_osmo_c.p_osmo)


Cout_perm=(Comp1.Cpermi+ Comp2.Cpermi+Comp3.Cpermi+Comp4.Cpermi+Comp5.Cpermi+Comp6.Cpermi)
print("Cout_perm is " + str(Cout_perm))
Cout_conc=(Comp1.Cconci+ Comp2.Cconci+Comp3.Cconci+Comp4.Cconci+Comp5.Cconci+Comp6.Cconci)
print("Cout_conc is " + str(Cout_conc))
Ci_out_c=[Comp1.Cconci, Comp2.Cconci, Comp3.Cconci, Comp4.Cconci, Comp5.Cconci, Comp6.Cconci]



Qconc=qc_out #kg/hr
Ci_out_c=Cconc_out
d_c=density_calc(25, sum(Ci_out_c))
print("Qconc is "+ str(Qconc))
Qperm=Comp1.Qperm #kg/hr
print("Qperm is "+ str(Qperm))
Ci_out_p=[Comp1.Cpermi, Comp2.Cpermi,Comp3.Cpermi, Comp4.Cpermi, Comp5.Cpermi, Comp6.Cpermi]
d_p=density_calc(25, sum(Ci_out_p))
Papplied_2=(P_osmo_c.p_osmo+P_osmo_f.p_osmo)/2-P_osmo_p.p_osmo+dp #print("applied pump pressure 2 pass " +str(Papplied_2))

Ppump_2=Papplied_2*Comp1.Qperm/d_p*1e5/3600#*1e5 #print("pump pressure 2 pass " +str(Ppump_2))

E_el_nf=Ppump_1/1000/0.8+Ppump_2/1000/0.8
print("E_el_nf"+str(E_el_nf)+"kw")
print("ci_out_p is "+ str(Ci_out_p))
print("ci_out_c is "+ str(Cconc))
spec=E_el_nf/(Qperm/d_p)
print("spec is "+str(spec))
SEC_el_feed=E_el_nf/(Qf/d_in) #kwh/m3
print("SEC_el is " +str(SEC_el_feed)+ " kWh/m3 feed")
rec=(Qperm/d_p)/(Qsw/d_in)

#chemical consumption 
Qhcl=0.27 #unit ml/h + convert it to kg 
Qantsc=5 #unit ml/hr+ convert it to kg 
#scale-up
Qantsc=5*(Qf/d_in)/2.5/1000
#mass balances 
for i in range(len(Cconc)):
    bal_i=Ci_in[i]*Qsw/d_in-(Qperm/d_p*Ci_out_p[i]+Qconc/d_c*Cconc[i])
    print(" for "+ str(i)+ " mass balance equals: " + str(bal_i))

dsw=1.02730778357126
lis=[]
bi=[]
#for i in range(len(Cconc)):
for i in range(2):
    lis.append([1/(Qperm*Ci_out_p[i]), 1/(Qconc*Ci_out_c[i])])
    bi.append([(Ci_in[i]*Qsw)/dsw])
print(lis)
A=np.array(lis)
B=np.array(bi)
print(B)
x=np.linalg.solve(A,B)
print("x is : " +str(x))
#d=1/x
#print(d)
print("Qconc_to_st2 is "+str(Qconc_to_st2))
balance=Qsw/dsw-Qperm/d_p-Qconc/d_c
balance=Qsw-Qperm-Qconc
print("balance is "+str(balance))
