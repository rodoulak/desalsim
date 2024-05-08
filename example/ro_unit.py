import math
from Desalsim.density_calc import density_calc 
from Desalsim import constants 
#%%
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
#        print(" ion concenctration of " + self.comp+ " in permeate stream is"+" "+ str(self.Cpermi))
#        print(" ion concenctration of " + self.comp+ " in concentrate stream is"+" "+str(self.Cconci))

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
        
    # def osmotic_pressure_calculation(self):
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
            
#%%Energy consumption 
class NfEnergy:   
    def __init__(self, P_osmo_c, P_osmo_f, P_osmo_p, dp, d_p, Qperm, Qf, d_in):
        self.P_osmo_c = P_osmo_c # Osmotic pressure of concentrate stream (units:bar)        
        self.P_osmo_f = P_osmo_f # Osmotic pressure of feed stream (units:bar)
        self.P_osmo_p = P_osmo_p # Osmotic pressure of permeate stream (units:bar)
        self.dp = dp #pressure drop (units:bar)
        self.d_p = d_p #Permeate stream density 
        self.Qperm = Qperm # Permeate flow rate
        self.Qf = Qf # Concentrate flow rate
        self.d_in = d_in #Feed stream density 

    def calculate_energy_consumption(self):
        Papplied = (self.P_osmo_c + self.P_osmo_f) / 2 - self.P_osmo_p + self.dp #Applied pressure (units: bar)
        Ppump = Papplied * self.Qperm / self.d_p * 1e5 / 3600 #Power for applied pressure (units:W)
        E_el_nf = (Ppump / 1000 / 0.8 ) #Electricity consumption (units:KW)
        spec = E_el_nf / (self.Qperm / self.d_p) #Specific energy consumption Kwh/m3 of permeate
        SEC_el_feed = E_el_nf / (self.Qf / self.d_in) #Specific energy consumption Kwh/m3 of feed
        
        return {
    "Applied pressure": Papplied,
    "Power for pump": Ppump,
    "E_el_nf": E_el_nf,
    "Spec": spec,
    "SEC_el_feed": SEC_el_feed
}

#%%
dp=2
        
Na=molarity(constants.MW_Na, 1, 11.9)

Ci_in=[12.33, 21.67, 0.45, 1.39, 0.45, 3.28]   
comp=['Na', 'Cl', 'K', 'Mg', 'Ca', 'SO4']
mg_in=sum(Ci_in)
print("mg_in is " + str(mg_in))  
d_in=density_calc.density_calc(25, mg_in)
print("d_in is " + str(d_in)) #kg/m3 
Qsw = 3000/24*d_in #149612.002308333#2290.858979 #kg/hr
Qf=Qsw 


rjr=[0.85, 0.85,0.85,0.88,0.78,0.96]
Wrec=0.4
#input data and call function
##(self, comp, Cfeedi, rjr, Wrec, Qf):        
 
Comp1=nfmass('Na', Ci_in[0], rjr[0], Wrec, Qf)    
Comp2=nfmass('Cl', Ci_in[1], rjr[1], Wrec, Qf) 
Comp3=nfmass('K', Ci_in[2], rjr[2], Wrec, Qf) 
Comp4=nfmass('Mg', Ci_in[3],rjr[3], Wrec, Qf) 
Comp5=nfmass('Ca', Ci_in[4], rjr[4], Wrec, Qf) 
Comp6=nfmass('SO4', Ci_in[5], rjr[5], Wrec, Qf) 


Comp1.permcal()
Comp2.permcal()
Comp3.permcal()
Comp4.permcal()
Comp5.permcal()
Comp6.permcal()

#output
Cconc_ro=[Comp1.Cconci, Comp2.Cconci, Comp3.Cconci, Comp4.Cconci, Comp5.Cconci, Comp6.Cconci]

Cin_conc_ro=(Comp1.Cfeedi+ Comp2.Cfeedi+Comp3.Cfeedi+Comp4.Cfeedi+Comp5.Cfeedi+Comp6.Cfeedi)
print("Cin is " + str(Cin_conc_ro))
Cperm_ro=[Comp1.Cpermi, Comp2.Cpermi, Comp3.Cpermi, Comp4.Cpermi, Comp5.Cpermi, Comp6.Cpermi]
Qperm_ro=Comp1.Qperm
Qconc_ro=Comp1.Qconc
print("Qperm_ro is "+str(Qperm_ro))
print("Qconc_ro is "+str(Qconc_ro))
Qantsc=60*(Qf/d_in)/1000

P_osmo_f=Osmotic_pressure(Comp1.Cfeedi, 1,Comp2.Cfeedi, -1, Comp3.Cfeedi,1, Comp4.Cfeedi,2, Comp5.Cfeedi, 2, Comp6.Cfeedi, -2)
print("P_osmo_f is "+str(P_osmo_f))
P_osmo_p=Osmotic_pressure(Comp1.Cpermi,1, Comp2.Cpermi, -1,Comp3.Cpermi,1, Comp4.Cpermi,2, Comp5.Cpermi, 2,Comp6.Cpermi, -2) 
P_osmo_c=Osmotic_pressure(Comp1.Cconci,1, Comp2.Cconci,-1, Comp3.Cconci,1, Comp4.Cconci,2, Comp5.Cconci, 2,Comp6.Cconci, -2)
print("P_osmo_c"+str(P_osmo_c.p_osmo))

#mass balances 
for i in range(len(Cconc_ro)):
    bal_i=Ci_in[i]*Qsw-(Qperm_ro*Cperm_ro[i]+Qconc_ro*Cconc_ro[i])
    print(" for "+ str(i)+ " mass balance equals: " + str(bal_i))

#electircal consumption 
Papplied_ro=(P_osmo_c.p_osmo+P_osmo_f.p_osmo)/2-P_osmo_p.p_osmo+dp #print("applied pump pressure 2 pass " +str(Papplied_2))
print("Papplied_ro is "+str(Papplied_ro))
Ppump_ro=Papplied_ro*Comp1.Qperm*100000/1000/3600#*1e5 #(W)
print("Ppump_ro is "+str(Ppump_ro))
#Ppump_ro=Papplied_ro*0.5*0.8*Comp1.Qperm/1000 #print("pump pressure 2 pass " +str(Ppump_2))
d_p=1

elec=NfEnergy(P_osmo_c.p_osmo, P_osmo_f.p_osmo, P_osmo_p.p_osmo, 2, d_p, Qperm_ro, Qf, d_in)
result=elec.calculate_energy_consumption()
for key, value in result.items():
        print(f"{key}: {value}")
