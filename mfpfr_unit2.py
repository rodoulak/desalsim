
import math
import nanofiltration_unit
import constants 
from density_calc import density_calc
import scaleup
Cin=nanofiltration_unit.Cconc
#from nanofiltration_unit import *
#%%
#conditions
T=20+273.15

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

#constants
kps_MgOH=5.61*0.000000000001
kps_CaOH=5.5*0.000001
d_mgoh_2=2.34 #kg/l
d_caoh_2=2.211 #kg/l
#assumptions
dp=0.5
npump=0.8

#%%Input parameters 
class inputpar:
    def __init__(self, Qin, Cc1, Cc2, Cc3, Cc4, Cc5, Cc6, Cc7, C_NaOH_1, C_NaOH_2, conv_1, conv_2):
        self.Qin=Qin
        self.CNa_in=Cc1/MW_Na #Na concentration (g/l) to (mol/l)
        self.CCl_in=Cc2/MW_cl #Cl concentration (g/l) to (mol/l)
        self.CK_in=Cc3/MW_K #K concentration (g/l) to (mol/l)
        self.CMg_in=Cc4/MW_Mg #Mg concentration (g/l) to (mol/l)
        self.CCa_in=Cc5/MW_Ca #Ca concentration (g/l) to (mol/l)
        self.CSO4_in=Cc6/MW_so4 #SO4 concentration (g/l) to (mol/l)
        self.CHCO3_in=Cc7/MW_HCO3 #HCO3 concentration (g/l) to (mol/l)
        self.C_NaOH_1=C_NaOH_1 #Concentration of NaOH solution for first step (mol/l)
        self.C_NaOH_2=C_NaOH_2 #Concentration of NaOH solution for second step (mol/l)
        self.conv_1=conv_1 #conversion step 1 
        self.conv_2=conv_2 #conversion step 2 
    
    def calc_step1(self):
        self.QMg_in=self.Qin*self.CMg_in #MOL/H
        self.QNaOH_1=(self.Qin*self.CMg_in*(self.conv_1/100)*2)/self.C_NaOH_1 #L/hr
        self.M_MgOH2_1=(self.Qin*self.CMg_in*(self.conv_1/100)*MW_MgOH)/1000 #kg/hr
        self.Qtot_out_1=self.Qin+self.QNaOH_1
        self.Mtot_out_1=self.Qin*d_in+self.QNaOH_1*1.04-self.M_MgOH2_1 #KG/hr   
        self.magma_d_1=self.M_MgOH2_1/self.Qtot_out_1 #The magma density: the quantity of solids produced per volume of slurry [kg/l]
        self.ph_1=14+math.log10(2*(kps_MgOH/4)**(1/3))
        self.Qtot_out_1=self.Qin+self.QNaOH_1-self.M_MgOH2_1/2.34 #l/hr #mgoh2 density is 2.34 kg/l
        self.CMg_out_1=(self.QMg_in*(1-self.conv_1/100))/self.Qtot_out_1 #mol/l
        self.CNa_out_1=(self.Qin*self.CNa_in+self.QNaOH_1*self.C_NaOH_1)/self.Qtot_out_1
        self.CCl_out_1=(self.Qin*self.CCl_in)/self.Qtot_out_1
        self.CK_out_1=(self.Qin*self.CK_in)/self.Qtot_out_1
        self.CCa_out_1=(self.Qin*self.CCa_in)/self.Qtot_out_1
        self.CSO4_out_1=(self.Qin*self.CSO4_in)/self.Qtot_out_1
        self.CHCO3_out_1=(self.Qin*self.CHCO3_in)/self.Qtot_out_1
        
    def calc_step2(self):
       self.QCa_in_2=self.Qtot_out_1*self.CCa_out_1 #mol/h
       self.QNaOH_2_st=(self.Qtot_out_1*(self.CCa_out_1*(self.conv_2/100)+self.CMg_out_1*1)*2)/self.C_NaOH_2 # the stoichiometric volumetric flow rate of sodium hydroxide for the second step L/hr
       self.COH_ph13=0.1 #mol/l
       self.QNaOH_2_add=((0.0216-self.COH_ph13)*(self.QNaOH_2_st+self.Qtot_out_1))/(self.COH_ph13-self.C_NaOH_2)
       self.Qtot_out_2=self.Qtot_out_1+self.QNaOH_2_st+self.QNaOH_2_add
       self.M_CaOH2_2=(self.Qtot_out_1*self.CCa_out_1*(self.conv_2/100)*MW_CaOH)/1000 #kg/hr
       self.M_MgOH2_2=self.Qtot_out_1*self.CMg_out_1*MW_MgOH/1000 #kg/hr
       self.magma_d_2=(self.M_MgOH2_2+self.M_CaOH2_2)/self.Qtot_out_2 #The magma density: the quantity of solids produced per volume of slurry [kg/l]
       self.Qout_2=self.Qtot_out_2-self.M_CaOH2_2/2.211-self.M_MgOH2_2/2.34
       self.CNa_out_2=(self.Qtot_out_1*self.CNa_out_1+(self.QNaOH_2_st+self.QNaOH_2_add)*self.C_NaOH_2)/self.Qtot_out_2       
       self.CCa_out_2=(self.QCa_in_2*(1-self.conv_2/100))/self.Qtot_out_2
       self.CCl_out_2=(self.Qtot_out_1*self.CCl_out_1)/self.Qtot_out_2
       self.CK_out_2=(self.Qtot_out_1*self.CK_out_1)/self.Qtot_out_2
       self.CMg_out_2=0
       self.CSO4_out_2=(self.Qtot_out_1*self.CSO4_out_1)/self.Qtot_out_2
       self.CHCO3_out_2=(self.Qtot_out_1*self.CHCO3_out_1)/self.Qtot_out_2

       self.ph_2=14+math.log10(0.1)
       
       self.Epump_1=(self.Qin*dp+self.QNaOH_1*dp)*1e5/3600/(1000*npump) #(W)
       self.Epump_2=(self.Qtot_out_1*dp+(self.QNaOH_2_add+self.QNaOH_2_st)*dp)*1e5/3600/(1000*npump)       #(W)
       
#%%Energy consumption 
class energycons:
    def energycalc(self, Qtot, QNaOH):
        self.Qtot=Qtot
        self.QNaOH=QNaOH
        self.Epump=(self.Qtot*dp+self.QNaOH*dp)/(1000*npump)
        print(self.Epump)
           
#%%mf-pfr
C_NaOH_1=1 #Concentration of NaOH solution first step MOL/L  
C_NaOH_2=1 #Concentration of NaOH solution first step MOL/L       
conv_1=95 # concertion rate step 1
conv_2=93 # concertion rate step 2
Cin_mfpfr=nanofiltration_unit.Cconc      
d_in=density_calc(25, sum(Cin_mfpfr))/1000
print("d is "+str(d_in)) 
Qin_mfpfr=nanofiltration_unit.Qconc/d_in #l/hr
print("nanofiltration_unit.Qconc"+str(nanofiltration_unit.Qconc))
#def __init__(self, Qin, Cc1, Cc2, Cc3, Cc4, Cc5, Cc6, Cc7, C_NaOH_1, C_NaOH_2, conv_1, conv_2):
mfpfr_dat=inputpar(Qin_mfpfr, Cin_mfpfr[0], Cin_mfpfr[1], Cin_mfpfr[2], Cin_mfpfr[3], Cin_mfpfr[4], Cin_mfpfr[5],  0, C_NaOH_1, C_NaOH_2, conv_1, conv_2)
#Cin_mfpf=Cout_conc        

print(mfpfr_dat.CNa_in)
mfpfr_dat.calc_step1()
mfpfr_dat.calc_step2()

Cour_mfpfr=mfpfr_dat.CNa_out_2+mfpfr_dat.CCa_out_2+mfpfr_dat.CCl_out_2+ mfpfr_dat.CK_out_2+mfpfr_dat.CMg_out_2+mfpfr_dat.CSO4_out_2#mol/l
CSO4_out_2=mfpfr_dat.CSO4_out_2 #mol/l
CNa_out_2=mfpfr_dat.CNa_out_2
Cnacl_out=CNa_out_2-2*CSO4_out_2
M_MgOH2_1=mfpfr_dat.M_MgOH2_1
M_CaOH2=mfpfr_dat.M_CaOH2_2
M_MgOH2=mfpfr_dat.M_MgOH2_2
print("na_out_2 is "+ str(CNa_out_2))
print("cout is "+ str(CSO4_out_2))
print("nacl_out_2 is "+ str(Cnacl_out))
print(Cour_mfpfr)
Cout_mfpfr_g=mfpfr_dat.CNa_out_2*MW_Na+mfpfr_dat.CCa_out_2*MW_Ca+mfpfr_dat.CCl_out_2*MW_cl+ mfpfr_dat.CK_out_2*MW_K+mfpfr_dat.CMg_out_2*MW_Mg+mfpfr_dat.CSO4_out_2*MW_so4 #g/l

print("total concentration in effluent is " + str(Cout_mfpfr_g))
Cout_all_m=[mfpfr_dat.CNa_out_2, mfpfr_dat.CCl_out_2, mfpfr_dat.CK_out_2, mfpfr_dat.CMg_out_2, mfpfr_dat.CCa_out_2, mfpfr_dat.CSO4_out_2] #mol/l
Cout_mfpfr_g=[mfpfr_dat.CNa_out_2*MW_Na,mfpfr_dat.CCl_out_2*MW_cl, mfpfr_dat.CK_out_2*MW_K, mfpfr_dat.CMg_out_2*MW_Mg, mfpfr_dat.CCa_out_2*MW_Ca,mfpfr_dat.CSO4_out_2*MW_so4]
Qout_2=mfpfr_dat.Qout_2

#chemical consumption 
QNAOH=mfpfr_dat.QNaOH_1+mfpfr_dat.QNaOH_2_add+mfpfr_dat.QNaOH_2_st #convert to kg 
print("QNAOH is "+str(QNAOH))
print("QNaOH_1 is "+str(mfpfr_dat.QNaOH_1)) 
print("mfpfr_dat.QNaOH_2_add is "+str(mfpfr_dat.QNaOH_2_add*24*300/1000)) #m3/year
print("mfpfr_dat.QNaOH_2_st is "+str(mfpfr_dat.QNaOH_2_st*24*300/1000)) #m3/year
#hcl to decrease ph to 7
#QHCl=((math.pow(10,-7)-math.pow(10,-mfpfr_dat.ph_2))*mfpfr_dat.Qout_2)/(1-math.pow(10,-7))

HCl_conc=1 #l of HCl 1M
OH_initial=math.pow(10, -14)/math.pow(10, -mfpfr_dat.ph_2)
OH_final=math.pow(10, -14)/math.pow(10, -7)
QHCl=(Qout_2*OH_initial-Qout_2*OH_final)/(OH_final+HCl_conc) #l of HCl 1M
Qout_2=Qout_2+QHCl
C_cl_out=(Cout_all_m[1]*mfpfr_dat.Qout_2+QHCl*HCl_conc)/Qout_2 #mol/l 
Cout_all_m[1]=C_cl_out
C_cl_out=C_cl_out*MW_cl
Cout_mfpfr_g[1]=C_cl_out
for i in range(2,6):
    Cout_mfpfr_g[i]=Cout_mfpfr_g[i]*mfpfr_dat.Qout_2/Qout_2
Cout_mfpfr_g[0]=Cout_mfpfr_g[0]*mfpfr_dat.Qout_2/Qout_2
d_out_s=density_calc(25,sum(Cout_mfpfr_g))/1000 #kg/m3
Mout_2=Qout_2*d_out_s
print("C_cl_out is "+ str(C_cl_out)+ " mol/h")
print("HCl flow rate is "+ str(QHCl)+ " l/h")
#electicity consumption 
E_el_mfpf=(mfpfr_dat.Epump_1+mfpfr_dat.Epump_2+(QHCl*0.3)*1e5/3600/(1000*npump))/1000
print("Electricity energy consumption is "+str(E_el_mfpf)+ " kwh")

E_fil=scaleup.scaleup(0.5, 0.3*1000, Mout_2)
E_el_mfpf=E_el_mfpf+E_fil
SEC_el_prod=(E_el_mfpf)/(M_MgOH2)
print("SEC_el_prod is "+str(SEC_el_prod)+" KWh/kg product ")
SEC_el_feed=(E_el_mfpf)/(Qin_mfpfr/1000)
print("SEC_el_feed is "+str(SEC_el_feed)+" KWh/m3 of feed ")
#mass balance:
bal=nanofiltration_unit.Qconc-M_MgOH2_1-M_CaOH2-M_MgOH2-Mout_2+QNAOH*1.04+QHCl*1.0585
print("bal is "+str(bal))

print("magma_d_1 "+str(mfpfr_dat.magma_d_1))

print("mfpfr_dat.magma_d_2"+str(mfpfr_dat.magma_d_2))
print("qout1 IS  "+str(mfpfr_dat.Qtot_out_1))
