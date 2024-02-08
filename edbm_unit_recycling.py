import math 
import numpy as np
import efc_unit
import scaleup 
import density_calc
import pandas as pd
import constants

#%%constants

I_d=400 # The electricl current desnity Am2

# A=0.16 #active area of the membrane across which ion permeation occurs (A)


F=96485.3 #Coulombs/mol
R_const=8.314462618 #kg⋅m2⋅s−2⋅K−1⋅mol−1
R_int=0.28#110#45 #ohm cm2
z=1
npump=0.8
dp=1.5
#Membrane characteristics

Cm_bp_H= 0.0000001 #mol/l 
Cm_bp_OH= 0.0000001 #mol/l 

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
MW_H=1.00784
MW_OH=17.008
MW_H2O=MW_H+MW_OH

T=20+273.15

#input data
r=(5-0.4)/(5) #recycling rate 
rs=(5-1)/5
print("r0 is "+str(r))
ph_s=4.71
ph_b=7
ph_a=7
d_a=1 #kg/l
d_s=1 #kg/l
d_b=1 #kg/l
#%%Calculations 
class edbmcacl:
    def __init__(self, Qin, Cc1, Cc2, Cc3, Cc4, Cc5, Cc6, Cc7, CH, COH, N,Cb1, Cb2,Cb3,Cb4,Cb5,Cb6,Cb7,Cb8,Cb9,Ca1,Ca2,Ca3,Ca4,Ca5,Ca6,Ca7,Ca8,Ca9):
        self.Qin=Qin
        self.CNa_s_in=Cc1/MW_Na #Na concentration (g/l) to (mol/l)
        self.CCl_s_in=Cc2/MW_cl #Cl concentration (g/l) to (mol/l)
        self.CK_s_in=Cc3/MW_K #K concentration (g/l) to (mol/l)
        self.CMg_s_in=Cc4/MW_Mg #Mg concentration (g/l) to (mol/l)
        self.CCa_s_in=Cc5/MW_Ca #Ca concentration (g/l) to (mol/l)
        self.CSO4_s_in=Cc6/MW_so4 #SO4 concentration (g/l) to (mol/l)
        self.CHCO3_s_in=Cc6/MW_HCO3 #HCO3 concentration (g/l) to (mol/l)      
        self.COH_s_in=COH/MW_OH
        self.PM_i=[MW_Na, MW_cl, MW_K, MW_Mg, MW_Ca, MW_so4, MW_HCO3, MW_H, MW_OH]
        self.N_trip=N #Number of triplets of channel (N_trip) 
        self.CH_s_in=CH/MW_H
        self.Ci_s_in=[self.CNa_s_in,self.CCl_s_in,self.CK_s_in,self.CMg_s_in,self.CCa_s_in,self.CSO4_s_in,self.CHCO3_s_in, self.CH_s_in,self.COH_s_in]
        self.Ci_b_in=[Cb1/MW_Na, Cb2/ MW_cl,Cb3/MW_K,Cb4/MW_Mg,Cb5/MW_Ca,Cb6/MW_so4,Cb7/MW_HCO3,1e-7,1e-7]
        self.Ci_a_in=[Ca1/MW_Na,Ca2/MW_cl,Ca3/MW_K,Ca4/MW_Mg,Ca5/MW_Ca,Ca6/MW_so4,Ca7/MW_HCO3,1e-7,1e-7]
#save initial concentration before recirculation: 
        self.Ci_s_in_0=self.Ci_s_in
        self.Ci_b_in_0=self.Ci_b_in
        self.Ci_a_in_0=self.Ci_a_in
        
    def flowrate(self):
        self.Q1_s_in=self.Qin/self.N_trip #units: l/h
        self.Q1_b_in=self.Qin/self.N_trip #units: l/h
        self.Q1_a_in=self.Qin/self.N_trip #units: l/h

    def in_mass_flow_rates(self):
        self.M_s_in=[[],[],[],[],[],[],[],[],[]]
        self.M_b_in=[[],[],[],[],[],[],[],[],[]]
        self.M_a_in=[[],[],[],[],[],[],[],[],[]]
                
        for i in range(0,9):
            self.M_s_in[i]=self.Q1_s_in*self.Ci_s_in[i]*self.PM_i[i]/1000 #units: kg/h
            self.M_b_in[i]=self.Q1_b_in*self.Ci_b_in[i]*self.PM_i[i]/1000 #units: kg/h
            self.M_a_in[i]=self.Q1_a_in*self.Ci_a_in[i]*self.PM_i[i]/1000 #units: kg/h

    
        self.M_h2o_a_in=self.Q1_a_in*d_a-sum(self.M_a_in)          
        self.M_h2o_b_in=self.Q1_b_in*d_b-sum(self.M_b_in) 
        self.M_h2o_s_in=self.Q1_s_in*d_s-sum(self.M_s_in) 
        self.KW_s_in=self.Ci_s_in[7]*self.Ci_s_in[8]
        self.KW_a_in=self.Ci_a_in[7]*self.Ci_a_in[8]
        self.KW_b_in=self.Ci_b_in[7]*self.Ci_b_in[8]
    
        if ph_s<7:
            self.M_s_in[1]=self.Q1_s_in*self.Ci_s_in[1]*self.PM_i[1]/1000+(self.M_s_in[7]/self.PM_i[7]-self.M_s_in[8]/self.PM_i[8])*(self.PM_i[1])
        elif ph_s>7:
            self.M_s_in[1]=self.Q1_s_in*self.Ci_s_in[1]*self.PM_i[1]/1000+(self.M_s_in[8]/self.PM_i[8]-self.M_s_in[7]/self.PM_i[7])*self.PM_i[1]
        

    def acid_channel(self):
        self.PM_i=[MW_Na, MW_cl, MW_K, MW_Mg, MW_Ca, MW_so4, MW_HCO3, MW_H, MW_OH]
        self.M_a_out=[[],[],[],[],[],[],[],[],[]]
        self.Ci_a_out=[[],[],[],[],[],[],[],[],[]]
        self.I_ext=A*I_d
        self.JA=3.6*self.I_ext/F/1000
        print("JA is "+str(self.JA))
        self.M_a_out[7]=self.M_a_in[7]+self.JA*self.PM_i[7]
        #print("self.M_a_out_h" + str(self.M_a_out[7]))
        self.M_a_out[1]=self.M_a_in[1]+self.JA*self.PM_i[1]
        #print("self.M_a_out_cl" + str(self.M_a_out[1]))
        self.M_a_out[8]=self.M_a_in[8]
        #print("self.M_a_out_oh" + str(self.M_a_out[8]))
        self.M_a_out[0]=self.M_a_in[0]
        #print("self.M_a_out_na" + str(self.M_a_out[0]))
        self.M_h2o_a_out=self.M_h2o_a_in-0.5*self.JA*MW_H2O
        #print("self.M_h2o_a_out" + str(self.M_h2o_a_out))
        self.M_a_out_t=self.M_h2o_a_out+self.M_a_out[7]+self.M_a_out[1]+self.M_a_out[8]+self.M_a_out[0]
        
        for i in range(2,7):
            self.M_a_out[i]=self.M_a_in[i]
            self.M_a_out_t=self.M_a_out_t+self.M_a_out[i]
        #print("self.M_a_out" + str(self.M_a_out))
        self.Q1_a_out=self.M_a_out_t/1
        #print("self.Q1_a_out" + str(self.Q1_a_out))
        for i in range(0,9):
            self.Ci_a_out[i]=self.M_a_out[i]/(self.Q1_a_out*self.PM_i[i]/1000)
        #print("self.Ci_a_out" + str(self.Ci_a_out))
        
    def base_channel(self):
        self.M_b_out=[[],[],[],[],[],[],[],[],[]]
        self.Ci_b_out=[[],[],[],[],[],[],[],[],[]]
        self.M_b_out[0]=self.M_b_in[0]+self.JA*self.PM_i[0]
        self.M_b_out[8]=self.M_b_in[8]+self.JA*self.PM_i[8]
        self.M_b_out[1]=self.M_b_in[1]
        self.M_b_out[7]=self.M_b_in[7]
        self.M_h2o_b_out=self.M_h2o_b_in+0.5*self.JA*MW_H2O
        self.M_b_out_t=self.M_h2o_b_out+self.M_b_out[7]+self.M_b_out[1]+self.M_b_out[8]+self.M_b_out[0]
    
        for i in range(2,7):
            self.M_b_out[i]=self.M_b_in[i]  
            self.M_b_out_t=self.M_b_out_t+self.M_b_out[i]
  
        self.Q1_b_out=self.M_b_out_t/1 #assume density is 1kg/l
        #print("self.Q1_b_out" + str(self.Q1_b_out))
        for i in range(0,9):
            self.Ci_b_out[i]=self.M_b_out[i]/(self.Q1_b_out*self.PM_i[i]/1000)
        print("self.Ci_b_out is m " + str(self.Ci_b_out))
    
    def salt_channel(self):
        self.M_s_out=[[],[],[],[],[],[],[],[],[]]
        self.Ci_s_out=[[],[],[],[],[],[],[],[],[]]
        self.M_s_out[0]=self.M_s_in[0]-(self.M_b_out[0]-self.M_b_in[0])
        print("self.M_s_in na is " + str(self.M_s_in[0]))
        print("self.M_b_out na is "+str(self.M_b_out[0]))
        print("self.M_b_in na is "+str(self.M_b_in[0]))
        self.M_s_out[1]=self.M_s_in[1]-(self.M_a_out[1]-self.M_a_in[1])
        self.M_s_out[8]=self.M_s_in[8]
        self.M_s_out[7]=self.M_s_in[7]
        self.M_h2o_s_out=self.M_h2o_s_in
        self.M_s_out_t=self.M_s_out[0]+self.M_s_out[1]+self.M_s_out[8]+self.M_s_out[7]+self.M_h2o_s_out
        for i in range(2,7):
            self.M_s_out[i]=self.M_s_in[i]
            self.M_s_out_t=self.M_s_out_t+self.M_s_out[i]
        self.Q1_s_out=self.M_s_out_t/1.1
        for i in range(0,9):
            self.Ci_s_out[i]=self.M_s_out[i]/(self.Q1_s_out*self.PM_i[i]/1000)     
        print("self.M_s_out of na is "+str(self.M_s_out[0]))
        print("self.Ci_s_out is m " + str(self.Ci_s_out))
        
        self.c1_emf=math.log(((self.Ci_a_in[7]+self.Ci_a_out[7])/2)/Cm_bp_H)        
        self.c2_emf=math.log(Cm_bp_OH/((self.Ci_b_in[8]+self.Ci_b_out[8])/2))        
        self.c3_emf=math.log(((self.Ci_s_in[1]+self.Ci_s_out[1])/2)/((self.Ci_a_in[1]+self.Ci_a_out[1])/2))
        print("self.c3_emf   " + str(self.c3_emf))
        self.c4_emf=math.log(((self.Ci_b_in[0]+self.Ci_b_out[0])/2)/((self.Ci_s_in[0]+self.Ci_s_out[0])/2))
        print("self.c4_emf   " + str(self.c4_emf))        
        self.EMF=((R_const*T/z/F)*(self.c1_emf+self.c2_emf+self.c3_emf+self.c4_emf))*self.N_trip
       # self.EMF=(((R_const*T)/z/F)*(math.log(((self.Ci_a_in[7]+self.Ci_a_out[7])/2)/Cm_bp_H)+math.log(Cm_bp_OH/(self.Ci_b_in[8]+self.Ci_b_out[8])/2)+math.log(((self.Ci_s_in[1]+self.Ci_s_out[1])/2)/((self.Ci_a_in[1]+self.Ci_a_out[1])/2))+ math.log(((self.Ci_b_in[0]+self.Ci_b_out[0])/2)/((self.Ci_s_in[0]+self.Ci_s_out[0])/2))))*self.N_trip#incomplete       
        self.V_ext=self.EMF+((self.I_ext*R_int)/(A*10000))*self.N_trip
        self.P=self.V_ext*self.I_ext 
        
    def recycl_below1M(self):
        self.M_s_in=self.M_s_out
        self.M_b_in=self.M_b_out
        self.M_a_in=self.M_a_out
        self.M_s_in_t=self.M_s_out_t
        self.M_a_in_t=self.M_a_out_t
        self.M_b_in_t=self.M_b_out_t
        self.Ci_s_in=self.Ci_s_out
        self.Ci_b_in=self.Ci_b_out
        self.Ci_a_in=self.Ci_a_out
        self.d_a_out=density_calc.density_calc(25, sum(self.Ci_a_out))/1000
        self.d_b_out=density_calc.density_calc(25, sum(self.Ci_b_out))/1000
        self.d_s_out=density_calc.density_calc(25, sum(self.Ci_s_out))/1000
        print("self.M_s_in" + str(self.M_s_in))
        print("self.Ci_s_in" + str(self.Ci_s_in))
        print("self.Ci_b_in" + str(self.Ci_b_in))
        print("self.M_b_in_t   " + str(self.M_b_in_t))
        

    def recycl_1M(self):
        self.Q_s_r=rs*self.M_s_out_t/self.d_s_out#r*self.Q1_s_in/(1-r)
        self.Q_s_t_out=self.Q1_s_in/(1-rs)
        self.Q_s_out=self.Q1_s_in
        self.Q_s_in_t=self.Q1_s_in+self.Q_s_r
        self.M_s_r=self.Q_s_r/self.d_s_out
        print("self.M_s_r "+str(self.M_s_r))
        self.M_s_t_out=self.Q_s_t_out/self.d_s_out
        #self.M_s_out=self.Q_s_out/self.d_s_out
        self.M_s_in_t=self.Q1_s_in/d_in+self.Q_s_r/self.d_s_out
        #self.M_s_in=self.M_s_in_t
        #acid
        self.Q_a_r=r*self.M_a_out_t/self.d_a_out#r*self.Q1_a_in/(1-r)
        self.Q_a_t_out=self.Q_a_r/r#self.Q1_a_in/(1-r)
        self.Q_a_out=self.Q1_a_in
        self.Q_a_in_t=self.Q1_a_in+self.Q_a_r
        self.M_a_r=self.Q_a_r/self.d_a_out
        self.M_a_t_out=self.Q_a_t_out/self.d_a_out
        #self.M_a_out=self.Q_a_out/self.d_a_out
        self.M_a_in_t=self.Q1_a_in/d_in+self.Q_a_r/self.d_a_out
        #self.M_a_in=self.M_a_in_t
        #base
        self.Q_b_r=r*self.M_b_out_t/self.d_b_out#r*self.Q1_b_in/(1-r)
        self.Q_b_t_out=self.Q1_b_in/(1-r)
        self.Q_b_out=self.Q1_b_in
        self.Q_b_in_t=self.Q1_b_in+self.Q_b_r
        self.M_b_r=self.Q_b_r/self.d_b_out
        self.M_b_t_out=self.Q_b_t_out/self.d_b_out
        #self.M_b_out=self.Q_b_out/self.d_b_out
        self.M_b_in_t=self.Q1_b_in/d_in+self.Q_b_r/self.d_b_out       
        #self.M_b_in=self.M_b_in_t
        print("Q1_a_in is "+str(self.Q1_s_in))
        print("Q_a_r is "+str(self.Q_s_r))
        print("Q_a_t_out is "+str(self.Q_s_t_out))
        print("Q_a_out is "+str(self.Q_s_out))
        print("Q_a_t_in is "+str(self.Q_s_in_t))
            
        
        for i in range(len(self.Ci_a_in)):
            self.Ci_s_in[i]=self.Ci_s_out[i]*self.M_s_r/self.M_s_in_t+self.Ci_s_in_0[i]*(self.M_s_in_t-self.M_s_r)/(self.M_s_in_t)      
            self.M_s_in[i]=self.Ci_s_in[i]*self.Q_s_in_t
            self.Ci_b_in[i]=self.Ci_b_out[i]*self.M_b_r/self.M_b_in_t+self.Ci_b_in_0[i]*(self.M_b_in_t-self.M_b_r)/(self.M_b_in_t) 
            self.M_b_in[i]=self.Ci_a_in[i]*self.Q_b_in_t
            self.Ci_a_in[i]=self.Ci_a_out[i]*self.Q_a_r/self.Q_a_in_t+self.Ci_a_in_0[i]*(self.Q1_a_in)/(self.Q_a_in_t) 
            self.M_a_in[i]=self.Ci_a_in[i]*self.Q_a_in_t
        if self.Ci_a_in[1]>=1:
            print("Ms recycling is "+str( self.M_s_r))
            print("Ma recycling is "+str( self.M_a_r))
            print("Mb recycling is "+str( self.M_b_r))
            print("recycling rate is " +str(r))
        print("Cs_new_in_i is "+str(self.Ci_s_in))
        print("Ca_new_in_i is "+str(self.Ci_a_in))
        print("Cb_new_in_i is "+str(self.Ci_b_in))
        # self.d_b=density_calc.density_calc(25,sum(self.Ci_b_out))/1000
        # self.Q_b_r=self.M_b_r/self.d_b
        # self.Q_b_out_f=self.Q1_b_out-self.Q_b_r
        self.M_b_out_f=self.Q_b_out*self.d_b_out
        # self.d_a=density_calc.density_calc(25,sum(self.Ci_a_out))/1000
        # self.Q_a_r=self.M_a_r/self.d_a
        # self.Q_a_out_f=self.Q1_a_out-self.Q_a_r
        self.M_a_out_f=self.Q_a_out*self.d_a_out
        # self.d_s=density_calc.density_calc(25,sum(self.Ci_s_out))/1000
        # self.Q_s_r=self.M_s_r/self.d_s
        # self.Q_s_out_f=self.Q1_s_out-self.Q_s_r
        self.M_s_out_f=self.Q_s_out*self.d_s_out
        
        self.M_s_in=[[],[],[],[],[],[],[],[],[]]
        self.M_b_in=[[],[],[],[],[],[],[],[],[]]
        self.M_a_in=[[],[],[],[],[],[],[],[],[]]
                
        for i in range(0,9):
            self.M_s_in[i]=self.Q_s_in_t*self.Ci_s_in[i]*self.PM_i[i]/1000 #units: kg/h
            self.M_b_in[i]=self.Q_b_in_t*self.Ci_b_in[i]*self.PM_i[i]/1000 #units: kg/h
            self.M_a_in[i]=self.Q_a_in_t*self.Ci_a_in[i]*self.PM_i[i]/1000 #units: kg/h
        print("self.M_a_in[7] is "+str(self.M_a_in[7]))
    
        self.M_h2o_a_in=self.Q_a_in_t*d_a-sum(self.M_a_in)          
        self.M_h2o_b_in=self.Q_b_in_t*d_b-sum(self.M_b_in) 
        self.M_h2o_s_in=self.Q_s_in_t*d_s-sum(self.M_s_in) 
        self.KW_s_in=self.Ci_s_in[7]*self.Ci_s_in[8]
        self.KW_a_in=self.Ci_a_in[7]*self.Ci_a_in[8]
        self.KW_b_in=self.Ci_b_in[7]*self.Ci_b_in[8]
    
        if ph_s<7:
            self.M_s_in[1]=self.Q_s_in_t*self.Ci_s_in[1]*self.PM_i[1]/1000+(self.M_s_in[7]/self.PM_i[7]-self.M_s_in[8]/self.PM_i[8])*(self.PM_i[1])
        elif ph_s>7:
            self.M_s_in[1]=self.Q_s_in_t*self.Ci_s_in[1]*self.PM_i[1]/1000+(self.M_s_in[8]/self.PM_i[8]-self.M_s_in[7]/self.PM_i[7])*self.PM_i[1]
            
        print("self.M_a_in[i] is "+str(self.M_a_in))
        
#%%Main      
 
Cna_s=[]
#Qin, Cc1, Cc2, Cc3, Cc4, Cc5, Cc6, Cc7, CH, COH, N,Cb1, Cb2,Cb3,Cb4,Cb5,Cb6,Cb7,Cb8,Cb9,Ca1,Ca2,Ca3,Ca4,Ca5,Ca6,Ca7,Ca8,Ca9

# Cin_edbm=[13.439897042995193, 20.719916965189295, 1.1353903085911996, 0, 0, 0.1812370596056371]#[23.85609967, 38.70863206, 0.802553685, 0, 1.458379022, 2.07087518]
Cin_edbm=efc_unit.Cout_efc_g
d_in=density_calc.density_calc(25,sum(Cin_edbm))/1000
Q_in_edbm=efc_unit.Mout/d_in#l/hr #*(1+r)#scaleup.scaleup(226.6711827,377.7853045,efc_conc_step.Qconc)/d_in#229
print("Q_in_edbm: 1+r "+str(Q_in_edbm))
Mw_in=Q_in_edbm/d_in 
N_trip=20
A=scaleup.scaleup(19.2, 300, Q_in_edbm)/3 #19.2*Q_in_edbm/300
#N_trip=round(scaleup.scaleup(20, 300, Q_in_edbm),0)       
# A=0.4#19/40
# A_tot=19*Q_in_edbm/60 #scaleup.scaleup(0.19, 100, Q_in_edbm)#19*Q_in_edbm/60/3#
# print("A is "+str(A))
# N_trip=round(A_tot/A,0)/10
# print("N_trip is "+str(N_trip))
# print("A_tot is "+str(A_tot))

edbm_dat=edbmcacl(Q_in_edbm, Cin_edbm[0], Cin_edbm[1], Cin_edbm[2], Cin_edbm[3], Cin_edbm[4], Cin_edbm[5], 0,  10**(-ph_s), 3.01551E-11, round(N_trip,0), 0,0,0,0,0,0,0,10**(-ph_b), 10**(-(14-ph_b)), 0,0,0,0,0,0,0,10**(-ph_a), 10**(-(14-ph_a)))

edbm_dat.flowrate()
edbm_dat.in_mass_flow_rates()
#edbm_dat.edbm_calc()
edbm_dat.acid_channel()
edbm_dat.base_channel()
edbm_dat.salt_channel()
Cna_s.append(edbm_dat.Ci_s_out[0])

print("Ci_s_out: "+ str(edbm_dat.Ci_s_out))

#
while (edbm_dat.Ci_b_out[0]<1) and (edbm_dat.Ci_s_out[0]>=0.05):
    edbm_dat.recycl_below1M()
    edbm_dat.acid_channel()
    edbm_dat.base_channel()
    edbm_dat.salt_channel()
    Cna_s.append(edbm_dat.Ci_s_out[0])
    if edbm_dat.Ci_s_out[0]<0:
#        edbm_dat.Ci_s_out[0]=Cna_s[len(Cna_s)-2]
#        print("len(Cna_s)is " +str(len(Cna_s)))
#        print("Cna_s[len(Cna_s)-1] is " +str(Cna_s[len(Cna_s)-2]))
        break
    print("Ci_s_out: "+ str(edbm_dat.Ci_s_out))
    print("Ci_a_out: "+ str(edbm_dat.Ci_a_out))
    print("Ci_b_out: "+ str(edbm_dat.Ci_b_out))
print(" conc >=1m")    

#print("ms in is " + str(edbm_dat.M_s_out_t))
#print("Ci_s_out: "+ str(edbm_dat.Ci_s_out))
#print("Ci_a_out: "+ str(edbm_dat.Ci_a_out))
#print("Ci_b_out: "+ str(edbm_dat.Ci_b_out))

for i in range(1): 
    edbm_dat.recycl_1M()
    while (edbm_dat.Ci_b_out[0]<1) and (edbm_dat.Ci_s_out[0]>=0.05):
        edbm_dat.recycl_below1M()
        edbm_dat.acid_channel()
        edbm_dat.base_channel()
        edbm_dat.salt_channel()
        Cna_s.append(edbm_dat.Ci_s_out[0])
        
        if edbm_dat.Ci_s_out[0]<0:
    # #        edbm_dat.Ci_s_out[0]=Cna_s[len(Cna_s)-2]
    # #        print("len(Cna_s)is " +str(len(Cna_s)))
    # #        print("Cna_s[len(Cna_s)-1] is " +str(Cna_s[len(Cna_s)-2]))
              break
        print("Ci_s_out: "+ str(edbm_dat.Ci_s_out))
        print("Ci_a_out: "+ str(edbm_dat.Ci_a_out))
        print("Ci_b_out: "+ str(edbm_dat.Ci_b_out))
    print(" conc >=1m")    
    # edbm_dat.acid_channel()
    # edbm_dat.base_channel()
    # edbm_dat.salt_channel()
    # Cna_s.append(edbm_dat.Ci_s_out[0])
    # print("is working !!!!")
    print("Ci_s_out: "+ str(edbm_dat.Ci_s_out))
    print("Ci_a_out: "+ str(edbm_dat.Ci_a_out))
    print("Ci_b_out: "+ str(edbm_dat.Ci_b_out))

Cbrine_out_t=sum(edbm_dat.Ci_s_out)
Cb_out=edbm_dat.Ci_b_out[0:6]
Cb_out_g=[Cb_out[0]*MW_Na, Cb_out[1]*MW_cl, Cb_out[2]*MW_K, Cb_out[3]*MW_Mg, Cb_out[4]*MW_Ca, Cb_out[5]*MW_so4]

print("Ci_b_out: "+ str(Cb_out))
Ca_out=edbm_dat.Ci_a_out
Ca_out=edbm_dat.Ci_a_out[0:6]
Ca_out_g=[Ca_out[0]*MW_Na, Ca_out[1]*MW_cl, Ca_out[2]*MW_K, Ca_out[3]*MW_Mg, Ca_out[4]*MW_Ca, Ca_out[5]*MW_so4]
Cbrine_out=edbm_dat.Ci_s_out #mol/l
Cbrine_out_g=[Cbrine_out[0]*MW_Na, Cbrine_out[1]*MW_cl, Cbrine_out[2]*MW_K, Cbrine_out[3]*MW_Mg, Cbrine_out[4]*MW_Ca, Cbrine_out[5]*MW_so4]
# print("emf is " +str(edbm_dat.EMF))
print("M_h2o_s_out is "+ str(edbm_dat.M_h2o_s_out))
print("v_ext "+ str(edbm_dat.V_ext))
print("p is " +str(edbm_dat.P))

M_a_out=edbm_dat.M_a_out_f*N_trip
print("M_a_out is "+str(M_a_out)+"kg/hr")
M_b_out=edbm_dat.M_b_out_f*N_trip
print("M_b_out is "+str(M_b_out)+"kg/hr")
M_s_out=edbm_dat.M_s_out_f*N_trip
print("M_s_out is "+str(M_s_out)+"kg/hr")
M_s_out_r=edbm_dat.M_s_r*N_trip
print("M_s_out recycling is "+str(M_s_out_r)+"kg/hr")
M_a_out_r=edbm_dat.M_a_r*N_trip
print("M_a_out recycling is "+str(M_a_out_r)+"kg/hr")
M_b_out_r=edbm_dat.M_b_r*N_trip
print("M_b_out recycling is "+str(M_b_out_r)+"kg/hr")

print("edbm_dat.M_s_in_t is "+str(edbm_dat.M_s_in_t))
print("edbm_dat.M_a_in_t is "+str(edbm_dat.M_a_in_t))
print("edbm_dat.M_b_in_t is "+str(edbm_dat.M_b_in_t))
print("edbm_dat.M_a_r is "+str(edbm_dat.M_a_r))
print("edbm_dat.M_b_r is "+str(edbm_dat.M_b_r))
print("edbm_dat.M_s_r is "+str(edbm_dat.M_s_r))
Q_a_out=edbm_dat.Q_a_out*N_trip
print("Q_a_out is "+str(Q_a_out)+"l/hr")
M_HCl_out=Q_a_out*constants.MW_HCl/1000 #kg/hr
Q_b_out=edbm_dat.Q_b_out*N_trip
print("Q_b_out is "+str(Q_b_out)+"l/hr")
M_NaOH_out=Q_b_out*edbm_dat.Ci_b_out[0]*constants.MW_NaOH/1000 #kg/hr
Q_s_out=edbm_dat.Q_s_out*N_trip
print("Q_s_out is "+str(Q_s_out)+"l/hr")
#assumptions for pressure drop: dp=1.2 for recirculation, dp=0.1-0.2bar for feed, dp=1 for extracting products 
Ppump=(3*Q_in_edbm*0.1+Q_s_out*1+Q_b_out*1+Q_a_out*1+edbm_dat.Q_s_r*N_trip*1.2+edbm_dat.Q_a_r*N_trip*1.2+edbm_dat.Q_b_r*N_trip*1.2)/1000/3600*1e5/npump #units: W -> l to m3 so /1000; bar to J 1e5N/m2*1J/m ; hr to 3660s  
print("Ppump is "+str(Ppump)+" W")
R_int=45 #ohm cm2
P_t=edbm_dat.P/1000+Ppump/1000
print("Total electrical consumption for EDBM is " + str(P_t)+ " KW")
E_el_Edbm=P_t
print("Q_in_edbm is "+str(Q_in_edbm))
print("Ntrip is "+str(N_trip))
SEC=(48*I_d*A)/(Q_b_out*(edbm_dat.Ci_b_out[0]-edbm_dat.Ci_b_in[0])*constants.MW_NaOH)
E_el_Edbm=48*I_d*A/1000+Ppump/1000
print("Total electrical consumption for EDBM is " + str(E_el_Edbm)+ " KW")
#E_el_Edbm=3*M_NaOH_out
print("sec is "+str(SEC)+"kwh/kg naoh")
Q_w_in=2*Q_in_edbm
bal=(Q_in_edbm*d_in+2*Q_in_edbm)-M_a_out-M_s_out-M_b_out#edbm_dat.M_a_t_out*N_trip-edbm_dat.M_b_t_out*N_trip-edbm_dat.M_s_t_out*N_trip+M_s_out_r+M_a_out_r+M_b_out_r
print("balance edbm is "+str(bal))

#%% sec based on calogero paper
#SEC= U*Id*A/(60*Qp*(Cp,out-Cp,in)*Mp) where -> U: electrical potential, Id: current density, Cp: product concentration, Mp: molar mass of desired product (g/mol)
SEC2=(48*I_d*A)/(Q_b_out*(edbm_dat.Ci_b_out[0]-edbm_dat.Ci_b_in[0])*constants.MW_NaOH)
print("sec2 is " +str(SEC2)) #*0.663
SEC3=F*48/(36*constants.MW_NaOH*66.3*N_trip)
print("sec3 is " +str(SEC3)) #*0.663
#%%
sum_table_C={'F1: Cs_in ->efc effluent': Cin_edbm,
           # 'F2: Cs_in_t': , 
           # 'F3: Cs_out_t': ,
           'F4: Cs_out': Cbrine_out_g, 
           # 'F5: Cs_r': ,
           'F6: Ca_in water': [0,0,0,0,0,0],
           # 'F7: Ca_in_t': , 
           # 'F8: Ca_out_t': ,
           'F9: Ca_out': Ca_out_g, 
           # 'F10: Ca_r': ,
           'F11: Cb_in water': [0,0,0,0,0,0],
           # 'F12: Cb_in_t': , 
           # 'F13: Cb_out_t': ,
           'F14: Cb_out': Cb_out_g, 
           # 'F15: Cb_r': ,

           }
sum_table_C=pd.DataFrame(sum_table_C)
print(sum_table_C)

sum_table_f={'F1: Qs_in->efc effluent': Q_in_edbm,
           # 'F2: Qs_in_t': , 
           # 'F3: Qs_out_t': ,
           'F4: Qs_out': Q_s_out, 
           'F5: Qs_r': edbm_dat.Q_s_r*N_trip,
           'F6: Qa_in water': Q_in_edbm,
           # 'F7: Qa_in_t': , 
           # 'F8: Qa_out_t': ,
           'F9: Qa_out': Q_a_out, 
            'F10: Qa_r': edbm_dat.Q_s_r*N_trip,
           'F11: Qb_in water': Q_in_edbm ,
           # 'F12: Qb_in_t': , 
           # 'F13: Qb_out_t': ,
           'F14: Qb_out': Q_b_out, 
           'F15: Qb_r': edbm_dat.Q_s_r*N_trip,
           }
sum_table_f=pd.DataFrame(sum_table_f, index=['flow rate'])
print(sum_table_f)
sum_table_M={'F1: Ms_in->efc effluent': efc_unit.Mout,
           #F2: Ms_in_t': , 
           #'F3: Ms_out_t': ,
           'F4: Ms_out': M_s_out, 
           'F5: Ms_r': M_s_out_r,
           'F6: Ma_in water': Q_in_edbm,
           #'F7: Ma_in_t': , 
           #'F8: Ma_out_t': ,
           'F9: Ma_out': M_a_out, 
           'F10: Ma_r': M_a_out_r,
           'F11: Mb_in water': Q_in_edbm,
           #'F12: Mb_in_t': , 
           #'F13: Mb_out_t': ,
           'F14: Mb_out': M_b_out, 
           'F15: Mb_r': M_b_out_r,
           }
sum_table_M=pd.DataFrame(sum_table_M, index=['Mass flow rate'])
print(sum_table_M)
with pd.ExcelWriter('C:\\Users\\rodoulaktori\\surfdrive\\PhD\\Process_modelling\\python\\edbm_results.xlsx') as writer:
    sum_table_f.to_excel(writer,sheet_name="water_mining_scenario")
    sum_table_M.to_excel(writer,sheet_name="water_mining_scenario", startcol=0, startrow=3, header=True)
    sum_table_C.to_excel(writer,sheet_name="water_mining_scenario", startcol=0, startrow=10, header=True)
