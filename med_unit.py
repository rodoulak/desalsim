import math
import nanofiltration_unit
import efc_conc_step
import density_calc
#%%
#conditions
# T=20+273.15
#d=1018.494
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


#assumptions:

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
    
#%%calculations 
class med_calc:
    def __init__(self, Qf, Cc1, Cc2, Cc3, Cc4, Cc5, Cc6):
        self.Qf=Qf
        self.CNa_in=Cc1 #Na concentration (g/l) 
        self.CCl_in=Cc2 #Cl concentration (g/l) 
        self.CK_in=Cc3 #K concentration (g/l) 
        self.CMg_in=Cc4 #Mg concentration (g/l) 
        self.CCa_in=Cc5 #Ca concentration (g/l) 
        self.CSO4_in=Cc6 #SO4 concentration (g/l) 
        #self.CHCO3_in=Cc6 #HCO3 concentration (g/l) 
        self.cons=[self.CNa_in,self.CCl_in,self.CK_in,self.CMg_in,self.CCa_in,self.CSO4_in]
        print("in conc is "+str(self.cons))

    def sal_calc(self):
        self.salinity_in=sum(self.cons)
        self.xf=self.salinity_in/d*1000
        print("xf is "+str(self.xf))
        self.Mf=self.Qf/3600 #kg/s
        
        
    def mass_bal_med(self):
        self.xn=200/d_b*1000 #g/l
        self.conc_f=self.xn/self.xf #[-]
        print("concentration factor is "+str(self.conc_f))
        self.Bn=self.xf*self.Mf/self.xn #kg/s
        self.Qb=self.Bn*3600 #kg/hr
        self.Mdist=self.Mf-self.Bn #kg/s
        self.Qdist=self.Mdist*3600 #kg/hr
        self.D1=self.Mdist/(1+self.lhv1/self.lhv2) #kg/s
        self.D2=self.D1*self.lhv1/self.lhv2 #kg/s
        print("Mf: "+str(self.Mf), "xn: "+ str(self.xn), "Bn: "+ str(self.Bn), "Mdist: "+ str(self.Mdist), "d1: "+ str(self.D1), "d2: "+ str(self.D2) )
        
    def temp_calc(self):
        self.DT_t=T_s-T_N #oC
        self.U1= 1.9695+0.012057*T3-0.000085989*T3**2+0.00000025651*T3**3
        self.U2=0.95*self.U1
        self.DT1=self.DT_t/(self.U1*(1/self.U1+1/self.U2)) #oC
        self.DT2=self.DT1*self.U1/self.U2 #oC
        self.T1=T_s-self.DT1 #oC
        self.T2=self.T1-self.DT2 #oC
        self.lhv1=2499.5698-0.9156*(self.T1-DT_loss)-0.048343*(self.T1-DT_loss)**2
        self.lhv2=2626.1
    
    def performance_par(self):
        self.Ms=self.D1*self.lhv1/lh_s #kg/s
        self.GOR=self.Mdist/self.Ms
        self.Qc=self.D2*self.lhv2 #Condenser thermal load (Qc) kj/s ot kw
        self.Q1= self.Ms*lh_s#thermal load in the first effect (Q1) kw
        print("Q1 is "+str(self.Q1))
        self.QCW= self.D2*self.lhv2*1000/(Cp_w*(T_in-T_cw_in))-self.Mf#Cooling water flow rate kg/s
        self.PR= self.Mdist/self.Mf#Performance ratio (PR) 
        self.Qsen=self.Mf*cp_sol*(T_in-25)
        self.Q_Tot=self.Q1+self.Qsen/1000

    def out_conc(self):
        self.CNa_out=self.CNa_in*self.Mf/(d/1000)/(self.Bn/(d_b/1000)) #Na concentration (g/l) 
        print("self.CNa_out is "+str(self.CNa_out))
        self.CCl_out=self.CCl_in*self.Mf/(d/1000)/(self.Bn/(d_b/1000))  #Cl concentration (g/l) 
        self.CK_out=self.CK_in*self.Mf/(d/1000)/(self.Bn/(d_b/1000))  #K concentration (g/l) 
        self.CMg_out=self.CMg_in*self.Mf/(d/1000)/(self.Bn/(d_b/1000))  #Mg concentration (g/l) 
        self.CCa_out=self.CCa_in*self.Mf/(d/1000)/(self.Bn/(d_b/1000)) #Ca concentration (g/l) 
        self.CSO4_out=self.CSO4_in*self.Mf/(d/1000)/(self.Bn/(d_b/1000))  #SO4 concentration (g/l)
        self.Na_w=(self.CNa_in*self.Mf/(d/1000)-self.CNa_out*self.Bn/(d_b/1000))
        self.Cl_w=(self.CCl_in*self.Mf/(d/1000)-self.CCl_out*self.Bn/(d_b/1000))
        self.K_w=(self.CK_in*self.Mf/(d/1000)-self.CK_out*self.Bn/(d_b/1000))
        print("self.Na_w is "+str(self.Na_w))
        print("self.cl_w is "+str(self.Cl_w))        
        self.Mg_w=(self.CMg_in*self.Mf/(d/1000)-self.CMg_out*self.Bn/(d_b/1000))
        self.Ca_w=(self.CCa_in*self.Mf/(d/1000)-self.CCa_out*self.Bn/(d_b/1000))
        print("self.mgw is "+str(self.Mg_w))
        print("self.ca_w is "+str(self.Ca_w))
        self.so4_w=(self.CSO4_in*self.Mf/(d/1000)-self.CSO4_out*self.Bn/(d_b/1000))
        print("self.mgw is "+str(self.so4_w))
#%%main 
#input variables:
Qf_med=nanofiltration_unit.Qperm+efc_conc_step.Qperm
Cin_med=[]
#Cin_med=sum([Comp1.Cpermi, Comp2.Cpermi,Comp3.Cpermi, Comp4.Cpermi, Comp5.Cpermi, Comp6.Cpermi])
for i in range(len(nanofiltration_unit.Ci_out_p)):
    Cin_med.append(nanofiltration_unit.Ci_out_p[i]*nanofiltration_unit.Qperm/Qf_med+efc_conc_step.Cperm[i]*efc_conc_step.Qperm/Qf_med)


Cin_med_t=sum(Cin_med)
d=density_calc.density_calc(20, Cin_med_t)
d_b=density_calc.density_calc(45, 200)
#Qf_med=Qf_med*d
print("Cin_ med is " + str(Cin_med))
med_dat=med_calc(Qf_med,Cin_med[0],Cin_med[1],Cin_med[2],Cin_med[3],Cin_med[4],Cin_med[5])
#med_dat=med_calc(Qf_med, Comp1.Cpermi, Comp2.Cpermi,Comp3.Cpermi, Comp4.Cpermi, Comp5.Cpermi, Comp6.Cpermi)
med_dat.sal_calc()
med_dat.temp_calc()
med_dat.mass_bal_med()
med_dat.performance_par()
med_dat.out_conc()
#sum results 
Cconc_med=[med_dat.CNa_out, med_dat.CCl_out, med_dat.CK_out,med_dat.CMg_out,med_dat.CCa_out,med_dat.CSO4_out]
d_b=density_calc.density_calc(45, sum(Cconc_med))
print(sum(Cconc_med))
Cconc_med_t=sum(Cconc_med)
Qout_med=med_dat.Qb
Qprod_med=med_dat.Qdist
Qr=5.5*Qf_med
Cw_out=[med_dat.Na_w,med_dat.Cl_w,med_dat.K_w,med_dat.Mg_w,med_dat.Ca_w,med_dat.so4_w]

#mass balances 
for i in range(len(Cconc_med)):
    bal_i=Cin_med[i]*Qf_med/(d/1000)-(Qout_med/(d_b/1000)*Cconc_med[i]+Qprod_med/1*Cw_out[i])
    print(" for "+ str(i)+ " mass balance equals: " + str(bal_i))

#energy
#assumptions for dp: dp=2 bar for cw 
#dp for feed is the same as the recirculation: 3.5bar 
#dp for distillate water 1bar 
E_el_med=((Qf_med*3.5+med_dat.QCW*3600*2+(Qr+med_dat.Qb)*3.5+med_dat.Qdist*1)/(1000*npump))*1e5/3600/1000 #kwh
print("E_el_med is " +str(E_el_med)+ " kWh")
SEC_el=E_el_med/(Qf_med/d) #kwh/m3
print("SEC_el is " +str(SEC_el)+ " kWh/m3 feed")
SEC_el_prod=E_el_med/(Qprod_med/1000) #kwh/m3
print("SEC_el per product is " +str(SEC_el_prod)+ " kWh/m3 dist water")
Qcw=med_dat.QCW*3600
print("Cooling water flow rate is "+str(Qcw)+" kg/hr")
E_th_med=med_dat.Q_Tot
print("Total thermal energy consumption is "+str(E_th_med)+" KW")
SEC_th=E_th_med/(Qf_med/d) #kwh/m3
print("SEC_th is " +str(SEC_th)+ " kWh_th/m3")
Cchem=0

print(E_el_med/(Qprod_med/1000))