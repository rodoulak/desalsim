import math
import density_calc
import med_unit
import scaleup
#%%constants
R=8.31446261815324 #gas constant
MW_Na=22.99
MW_cl=35.453
MW_so4=32.064+4*15.999
MW_K=39.102
MW_Ca=40.08
MW_Mg=24.31
MW_HCO3=1.008+12.011+3*15.999
MW_NaCl=58.442729
MW_KCl=MW_K+MW_cl
MW_Mgso4=MW_Mg+MW_so4
MW_Caso4=MW_Ca+MW_so4
MW_Na2so4=2*MW_Na+MW_so4
MW_MgOH=58.3197

CaSO4_sat=2800  #ppm
Cp_f=3.14 # units: KJ* Kg*oC
CP_cw=4.187 # units: KJ* Kg*oC
UA=45990
salt_mois=20
LHV_v=2107.92 # kj/kg (gathered from steam table)
LHV_s=2357.69 # kj/kg (gathered from steam table)
dp=0.1
dp_slurry=1
npump=0.8
#assumptions 
#%%Calculations
class thermal_calc:
    def __init__(self, T_op, Qf, Cf_s, Cf_caso4):
        #mass balance 
        self.T_op= T_op #oC
        self.Qf=Qf #kg/h
        self.Cf_s=Cf_s #g/l
        self.Cf_caso4=Cf_caso4
        
    def mass_bal_cryst(self):
        self.nacl_sat_wt=0.136*(self.T_op*9/5*32)*25.009/100 #wt%
        self.nacl_sol=self.Qf*Cf_in[0]/(d_sol)*MW_NaCl/MW_Na/1000 #kg/hr
        self.salt_solids=self.Qf*self.Cf_s/d_sol/1000 #kg/h
        self.caso4_sol=self.Qf*self.Cf_caso4/d_sol/1000 #kg/h
        self.solid_mass=(self.salt_solids)/(1-salt_mois/100) #kg/h
        self.ev_mass=self.Qf-self.solid_mass #kg/h

    def heat_bal_cryst(self):
        #self.BPE=0.0127*(self.nacl_sat_wt/10)**2+0.1743*(self.nacl_sat_wt/10)+(0.3443)*5/9
        self.BPE=0.0127*(self.nacl_sat_wt/10)**2+0.1743*(self.nacl_sat_wt/10+0.3443)*5/9 #oC
        self.heat_req=LHV_v*self.ev_mass+self.Qf*Cp_f*(T_op-T_in) #kj/h
        self.Qt=self.heat_req/3600
        self.T_s=self.heat_req/UA+self.T_op #steam temperature oC
        self.steam_mass=self.heat_req/LHV_s #kg/hr
        self.cw_mass=(self.ev_mass*LHV_v)/(CP_cw*(T_cw_o-T_cw_f))

class conc_cal:
    def __init__(self, Qf, solid_mass, C1, Cc1, C2, Cc2, C3, Cc3, C4, Cc4, C5, Cc5, C6, Cc6):
        self.Qf=Qf
        self.solid_mass=solid_mass
        self.C1=C1
        self.C2=C2
        self.C3=C3
        self.C4=C4
        self.C5=C5
        self.C6=C6
        self.Cc1=Cc1
        self.Cc2=Cc2
        self.Cc3=Cc3
        self.Cc4=Cc4
        self.Cc5=Cc5
        self.Cc6=Cc6
    
    def molarity(self):        
        self.Cc1_mol=self.Cc1/MW_Na
        self.Cc2_mol=self.Cc2/MW_cl
        self.Cc3_mol=self.Cc3/MW_K
        self.Cc4_mol=self.Cc4/MW_Mg
        self.Cc5_mol=self.Cc5/MW_Ca
        self.Cc6_mol=self.Cc6/MW_so4
        print("mol cl-na-k"+str(self.Cc2_mol-self.Cc1_mol-self.Cc3_mol))
        print(self.Cc1_mol,self.Cc2_mol,self.Cc3_mol,self.Cc4_mol,self.Cc5_mol,self.Cc6_mol)
        if self.Cc1_mol>self.Cc2_mol: 
            self.CNacl_mol=self.Cc2_mol-self.Cc3_mol
            self.rem_na_mol=self.Cc1_mol-self.Cc2_mol
        else:
            self.CNacl_mol=self.Cc1_mol
            
        print("ncl mol is "+str(self.CNacl_mol))
        self.CKCl_mol=self.Cc3_mol
        #self.CMgS04_mol=self.Cc4
        self.CCaso4_mol=self.Cc6_mol
        #self.CCaoh_mol=self.Cc6_mol-self.Cc5_mol
    
    def conc_salt_comp(self):
        self.CNacl=self.CNacl_mol*MW_NaCl #g/l
        print("self.CNacl in g/l is "+str(self.CNacl))
        self.CKCl=self.CKCl_mol*MW_KCl
        #self.CMgSO4=self.CMgS04_mol*MW_Mgso4
        self.CCaso4=self.CCaso4_mol*MW_Caso4
                
    def salt_conc(self):
        self.CNacl_out=self.CNacl*self.Qf/d_sol/(self.solid_mass/(1974.46459880993/1000)) #salt density is 1974.46459880993kg/m3
        #self.CNacl_out=self.CNacl*self.Qf/self.solid_mass
#        self.CKCl_out=self.CKCl*self.Qf/self.solid_mass
#        self.CMgoh_out=self.CMgoh*self.Qf/self.solid_mass
#        self.CCaso4_out=self.CCaso4*self.Qf/self.solid_mass
        self.CKCl_out=self.CKCl*self.Qf/d_sol/(self.solid_mass*1000/1974.46459880993)
        #self.CMgSO4_out=self.CMgSO4*self.Qf/d_sol/(self.solid_mass*1000/1974.46459880993)
        self.CCaso4_out=self.CCaso4*self.Qf/d_sol/(self.solid_mass*1000/1974.46459880993)
        print('CNacl_out_out is '+str(self.CNacl_out))
        #print('CKCl_out is '+str(self.CKCl_out))
        #print('CMgoh_out is '+str(self.CMgSO4_out))
        print('CCaso4_outt is '+str(self.CCaso4_out))
        self.CNa=self.Cc1*self.Qf/d_sol/(self.solid_mass/1.9745)#self.CNacl_out/MW_NaCl*MW_Na+self.rem_na_mol*MW_Na
        self.CCl=self.CNacl_out/MW_NaCl*MW_cl+self.CKCl_out/MW_KCl*MW_cl
        print("self.CCl is "+str(self.CCl))
        self.CK=self.CKCl_out/MW_KCl*MW_K
        self.CCa=self.Cc5*self.Qf/d_sol/(self.solid_mass/1.9745)#self.CCaso4/MW_Caso4*MW_Ca
        self.CSO4=self.Cc6*self.Qf/d_sol/(self.solid_mass/1.9745)#self.CCaso4/MW_Caso4*MW_so4
        self.CMg=self.Cc4*self.Qf/d_sol/(self.solid_mass/1.9745) #density_salt=1974.5 kg/m3
        
#%%Main

#Input parameters
T_in=40 #oC
T_top=60 #oC
T_cw_f=25 #oC
T_cw_o=40 #oC
T_op=60
Qf=med_unit.Qout_med #kg/h
Cf_in=med_unit.Cconc_med
print("Cf_in is " + str(Cf_in))
Cf_s=med_unit.Cconc_med_t
print("cf_s is "+str(Cf_s))

Cf_s=med_unit.Cconc_med_t
d_sol=density_calc.density_calc(T_in, Cf_s)/1000
Cf_caso4=0.075132776
th_cryst_dat=thermal_calc(T_op, Qf, Cf_s, Cf_caso4)
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
M_Nacl=th_cryst_dat.solid_mass
th_cryst_dat_2=conc_cal(Qf, M_Nacl , 'Na',Cf_in[0], 'cl',Cf_in[1],'k', Cf_in[2], 'mg', Cf_in[3], 'ca', Cf_in[4], 'so4', Cf_in[5])
th_cryst_dat_2.molarity()
th_cryst_dat_2.conc_salt_comp()
th_cryst_dat_2.salt_conc()
Qcw=th_cryst_dat.cw_mass
Q_evap_mass=th_cryst_dat.ev_mass
print("Q_evap_mass is "+str(Q_evap_mass))
print("M_Nacl "+str(M_Nacl))
Csalt_out=[th_cryst_dat_2.CNa,th_cryst_dat_2.CCl, th_cryst_dat_2.CK, th_cryst_dat_2.CMg, th_cryst_dat_2.CCa, th_cryst_dat_2.CSO4]
E_el_th_Cr=((Qf*3.5+th_cryst_dat.ev_mass*1+th_cryst_dat.nacl_sat_wt*3.5+2*Qcw*2)*1e5/3600/(1000*npump))/1000 #kwh
print("E_el_th_Cr is "+str(E_el_th_Cr))
E_th_th_Cr=th_cryst_dat.heat_req/3600 #kwh
print("E_th_th_Cr is "+str(E_th_th_Cr))

E_fil=scaleup.scaleup(2.4, 4800, Qf)
E_el_th_Cr=E_el_th_Cr+E_fil
SEC_el=(E_fil+E_el_th_Cr)/(Qf*d_sol*1000)
SEC_el_prod=(E_fil+E_el_th_Cr)/(Q_evap_mass/1000)
print("SEC_el_prod is "+str(SEC_el_prod)+" KWh/m3 product ")
SEC_el_prod2=(E_fil+E_el_th_Cr)/(M_Nacl)
print("SEC_el_prod2 is "+str(SEC_el_prod2)+" KWh/kg product ")
Qchem=0
xp=E_th_th_Cr/Qf
print("xp is "+str(xp))

#mass balance 
for i in range(len(Cf_in)):
    bal_i=Cf_in[i]*Qf/d_sol-(M_Nacl/1.974465*Csalt_out[i])
    print(" for "+ str(i)+ " mass balance equals: " + str(bal_i))
print("Csalt_out na is "+str(Csalt_out[0]))
print("Csalt_in na is "+str(Cf_in[0]))
print("Qf is "+str(Qf))
    
print("nacl in sol "+str(th_cryst_dat.nacl_sol))
print("salt solid in sol "+str(th_cryst_dat.solid_mass))
over_bal=Qf-M_Nacl-Q_evap_mass
print(over_bal)