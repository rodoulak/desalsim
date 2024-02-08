#Pitzer_Model_Na_K_Cl_NO3_SO4

import math 
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
#import mfpfr_unit2
import efc_conc_step
import constants
import density_calc
import pdb
#lists
line_ice=[[],[],[],[]]
line_na2so4=[[],[],[],[]]
line_nacl=[[],[],[]]

#%%
#Definition of all salt ions and their charges

#ION_DATA =["Type", "Name", "Charge", "Cation", "Na", 1,  "Cation", "K",1, "Anion","Cl",-1,"Anion","NO3",-1,"Anion","SO4",-2]

#
# // Ion interaction data of ion pairs
ion_pair_data=[["Cation", "Anion", 1,-1, "Na", "Cl","<a0>(B0)","<a1>", "<a2>","<a3>", "<a4>","<a5>", 2.569010e-01, 0.0, -4.927911e-01, 0.0, 0.0, 0.0 ,"<a0>(B1)","<a1>", "<a2>","<a3>", "<a4>","<a5>", 5.657432e-01, -3.358675e+02, -1.878387e+00, 0.0, 0.0, 0.0,  "<a0>(B2)","<a1>", "<a2>","<a3>", "<a4>","<a5>", -1.091297e+00, 0.0, 9.013925e-01, 0.0, 0.0, 1.489073e-01,  "<a0>(B3)","<a1>", "<a2>","<a3>", "<a4>","<a5>",0.0, 0.0, 0.0, 0.0, 0.0, 0.0,"<a0>(C)","<a1>", "<a2>","<a3>", "<a4>","<a5>",  -4.726933e-03, 0.0 ,8.854208e-03, 0.0, 0.0, 0.0, "<alpha1>","<alpha2>","<alpha3>",1.4, 0.5, 0.0],
                ["Cation", "Anion", 1,-1, "Na", "NO3","<a0>(B0)","<a1>", "<a2>","<a3>", "<a4>","<a5>",4.227784e-01, -2.099681e+02, -3.306167e-01, 0.0, 0.0, -9.754582e-02, "<a0>(B1)","<a1>", "<a2>","<a3>", "<a4>","<a5>",7.254621e+00, -1.710939e+04, -8.730723e+01, 1.224729e-01, 0.0, -1.644950e+00,"<a0>(B2)","<a1>", "<a2>","<a3>", "<a4>","<a5>",  -3.950260e+00, 1.596136e+03, 2.670817e+00, 0.0, 0.0, 9.176415e-01,"<a0>(B3)","<a1>", "<a2>","<a3>", "<a4>","<a5>",  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,"<a0>(C)","<a1>", "<a2>","<a3>", "<a4>","<a5>", -5.566523e-03, 3.177426e+00, 5.569587e-03, 0.0, 0.0, 1.284590e-03, "<alpha1>","<alpha2>","<alpha3>",1.7, 0.8, 0.0],
                ["Cation", "Anion", 1,-2, "Na", "SO4","<a0>(B0)","<a1>", "<a2>","<a3>", "<a4>","<a5>",5.002718e-01, 7.362932e+03, 4.699443e+01, -7.941344e-02, 0.0, 0.0,"<a0>(B1)","<a1>", "<a2>","<a3>", "<a4>","<a5>",1.021684e+00, -4.221184e+03, -1.325564e+01, 0.0, 0.0, 0.0,"<a0>(B2)","<a1>", "<a2>","<a3>", "<a4>","<a5>",-2.118405e+00, -1.042977e+04, -6.910828e+01, 1.184556e-01, 0.0, 3.144602e-01,"<a0>(B3)","<a1>", "<a2>","<a3>", "<a4>","<a5>",  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,"<a0>(C)","<a1>", "<a2>","<a3>", "<a4>","<a5>",-5.061284e-03, 0.0, -4.425431e-02, 1.847906e-04, 0.0, 0.0,"<alpha1>","<alpha2>","<alpha3>",1.2, 0.2, 0.0],
                ["Cation", "Anion", 1,-1, "K", "Cl","<a0>(B0)","<a1>", "<a2>","<a3>", "<a4>","<a5>",5.087633e-02, -7.267703e+02, -3.572020e+00, 4.444364e-03, 0.0, 0.0,"<a0>(B1)","<a1>", "<a2>","<a3>", "<a4>","<a5>",-7.100622e+00, 2.834008e+04, 1.987679e+02, -4.908723e-01, 2.026021e-04, 1.701432e+00,"<a0>(B2)","<a1>", "<a2>","<a3>", "<a4>","<a5>",  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,"<a0>(B3)","<a1>", "<a2>","<a3>", "<a4>","<a5>",  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  "<a0>(C)","<a1>", "<a2>","<a3>", "<a4>","<a5>",-6.456211e-04, 5.774723e+01, 3.021599e-01, -4.001186e-04, 0.0, 0.0, "<alpha1>","<alpha2>","<alpha3>",2.0, 0.0, 0.0],
                ["Cation", "Anion", 1,-1, "K", "NO3",    "<a0>(B0)","<a1>", "<a2>","<a3>", "<a4>","<a5>",
-3.837176e-03, 0.0, 0.0 ,0.0, 0.0, 0.0,"<a0>(B1)","<a1>", "<a2>","<a3>", "<a4>","<a5>",
-2.893107e+00, 0.0 ,0.0, -5.143481e-03, 0.0, 7.044389e-01,"<a0>(B2)","<a1>", "<a2>","<a3>", "<a4>","<a5>",
-1.575356e-01, 0.0 ,5.093000e-01, 0.0 ,0.0, 0.0,"<a0>(B3)","<a1>", "<a2>","<a3>", "<a4>","<a5>",
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  "<a0>(C)","<a1>", "<a2>","<a3>", "<a4>","<a5>",
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,"<alpha1>","<alpha2>","<alpha3>",1.4 ,0.5, 0.0],
    ["Cation", "Anion",1,-2, "K", "SO4","<a0>(B0)","<a1>", "<a2>","<a3>", "<a4>","<a5>",
0.0, -6.080135e+02, -3.088000e+00, 4.265000e-03, 0.0,0.0,"<a0>(B1)","<a1>", "<a2>","<a3>", "<a4>","<a5>",
-2.358404e+01, 1.143312e+05, 8.818416e+02, -2.397456e+00, 1.099367e-03, 5.638375e+00,
"<a0>(B2)","<a1>", "<a2>","<a3>", "<a4>","<a5>",3.236876e-03, 0.0, 5.408952e-02, -2.342338e-04, 8.868560e-08, 0.0,
"<a0>(B3)","<a1>", "<a2>","<a3>", "<a4>","<a5>",  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  "<a0>(C)","<a1>", "<a2>","<a3>", "<a4>","<a5>",  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
"<alpha1>","<alpha2>","<alpha3>",1.4 ,0.0, 0.0]]
#
## 
## // Parameters for ion pairs of same charge type
                  #"Charge1","Charge2", "Cation/Anion","Cation/Anion", "<a0>(C)","<a1>", "<a2>","<a3>", "<a4>","<a5>",
#SAME_KIND_DATA 
same_kind_data=[[1, 1, "Na", "K", -2.161518e-02 , 0 ,2.581155e-01, -5.736134e-04,  0 , 0],[-1,-1, "Cl",  "NO3", -1.874000e-02, 0,0  ,0 , 0  ,0],[-1,-2, "Cl",  "SO4",  1.982079e-02, 5.905309e+01,  0 ,0,0 , 0],[-1,-2, "NO3", "SO4",  1.200000e-01, 0, 0,0 ,0,0]]
## 
## // Parameters for ion triples
#"CIApar", "<cation name>"," <middle ion name>"," <anion name>", " <a0>(?)", "<a1>", " <a2>"," <a3>", " <a4>",  "<a5>",
ion_triple_data =[["CIApar" ,"Na","K ", " Cl"  , 5.406252e-04  , 2.296865e+01 ,  6.346372e-02  , 0  ,0 ,  0],
                  ["CIApar", "Na" ,"K ",  "NO3 ", 1.226880e-03,   0           ,  -4.352200e-02,   1.154705e-04,   0,   0],
                  ["CIApar"," Na" ,"K"  , "SO4" , 5.277144e-03,   1.057085e+01,   0           ,   0           ,   0 ,  0],
                  [" CIApar", "Na" ,"Cl"  ,"NO3" , 4.046912e-04,  -1.229644e+01,  -6.148374e-02 ,  8.806542e-05,   0  , 0],
                  [" CIApar", "Na" ,"Cl" , "SO4" ,-2.409734e-03,  -1.067557e+02,  -2.806344e-01 ,  0           ,   0 ,  0],
                  [ "CIApar"," Na" ,"NO3" ,"SO4" ,-1.257539e-03,  -1.811048e+03,  -1.476757e+01 ,  3.924687e-02,  -1.698107e-05, 0],
                  ["CIApar"," K"  ,"Cl" , "NO3" , 1.758042e-03,   2.579013e+01,   7.915893e-02 ,  0           ,   0 ,  0],
                  ["CIApar", "K" , "Cl" , "SO4" , 6.185236e-03,  -1.123802e+03,  -6.475615e+00 ,  9.349162e-03,   0 ,  0],
                  ["CIApar", "K", " NO3" ,"SO4", -2.043612e-02,   4.438684e+03,   4.418141e+01 , -1.439547e-01,   7.715206e-05, 0]]
#
## // Parameters for the salts and hydrate phases
SALT_DATA_2comp =[["Ice", 1,"H", 2,"O",1,  "<a0>(K0)"," <a1>"," <a2>"," <a3>"," <a4>"," <a5>",  0.2385066,	1.498147e03,	1.165961e01,	-1.290866e-02,	0.0,		0.0],
                   ["NaCl", 0 ,"Na", 1," Cl", 1,  "<a0>(K0)"," <a1>"," <a2>"," <a3>"," <a4>"," <a5>",  -3.194744e-1,	0.0,	4.057681,	-2.085546e-2,	0.0,	9.243879e-1],
                   ["NaCl.2H2O", 2, "Na", 1 ,"Cl", 1,  "<a0>(K0)"," <a1>"," <a2>"," <a3>"," <a4>"," <a5>",  3.826504e00,		0.0,		5.798855e02,	-4.167374e00,	3.805646e-03, 0.0],
                   ["NaNO3", 0, "Na", 1, "NO3" ,1,  "<a0>(K0)"," <a1>"," <a2>"," <a3>"," <a4>"," <a5>",  2.503743e00,	-2.045335e04,	-1.050257e02,	1.495354e-01, 0.0, 0.0],
                   ["Na2SO4", 0, "Na", 2, "SO4", 1,  "<a0>(K0)"," <a1>"," <a2>"," <a3>"," <a4>"," <a5>",  -7.072836e-01,	-9.860330e03,	-3.397271e01,	0.0,	0.0	,0.0],
                   ["Na2SO4.10H2O", 10," Na", 2, "SO4", 1,  "<a0>(K0)"," <a1>"," <a2>"," <a3>"," <a4>"," <a5>",  -2.845545e00	,	0.0	,7.977607e01,	-1.610020e-01,	0.0,	0.0],
                   ["KCl", 0,", K", 1, "Cl", 1,  "<a0>(K0)"," <a1>"," <a2>"," <a3>"," <a4>"," <a5>",  2.081763e00 ,	7.552778e03,	1.830581e01, 0.0 ,0.0, 0.0],
                   ["KNO3", 0," K", 1, "NO3", 1,   "<a0>(K0)"," <a1>"," <a2>"," <a3>"," <a4>"," <a5>",-2.567651e-01,	-1.369397e05,	-1.150760e03,	3.320581e00,	-1.602324e-03, 0.0],
                   ["K2SO4", 0,"K", 2," SO4", 1,  "<a0>(K0)"," <a1>"," <a2>"," <a3>"," <a4>"," <a5>", -4.125062,	-2.495771e4,	-1.010588e2,	9.071389e-2,	0.0,	0.0]]

#" <full salt name>"," <nH2O>", " <ion name 1>"," <n1>"," <ion name 2>"," <n2>","<ion name 3>"," <n3>" ,
SALT_DATA_3comp =[["Na2K6(SO4)4", 0 ,"Na", 2, "K", 6, "SO4", 4,"<a0>(K0)"," <a1>"," <a2>"," <a3>"," <a4>"," <a5>",",",",", -1.736574e+01,	-6.378531e+04,	-1.798858e+02,	0.0,	0.0	,0.0,",",","],
                   ["Na3SO4NO3.H2O", 1," Na", 3," SO4", 1, "NO3", 1,"<a0>(K0)"," <a1>"," <a2>"," <a3>"," <a4>"," <a5>",",",",", 9.445820e-01,-1.427731e+04,	-3.845114e+01,	0.0,	0.0,	0.0,",",","]]

##Parameters ak
ak_I=[1.925154014814667,
-0.060076477753119,
-0.029779077456514,
-0.007299499690937,
+0.000388260636404,
+0.000636874599598,
+0.000036583601823,
-0.000045036975204,
-0.000004537895710,
+0.000002937706971,
+0.000000396566462,
-0.000000202099617,
-0.000000025267769,
+0.000000001352261,
+0.000000001229405,
-0.000000000821969,
-0.000000000050847,
+0.000000000046333,
+0.000000000001943,
-0.000000000002563,
-0.000000000010991]
ak_II=[+0.628023320520852,
+0.462762985338493,
+0.150044637187895,
-0.028796057604906,
-0.036552745910311,
-0.001668087945272,
+0.006519840398744,
+0.001130378079086,
-0.000887171310131,
-0.000242107641309,
+0.000087294451594,
+0.000034682122751,
-0.000004583768938,
-0.000003548684306,
-0.000000025045388,
+0.000000216991779,
+0.000000008077957,
+0.000000004558555,
-0.000000006944757,
-0.000000002849257,
+0.000000000237816]

#Equations
# Reference temperature for temperature polynomial
REF_TEMP_POLYNOM = 298.15
# Gas constant [J/mol/K]
GASCONST_R = 8.31451 
# Empirical parameter. value=1.2 (kg/mol)**1/2 at 25C
CONST_B = 1.2 
#%%
#Inputs
#Composition of the solution at 298.15 K:
# C_t=mfpfr_unit2.Cout_all
# print("C_t"+str(C_t))
# M_out=mfpfr_unit2.Mout_2

# Na    = [C_t[0]/mfpfr_unit2.d_out_s*1000,1]     # mol/kg
# K     = [C_t[2]/mfpfr_unit2.d_out_s*1000,1]   # mol/kg
# Cl    = [C_t[1]/mfpfr_unit2.d_out_s*1000,-1]    # mol/kg
# NO3   = [0,-1]          # mol/kg
# SO4   = [C_t[5]/mfpfr_unit2.d_out_s*1000,-2]   # mol/kg
na_conc_list=[]
cl_conc_list=[]
so4_conc=[]
M_H2O_list=[]
naCL_wt=[]
NA2so4_wt=[]
M_H2O_wt=[]
U_salt_all=[]
w_sol= 100#efc_conc_step.Qconc  #kg 
for cw in range(1,101):
    wt_na2_so4=cw #wt%
    Cna2so4=wt_na2_so4*w_sol/100 #kg
    c_na_1= Cna2so4*2*1000/constants.MW_Na2SO4 #mol
    c_so4= Cna2so4*1000/constants.MW_Na2SO4  #mol
    for cl_w in range(1,101-cw):
        wt_nacl=cl_w #wt%
        Cnacl=wt_nacl*w_sol/100 #kg
        c_na_2= Cnacl*1*1000/constants.MW_NaCl
        c_cl= Cnacl*1000/constants.MW_NaCl
        c_k=0
        M_H2O = 0.01810534#(100-wt_na2_so4-wt_nacl)*w_sol*1000/(100*18.01528*w_sol) #mol/kg
        Na    = [(c_na_1+c_na_2)/w_sol,1]     # mol/kg
        K     = [c_k,1]   # mol/kg
        Cl    = [c_cl/w_sol,-1]    # mol/kg
        NO3   = [0,-1]          # mol/kg
        SO4   = [c_so4/w_sol,-2]   # mol/kg
        so4_conc.append(SO4[0])
        na_conc_list.append(Na[0])
        cl_conc_list.append(Cl[0])
        M_H2O_list.append(M_H2O)
        naCL_wt.append(wt_nacl)
        NA2so4_wt.append(wt_na2_so4)
        M_H2O_wt.append(100-wt_na2_so4-wt_nacl)
        ION_DATA =[[],[],[],[],[]]
        ION_DATA =[["Cation", "Na", Na[0], 1],["Anion","Cl", Cl[0], -1],["Cation", "K",K[0], 1,],["Anion","NO3",NO3[0],-1],["Anion","SO4",SO4[0],-2]]
        T_ini=20
        T_final=-20
        T_interval=-2
        
        #%%
        CONST_B = 1.2 
        o=0
        U_salt=[[],[],[],[],[]]
        for T_C in range(T_ini,T_final,T_interval):
            o=o+1
            T=273.15+T_C
            T2=T**2
            T3=T**3
            sqT3=math.sqrt(T3)    
            ## Calculates the Debye-Hckel-Term A_phi (Equation 10)for equation 22 
            APhi=-0.817653 - 0.8685276/(T-222) + 19251/T2 + 5.251284e-3*T - 7.149397e-6*T2 + 9.338559e-12*T2*T2 
            ##Calculates the Debye-Hckel-Term A_L (Equation 23)
            A_L = 1297752 + 234787.2/(T - 222) - 2.700596E10/T2 + 6.1048097E7/T - 5.2177636E7/math.sqrt(T)- 0.01371691*math.sqrt(T2*T2*T) + 5.337301E-18*T2*T2*T2*T2    
            ##Calculates the Debye-Hckel-Term A_J (Equation 24)
            A_J=- 0.2347872E6/((T-222)*(T-222)) - 6.1048097E8/(T*T) + 5.401193E10/(T3) - 2.608882E7/sqT3- 3.429228E-2*sqT3 + 4.269841E-17*T3*T3*T
            Salt_list=[ [], [], [], [], [], [], [], [], []]
            ##Calculates the temperature polynomial of the ion parameters (Equation 25)
            for i in range(0,9):
                a0=np.array(SALT_DATA_2comp[i][12])
                a1=np.array(SALT_DATA_2comp[i][13])
                a2=np.array(SALT_DATA_2comp[i][14])
                a3=np.array(SALT_DATA_2comp[i][15])
                a4=np.array(SALT_DATA_2comp[i][16])
                a5=np.array(SALT_DATA_2comp[i][17])
                T_polynom= a0 + a1*(1/T - 1/REF_TEMP_POLYNOM)+ a2*math.log(T/REF_TEMP_POLYNOM)+ a3*(T - REF_TEMP_POLYNOM)+ a4*(T*T - REF_TEMP_POLYNOM*REF_TEMP_POLYNOM) + a5*math.log(T - 225)
        
                Salt_list[0].append(SALT_DATA_2comp[i][0])
                Salt_list[1].append(T_polynom) #activity
       
       
            for i in range(0,2):
                Salt_list[0].append(SALT_DATA_3comp[i][0])
                a0=np.array(SALT_DATA_3comp[i][16])
                a1=np.array(SALT_DATA_3comp[i][17])
                a2=np.array(SALT_DATA_3comp[i][18])
                a3=np.array(SALT_DATA_3comp[i][19])
                a4=np.array(SALT_DATA_3comp[i][20])
                a5=np.array(SALT_DATA_3comp[i][21])
                T_polynom= a0 + a1*(1/T - 1/REF_TEMP_POLYNOM)+ a2*math.log(T/REF_TEMP_POLYNOM) + a3*(T - REF_TEMP_POLYNOM)+ a4*(T*T - REF_TEMP_POLYNOM*REF_TEMP_POLYNOM) + a5*math.log(T - 225)
                Salt_list[1].append(T_polynom)
            #end
            b=1.2 # (Kg/mol)**1/2
            Z_=(Na[0]*abs(Na[1]))+(K[0]*abs(K[1]))+(Cl[0]*abs(Cl[1])) +(SO4[0]*abs(SO4[1]))
            I_=((Na[0]*abs(Na[1])**2)+(K[0]*abs(K[1])**2)+(Cl[0]*abs(Cl[1])**2) +(SO4[0]*abs(SO4[1])**2))/2    
            #ion pair data
            ion_pair=[[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
            for i in range (0,6):
                ion_pair[0].append("Ion pair"+"_"+ion_pair_data[i][4]+ion_pair_data[i][5])
                a0C=np.array(ion_pair_data[i][60])
                a1C=np.array(ion_pair_data[i][61])
                a2C=np.array(ion_pair_data[i][62])
                a3C=np.array(ion_pair_data[i][63])
                a4C=np.array(ion_pair_data[i][64])
                a5C=np.array(ion_pair_data[i][65])                                  
                para_C= a0C + (a1C*((1/T) - (1/REF_TEMP_POLYNOM)))+ (a2C*math.log(T/REF_TEMP_POLYNOM)) + (a3C*(T - REF_TEMP_POLYNOM))+ (a4C*((T*T) - (REF_TEMP_POLYNOM*REF_TEMP_POLYNOM)))+ (a5C*math.log(T - 225))                  
                ion_pair[1].append(para_C)
                  
                a0B0=np.array(ion_pair_data[i][12])
                a1B0=np.array(ion_pair_data[i][13])
                a2B0=np.array(ion_pair_data[i][14])
                a3B0=np.array(ion_pair_data[i][15])
                a4B0=np.array(ion_pair_data[i][16])
                a5B0=np.array(ion_pair_data[i][17])                                           
                beta_0= a0B0 + (a1B0*((1/T) - (1/REF_TEMP_POLYNOM)))+ (a2B0*math.log(T/REF_TEMP_POLYNOM))+ (a3B0*(T - REF_TEMP_POLYNOM))+ (a4B0*((T*T) - (REF_TEMP_POLYNOM*REF_TEMP_POLYNOM)))+ (a5B0*math.log(T - 225))
                ion_pair[2].append(beta_0)
                  
                a0B1=np.array(ion_pair_data[i][24])
                a1B1=np.array(ion_pair_data[i][25])
                a2B1=np.array(ion_pair_data[i][26])
                a3B1=np.array(ion_pair_data[i][27])
                a4B1=np.array(ion_pair_data[i][28])
                a5B1=np.array(ion_pair_data[i][29])
                                                  
                beta_1= a0B1 + (a1B1*((1/T) - (1/REF_TEMP_POLYNOM)))+ (a2B1*math.log(T/REF_TEMP_POLYNOM)) + (a3B1*(T - REF_TEMP_POLYNOM))+ (a4B1*((T*T) - (REF_TEMP_POLYNOM*REF_TEMP_POLYNOM)))+ (a5B1*math.log(T - 225))
                ion_pair[3].append(beta_1)
            
                a0B2=np.array(ion_pair_data[i][36])
                a1B2=np.array(ion_pair_data[i][37])
                a2B2=np.array(ion_pair_data[i][38])
                a3B2=np.array(ion_pair_data[i][39])
                a4B2=np.array(ion_pair_data[i][40])
                a5B2=np.array(ion_pair_data[i][41])
                            
                beta_2= a0B2 + (a1B2*((1/T) - (1/REF_TEMP_POLYNOM)))+ (a2B2*math.log(T/REF_TEMP_POLYNOM)) + (a3B2*(T - REF_TEMP_POLYNOM))+ (a4B2*((T*T) - (REF_TEMP_POLYNOM*REF_TEMP_POLYNOM)))+ (a5B2*math.log(T - 225))
                ion_pair[4].append(beta_2)         
            
                a0B3=np.array(ion_pair_data[i][48])
                a1B3=np.array(ion_pair_data[i][49])
                a2B3=np.array(ion_pair_data[i][50])
                a3B3=np.array(ion_pair_data[i][51])
                a4B3=np.array(ion_pair_data[i][52])
                a5B3=np.array(ion_pair_data[i][53])
                                                 
                beta_3= a0B3 + (a1B3*((1/T) - (1/REF_TEMP_POLYNOM)))+ (a2B3*math.log(T/REF_TEMP_POLYNOM))+ (a3B3*(T - REF_TEMP_POLYNOM))+ (a4B3*((T*T) - (REF_TEMP_POLYNOM*REF_TEMP_POLYNOM))) + (a5B3*math.log(T - 225))
             
                     
                ion_pair[5].append(beta_3)
                  
                  
                alpha_0=np.array(ion_pair_data[i][69])
                alpha_1=np.array(ion_pair_data[i][70])
                alpha_2=np.array(ion_pair_data[i][71])
                
                ion_pair[6].append(alpha_0)
                ion_pair[7].append(alpha_1)
                ion_pair[8].append(alpha_2)
        
        
                 
                z1=np.array(ion_pair_data[i][2])
                z2=np.array(ion_pair_data[i][3])
                          
                x1=float(ion_pair_data[i][69])*float(math.sqrt(I_))
                
                if x1==0.0:
                    g1=0
                else:
                    g1=(2/x1**2)*(1-((1+x1)*math.exp(-x1)))
           
                x2=float(ion_pair_data[i][70])*math.sqrt(I_) 
            
                if x2==0.0:
                    g2=0
                else:
                    g2=(2/x2**2)*(1-((1+x2)*math.exp(-x2)))    
                    
                x3=float(ion_pair_data[i][71])*math.sqrt(I_) 
            
                if x3==0.0: 
                    g3=0
                else:
                    g3=(2/x3**2)*(1-((1+x3)*math.exp(-x3)))
            
            
                if x1==0.0:
                    dg1=0
                else:
                    dg1=-2*(1-((1+x1+(x1**2)/2)*math.exp(-x1)))/(x1**2)
            
                if x2==0.0: 
                
                    dg2=0
                else:
                    dg2=-2*(1-((1+x2+(x2**2)/2)*math.exp(-x2)))/(x2**2)
                
                if x3==0.0:
                
                    dg3=0
                else:
                    dg3=-2*(1-((1+x3+(x3**2)/2)*math.exp(-x3)))/(x3**2)
        
            
                para_B_Phi=beta_0+(beta_1*math.exp(-alpha_0*math.sqrt(I_)))+ (beta_2*math.exp(-alpha_1*math.sqrt(I_)))+ (beta_3*math.exp(-alpha_2*math.sqrt(I_)))
            
                para_B= beta_0 +(beta_1*g1)+(beta_2*g2)+(beta_3*g3)
                                 
                para_dB = ((beta_1*dg1) +(beta_2*dg2)+(beta_3*dg3))/I_
                     
                ion_pair[9].append(para_B_Phi)
                ion_pair[10].append(para_B)
                ion_pair[11].append(para_dB)
                
            #end 
        
            #NaCl
            ion_pair[12].append(Na[0]) #cationmol
            ion_pair[13].append(Cl[0]) #anionmol
              #NaNO3
            ion_pair[12].append(Na[0]) 
            ion_pair[13].append(NO3[0])
              #NaSO4
            ion_pair[12].append(Na[0]) 
            ion_pair[13].append(SO4[0])
              #KCl
            ion_pair[12].append(K[0]) 
            ion_pair[13].append(Cl[0])
              #KNO3
            ion_pair[12].append(K[0]) 
            ion_pair[13].append(NO3[0])
              #KSO4
            ion_pair[12].append(K[0]) 
            ion_pair[13].append(SO4[0])
        
          #ion triple caculates Psi
            ion_triple=[[],[],[],[],[],[],[]]
            for i in range(0,9):
                ion_triple[0].append("Ion triple"+"_"+ion_triple_data[i][1]+ion_triple_data[i][2]+ion_triple_data[i][3])
                a0=np.array(ion_triple_data[i][4])
                a1=np.array(ion_triple_data[i][5])
                a2=np.array(ion_triple_data[i][6])
                a3=np.array(ion_triple_data[i][7])
                a4=np.array(ion_triple_data[i][8])
                a5=np.array(ion_triple_data[i][9])          
                para_Psi= a0+(a1*((1/T)-(1/REF_TEMP_POLYNOM)))+(a2*math.log(T/REF_TEMP_POLYNOM))+ (a3*(T - REF_TEMP_POLYNOM))+ (a4*((T*T) - (REF_TEMP_POLYNOM*REF_TEMP_POLYNOM)))+ (a5*math.log(T - 225))     
                ion_triple[1].append(para_Psi)
        
        #ion pair same kind
            ion_pair_same_kind=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
            for i in range(0,4):
                ion_pair_same_kind[0].append("Ion same kind"+"_"+same_kind_data[i][2]+same_kind_data[i][3])
                a0=np.array(same_kind_data[i][4])
                a1=np.array(same_kind_data[i][5])
                a2=np.array(same_kind_data[i][6])
                a3=np.array(same_kind_data[i][7])
                a4=np.array(same_kind_data[i][8])
                a5=np.array(same_kind_data[i][9])
                
                z1=np.array(same_kind_data[i][0])
                z2=np.array(same_kind_data[i][1])
                
                x = 6*APhi*math.sqrt(I_)
                x_ij = z1*z2*x
                x_ii = z1*z1*x
                x_jj = z2*z2*x
                xc=[x_ij,x_ii,x_jj]
               
                para_theta= a0+(a1*((1/T)-(1/REF_TEMP_POLYNOM)))+(a2*math.log(T/REF_TEMP_POLYNOM))+ (a3*(T - REF_TEMP_POLYNOM))+ (a4*((T*T) - (REF_TEMP_POLYNOM*REF_TEMP_POLYNOM)))+ (a5*math.log(T - 225))
         
                if np.array(same_kind_data[i][0])==np.array(same_kind_data[i][1]):
                    ion_pair_same_kind[1].append(para_theta)
                    ion_pair_same_kind[2].append(0) #E_theta=0
                    ion_pair_same_kind[3].append(0) #E_d_theta=0
                    ion_pair_same_kind[4].append(para_theta) #PhiPhi= para_theta
                    ion_pair_same_kind[5].append(para_theta) #Phi= para_theta
                    ion_pair_same_kind[6].append(0) #dPhi=0
                    ion_pair_same_kind[7].append([])
                    ion_pair_same_kind[8].append([])
                    ion_pair_same_kind[9].append([])
                else:
                    b_k = 0 
                    b_k1 = 0 
                    b_k2=0
                    d_k = 0 
                    d_k1 = 0 
                    d_k2=0
                   
                    if x_ij<1:
                        z = (4*x_ij**0.2)-2
                        dzdx = (4*x_ij**(-0.8))/5
                        for k in range(0, 21,-1):##=21:-1:1
                            ak_Ix=ak_I(k)
                            ak_IIx=ak_II(k)
                       
                        #transfer variables from previous iteration
                            b_k2 = b_k1 
                            b_k1 = b_k    
                            d_k2 = d_k1 
                            d_k1 = d_k  
            #            #calculate new iteration step
            
                            b_k = z*b_k1 - b_k2 + ak_I(k)
                            d_k = b_k1 + z*d_k1 - d_k2
                    else:
                        z = (((40*x_ij**(-0.1))-22))/9
                        dzdx = (-4*x_ij**(-1.1))/9
                       
                        for k in range(0, 21,-1):##=21:-1:1
                            ak_Ix=ak_I(k)
                            ak_IIx=ak_II(k)
                                     
                        #transfer variables from previous iteration
                            b_k2 = b_k1 
                            b_k1 = b_k    
                            d_k2 = d_k1 
                            d_k1 = d_k  
                        # calculate new iteration step
                            b_k = (z*b_k1)-b_k2+ak_II(k)
                            d_k = b_k1+(z*d_k1)-d_k2
                        
                        J_ij  = (x_ij/4)-1+((b_k-b_k2))/2
                        dJ_ij = 0.25+dzdx*(d_k-d_k2)              
                        ion_pair_same_kind[7].append([x_ij,J_ij, dJ_ij])
            
                          
                    b_k = 0 
                    b_k1 = 0 
                    b_k2=0
                    d_k = 0 
                    d_k1 = 0 
                    d_k2=0
                   
                    if x_ii<1:
                        z = (4*x_ii**0.2)-2
                        dzdx = (4*x_ii**(-0.8))/5
                        for k in range(0, 21,-1):
                            ak_Ix=ak_I(k)
                            ak_IIx=ak_II(k)
              # transfer variables from previous iteration
                            b_k2 = b_k1 
                            b_k1 = b_k    
                            d_k2 = d_k1 
                            d_k1 = d_k 
                        #calculate new iteration step
                            b_k = z*b_k1 - b_k2 + ak_I(k)
                            d_k = b_k1 + z*d_k1 - d_k2 
                    else:
                        for k in range(0, 21,-1):
                            ak_Ix=ak_I(k)
                            ak_IIx=ak_II(k)
                            z = ((40*x_ii**(-0.1)-22))/9
                            dzdx = (-4*x_ii**(-1.1))/9      
                          #transfer variables from previous iteration
                            b_k2 = b_k1 
                            b_k1 = b_k    
                            d_k2 = d_k1 
                            d_k1 = d_k 
                                # calculate new iteration step
                            b_k = z*b_k1 - b_k2 + ak_II(k)
                            d_k = b_k1 + z*d_k1 - d_k2
                      
                        J_ii  = (x_ii/4)-1+((b_k-b_k2))/2
                        dJ_ii = 0.25+dzdx*(d_k-d_k2)
                   
                        ion_pair_same_kind[8].append([x_ii,J_ii, dJ_ii])
                   
                    b_k = 0 
                    b_k1 = 0 
                    b_k2=0
                    d_k = 0 
                    d_k1 = 0 
                    d_k2=0
                        
                    if x_jj<1:
                        z = (4*x_jj**0.2)-2
                        dzdx = (4*x_jj**(-0.8))/5
                        for k in range(0, 21,-1):##=21:-1:1
                            ak_Ix=ak_I(k)
                            ak_IIx=ak_II(k)
        # transfer variables from previous iteration
                            b_k2 = b_k1 
                            b_k1 = b_k    
                            d_k2 = d_k1 
                            d_k1 = d_k 
                    #calculate new iteration step
                            b_k = z*b_k1 - b_k2 + ak_I(k)
                            d_k = b_k1 + z*d_k1 - d_k2 
        
                    else:
                        for k in range(0, 21,-1):##=21:-1:1
                            ak_Ix=ak_I(k)
                            ak_IIx=ak_II(k)
                            z = ((40*x_jj**(-0.1)-22))/9
                            dzdx = (-4*x_jj**(-1.1))/9        
                    #transfer variables from previous iteration
                            b_k2 = b_k1 
                            b_k1 = b_k    
                            d_k2 = d_k1 
                            d_k1 = d_k 
                    # calculate new iteration step
                            b_k = z*b_k1 - b_k2 + ak_II(k)
                            d_k = b_k1 + z*d_k1 - d_k2
                    
                        J_jj  = (x_jj/4)-1+((b_k-b_k2))/2
                        dJ_jj = 0.25+dzdx*(d_k-d_k2)
                   
                        ion_pair_same_kind[9].append([x_jj,J_jj, dJ_jj])
    
                    para_E_theta = z1*z2/4.0/I_*(J_ij - J_ii/2 - J_jj/2)
                    #print('para_E_theta is ' + str(para_E_theta))
             
                    para_E_d_theta = -para_E_theta/I_ + z1*z2/8.0/(I_*I_)*(x_ij*dJ_ij - x_ii*dJ_ii/2 - x_jj*dJ_jj/2)
                     
                    para_PhiPhi = para_theta + para_E_theta + I_*para_E_d_theta
            
                    para_Phi=para_theta + para_E_theta
             
                    para_dPhi=para_E_d_theta
    
                    ion_pair_same_kind[1].append(para_theta)
                    ion_pair_same_kind[2].append(para_E_theta) #E_theta=para_E_theta
                    ion_pair_same_kind[3].append(para_E_d_theta) #E_d_theta=para_E_d_theta
                    ion_pair_same_kind[4].append(para_PhiPhi) #PhiPhi= para_PhiPhi
                    ion_pair_same_kind[5].append(para_Phi) #Phi= para_Phi 
                    ion_pair_same_kind[6].append(para_dPhi) #dPhi=para_dPhii
        
        
            ion_pair_same_kind[10].append(Na[0]*K[0]) #molprod  
            ion_pair_same_kind[11].append(Na[0]) #molcation1
            ion_pair_same_kind[12].append(K[0]) #molcation2
            ion_pair_same_kind[13].append([])
            ion_pair_same_kind[14].append([])
            
            ion_pair_same_kind[10].append(Cl[0]*NO3[0])
            ion_pair_same_kind[11].append([])
            ion_pair_same_kind[12].append([])
            ion_pair_same_kind[13].append(Cl[0])
            ion_pair_same_kind[14].append(NO3[0])
            
            ion_pair_same_kind[10].append(Cl[0]*SO4[0])
            ion_pair_same_kind[11].append([])
            ion_pair_same_kind[12].append([])
            ion_pair_same_kind[13].append(Cl[0])
            ion_pair_same_kind[14].append(SO4[0])
            
            ion_pair_same_kind[10].append(NO3[0]*SO4[0])
            ion_pair_same_kind[11].append([])
            ion_pair_same_kind[12].append([])
            ion_pair_same_kind[13].append(NO3[0])
            ion_pair_same_kind[14].append(SO4[0])
        
        ##Term1 constant
        #
            term1 = -APhi*math.sqrt(I_)*math.sqrt(I_)*math.sqrt(I_)/(1 + CONST_B*math.sqrt(I_))
        #
        ##Term2 Pairs
        
            term2f=0
            Osmoticterms=[[],[],[],[],[],[],[]]  
            for i in range(0,6):
                term2= (ion_pair[12][i]*ion_pair[13][i])*(ion_pair[9][i]+(Z_*ion_pair[1][i]))
                term2f=term2f+term2
        
            Osmoticterms[0].append(term1)
            Osmoticterms[1].append(term2f) #!!!!!!!!!!Change!!!!!!!!0.311447#
        
        ##Term3 triple ions
        #
        ##NaKCl
            ion_triple[2].append(Na[0]*K[0])   #cationcationmolprod
            ion_triple[3].append(ion_pair_same_kind[4][0]) #paraPhiPhi
            ion_triple[4].append(Cl[0]) #anion
            Osmoticterms[2].append((ion_triple[4][0])*ion_triple[1][0]) #anionmolPsiprod
            Osmoticterms[3].append([])
            ion_triple[6].append([])
            ion_triple[5].append([])   
        #
        ##NaKNO3
            ion_triple[2].append(Na[0]*K[0])
            ion_triple[3].append(ion_pair_same_kind[4][0])
            ion_triple[4].append(NO3[0])
            Osmoticterms[2].append((ion_triple[4][1])*ion_triple[1][1]) #Osmoticterms(2).anionmolPsiprod=ion_triple(2).anion*ion_triple(2).Psi
            Osmoticterms[3].append([])
            ion_triple[6].append([])
            ion_triple[5].append([])    
        ##NaKSO4
            ion_triple[2].append(Na[0]*K[0])
            ion_triple[3].append(ion_pair_same_kind[4][0])
            ion_triple[4].append(SO4[0])
            Osmoticterms[2].append((ion_triple[4][2])*ion_triple[1][2])
            Osmoticterms[3].append([])
            ion_triple[6].append([])
            ion_triple[5].append([])
        ##NaClNO3
            ion_triple[5].append(Cl[0]*NO3[0])
            ion_triple[3].append(ion_pair_same_kind[4][1]) #paraPhiPhi append ion_pair_same_kind(2).PhiPhi
            ion_triple[6].append(Na[0])
            Osmoticterms[3].append(ion_triple[6][3]*ion_triple[1][3])
        ##NaClSO4
            ion_triple[5].append(Cl[0]*SO4[0])
            ion_triple[3].append(ion_pair_same_kind[4][2])
            ion_triple[6].append(Na[0])
            Osmoticterms[3].append(ion_triple[6][4]*ion_triple[1][4])
        ##NaNO3SO4
            ion_triple[5].append(SO4[0]*NO3[0])
            ion_triple[3].append(ion_pair_same_kind[4][3])
            ion_triple[6].append(Na[0])
            Osmoticterms[3].append(ion_triple[6][5]*ion_triple[1][5])
        ##KClNO3
            ion_triple[5].append(Cl[0]*NO3[0])
            ion_triple[3].append(ion_pair_same_kind[4][1])
            ion_triple[6].append(K[0])
            Osmoticterms[3].append(ion_triple[6][6]*ion_triple[1][6])
        ##KClSO4
            ion_triple[5].append(Cl[0]*SO4[0])
            ion_triple[3].append(ion_pair_same_kind[4][2])
            ion_triple[6].append(K[0])
            Osmoticterms[3].append(ion_triple[6][7]*ion_triple[1][7])
        ##KNO3SO4
            ion_triple[5].append(SO4[0]*NO3[0])
            ion_triple[3].append(ion_pair_same_kind[4][3])
            ion_triple[6].append(K[0])
            Osmoticterms[3].append(ion_triple[6][8]*ion_triple[1][8])
        ##Term3
            sum_psi=0
            for i in range (0,3): #i=1:3
                j=Osmoticterms[2][i]
                sum_psi=j+sum_psi
            term3=ion_triple[2][0]*(ion_pair_same_kind[4][0]+sum_psi)
            Osmoticterms[4].append(term3)
        
        ##Term4
            term4=0
            for i in range (3,6): #i=4:6
                j=ion_triple[5][i]*(ion_triple[3][i]+ion_triple[6][i]*ion_triple[1][i]+ion_triple[6][i+3]*ion_triple[1][i+3])
                term4=j+term4
            
            Osmoticterms[5].append(term4)
        
        
            H=Osmoticterms[0][0]+Osmoticterms[1][0]+Osmoticterms[4][0]+Osmoticterms[5][0]
            osmot_coeff=2*H
            Osmoticterms[6].append(osmot_coeff)
            
            ion_mol=[[],[],[],[]]
          #clear
            ion_mol[0].append("Na") #name
            ion_mol[1].append(Na[0]) #mol
            ion_mol[2].append(Na[1]) #charge
            
            ion_mol[0].append(("K"))
            ion_mol[1].append(K[0])
            ion_mol[2].append(K[1])
            
            ion_mol[0].append(("Cl"))
            ion_mol[1].append(Cl[0])
            ion_mol[2].append(Cl[1])    
              
            ion_mol[0].append(("NO3"))
            ion_mol[1].append(NO3[0])
            ion_mol[2].append(NO3[1])
               
            ion_mol[0].append(("SO4"))
            ion_mol[1].append(SO4[0])
            ion_mol[2].append(SO4[1])
        #  
            sum_ions=0
            for i in range(0,5):#i=1:5
                j=ion_mol[1][i]
                sum_ions=j+sum_ions
        
            lna_H2O=-M_H2O*(Osmoticterms[6][0]+(sum_ions))
        
        ##Function F
            Function_terms=[]
        ##Term1 
            term1F = -APhi*(((math.sqrt(I_))/(1+CONST_B*math.sqrt(I_)))+(2/CONST_B)*math.log(1+CONST_B*math.sqrt(I_))) 
            Function_terms.append(term1F) 
            term2F = 0.0
        
            for i in range(0,6): #i=1:6
                term2= (ion_pair[12][i]*ion_pair[13][i])*ion_pair[11][i]
                term2F=term2F+term2
            Function_terms.append(term2F)
            
            term3F=0        
            for i in range(0,4): #i=1:4
                term3= ion_pair_same_kind[10][i]*ion_pair_same_kind[6][i]
                term3F=term3F+term3
            #end
            Function_terms.append(term3F)
        
            F=sum(Function_terms)
            Function_terms.append(F)
             
        ##activities
            for i in range (0,5): #i=1:5
                if ion_mol[1][i]==0:
                    ion_mol[3].append(0)
                else: 
                    ion_mol[3].append(math.log(ion_mol[1][i]))
            #end
            ln_gamma_M=[[],[]]
            for i in range(0,2): #i=1:2
                ln_gamma_M[i].append(ion_mol[0][i])
                ln_gamma_M[i].append((ion_mol[2][i])**2*F)    
            #end
        
        ##Term2 for Na
            term2ln=0
        
            for i in range(0,3): #i=1:3    
                term2= (ion_pair[13][i])*(2*ion_pair[10][i]+(Z_*ion_pair[1][i]))
                term2ln=term2ln+term2  
            #end
        
            ln_gamma_M[0].append(term2ln)
         
        # #Term3 for Na
            term3ln=0
        
            for i in range(0,6): #i=1:6
                term3= (ion_pair[13][i])*(ion_pair[12][i])*ion_pair[1][i]
                term3ln=term3ln+term3  
            #end
        
            ln_gamma_M[0].append(abs(ion_mol[2][0])*term3ln)
          #Term4 for Na
            term4ln=0
        
            for i in range(0,3): #i=1:3
                term4= (ion_triple[4][i]*ion_triple[1][i])
                term4ln=term4ln+term4  
            #end
            ln_gamma_M[0].append((ion_pair_same_kind[12][0])*((2*ion_pair_same_kind[5][0])+term4ln))
         
          #Term5 for Na
            term5ln=0
        
            for i in range(3,6): #i=4:6
                term5= (ion_triple[6][i]*ion_triple[1][i])     
                term5ln=term5ln+term5  
            #end
            ln_gamma_M[0].append(term5ln)
        
        ##Term2 for K
            term2ln=0
            for i in range(3,6): #i=4:6
                term2= (ion_pair[13][i]*(2*ion_pair[10][i]+(Z_*ion_pair[1][i])))
                term2ln=term2ln+term2   
            #end
            ln_gamma_M[1].append(term2ln)
            ln_gamma_M[1].append(abs(ion_mol[2][1])*term3ln)
        
        ##Term4 for K
        
            term4ln=0
            
            for i in range(0,3): #i=1:3
                term4= (ion_triple[4][i]*ion_triple[1][i])
                term4ln=term4ln+term4   
            #end
            ln_gamma_M[1].append((ion_pair_same_kind[11][0])*((2*ion_pair_same_kind[6][0])+term4ln))
         
            term5ln=0
            for i in range (6,9):
                term5=(ion_triple[5][i]*ion_triple[1][i])    
                term5ln=term5ln+term5   
            #end
            ln_gamma_M[1].append(term5ln)
        
        ##ln_gamma_X
            ln_gamma_X=[[],[],[]]
            for i in range (2,5):
                j=i-2
                ln_gamma_X[j].append(ion_mol[0][i])
                ln_gamma_X[j].append((ion_mol[2][i])**2*F)   
                #end
        ##Term2 for Cl
            term2ln=0
            term2= (ion_pair[12][0])*(2*ion_pair[10][0]+(Z_*ion_pair[1][0]))+(ion_pair[12][3])*(2*ion_pair[10][3]+(Z_*ion_pair[1][3]))
            term2ln=term2ln+term2
            ln_gamma_X[0].append(term2ln)
            ln_gamma_X[0].append(abs(ion_mol[2][2])*term3ln)
        ##Term4 for Cl
            term4ln=0
            #4,5,7,8:  
            term4= (ion_triple[6][3]*ion_triple[1][3])+ (ion_triple[6][4]*ion_triple[1][4])+ (ion_triple[6][6]*ion_triple[1][6])+ (ion_triple[6][7]*ion_triple[1][7])
            term4ln=term4ln+term4
            term4psk=0
        ##i=3,2
            term4=(ion_pair_same_kind[14][1])*((2*ion_pair_same_kind[5][1])+term4ln)+(ion_pair_same_kind[14][2])*((2*ion_pair_same_kind[5][2])+term4ln)
            term4psk=term4psk+term4
            ln_gamma_X[0].append(term4psk)
            ln_gamma_X[0].append((ion_triple[2][0]*ion_triple[1][0]))
        ##Term2 for NO3
            term2ln=0
        ## i=2,5
            term2= (ion_pair[12][1])*(2*ion_pair[10][1]+(Z_*ion_pair[1][1]))+(ion_pair[12][4])*(2*ion_pair[10][4]+(Z_*ion_pair[1][4]))   
            term2ln=term2ln+term2
            ln_gamma_X[1].append(term2ln)
            ln_gamma_X[1].append(abs(ion_mol[2][3])*term3ln)
        ##Term4 for NO3
        
            term4ln=0
        ## i=4,6,7,9
            term4=(ion_triple[6][3]*ion_triple[1][3])+(ion_triple[6][5]*ion_triple[1][5])+(ion_triple[6][6]*ion_triple[1][6])+(ion_triple[6][8]*ion_triple[1][8])
            term4ln=term4ln+term4
            term4psk=0
            #term4=(ion_pair_same_kind[13][1])*((2*ion_pair_same_kind[5][1])+term4ln)+(ion_pair_same_kind[13][3])*((2*ion_pair_same_kind[5][3])+term4ln)
            term4psk=term4psk+term4
        
            ln_gamma_X[1].append(term4psk)
            ln_gamma_X[1].append((ion_triple[2][1]*ion_triple[1][1]))
        
        ##Term2 for SO4
            term2ln=0
        ## i=3,5
            term2= (ion_pair[12][2])*(2*ion_pair[10][2]+(Z_*ion_pair[1][2])) +(ion_pair[12][4])*(2*ion_pair[10][4]+(Z_*ion_pair[1][4]))    
            term2ln=term2ln+term2
            
            ln_gamma_X[2].append(term2ln)
            ln_gamma_X[2].append(abs(ion_mol[2][4])*term3ln)
        ##Term4 for SO4
        ## i=5,6,8,9
            term4= (ion_triple[6][4]*ion_triple[1][4])+(ion_triple[6][5]*ion_triple[1][5])+(ion_triple[6][7]*ion_triple[1][7]) +(ion_triple[6][8]*ion_triple[1][8])     
            term4ln=term4ln+term4
        
            term4psk=0
        ##i=2,3
            #term4=(ion_pair_same_kind[13][3])*((2*ion_pair_same_kind[5][3])+term4ln) +(ion_pair_same_kind[14][2])*((2*ion_pair_same_kind[5][2])+term4ln)
            term4psk=term4psk+term4
            ln_gamma_X[2].append(term4psk)
            ln_gamma_X[2].append((ion_triple[2][2]*ion_triple[1][2]))
        
        ##ln_gamma_M
            ln_a=[[],[],[],[],[]]
            for i in range(0,2):
                ln_gamma_M[i].append(ln_gamma_M[i][1]+ln_gamma_M[i][2]+ln_gamma_M[i][3]+ln_gamma_M[i][4]+ln_gamma_M[i][5])#(i).term1+ln_gamma_M(i).term2+ln_gamma_M(i).term3+ln_gamma_M(i).term4+ln_gamma_M(i).term5)
            #end
            for i in range(0,3):
                ln_gamma_X[i].append(ln_gamma_X[i][1]+ln_gamma_X[i][2]+ln_gamma_X[i][3]+ln_gamma_X[i][4]+ln_gamma_X[i][5])
            #end
            for j in range(0,5):
                if j<2:
                    ln_a[j].append(ion_mol[0][j])
                    ln_a[j].append(ln_gamma_M[j][6]+math.log1p(ion_mol[1][j]))
                else:
                    ln_a[j].append(ion_mol[0][j])
                    if ion_mol[1][j]!=0:
                        ln_a[j].append(ln_gamma_X[j-2][6]+math.log(ion_mol[1][j]))
                    else:
                        # ln_a[j].append("inf")
                        ln_a[j].append(0)
                #end
            #end
        ##lnK
        ##Ice
            Salt_list[2].append(lna_H2O)
            Salt_list[8].append(1) #num h20
            for i in range(3,8):
                Salt_list[i].append([])    
        # ##NaCl
            Salt_list[8].append(0) #num h20
            Salt_list[3].append(1) #num_Na
            Salt_list[4].append(0) #num_K
            Salt_list[5].append(1) #num_Cl
            Salt_list[6].append(0) #num_NO3
            Salt_list[7].append(0) #num_SO4    
            Salt_list[8].append(0) #num_H2O
            Salt_list[2].append(Salt_list[3][1]*ln_a[0][1]+ Salt_list[5][1]*ln_a[2][1] +Salt_list[8][1]*lna_H2O) #lnk
        ##NaCl.2H2O
            Salt_list[8].append(2) #num h20
            Salt_list[3].append(1) #num_Na
            Salt_list[4].append(0) #num_K
            Salt_list[5].append(1) #num_Cl
            Salt_list[6].append(0) #num_NO3
            Salt_list[7].append(0) #num_SO4    
            Salt_list[8].append(2) #num_H2O
            Salt_list[2].append(Salt_list[3][2]*ln_a[0][1]+ Salt_list[5][2]*ln_a[2][1] +Salt_list[8][1]*lna_H2O)
        ##NaNO3
            Salt_list[8].append(0) #num h20
            Salt_list[3].append(1) #num_Na
            Salt_list[4].append(0) #num_K
            Salt_list[5].append(0) #num_Cl
            Salt_list[6].append(1) #num_NO3
            Salt_list[7].append(0) #num_SO4    
            Salt_list[8].append(0) #num_H2O
            Salt_list[2].append(0)#Salt_list[3][3]*ln_a[0][1]+Salt_list[6][3]*0+Salt_list[8][3]*lna_H2O) #lna_activitiy for no3 =-inf ->0
        ##Na2SO4
            Salt_list[8].append(0) #num h20
            Salt_list[3].append(2) #num_Na
            Salt_list[4].append(0) #num_K
            Salt_list[5].append(0) #num_Cl
            Salt_list[6].append(0) #num_NO3
            Salt_list[7].append(1) #num_SO4    
            Salt_list[8].append(0) #num_H2O
            Salt_list[2].append(Salt_list[3][4]*ln_a[0][1]+ Salt_list[7][4]*ln_a[4][1]+Salt_list[8][4]*lna_H2O)
        ##Na2SO4.10H2O
            Salt_list[8].append(10) #num h20
            Salt_list[3].append(2) #num_Na
            Salt_list[4].append(0) #num_K
            Salt_list[5].append(0) #num_Cl
            Salt_list[6].append(0) #num_NO3
            Salt_list[7].append(1) #num_SO4    
            Salt_list[8].append(10) #num_H2O
            Salt_list[2].append(Salt_list[3][5]*ln_a[0][1]+ Salt_list[7][5]*ln_a[4][1]+Salt_list[8][5]*lna_H2O)
        ##KCl
            Salt_list[8].append(0) #num h20
            Salt_list[3].append(0) #num_Na
            Salt_list[4].append(1) #num_K
            Salt_list[5].append(1) #num_Cl
            Salt_list[6].append(0) #num_NO3
            Salt_list[7].append(0) #num_SO4    
            Salt_list[8].append(0) #num_H2O
            Salt_list[2].append(Salt_list[4][6]*ln_a[1][1]+ Salt_list[5][6]*ln_a[2][1]+Salt_list[8][6]*lna_H2O)
        ##KNO3
            Salt_list[8].append(0) #num h20
            Salt_list[3].append(0) #num_Na
            Salt_list[4].append(1) #num_K
            Salt_list[5].append(0) #num_Cl
            Salt_list[6].append(1) #num_NO3
            Salt_list[7].append(0) #num_SO4    
            Salt_list[8].append(0) #num_H2O
            Salt_list[2].append(Salt_list[4][7]*ln_a[1][1]+Salt_list[6][7]*0+Salt_list[8][7]*lna_H2O)
        ##K2SO4
            Salt_list[8].append(0) #num h20
            Salt_list[3].append(0) #num_Na
            Salt_list[4].append(2) #num_K
            Salt_list[5].append(0) #num_Cl
            Salt_list[6].append(0) #num_NO3
            Salt_list[7].append(1) #num_SO4    
            Salt_list[8].append(0) #num_H2O
            Salt_list[2].append(Salt_list[4][8]*ln_a[1][1]+ Salt_list[7][8]*ln_a[4][1]+Salt_list[8][8]*lna_H2O)
        ##Na2K6(SO4)4
            Salt_list[8].append(0) #num h20
            Salt_list[3].append(2) #num_Na
            Salt_list[4].append(6) #num_K
            Salt_list[5].append(0) #num_Cl
            Salt_list[6].append(0) #num_NO3
            Salt_list[7].append(4) #num_SO4    
            Salt_list[8].append(0) #num_H2O
            Salt_list[2].append(Salt_list[3][9]*ln_a[0][1]+ Salt_list[4][9]*ln_a[1][1]+Salt_list[7][9]*ln_a[4][1]+Salt_list[8][9]*lna_H2O)
        ##Na3SO4NO3.H2O
            Salt_list[8].append(1) #num h20
            Salt_list[3].append(3) #num_Na
            Salt_list[4].append(0) #num_K
            Salt_list[5].append(0) #num_Cl
            Salt_list[6].append(1) #num_NO3
            Salt_list[7].append(1) #num_SO4    
            Salt_list[8].append(1) #num_H2O
            Salt_list[2].append(Salt_list[3][10]*ln_a[0][1]+ Salt_list[6][10]*0+Salt_list[7][10]*ln_a[4][1]+Salt_list[8][10]*lna_H2O)
            
            ln_K0=[]
            ln_K=[]
            n_step=int((T_C-T_ini)/T_interval)
            # print("n_step is "+str(n_step))
            for j in range(0,10):
                ln_K0.append(Salt_list[1][j])
                ln_K.append(Salt_list[2][j])
 
                U_salt[0].append(Salt_list[0][j])
                U_salt[1].append(T_C) 
                U_salt[2].append(ln_K0[j])
                U_salt[3].append(ln_K[j])
                try:
                    #U_salt[4].append(math.exp(round(ln_K[j],2)-round(ln_K0[j],2)))
                    U_salt[4].append(math.exp((ln_K[j])-(ln_K0[j])))
                except OverflowError:
                    # U_salt[4].append(float('inf'))
                    U_salt[4].append(float(0))
            for j in range(0,len(U_salt[0])):
                # print(j)
                if U_salt[4][j]>=1.0 and U_salt[4][j]<1.5:
                    if (U_salt[0][j]=='Ice') and not any(T_C == temp for temp in line_ice[0])  and not any((SO4[0]) == c for c in line_ice[2]):
                        line_ice[0].append(T_C)
                        line_ice[1].append(U_salt[4][j])
                        line_ice[2].append(SO4[0])
                        #print("sat point ice"+str(U_salt[4][j]))
                        line_ice[3].append(j)
                    if (U_salt[0][j]=='NaCl.2H2O') and not any(T == temp for temp in line_nacl[0]) and not any((c_cl) == c for c in line_nacl[2]):
                        line_nacl[0].append(T)
                        line_nacl[1].append(U_salt[4][j])
                        line_nacl[2].append(Cl[0])
                    if (U_salt[0][j]=='Na2SO4.10H2O') and not any(T == temp for temp in line_na2so4[0]) and not any(SO4[0] == c for c in line_na2so4[2]):
                        line_na2so4[0].append(T)
                        line_na2so4[1].append(U_salt[4][j])
                        line_na2so4[2].append(SO4[0])       #(c_na_1+c_na_2)
                        line_na2so4[3].append(U_salt[4][0])
                        
                     
            # U_salt_all.append(U_salt)            
        #end
    
    #end
    
print("na conc"+str(Na[0]))
print("so4 conc "+str(SO4[0]))
T_Na2SO4_sol_line=line_na2so4[0]
C_Na2SO4_sol_line=line_na2so4[2]
C_ice_sol_line=line_ice[2]
T_ice_sol_line=line_ice[0]

df=pd.DataFrame(np.array(T_Na2SO4_sol_line), np.array(C_Na2SO4_sol_line))
df2=pd.DataFrame(np.array(T_ice_sol_line), np.array(C_ice_sol_line))
with pd.ExcelWriter('C:\\Users\\rodoulaktori\\surfdrive\\PhD\\Process_modelling\\python\\efc_results.xlsx') as writer:
    df.to_excel(writer,sheet_name="C_Na2SO4_sol_line")
    df2.to_excel(writer,sheet_name="ICE_sol_line")


# #figure (1): x: temperature (U_Salt), y:saturation (U_Salt)
# h=0
# l=[]
# for j in [0,1,5]:#[1:2,6]
#     h=h+1
#     for f in range (0,int(g)): #=1:i
#         plt.plot([U_salt[1][j]],[U_salt[4][j]]) #, label=U_salt[0][j]
#         #plt.plot(U_salt[1][j+((f-1)*11)],U_salt[4][j+((f-1)*11)], label=U_salt[0][j])
#         #plt.legend()
#     l.append(U_salt[0][j])
# k=0
# x1=[]
# y1=[]
# t=[]
# for k in range(0,210,10):
#     x1.append(U_salt[1][k])
#     y1.append(U_salt[4][k])
# plt.plot(x1,y1,label='ice')
# x2=[]
# y2=[]
# k=2
# for k in range(2,210,11):
#     x2.append(U_salt[1][k])
#     y2.append(U_salt[4][k])
#     t.append(1)
# plt.plot(x2,y2,label='nacl2h2o')
# x3=[]
# y3=[]
# k=5
# for k in range(5,210,11):
#     x3.append(U_salt[1][k])
#     y3.append(U_salt[4][k])
# plt.plot(x3,y3,label='na2so4.10h2o')
# plt.plot(x3,t, label='saturation line')
# plt.title('saturation vs temperature')
# plt.legend()
# plt.show()
# plt.savefig('saturation vs temperature.png')

# df=pd.DataFrame(data=ion_pair_same_kind)

#     end
# end


