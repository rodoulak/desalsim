#Crystal Size Disribution in a Eutectic System

'''Model for the estimation of the distribution of crystal sizes 
during the cooling and eutectic freeze crystallization of an
aqueous solution of Compound 1 (Na2SO4/Sodium sulfate), that initially includes traces of Compound 2 (considered as impurity: Extracellular Polymeric Substances/Sodium chloride), in a batch/continuous reactor'''

#%%
#constants 
dp=0.1
dp_slurry=1
npump=0.8

#Import the libraries that will be used
import math
import pandas as pd
pi=math.pi
exp=math.exp
log=math.log

import matplotlib.pyplot as plt
fig1=plt.figure(1)
fig2=plt.figure(2)
fig3=plt.figure(3)
fig4=plt.figure(4)
fig5=plt.figure(5)
fig6=plt.figure(6)

import time
import numpy as np
from scipy import interpolate
import scipy.interpolate as interpolate
from sklearn.linear_model import LinearRegression
from collections import namedtuple
from desalsim import scaleup
from desalsim.density_calc import density_calc 
# import pitzer_model_3phases_w
from desalsim import constants 
#Set initial time to calculate the elapsed time of the calculations
time0=time.time()
MW_so4=32.064+4*15.999

#input data
Qf=100 #m3/hr
C_i_in=[20.16, 27.87, 0.50, 0.0, 0.12, 14.66]

T_initial =273+25
d_in=density_calc(25, sum(C_i_in))
#%%

#Declaration of constant process parameters - Part 1: Process Parameters

Vreactor=Qf/0.7            #Volume of the reactor : m^3

#Influx
Fin_water=0                                                                                                         #Influx of water (in mother liquor) in the reactor : kg/s
Fin_compound1=0                                                                                                     #Influx of Compound 1 (in mother liquor) in the reactor : kg/s
Fin_compound2=0                                                                                                     #Influx of Compound 2 (in mother liquor) in the reactor : kg/s
Fin_mother_liquor=Fin_water+Fin_compound1+Fin_compound2                                                             #Influx of mother liquor in the reactor : kg/s
Fin_liquid_tot=Fin_mother_liquor                                                                                    #Total influx of liquid in the reactor : kg/s

Fin_ice=0                                                                                                           #Influx of ice crystals in the reactor : kg/s
Fin_compound1_cr=0                                                                                                  #Influx of Compound 1 crystals in the reactor : kg/s
Fin_solids_tot=Fin_ice+Fin_compound1_cr                                                                             #Total influx of solid crystals in the reactor : kg/s

#Outflux
Fout_water=0                                                                                                        #Outflux of water (in mother liquor) from the reactor : kg/s
Fout_compound1=0                                                                                                    #Outflux of Compound 1 (in mother liquor) from the reactor : kg/s
Fout_compound2=0                                                                                                    #Outflux of Compound 2 (in mother liquor) from the reactor : kg/s
Fout_mother_liquor=Fin_water+Fin_compound1+Fin_compound2                                                            #Outflux of mother liquor from the reactor : kg/s
Fout_liquid_tot=Fout_mother_liquor                                                                                  #Total outflux of liquid from the reactor : kg/s

Fout_ice=0                                                                                                          #Outflux of ice crystals from the reactor : kg/s
Fout_compound1_cr=0                                                                                                 #Outflux of Compound 1 crystals from the reactor : kg/s
Fout_solids_tot=Fout_ice+Fout_compound1_cr                                                                          #Total outflux of solid crystals from the reactor : kg/s

Fout_top=0                                                                                                          #Outflux from the top (ice and mother liquor) : kg/s
Fout_bottom=0                                                                                                       #Outflux from the bottom (Compound 1 and mother liquor) : kg/s

#Net Heat flux 
Cooling_rate=2500000                                                                                                 #Assumed NET cooling rate to cool down : W = J/s :Assumption !!! Includes the heat removed by the heat exchanger and the heat provided by the environment and the system
    
#General Parameters
classes=300                                                                                                         #Number of crystal size classes
maxl=1e-3                                                                                                           #Maximum crystal length : m
deltal=maxl/classes                                                                                                 #Length interval per class : m
deltat=5                                                                                                           #Time interval : sec  

class_length=np.array([])                                                                                           #Array of the lengths for the crystal size classes used for plotting
for i in range(0,classes):
       class_length=np.append(class_length,((i*deltal)+deltal))
       
#Storage lists for different values
t_store=np.array([])
sigma_compound1_store=np.array([])
DT_ice_store=np.array([])
M_ice_store=np.array([])
M_compound1_cr_store=np.array([])
Treactor_store=np.array([])


#%%

#Declaration of constant process parameters - Part 2: WATER Properties library

#Physical properties
d_water=1*1000                                                                                                      #Density of water (25oC) : kg/m^3 (source:Wikipedia) 
Cp_water=4.1813*1000                                                                                                #Heat capacity of water (25 oC) : J/(kg * K) (source:Wikipedia)
#%%

#Declaration of constant process parameters - Part 3: ICE Properties 

#Kinetic and Shape Parameters

#Ice
kg_ice=5e-8                                                                                                         #Growth rate constant for ice crystals : m/(K*s) (source: Characterization and Population....)
kb_ice=5e4                                                                                                          #Nucleation rate constant for ice crystals : 1/(m^2*s*K^2) (source: Characterization and Population....)
kd_ice=1e-10                                                                                                        #ASSUMED Dissolution rate constant for ice crystals : m/(K*s)

r_ice=10                                                                                                            #L/h aspect ratio of ice crystals, assuming a disk shape where length=10*height : -
kv_ice=pi/(4*r_ice)                                                                                                 #Shape factor of ice crystals, where V=kv*l^3 (assuming disk shape and l=longest side) : -
ka_ice=pi*(1/r_ice+1)                                                                                               #Surface factor of ice crystals, where A=ka*l^2 (assuming disk shape and l=longest side) : -

#Physical Properties
d_ice=0.917*1000                                                                                                    #Density of ice (25oC) : kg/m^3 (source:https://hypertextbook.com/facts/2000/AlexDallas.shtml)
Cp_ice=2.05*1000                                                                                                    #Heat capacity of ice (0 oC) : J/(kg * K) (source:https://www.engineeringtoolbox.com/ice-thermal-properties-d_576.html)

#Thermodynamic Properties
DHfus_ice=333.55*1000                                                                                               #Enthalpy of fusion of ice : J/kg (source:Wikipedia)
#%%
#Declaration of constant process parameters - Part 5: Sodium Sulfate Properties 

#Physical Properties
MW_Na2SO4=142.04                
MW_Na2SO4_10H2O=322.192258                                                                         #Molecular weight of sodium sulfate : g/mol
d_Na2SO4=2.66*1000                                                                                                  #Density of sodium sulfate (25oC) : kg/m^3 (source:???) 
Cp_Na2SO4=(128.2/MW_Na2SO4)*1000                                                                                      #Heat capacity of sodium sulfate (25 oC) : J/(kg * K) (source:source: https://en.wikipedia.org/wiki/Sodium_sulfate_(data_page))

#Kinetic and Shape Parameters
kg_Na2SO4=(43960*exp(-24600/((8.314)*T_initial )))/(1000*3600)                                                             #Growth rate constant for sodium sulfate crystals : m/s : Mirabilite at 25oC(source: crystallization kinetics of sodium sulfate... #UBC TBC)
kb_Na2SO4=(0.064*exp(70900/((8.314)*T_initial )))*1000/(3600*d_Na2SO4)                                                     #Nucleation rate constant for sodium sulfate crystals : 1/(kg*s) : Mirabilite (source: crystallization kinetics of sodium sulfate... #UBC TBC)
kd_Na2SO4=1e-10                                                                                                     #ASSUMED Dissolution rate constant for sodium sulfate crystals : m/s

r_Na2SO4=1/10                                                                                                       #L/h aspect ratio of sodium sulfate crystals, assuming a cylindrical shape where length=height/10 : -
kv_Na2SO4=(pi*r_Na2SO4**2)/4                                                                                        #Shape factor of sodium sulfate crystals, where V=kv*h^3 (assuming cylindrical shape and h=longest side) : -
ka_Na2SO4=(pi*r_Na2SO4)*(1+r_Na2SO4/2)                                                                              #Surface factor of sodium sulfate crystals, where A=ka*h^2 (assuming cylindrical shape and h=longest side) : -

#Thermodynamic Properties
DHfus_Na2SO4=252*1000                                                                                               #Enthalpy of fusion of sodium sulfate : J/kg (source:Solubility of Na2SO4...)

#Create the library with the constant parameters for sodium sulfate
Na2SO4Parameters = namedtuple('Na2SO4Parameters', ['MW','d','Cp',
                                                   'DHfus','kg','kb','kd','r',
                                                   'kv','ka'])

Na2SO4Par = Na2SO4Parameters(MW=MW_Na2SO4,        
                             d=d_Na2SO4,        
                             Cp=Cp_Na2SO4,    
                             DHfus=DHfus_Na2SO4,  
                             kg=kg_Na2SO4,      
                             kb=kb_Na2SO4,        
                             kd=kd_Na2SO4,
                             r=r_Na2SO4,
                             kv=kv_Na2SO4,         
                             ka=ka_Na2SO4)  
#%%
#Declaration of constant process parameters - Part 7: Sodium Chloride Properties

#Physical Properties
MW_NaCl=58.4                                                                                                        #Molecular weight of NaCl : g/mol (ASSUMED??????) 
d_NaCl=2.16*1000                                                                                                    #Density of NaCl (25oC) : kg/m^3 (source:https://pubchem.ncbi.nlm.nih.gov/compound/sodium_chloride#section=Solubility) 
Cp_NaCl=(50.5/MW_NaCl)*1000                                                                                         #Heat capacity of NaCl (25 oC) : J/(kg * K) (source:http://webbook.nist.gov/cgi/cbook.cgi?ID=C7647145&Type=JANAFS&Table=on)

#Create the library with the constant parameters for sodium chloride
NaClParameters = namedtuple('NaClParameters', ['MW','d','Cp'])

NaClPar = NaClParameters(MW=MW_NaCl,        
                       d=d_NaCl,        
                       Cp=Cp_NaCl)
#%%

#Assign all the used (non-constant) variables to a class
 
class Variables:
    def __init__(self,t,N_ice,N_compound1_cr,V_ice,V_compound1_cr,V_solids_tot,M_ice,M_compound1_cr,M_solids_tot,Ntot_ice,Ntot_compound1_cr,Ntot,C_compound1,C_compound2,Treactor,V_liquid_tot,V_mother_liquor,V_water,V_compound1,V_compound2,M_water,M_compound1,M_compound2,M_mother_liquor,M_liquid_tot,Vtot,Mtot,solids_weight_percentage,V_ice_distr,V_compound1_cr_distr,V_ice_previous_loop,V_compound1_cr_previous_loop):
        self.t=t
        self.N_ice=N_ice
        self.N_compound1_cr=N_compound1_cr
        self.V_ice=V_ice
        self.V_compound1_cr=V_compound1_cr
        self.V_solids_tot=V_solids_tot
        self.M_ice=M_ice
        self.M_compound1_cr=M_compound1_cr
        self.M_solids_tot=M_solids_tot
        self.Ntot_ice=Ntot_ice
        self.Ntot_compound1_cr=Ntot_compound1_cr
        self.Ntot=Ntot
        self.C_compound1=C_compound1
        self.C_compound2=C_compound2
        self.Treactor=Treactor
        self.V_liquid_tot=V_liquid_tot
        self.V_mother_liquor=V_mother_liquor
        self.V_water=V_water
        self.V_compound1=V_compound1
        self.V_compound2=V_compound2
        self.M_water=M_water
        self.M_liquid_tot=M_liquid_tot
        self.M_mother_liquor=M_mother_liquor
        self.M_compound1=M_compound1
        self.M_compound2=M_compound2
        self.Mtot=Mtot
        self.Vtot=Vtot
        self.solids_weight_percentage=solids_weight_percentage
        self.V_ice_distr=V_ice_distr
        self.V_compound1_cr_distr=V_compound1_cr_distr
        self.V_ice_previous_loop=V_ice_previous_loop
        self.V_compound1_cr_previous_loop=V_compound1_cr_previous_loop
       
Var=Variables

#%%

#Function for the selection of the compounds that will be used for the simulation. 

def Compounds():
    global Compound1,Compound2,tck,lm,Linear_regr,C_compound1_initial,C_compound2_initial
    
    #Initial concentrations of the different compounds
    #Initial sodium chloride concentration in the mother liquor in the reactor : molal
    C_so4_in=(C_i_in[5]/constants.MW_so4)/(d_in/1000) #mol/kg
    print("C_so4_i is "+str(C_so4_in))
    C_nacl_in=(C_i_in[0]/constants.MW_Na-2*C_so4_in)/(d_in/1000) #mol/kg
    print("C_nacl_in" + str(C_nacl_in))
    C_Na2SO4_initial=np.array([C_so4_in]).reshape(-1,1)   #np.array([0.364]).reshape(-1,1)                                                                                             #Initial sodium sulfate concentration in the mother liquor in the reactor : molal
    C_NaCl_initial=np.array([C_nacl_in]).reshape(-1,1)                               # 0.001                                                        
    #The user is given a choice for Compound 1 (which is the main compound crystallizing), in this case sodium sulfate  
    
    Compound1=Na2SO4Par
       
    print('Compound 1 : Sodium Sulfate')
        
    #Initial concentration of Compound 1
    C_compound1_initial=C_Na2SO4_initial
    print("C_compound1_initial is "+str(C_compound1_initial))
    T_Na2SO4_sol_line=[]
    C_Na2SO4_sol_line=[]
    #Interpolation/Extrapolation of the Na2SO4 solubility line with a cubic spline 
    import os
    import pandas as pd

    # Get the directory where the script is located
    current_directory = os.path.dirname(__file__)

    #   Define the filename 
    filename = "efc_results.xlsx"

    # Combine the directory path with the filename
    file_path = os.path.join(current_directory, filename)

    # Read the Excel file
    df = pd.read_excel(file_path, sheet_name="C_Na2SO4_sol_line")

    #Sodium sulfate (Mirabilite) Solubility Line, Data from source: The solution ref on dropbox (1100-1200; P.1120)
    T_Na2SO4_sol_line= df.iloc[:, 1].tolist()
    C_Na2SO4_sol_line=df.iloc[:, 0].tolist()
    C_Na2SO4_sol_line=np.array(C_Na2SO4_sol_line)

    print("T_Na2SO4_sol_line "+str(T_Na2SO4_sol_line))

    print("C_Na2SO4_sol_line "+str(C_Na2SO4_sol_line))
                                                            

    tck= interpolate.splrep(T_Na2SO4_sol_line, C_Na2SO4_sol_line, k=1)
    print("tck is " +str(tck))

    #Interpolation/Extrapolation of the ice line using linear regression
    df2 = pd.read_excel(file_path, sheet_name="ICE_sol_line")

    Fr_point_depr_Na2SO4= np.array(df2.iloc[:, 1].tolist())
    C_Na2SO4_ice_line=df2.iloc[:, 0].tolist()    
    print("C_Na2SO4_ice_line is "+str(C_Na2SO4_ice_line))
    #Ice line in the Na2SO4-Water phase diagram : Source : CONCENTRATIVE PROPERTIES OF AQUEOUS SOLUTIONS:DENSITY, REFRACTIVE INDEX, FREEZING POINT DEPRESSION, AND VISCOSITY

    T_fr_Na2SO4=273+Fr_point_depr_Na2SO4    

    print("T_fr_na2so4 is " + str(T_fr_Na2SO4))
    lm = LinearRegression()
    Linear_regr=lm.fit(np.array(C_Na2SO4_ice_line).reshape(-1, 1), T_fr_Na2SO4)
    #The user is also given a choice for Compound 2 (initially considered as impurity), in this case sodium chloride      
    Compound2=NaClPar

    print('Compound 2 : Sodium Chloride')
        
    #Initial concentration of Compound 2
    C_compound2_initial=C_NaCl_initial
    return(C_Na2SO4_sol_line, tck)
   
#%%

#Function for the initialization of all the (non-constant) variables, that are assigned to the Variables class 
def Initialize(Var):
    global Compound1,Compound2,C_compound1_initial,C_compound2_initial
    #Initial values
    
    Var.t=0                                                                                                             #Time : sec
    
    Var.N_ice=np.array([])                                                                                              #Number of ice crystals in each class in the reactor : absolute number of ice crystals
    Var.N_compound1_cr=np.array([])                                                                                     #Number of compound 1 crystals in each class in the reactor : absolute number of Na2SO4 crystals
    
    Var.V_ice=0                                                                                                         #Total volume of ice crystals in the reactor: m^3
    Var.V_compound1_cr=0                                                                                                #Total volume of compound 1 crystals in the reactor: m^3
    Var.V_solids_tot=Var.V_ice+Var.V_compound1_cr                                                                       #Total volume of solids in the reactor : m^3
    
    Var.V_ice_distr=np.zeros(classes)                                                                                   #Volume distribution of ice crystals in each class : m^3
    Var.V_compound1_cr_distr=np.zeros(classes)                                                                          #Volume distribution of compound 1 crystals in each class : m^3
    
    Var.M_ice=0                                                                                                         #Total mass of ice in the reactor : kg
    Var.M_compound1_cr=0                                                                                                #Total mass of compound 1 crystals in the reactor : kg
    Var.M_solids_tot=Var.M_ice+Var.M_compound1_cr                                                                       #Total mass of solids in the reactor : kg
    
    Var.Ntot_ice=0                                                                                                      #Total number of ice crystals in the reactor
    Var.Ntot_compound1_cr=0                                                                                             #Total number of compound 1 crystals in the reactor
    Var.Ntot=Var.Ntot_ice+Var.Ntot_compound1_cr                                                                         #Total number of crystals in the reactor
    
    Var.C_compound1=C_compound1_initial                                                                                 #Initial compound 1 concentration in the mother liquor in the reactor : molal = mol compound 1/kg water
    Var.C_compound2=C_compound2_initial                                                                                 #Initial compound 2 concentration in the mother liquor in the reactor : molal = mol compound 2/kg water
    
    Var.Treactor=T_initial                                                                                                  #Initial temperature of the bulk liquid in the crystallizer : K
    
    Var.V_liquid_tot=0.7*Vreactor                                                                                       #Total volume of liquid in the reactor : m^3 (ASSSUMED to be 70 % of the reactor volume)
    Var.V_mother_liquor=Var.V_liquid_tot                                                                                #Volume of the mother liquor in the reactor = Vwater+V_Na2SO4+Vimp: m^3
    print("v mother liquid is "+str(Var.V_mother_liquor))
    
    Var.V_water=Var.V_mother_liquor/(1+Var.C_compound1*Compound1.MW*d_water/(1000*Compound1.d)+Var.C_compound2*Compound2.MW*d_water/(1000*Compound2.d))       #Volume of water in the reactor : m^3
    Var.V_compound1=(Var.C_compound1*Compound1.MW*Var.V_water*d_water)/(1000*Compound1.d)                                                                     #Volume of compound 1 in the reactor : m^3
    Var.V_compound2=Var.V_mother_liquor-Var.V_water-Var.V_compound1  
    print("Var.V_water is "+str(Var.V_water))

    Var.M_water=Var.V_water*d_water                                                                                     #Mass of water in the reactor : kg
    Var.M_compound1=Var.V_compound1*Compound1.d                                                                         #Mass of compound 1 in the reactor : kg
    print("Var.M_compound1 is "+str(Var.M_compound1))
    Var.M_compound2=Var.V_compound2*Compound2.d                                                                         #Mass of compound 2 in the reactor : kg
    Var.M_mother_liquor=Var.M_water+Var.M_compound1+Var.M_compound2                                                     #Mass of the mother liquor in the reactor :kg
    Var.M_liquid_tot=Var.M_mother_liquor                                                                                #Total mass of liquid in the reactor :kg
    
    Var.Vtot=Var.V_solids_tot+Var.V_liquid_tot                                                                          #Total volume of solids and liquid in the reactor : m^3
    Var.Mtot=Var.M_solids_tot+Var.M_liquid_tot                                                                          #Total mass of solids and liquid in the reactor : kg
    print("Mtot in loop is "+str(Var.Mtot))
    Var.V_ice_previous_loop=Var.V_ice                                                                                   #Total volume of ice crystals stored from the previous time interval : m^3
    Var.V_compound1_cr_previous_loop=Var.V_compound1_cr                                                                 #Total volume of compound 1 crystals stored from the previous time interval : m^3
    
    Var.solids_weight_percentage=(Var.M_solids_tot/Var.Mtot)*100                                                        #Solids weight percentage in the reactor : solids wt % = (total mass of solids in the reactor / total mass in the reactor)*100 % (neglecting the potential solid impurities)
    
    for i in range(0,classes):                                                                                          #Loop over the number of (classes-1) to initialize all the values (The range of the loop is from 0 to classes due to python syntax, the loop actually runs from 0 to classes-1=299 )
        Init_size_distr=2e16*deltal*(exp((-i*deltal)/2.5e-5))                                                           #Initialize the number of crystals in each with this standard distribution (source:???)
        
        Var.N_ice=np.append(Var.N_ice,Init_size_distr)
        Var.N_compound1_cr=np.append(Var.N_compound1_cr,Init_size_distr)
        
        Var.V_ice=Var.V_ice+Var.N_ice[i]*kv_ice*exp(3*log(deltal*(i+1)))
        Var.V_compound1_cr=Var.V_compound1_cr+Var.N_compound1_cr[i]*Compound1.kv*exp(3*log(deltal*(i+1)))
    
    Var.V_ice_distr[0]=Var.N_ice[0]*kv_ice*0.25*(deltal**4)/deltal
    Var.V_compound1_cr_distr[0]=Var.N_compound1_cr[0]*Compound1.kv*0.25*(deltal**4)/deltal
    
    for i in range (1,classes-1):
        
        Var.V_ice_distr[i]=Var.N_ice[i]*kv_ice*0.25*((((i+1)*deltal)**4)-((i*deltal)**4))/deltal
        Var.V_compound1_cr_distr[i]=Var.N_compound1_cr[i]*Compound1.kv*0.25*((((i+1)*deltal)**4)-((i*deltal)**4))/deltal
            
    Var.V_ice_distr[classes-1]=Var.N_ice[classes-1]*kv_ice*0.25*(((classes*deltal)**4)-(((classes-1)*deltal)**4))/deltal
    Var.V_compound1_cr_distr[classes-1]=Var.N_compound1_cr[classes-1]*Compound1.kv*0.25*(((classes*deltal)**4)-(((classes-1)*deltal)**4))/deltal
    
    Var.V_solids_tot=Var.V_ice+Var.V_compound1_cr
        
    Var.M_ice=Var.M_ice+Var.V_ice*d_ice
    Var.M_compound1_cr=Var.M_compound1_cr+Var.V_compound1_cr*Compound1.d

    Var.M_solids_tot=Var.M_ice+Var.M_compound1_cr 
       
    Var.Ntot_ice=Var.Ntot_ice+sum(Var.N_ice)
    Var.Ntot_compound1_cr=Var.Ntot_compound1_cr+sum(Var.N_compound1_cr)
    Var.Ntot=Var.Ntot_ice+Var.Ntot_compound1_cr
    
    Var.Vtot=Var.V_solids_tot+Var.V_liquid_tot
    Var.Mtot=Var.M_solids_tot+Var.M_liquid_tot
    
    Var.V_ice_previous_loop=Var.V_ice                                                               
    Var.V_compound1_cr_previous_loop=Var.V_compound1_cr                                                       
    
    Var.solids_weight_percentage=(Var.M_solids_tot/Var.Mtot)*100                                                                   
    
    return (Var)
#%%
    
 #Function for the estimation of crystal size distribution - Part 1 : Initialization of the local parameters and time interval 
    
def CrystalSizeDistribution():
    global sigma_compound1,DT_ice,V_ice_diff,V_compound1_cr_diff,M_ice_diff,M_compound1_cr_diff,M_ice_reactions,M_compound1_cr_reactions,Compound1,Compound2,tck,lm,Linear_regr        
        
    #Time interval
    Var.t=Var.t+deltat
    
    atot_ice=0                                                                                                      #Total surface area of ice crystals in the reactor in every loop: m^2
    atot_compound1_cr=0                                                                                             #Total surface area of compound 1 crystals in the reactor in every loop: m^2
    
    #Volume and mass difference(increase/decrease) of crystals in every time interval
    V_ice_diff=0                                                                                                    #Volume difference of ice crystals in every time interval : m^3
    V_compound1_cr_diff=0                                                                                           #Volume difference of compound 1 crystals in every time interval : m^3
    
    M_ice_diff=0                                                                                                    #Mass difference of ice crystals in every time interval : kg
    M_compound1_cr_diff=0                                                                                           #Mass difference of compound 1 crystals in every time interval : kg
    
    #Mass of crystals produced/consumed via reacting : growth and/or nucleation and/or dissolution in every time interval
    M_ice_reactions=0                                                                                               #Mass of ice crystals produced/consumed via reacting : growth and/or nucleation and/or dissolution in every time interval : kg
    M_compound1_cr_reactions=0                                                                                      #Mass of compound 1 crystals produced/consumed via reacting : growth and/or nucleation and/or dissolution in every time interval : kg
    
    #Number of crystals produced by nucleation in every time interval
    Nb_ice=0
    Nb_compound1_cr=0
    
    #Number of crystals that  grow in OR out of every class size in every time interval (ASSUMED POSITIVE for growth INTO a class and NEGATIVE for growth OUT of a class)
    Ng_ice=0
    Ng_compound1_cr=0   

    #Number of crystals dissolved in OR out of every class size in every time interval (ASSUMED POSITIVE for dissolution INTO a class and NEGATIVE for dissolution OUT of a class)
    Nd_ice=0
    Nd_compound1_cr=0
    
    #Calculation of the equilibrium values using interpolation/extrapolation of the solubility and ice lines
    C_compound1_eq=interpolate.splev(Var.Treactor, tck, der=0)
    
    C_compound1_eq=float(C_compound1_eq)
    print("C_compound1_eq is "+str(C_compound1_eq))
          
    T_ice_eq_array=lm.predict(Var.C_compound1)
    T_ice_eq=T_ice_eq_array[0]
    print("T_ice_eq is "+str(T_ice_eq))
    
    #Compound 1 Supersaturation  
    sigma_compound1=(Var.C_compound1/C_compound1_eq)-1
    
    #Ice undercooling
    DT_ice=T_ice_eq-Var.Treactor
   
    #%%
    
    #Function for the estimation of crystal size distribution - Part 2 : Compound 1 crystal size distribution / Mass balance
    if sigma_compound1<0:                                                                                             #If there is no supersaturation,the system is undersaturated and there is only dissolution of compound 1 crystals
        #Compound 1 crystals growth rate : m/s
        G_compound1_cr=0
        
        #Compound 1 crystals nucleation rate : crystals/s 
        B_compound1_cr=0
        
        #Compound 1 crystals dissolution rate : m/s
        D_compound1_cr=Compound1.kd*(-sigma_compound1)
        
    else:
        
        G_compound1_cr=Compound1.kg*(sigma_compound1**2)
        B_compound1_cr=Compound1.kb*Var.M_compound1_cr*(sigma_compound1**2)
        D_compound1_cr=0
    
    #0 compound 1 crystal class size due to secondary nucleation (neglecting primary nucleation : ASSUMPTION)
    Nb_compound1_cr=B_compound1_cr*deltat
    Ng_compound1_cr=-(G_compound1_cr*Var.N_compound1_cr[0]*deltat)/deltal
    Nd_compound1_cr=-(D_compound1_cr*Var.N_compound1_cr[0]*deltat)/deltal
   
    Var.N_compound1_cr[0]=Var.N_compound1_cr[0]+Nb_compound1_cr+Ng_compound1_cr+Nd_compound1_cr-(Fout_compound1_cr*deltat*Var.N_compound1_cr[0])/Var.M_compound1_cr
    Var.V_compound1_cr_distr[0]=Var.N_compound1_cr[0]*Compound1.kv*0.25*(deltal**4)/deltal
    atot_compound1_cr=atot_compound1_cr+Var.N_compound1_cr[0]*Compound1.ka*(deltal**2)*(0**2)
    
    #1-298 crystal class size                                                                                    #The range of the loop is from 1 to classes-1=299 due to python syntax, the loop actually runs from 1 to 298
    for i in range(1,(classes-1)):
        Ng_compound1_cr=(G_compound1_cr*deltat*(Var.N_compound1_cr[i-1]-Var.N_compound1_cr[i]))/deltal
        Nd_compound1_cr=(D_compound1_cr*deltat*(Var.N_compound1_cr[i+1]-Var.N_compound1_cr[i]))/deltal

        Var.N_compound1_cr[i]=Var.N_compound1_cr[i]+Ng_compound1_cr+Nd_compound1_cr-(Var.N_compound1_cr[i]*Fout_compound1_cr*deltat)/Var.M_compound1_cr
        Var.V_compound1_cr_distr[i]=Var.N_compound1_cr[i]*Compound1.kv*0.25*((((i+1)*deltal)**4)-((i*deltal)**4))/deltal
        atot_compound1_cr=atot_compound1_cr+Var.N_compound1_cr[i]*Compound1.ka*(deltal**2)*(i**2)     
        
    #299 crystal class size
    Ng_compound1_cr=(G_compound1_cr*deltat*Var.N_compound1_cr[classes-2])/deltal
    Nd_compound1_cr=-(D_compound1_cr*deltat*Var.N_compound1_cr[classes-1])/deltal
    
    Var.N_compound1_cr[classes-1]=Var.N_compound1_cr[classes-1]+Ng_compound1_cr+Nd_compound1_cr-(Var.N_compound1_cr[classes-1]*deltat*Fout_compound1_cr)/Var.M_compound1_cr
    Var.V_compound1_cr_distr[classes-1]=Var.N_compound1_cr[classes-1]*Compound1.kv*0.25*(((classes*deltal)**4)-(((classes-1)*deltal)**4))/deltal
    atot_compound1_cr=atot_compound1_cr+Var.N_compound1_cr[classes-1]*Compound1.ka*(deltal**2)*((classes-1)**2)  
        
    #Volume difference(increase/decrease) of compound 1 crystals in this interval 
    V_compound1_cr_diff=sum(Var.V_compound1_cr_distr)-Var.V_compound1_cr_previous_loop
    
    #Mass difference(increase/decrease) of compound 1 crystals in this interval 
    M_compound1_cr_diff=V_compound1_cr_diff*Compound1.d
    
    #Mass of compound 1 crystals produced/consumed via reacting : growth and/or nucleation and/or dissolution in this time interval
    M_compound1_cr_reactions=M_compound1_cr_diff-Fin_compound1_cr*deltat+Fout_compound1_cr*deltat
    
     #%%
    
     #Function for the estimation of crystal size distribution - Part 3 :Ice crystal size distribution / Mass balance
    if DT_ice<0:                                                                                                 #If there is no ice undercooling, there is only dissolution of ice crystals
        #Ice crystals growth rate : m/s
        G_ice=0
        
        #Ice crystals nucleation rate : crystals/s 
        B_ice=0
        
        #Ice crystals dissolution rate : m/s
        D_ice=kd_ice*(-DT_ice)
        
    else:
        
        G_ice=kg_ice*DT_ice
        B_ice=kb_ice*atot_ice*(DT_ice**2)
        D_ice=0
    
    #0 ice crystal class size due to secondary nucleation (neglecting primary nucleation : ASSUMPTION)
    Nb_ice=B_ice*deltat
    Ng_ice=-(G_ice*Var.N_ice[0]*deltat)/deltal
    Nd_ice=-(D_ice*Var.N_ice[0]*deltat)/deltal

    Var.N_ice[0]=Var.N_ice[0]+Nb_ice+Ng_ice+Nd_ice-(Fout_ice*deltat*Var.N_ice[0])/Var.M_ice
    Var.V_ice_distr[0]=Var.N_ice[0]*kv_ice*0.25*(deltal**4)/deltal
    atot_ice=atot_ice+Var.N_ice[0]*ka_ice*(deltal**2)*(0**2)
     
    #1-298 crystal class size                                                                                    #The range of the loop is from 1 to classes-1=299 due to python syntax, the loop actually runs from 1 to 298
    for i in range(1,(classes-1)):
        Ng_ice=(G_ice*deltat*(Var.N_ice[i-1]-Var.N_ice[i]))/deltal
        Nd_ice=(D_ice*deltat*(Var.N_ice[i+1]-Var.N_ice[i]))/deltal
        
        Var.N_ice[i]=Var.N_ice[i]+Ng_ice+Nd_ice-(Var.N_ice[i]*Fout_ice*deltat)/Var.M_ice
        Var.V_ice_distr[i]=Var.N_ice[i]*kv_ice*0.25*((((i+1)*deltal)**4)-((i*deltal)**4))/deltal
        atot_ice=atot_ice+Var.N_ice[i]*ka_ice*(deltal**2)*(i**2)
        
    #299 crystal class size
    Ng_ice=(G_ice*deltat*Var.N_ice[classes-2])/deltal
    Nd_ice=-(D_ice*deltat*Var.N_ice[classes-1])/deltal
      
    Var.N_ice[classes-1]=Var.N_ice[classes-1]+Ng_ice+Nd_ice-(Var.N_ice[classes-1]*deltat*Fout_ice)/Var.M_ice
    Var.V_ice_distr[classes-1]=Var.N_ice[classes-1]*kv_ice*0.25*(((classes*deltal)**4)-(((classes-1)*deltal)**4))/deltal
    atot_ice=atot_ice+Var.N_ice[classes-1]*ka_ice*(deltal**2)*((classes-1)**2)
    
    #Volume difference(increase/decrease) of ice crystals in this interval
    V_ice_diff=sum(Var.V_ice_distr)-Var.V_ice_previous_loop
    
    #Mass difference(increase/decrease) of ice crystals in this interval
    M_ice_diff=V_ice_diff*d_ice
    
    #Mass of ice crystals produced/consumed via reacting : growth and/or nucleation and/or dissolution in this time interval
    M_ice_reactions=M_ice_diff-Fin_ice*deltat+Fout_ice*deltat

    return (Var)
    #%%
    
    #Function for the calculation of the mass balances for the different compounds
    
def MassBalances():    
    
    global V_ice_diff,V_compound1_cr_diff,M_ice_diff,M_compound1_cr_diff,M_ice_reactions,M_compound1_cr_reactions
    
    #Calculation of new values based on mass balances
    
    #Calculation of the new mass for every compound in the reactor
    Var.M_ice=Var.M_ice+M_ice_diff
    Var.M_compound1_cr=Var.M_compound1_cr+M_compound1_cr_diff
    #print("M_compound1_cr is "+str(Var.M_compound1_cr))
    Var.M_solids_tot=Var.M_ice+Var.M_compound1_cr
    print("Var.M_compound1 is "+str(Var.M_compound1))  
    Var.M_water=Var.M_water+Fin_water*deltat-Fout_water*deltat-M_ice_reactions                                                                
    Var.M_compound1=Var.M_compound1+Fin_compound1*deltat-Fout_compound1*deltat-M_compound1_cr_reactions   

                             
    Var.M_compound2=Var.M_compound2+Fin_compound2*deltat-Fout_compound2*deltat                                       
    Var.M_mother_liquor=Var.M_water+Var.M_compound1+Var.M_compound2             
    Var.M_liquid_tot=Var.M_mother_liquor
    
    Var.Mtot=Var.M_solids_tot+Var.M_liquid_tot
    print("Mtot in the loop is "+str(Var.Mtot))
    
    #Calculation of the new volume for every compound in the reactor
    Var.V_ice=Var.V_ice+V_ice_diff
    Var.V_compound1_cr=Var.V_compound1_cr+V_compound1_cr_diff
    Var.V_solids_tot=Var.V_ice+Var.V_compound1_cr

    Var.V_water=Var.M_water/d_water                                                                
    Var.V_compound1=Var.M_compound1/Compound1.d                                  
    Var.V_compound2=Var.M_compound2/Compound2.d                                       
    Var.V_mother_liquor=Var.V_water+Var.V_compound1+Var.V_compound2             
    Var.V_liquid_tot=Var.V_mother_liquor
    
    Var.Vtot=Var.V_solids_tot+Var.V_liquid_tot
    
    #New concentration of Na2SO4 in the bulk liquid in the reactor
    Var.C_compound1=(Var.M_compound1*1000)/(Compound1.MW*Var.M_water)                                               # Units : molals = mol compound 1/kg water
    Var.C_compound1_mliq=Var.C_compound1
    
    #Store the new volumes of ice and compound 1 crystals for the calculation of V_diff in the next time interval
    Var.V_ice_previous_loop=Var.V_ice
    Var.V_compound1_cr_previous_loop=Var.V_compound1_cr
    
    #Calculation of the new solids contents in the reactor
    Var.solids_weight_percentage=(Var.M_solids_tot/Var.Mtot)*100 
    
    return (Var)
  #%%
 
  #Function for the calculation of the heat balance 
def HeatBalances():
    global M_ice_reactions,M_compound1_cr_reactions, Qtot
    
    #Heat balance - New temperature in the reactor (!!!Assumed heat flux into the system is POSITIVE)
    Qcryst_ice=(M_ice_reactions*DHfus_ice)/deltat                                                                     #Heat flux to the system due to ice crystallization : J/s
    Qcryst_compound1_cr=(M_compound1_cr_reactions*Compound1.DHfus)/deltat                                             #Heat flux to the system due to compound 1 crystallization : J/s
    
    Qtot=Qcryst_ice+Qcryst_compound1_cr-Cooling_rate                                                                  #Total heat flux to the system : J/s
    print("Qcryst_ice is "+str(Qcryst_ice))
    print("Qcryst_compound1_cr is "+str(Qcryst_compound1_cr))
    print("Qtot is "+str(Qtot))
    Na2SO4_weight_percentage=(Var.M_compound1/Var.M_liquid_tot)*100                                                     #Na2SO4 weight percentage in the liquid in the reactor for the estimation of the heat capacity of the solution : Na2SO4 wt % = (total mass of Na2SO4 in the liquid in the reactor / total mass of liquid in the reactor)*100 %
    print("Na2SO4_weight_percentage is "+str(Na2SO4_weight_percentage))
    print("Var.M_compound1 is "+str(Var.M_compound1))
    print("Var.M_liquid_tot is "+str(Var.M_liquid_tot))

    Cp_solution=(0.9988-0.006494*Na2SO4_weight_percentage+0.00003025*(Na2SO4_weight_percentage**2)-0.0000001286*(Na2SO4_weight_percentage**3))*4.184*1000     #Heat capacity of the solution (Calculated for 25 oC, Neglecting the impurities and assuming no change of heat capacity due to temperature) (source:Heat of solution, heat capacity...) :   J/(kg*K) 
    
    Var.Treactor=(Qtot*deltat)/(Var.M_liquid_tot*Cp_solution+Var.M_compound1_cr*Compound1.Cp+Var.M_ice*Cp_ice)+Var.Treactor
    print("Treactor is "+str(Var.Treactor))
    return (Var, Qtot, Cp_solution)
#%%
        
#Function for storing several values        
def Store():
    global t_store,sigma_compound1_store,DT_ice_store,M_ice_store,M_compound1_cr_store,Treactor_store   
    t_store=np.append(t_store,Var.t)
    sigma_compound1_store=np.append(sigma_compound1_store,sigma_compound1)
    DT_ice_store=np.append(DT_ice_store,DT_ice)
    M_ice_store=np.append(M_ice_store,Var.M_ice)
    M_compound1_cr_store=np.append(M_compound1_cr_store,Var.M_compound1_cr)
    Treactor_store=np.append(Treactor_store,Var.Treactor)

 
#%%
#Function for plotting the initial values of several variables    
    
def InitialPlot():
    global graph1,graph2,graph3,graph4

    graph1=fig1.add_subplot(311)
    graph1.set_title('Number of ice crystals in each class')
    graph1.set_xlabel('Class Length (m)')
    graph1.set_ylabel('Number of ice crystals')   
    graph1.plot(class_length, Var.N_ice, 'r-',label='Initial values')
    
    graph2=fig1.add_subplot(313)
    graph2.set_title('Number of Compound 1 crystals in each class')
    graph2.set_xlabel('Class Length (m)')
    graph2.set_ylabel('Number of Comp. 1 crystals')   
    graph2.plot(class_length, Var.N_compound1_cr, 'r-',label='Initial values')
    
    graph3=fig2.add_subplot(311)
    graph3.set_title('Volume of ice crystals in each class')
    graph3.set_xlabel('Class Length (m)')
    graph3.set_ylabel('Volume of ice crystals (m^3)')   
    graph3.plot(class_length, Var.V_ice_distr, 'r-',label='Initial values')
        
    graph4=fig2.add_subplot(313)
    graph4.set_title('Volume of Compound 1 crystals in each class')
    graph4.set_xlabel('Class Length (m)')
    graph4.set_ylabel('Volume of Comp. 1 crystals (m^3)')   
    graph4.plot(class_length, Var.V_compound1_cr_distr, 'r-',label='Initial values')   
 #%%
#Function for plotting final results
    
def Plot():
    global t_store,sigma_compound1_store,DT_ice_store,M_ice_store,M_compound1_cr_store,Treactor_store,graph1,graph2,graph3,graph4 
    
    graph1.plot(class_length,Var.N_ice,'b-',label='Final values')

    graph2.plot(class_length,Var.N_compound1_cr,'b-',label='Final values')
    fig1.legend(loc='center right')
    fig1.show
    
    graph3.plot(class_length,Var.V_ice_distr,'b-',label='Final values')

    graph4.plot(class_length,Var.V_compound1_cr_distr,'b-',label='Final values')
    fig2.legend(loc='center right')
    fig2.show

    graph5=fig3.add_subplot(111)
    graph5.plot(t_store,sigma_compound1_store,'b-')
    graph5.set_title('Supersaturation of Compound 1 vs time')
    graph5.set_xlabel('Time(sec)')
    graph5.set_ylabel('Superasturation of Compound 1 (-)')
    fig3.show
    
    graph6=fig4.add_subplot(111)
    graph6.plot(t_store,DT_ice_store,'b-')
    graph6.set_title('Undercooling of ice vs time')
    graph6.set_xlabel('Time(sec)')
    graph6.set_ylabel('Undercooling of ice (K)')
    fig4.show
    
    graph7=fig5.add_subplot(311)
    graph7.set_title('Mass of ice crystals vs time')
    graph7.set_xlabel('Time(sec)')
    graph7.set_ylabel('Mass of ice crystals (kg)')
    graph7.plot(t_store,M_ice_store,'b-')
    
    graph8=fig5.add_subplot(313)
    graph8.set_title('Mass of Compound 1 crystals vs time')
    graph8.set_xlabel('Time(sec)')
    graph8.set_ylabel('Mass of Comp. 1 crystals (kg)')
    graph8.plot(t_store,M_compound1_cr_store,'b-')
    fig5.show
    
    graph9=fig6.add_subplot(111)
    graph9.set_title('Temperature in the reactor vs time')
    graph9.set_xlabel('Time(sec)')
    graph9.set_ylabel('Temperature (K)')
    graph9.plot(t_store,Treactor_store,'b-')
 #%%
# #calculate final concentration for all ions 
class conc_f:
    def __init__(self, Cna_in, Ccl_in, Ck_in, Cmg_in, Cca_in, Cso4_in, Qin, Min, Mout, M_compound1_cr2):
        self.Cna_in=Cna_in
        self.Ccl_in=Ccl_in
        self.Ck_in=Ck_in
        self.Cmg_in=Cmg_in
        self.Cca_in=Cca_in
        self.Cso4_in=Cso4_in
        self.Qin=Qin
        self.Min=Min
        self.Mout=Mout
        self.M_compound1_cr2=M_compound1_cr2
        
        
    def calc_f_c(self):
        print("self.Min is "+str(self.Min))
        print("self.Mout is "+str(self.Mout))
        print("self.M_compound1_cr2 is "+str(self.M_compound1_cr2))
        self.Ccl_out=self.Ccl_in*self.Min/self.Mout
        self.Ck_out=self.Ck_in*self.Min/self.Mout
        self.Cmg_out=self.Cmg_in*self.Min/self.Mout
        self.Cca_out=self.Cca_in*self.Min/self.Mout
        self.Mna_out=self.Cna_in*self.Min/(1000*d_in/1000)-self.M_compound1_cr2*2*constants.MW_Na/MW_Na2SO4_10H2O        
        self.Cna_out=((self.Cna_in*self.Min/(1000*d_in/1000)-self.M_compound1_cr2*2*constants.MW_Na/MW_Na2SO4_10H2O)/self.Mout)/1000#g/l
        self.Cso4_out=((self.Cso4_in*self.Min/(1000*d_in/1000)-self.M_compound1_cr2*constants.MW_so4/MW_Na2SO4_10H2O)/self.Mout)/1000#g/l
        print("Cna_out is "+str(self.Cna_out))
        print("Cso4_out is "+str(self.Cso4_out))
    
#%%
    
#Main Procedure

# Compounds()
C_Na2SO4_sol_line, tck = Compounds()
print("tck is " + str(tck))

print(C_Na2SO4_sol_line)

Var=Initialize(Var)
  
InitialPlot()

while Var.t<=55000:
        if Var.solids_weight_percentage>=40:                      #Set a maximum limit based on the solids weight percentage (wt%)
            print('Time:',Var.t,'sec,','Maximum solids weight percentage in the reactor reached : 40 wt% solids')
            break
        else:
            Var=CrystalSizeDistribution()
            Var=MassBalances()
            Var,Qtot, Cp_solution=HeatBalances()
            if ((Var.t)%100)==0:
                Store()
            else:
                pass
Plot()
print("Var.t is " +str(Var.t))
print("Var.solids_weight_percentage is " +str(Var.solids_weight_percentage))

#Calculate and print the elapsed time after the calculations
elapsed=time.time()-time0
print('Time elapsed : ',elapsed)#Calculate and print the elapsed time after the calculations
Residence_time=Var.t/3600
res=Residence_time#Vreactor/(Qf)
print(res)
print('Residence Time : ', Residence_time,'hours')
print('Volume reactor : ', Vreactor*1000,'liters')
Flow_rate=Qf/(d_in)/Residence_time
print('Continous Flow Rate : ', Flow_rate,'liters/hour')
Min=Qf#Residence_time #kg/hr 
#mass flowrate of na2so4 
M_compound1_cr=Var.M_compound1_cr*MW_Na2SO4_10H2O/MW_Na2SO4/res#Residence_time
print("mass flowrate of na2so4 "+ str(M_compound1_cr)+"kg/hr")
#mass flow rate of ice 
M_ice= Var.M_ice/res#Residence_time
print("ice mass flow rate " + str(M_ice)+ "kg/hr")
#Mass flow rate for effluent

Mout = Var.M_liquid_tot.item() / res # Residence_time
print("Mass flow rate for effluent is "+ str(Mout)+"kg/hr")


conc_f_in=conc_f(C_i_in[0],C_i_in[1], C_i_in[2], C_i_in[3], C_i_in[4], C_i_in[5],Qf/d_in, Qf, Var.M_liquid_tot.item(), M_compound1_cr)

dens=1.02
conc_f_in.calc_f_c()
Cout= Var.C_compound1_mliq #out concentration for na2so4. only na, so4 change 
Cout_so4= Cout.item()
Cout_cl=Var.C_compound2
Cout_cl= Cout_cl.item()
Cout_na= conc_f_in.Cna_out
Cout_so4= conc_f_in.Cso4_out
Cout_efc_g=[Cout_na+conc_f_in.Ccl_out/constants.MW_cl*constants.MW_Na, conc_f_in.Ccl_out, conc_f_in.Ck_out, 0, 0, Cout_so4] #g/l 
print("Cout_efc is "+str(Cout_efc_g))


M_solids_tot=Var.M_solids_tot/res#Residence_time
print("M_solids_tot is "+str(M_solids_tot))
Mtot=Var.Mtot/res#Residence_time
print("Mtot is " + str(Mtot))
M_ice=Var.M_ice/res#Residence_time
print("M_ice is "+str(M_ice))

V_compound1=Var.V_compound1/res
print("V_compound1 is "+str(V_compound1))

#operating temperature 
Treactor=Var.Treactor-273
print("operating temperature is "+str(*Treactor[0])+"oC")

#required energy 
E_el=Qtot/1000*res # KW 
Q_t=abs(Qtot/1000)
print("Q_t"+str(Q_t))
E_p=(Min*dp+Mout*dp+(M_ice+M_compound1_cr)*dp_slurry)*1e5/3600/(1000*npump)/1000 #kwh
E_fil=scaleup.scaleup(0.5, 0.3*1000, Qf)
E_t=Q_t+E_p+E_fil
print("E_t"+str(E_t))
sec_f=E_t/Qf
sec_eff=E_t/(Mout/1020)
sec_product=E_t/M_compound1_cr
w_in_cryst=M_compound1_cr*10*18.01528/MW_Na2SO4_10H2O
print(w_in_cryst)
print("M_ice*res is "+str(M_ice*res))
print("Mout*res is "+str(Mout*res))
print("M_compound1_cr*res is "+str(M_compound1_cr*res))

M_efc_out=Var.M_liquid_tot.item()
M_compound1_cr=M_compound1_cr*res
M_ice=M_ice*res

