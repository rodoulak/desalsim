# Tutotrial: Thermal crystallizer  (TCr) unit

TCr is a thermal based process that is used to crystallize salt. In this work, TCr aims to recover high quality water from brine solutions and crystallize salts.
![tcr](https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/assets/150446818/1cd7f730-e4fa-4a0c-b332-822899e02433)


In **DesalSim** package, the TCr unit is used to model the operation of a Thermal crystallizer  technology. Upon simulation, it will generate the influent/effluent mass flows and their concentrations, the cooling water, and the energy requirements.
The TCr function consists of three classes: ['thermal_calc'](#use-thermal_calc-class), ['conc_cal'](#use-conc_cal-class), ['calculate_energy'](#use-calculate_energy-class).  
In this tutorial, we will focus on how to use the the classes and their methods. 

**Table 1** gives an overview of the main inputs and outputs for each process unit of Thermal crystallizer . 
| Process                                   | Input                                       | Output                                                |
|-------------------------------------------|---------------------------------------------|-------------------------------------------------------|
| Thermal crystallizer                  | Feed flow rate [m3/h]                        | Flow rate of water [kg/h]                                  |
|                                        | Ion concentration [g/L]                      | Flow rate of NaCl [kg/h]                                   |
|                                        | Feed temperature [°C]                        | Cooling water flow rate [m3/h]                             |
|                                        | Steam temperature [°C]                       | Electrical [kWhel] and thermal [kWhth] requirements        |

The mathematical description of Thermal Crystallizer technology is given in [Mathematical description](https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/tree/main/paper/Mathematical_description.pdf), see Section A.3. 

## 1. Getting started 
### 1.1. Import class
```python
import Desalsim
```
Then import the class:  
```python
from Desalsim.thermal_cryst_f import thermal_calc
from Desalsim.thermal_cryst_f import conc_cal
from Desalsim.thermal_cryst_f import calculate_energy
```
Additionally, function for calculating density (`density_calc.py`) need to be imported or constants (comparison.py) where user can add constant values like MW, prices etc, need to be imported.
```python
from Desalsim.density_calc import density_calc 
from Desalsim import constants 
from Desalsim import scaleup
```
### 1.2. Define feed characteristics
You can initialize the feed solution by setting the flow rate, specifying the focus components and their concentration. 
```python
  # Feed concentration
components = ['Na', 'Cl', 'K', 'Mg', 'Ca', 'SO4']
Cf_tcr_in=[80.42, 116.69, 2.29, 0.01, 0.04, 0.54]
# Total salt concentration (g/L)
Cf_s=sum(Cf_tcr_in)

# CaSO4 concentration in feed solution (g/L)
Cf_caso4=Cf_tcr_in[4]*MW_Caso4/MW_Ca

#Feed flow rate 
Qf=1000 #kg/h

# input conditions 
T_in=40 #oC
```
Note that if you want to add more components, you need to update the components list and include the concentration of the new component in the _Cf_tcr_in_. 

You can calculate the density of the feed solution:
```python
d_sol=density_calc(T_in, Cf_s)/1000 #density of feed

```
### 1.3. Set operating assumptions  
You need to set operating assumptions related to temperatures such as the temperature in the last effect, the intake/outake cooling water temperature. Additionally, the moinsture in salt stream needs to be set.
```python
# Assumptions:
T_in=40 # (oC)
T_cw_f=25 #intake cooling water temperature (oC)
T_cw_o=40 #out cooling water temperature (oC)
T_top=60 # Operating temperature (oC)

salt_mois=20 # % moinsture in salt stream (units: %)
```
#### Thermodyanamic data 
The latent heat of motive steam and vapor is calculated based on steam tables. 
```python
# 𝐻𝑒𝑎𝑡 𝑇𝑟𝑎𝑛𝑠𝑓𝑒𝑟 𝐶𝑜𝑛𝑠𝑡𝑎𝑛t
UA=45990
LHV_v=2199.7#2107.92 # kj/kg (gathered from steam table)
LHV_s=2357.69 # kj/kg (gathered from steam table)
```
Finally, you need to set assumptions related to pumping like pressure drop (_dp_) and pump efficiency (_npump_). 
```python
dp=0.1 # pressure drop (units: bar)
dp_slurry=3.5 #pressure drop for slurry streams (units: bar)
dp_f=3.5 #pressure drop for feed (units: bar)
dp_w=1 #pressure drop for qater streams (units: bar)
dp_cw=2 #pressure drop for cooling water (units: bar)

npump=0.8 #pump efficiency (units: -)
```
### 1.4. Set constants 
You need to set constant parameters like the specific heat capacity of water
```python
Cp_f=3.14 # Feed specific heat capacity (units: KJ* Kg*oC)
CP_cw=4.187 # Water specific heat capacity (units: KJ* Kg*oC)
```
After setting all the required inputs, then you can create the functions' objectives. 

## 2. Use thermal_calc class   
thermal_calc is a class used to represent Mass Balance for Thermal Crystallizer Unit. In particular, it calculates the evaporated mass and mass of recovered solid, the cooling water requirements and the energy requirements. 
thermal_calc takes as input the operating temperature (_T_op_), the inlet flow rate (_Qf_), the total ion concentration of solution (_Cf_s_), the Concentration of CaSO4 in the solution (_Cf_caso4_), the  inlet feed temperature (_T_in_), the list of concentration of ions in the solution (_Cf_in_), the salt moisture in slurry stream (_salt_mois_),  the latent heat of vapor (_LHV_v_),  the latent heat of steam (_LHV_s_), Temperature of intake/outake cooling water (_T_cw_o, T_cw_f_).
### 2.1. Overview 
The following attributes are available within the thermal_calc class:  
- `T_op`: Operating temperature (°C)
- `Qf`: Inlet flow rate (kg/hr)
- `Cf_s`: Total ion concentration of solution (g/L)
- `Cf_caso4`: Concentration of CaSO4 in the solution (g/L)
- `T_in`: Inlet feed temperature (°C)
- `d_sol`: Inlet feed density (kg/L)
- `Cf_in`: List of concentration of ions in the solution (g/L)
- `nacl_sat_wt`: Weight percent of NaCl at saturation
- `nacl_sol`: NaCl in the solution (kg/h)
- `salt_solids`: Salt solids (kg/h)
- `caso4_sol`: CaSO4 in the solution (kg/h)
- `solid_mass`: Solid mass (kg/h)
- `ev_mass`: Evaporated mass (kg/h)
- `BPE`: Boiling point elevation (°C)
- `heat_req`: Heat required (kJ/h)
- `T_s`: Steam temperature (°C)
- `steam_mass`: Steam mass (kg/h)
- `cw_mass`: Cooling water mass (kg/h)


The thermal_calc class provides the following two methods:
```python
# This method calculates the mass balance for thermal crystallization
mass_bal_cryst()
```
```python
# This method calculates the heat balance for thermal crystallization
heat_bal_cryst()
```
They are called upon initialization of the class instance.

### 2.2. Create thermal_calc objects
thermal_calc takes as input the operating temperature (_T_op_), the inlet flow rate (_Qf_), the total ion concentration of solution (_Cf_s_), the Concentration of CaSO4 in the solution (_Cf_caso4_), the  inlet feed temperature (_T_in_), the list of concentration of ions in the solution (_Cf_in_), the salt moisture in slurry stream (_salt_mois_),  the latent heat of vapor (_LHV_v_),  the latent heat of steam (_LHV_s_), Temperature of intake/outake cooling water (_T_cw_o, T_cw_f_).   
```python
# Create an instance of the thermal_calc class
th_cryst_dat=thermal_calc(T_op, Qf, Cf_s, Cf_caso4, T_in, Cf_tcr_in, salt_mois, LHV_v, LHV_s, T_cw_o, T_cw_f)
```

### 2.3. Use 'mass_bal_cryst' method
The method calculates the mass balance for thermal crystallization.   In particular, it calculates the recovered salt flow rate (_solid_mass_) and the evaporation mass flow rate (_ev_mass_). 
```python
th_cryst_dat.mass_bal_cryst()
```
It doesn't take additional inputs. 
### 2.3.1. Assigned the results to output parameters 
After the mass calculation, you can assigned the results to output parameters: 
```python
# Recovered salt flow rate (kg/h)
M_Nacl=th_cryst_dat.solid_mass

# Evaporation mass flow rate (water recovery), (kg/h)
Q_evap_mass=th_cryst_dat.ev_mass
```
### 2.4. Use ''heat_bal_cryst'' method 
The method calculates the energy balance for thermal crystallization. In particular, it calculates the steam_mass (_steam_mass_), and the total distillate flow rate (_cw_mass_) and required heat (_heat_req_). 

```python
th_cryst_dat.heat_bal_cryst()
```
It doesn't take additional inputs. 
### 2.4.1. Assigned the results to output parameters 
After the mass calculation, you can assigned the results to output parameters: 
```python
#Cooling water requirements (kg/h)
Qcw_tcr=th_cryst_dat.cw_mass

#Calculate energy consumption 
heat_req=th_cryst_dat.heat_req
```

### 2.5. Print results 
You can print results from mass and energy calculations. 
```python
print("Mass flowrate of recovered salt is "+str(round(M_Nacl,2))+"kg/hr")
print("Mass flowrate of recovered water is "+str(round(Q_evap_mass,2))+"kg/hr")

#Cooling water requirements (kg/h)
Qcw_tcr=th_cryst_dat.cw_mass
```
Mass flowrate of recovered salt is 0.22kg/hr  
Mass flowrate of recovered water is 999.78kg/hr  

Mass flowrate of required cooling water is 35016.66kg/hr

## 3. Use conc_cal class 
conc_cal is a class used to calculate the concentration of different ions in a solution. For this calculation, you need to use the solution flow rate (_Qf_), the mass of solid stream (_solid_mass_), the name of ion and its concentration in the stream (_C1, Cc1, C2, Cc2, C3, Cc3, C4, Cc4, C5, Cc5, C6, Cc6_).   
### 3.1. Oveview
The following attributes are available within the conc_cal class:  
- `Qf`: Flow rate (units: m³/h)
- `solid_mass`: Mass of solid stream (units: kg)
- `Cc1`, `Cc2`, `Cc3`, `Cc4`, `Cc5`, `Cc6`: Concentrations of six different ions (units: ppm)
- `d_sol`: Density of the solution at 40°C (units: kg/m³)


The conc_cal class provides the following method:
```python
conc_cal()
```
This class consists of three methods that calculate the molarity of different ions in the solution ('molarity'), the concentration of different salts in the solution ('conc_salt_comp') and the the concentration of different salts in the outlet stream ('salt_conc'). 

### 3.2. Create conc_cal objectives and call functions 

```python
 Calculation of the concentration of different ions in the solution
th_cryst_dat_2=conc_cal(Qf, M_Nacl , 'Na',Cf_tcr_in[0], 'cl',Cf_tcr_in[1],'k', Cf_tcr_in[2], 'mg', Cf_tcr_in[3], 'ca', Cf_tcr_in[4], 'so4', Cf_tcr_in[5], T_in)
th_cryst_dat_2.molarity()
th_cryst_dat_2.conc_salt_comp()
th_cryst_dat_2.salt_conc()
```
### 3.2. Assigned the results to output parameters
You can assigned the results to output parameters:
```python
#Salt stream concentration (g/L)
Csalt_out=[th_cryst_dat_2.CNa,th_cryst_dat_2.CCl, th_cryst_dat_2.CK, th_cryst_dat_2.CMg, th_cryst_dat_2.CCa, th_cryst_dat_2.CSO4]

```

## 4. Use 'calculate_energy' class
'calculate_energy' is a class used to represent the calculation of energy consumption and the specific energy consumption for Thermal Crystallizer Unit. For this calculation, the flow rate of the feed solution (_Qf_), the distillate flow rate (_`Q_evap_mass_), flow rate of cooling water (_Qcw_), the mass of recovered NaCl (_M_Nacl_), the required heat (_heat_req_) calculated in "thermal_calc", density of the feed solution (_d_sol_) are used. In addition, the 'calculate_energy' takes as input the expected pressure drop in cooling water stream, slurry stream and water stream (_dp_cw, dp_slurry, dp_w_) and the pump efficiency (_npump_).
### 4.1. Oveview 
The following attributes are available within the 'calculate_energy' class:  
- `Qf`: Flow rate of the feed solution (m³/h).
- `Q_evap_mass`: Distillate flow rate (m³/h).
- `Qcw`: Flow rate of cooling water (m³/h).
- `d_sol`: Density of the feed solution (kg/m³).
- `npump`: Pump efficiency.

The  'calculate_energy' class provides the following method:
```python
calculate_energy_consumption()
```
This method returns the the Specific energy consumption for electrical energy (_SEC_el_), Specific energy consumption for electrical energy per unit mass of evaporated water in kWh/kg (_SEC_el_prod_), Specific energy consumption for electrical energy per unit mass of NaCl produced in kWh/kg (_SEC_el_prod2_)

### 4.2. Create nf_energy objectives, calculate Energy Consumption and assign the results to output parameters 
The following objective is created for energy consumption. Assumptions for pressure drop and pump efficiency need to be made. 

```python
E_el_th_Cr, E_th_th_Cr, SEC_el_f, SEC_el_NaCl, SEC_el_w = calculate_energy(Qf, Q_evap_mass, Qcw_tcr, M_Nacl, heat_req, d_sol, dp_f, dp_w, dp_slurry, dp_cw, npump)

```
### 4.3. Print results 
You can print results from energy calculations. The specific energy consumption is also calculated so you can validate easier the results. 
```python
print(f"SEC_el_prod is {SEC_el_w} KWh/m3 product")
print(f"SEC_el_prod2 is {SEC_el_NaCl} KWh/kg of NaCl product")
print("SEC_th_prod is "+str( round((E_th_th_Cr/(Qf/1000)),2))+" KWh/m3 intake")
```
SEC_el_prod is 5.02 KWh/m3 product  
SEC_el_prod2 is 23.06 KWh/kg of NaCl product  
SEC_th_prod is 628.34 KWh/m3 intake  

## References 



