# Tutotrial: Multi-plug flow reactor (MF-PFR) unit

Multi-plug flow reactor (MF-PFR) is a chemical precipitation technology focusing on the precipitation of magnesium hydroxide and calcium hydroxide with the addition of an alkaline solution. It is a two-step technology, and a chemical reaction takes place between the hydroxyl ions present in the alkaline solution and magnesium ions in the brine, promoting the precipitation of magnesium and calcium hydroxide crystals.

![mfpfr](https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/assets/150446818/65a5fccb-e2ac-4d27-9c66-69dca9157656)

In DesalSim package, the MF-PFR unit is used to model the operation of a chemical precipitation technology. Upon simulation, it will generate the influent/effluent mass flows and their concentrations, the chemical requirements, and the energy requirements.
The MF-PFR function consists of three classes: [MFPFRCALC](#use-mfpfrcalc-class), [energycons](#use-energycons-class) and [HClAddition](#use-hcladdition-class).  In this tutorial, we will focus on how to use the three classes.   

**Table 2** gives an overview of the main inputs and outputs for each process unit of Multi-plug flow reactor. 
| Process                                   | Input                                       | Output                                                |
|-------------------------------------------|---------------------------------------------|-------------------------------------------------------|
| Multi-plug flow reactor                   | Feed flow rate [m³/h]                       | Alkaline solution flow rate [L/h]                    |
|                                           | Ion concentration [g/L]                     | Flow rate of Mg(OH)₂ [kg/h]                          |
|                                           | Concentration of the alkaline solution [M] | Flow rate of Ca(OH)₂ [kg/h]                          |
|                                           | Concentration of the acid solution [M]     | Acid solution flow rate [L/h]                        |
|                                           | Products characteristics e.g. solubility...| Effluent flow rate [m³/h] and composition [g/L]      |
|                                           |                                             | Electricity requirements [kWhel]                     |

The mathematical description of MF-PFR technology is given in [Mathematical description](https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/tree/main/paper/Mathematical_description.pdf), see Section A.4. 
## 1. Getting started 
### 1.1. Import classes 
```python
import Desalsim
```
Then import the three classes:  
```python
from Desalsim.mfpfr_unit_f import MFPFRCALC
from Desalsim.mfpfr_unit_f import HClAddition
from Desalsim.mfpfr_unit_f import energycons
```
Additionally, function for calculating density (`density_calc.py`), constants (`comparison.py`) where user can add constant values like MW, prices etc, scaleup ('scaleup.py') need to be imported. 
```python
from Desalsim.density_calc import density_calc 
from Desalsim import constants 
from Desalsim import scaleup
import math
```
### 1.2. Define feed characteristics
You can initialize the feed solution by setting the flow rate, specifying the focus components and their concentration. 
```python
    # Feed concentration
components = ['Na', 'Cl', 'K', 'Mg', 'Ca', 'SO4']
Cin_mfpfr = [17.3, 38.6, 0.6, 4.3, 1.2, 9.9]  
MW = [MW_Na, MW_Ca, MW_Cl, MW_K, MW_Mg, MW_SO4]

    # Feed flowrate
Qin_mfpfr = 1000  # Flow rate in l/hr
```
Note that if you want to add more components, you need to update the components list and include the concentration of the new component in the _Ci_in_

You can calculate the density of the feed solution:
```python
mg_in = sum(Ci_in)
T=20+273 #Operating temperature (units: K)

    # Feed flow density 
d_in = density_calc(T-273, mg_in)  # kg/m3
```
### 1.3. Set Porducts Characteristics  
You need to set product characteristics such as product solubility and density. 
```python
# Product solublity of Mg(OH)2
kps_MgOH=5.61*0.000000000001
# Product solublity of Ca(OH)2
kps_CaOH=5.5*0.000001

# Mg(OH)2 density (units: kg/l)
d_mgoh_2=2.34
# Ca(OH)2 density (units: kg/l)
d_caoh_2=2.211 
```
After setting all the required inputs, then you can create the functions' objectives. 
### 1.3. Set opearional characteristics 
You need to set operational characteristics such as conversion rate for every step (_conv_) based on literature or experimental work and the concentration of the alkaline solution (_C_NaOH_). 
```python
# Concentration of NaOH solution for step 1 and step 2 in MOL/L
C_NaOH = [1, 1]

# Conversion rate for step 1 and step 2
conv = [95, 93]  
```
Additionally, you need to set assumptions related to pumping like pressure drop (_dp_) and pump efficiency (_npump_). 
```python
dp=0.5 # pressure drop (units: bar)
dp_HCl=0.3# pressure drop HCl solution (units: bar)
npump=0.8 #pump efficiency (units: -)
```

## 2. Use MFPFRCALC class   
'MFPFRCALC' is a class used to represent Mass Balance for MF-PFR Unit. In particular, it calculates the flowrates, the concentration of the streams and the requirements of alkaline solution in the two steps. 
MFPFRCALC takes as input the feed flow rate (_Qin_mfpfr_), the ion concentration in the feed (_Cin_mfpfr_), the concentration of NaOH solution for step 1 and step 2 (_C_NaOH_), the Conversion rate (_conv_).  
### 2.1. Overview 
The following attributes are available within the NFMass class:  
- `MW`: (float) Molecular weight of the solute (g/mol)
- `Ci_in`: (float) Initial concentration of the solute (g/L)
- `conv_1`: (float) Conversion rate for Magnesium precipitation in the first step 
- `conv_2`: (float) Conversion rate for Calcium precipitation in the second step    
- `C_NaOH_1`: (float) Concentration of NaOH solution for the first step (mol/L)
- `Qin`: (float) Flow rate of input (unit unspecified)
- `QMg_in`: (float) Flow rate of magnesium input (unit unspecified)
- `QNaOH_1`: (float) Volumetric flow rate of sodium hydroxide in the first step (L/h)
- `M_MgOH2_1`: (float) Outlet mass flow rate of magnesium hydroxide produced in the first step (kg/h)
- `Qtot_out_1`: (float) Outlet volumetric flow rate in the first step (L/h)
- `Mtot_out_1`: (float) Step 1 outlet mass flow rate (kg/h)
- `magma_d_1`: (float) Magma density: the quantity of solids produced per volume of slurry (kg/L)
- `ph_1`: (float) pH of solution during the first step 
- `kps_MgOH`: (float) Product solubility of Mg(OH)2 
- `Ci_out_1`: (float) The outlet ion concentration from step 1 (mol/L)
- `QNaOH_2_st`: (float) The stoichiometric volumetric flow rate of sodium hydroxide for the second step (L/hr)
- `QNaOH_2_add`: (float) The added volumetric flow rate of sodium hydroxide needed to reach a pH of 13 (L/h)
- `M_CaOH2_2`: (float) The outlet mass flow rate of calcium hydroxide produced during the second step (kg/hr)
- `M_MgOH2_2`: (float) The outlet mass flow rate of magnesium hydroxide produced during the second step (kg/hr)
- `magma_d_2`: (float) Magma density: the quantity of solids produced per volume of slurry (kg/L)
- `Qtot_out_2`: (float) Total outlet volumetric flow rate for the second step (L/h)
- `Ci_out_2`: (float) The outlet ion concentration from step 2 (mol/L)
 

TheMFPFRCALC class provides the following methods:
```python
# Calculation for the 1st step 
calc_step1()

# Calculation for the 2nd step 
calc_step2():
```
The two methods calculate the flowrates, the concentration of the streams and the requirements of alkaline solution in the two steps. It is automatically called upon initialization of the class instance.

### 2.2. Create MFPFRCALC objects
'MFPFRCALC' takes as input the feed flow rate (_Qin_mfpfr_), the ion concentration in the feed (_Cin_mfpfr_), the concentration of NaOH solution for step 1 and step 2 (_C_NaOH_), the Conversion rate (_conv_).
```python
# Create an instance of the inputpar class with the defined parameters
mfpfr_dat = MFPFRCALC(Qin_mfpfr, Cin_mfpfr, *C_NaOH, *conv)
```
Then the first method for step 1 is called. It takes as input the product solublity of Mg(OH)2 (_kps_MgOH_) and the Mg(OH)2 density (_ d_mgoh_2_). It calculates the flowrates, the concentration of the streams and the requirements of alkaline solution in the 1st step. 
First, it calculate the molar flow rate of magnesium in the reactor during the 1° stepin and then based on the concentration it clculate the volumetric flow rate of sodium hydroxide. Then it calculates magma density: the quantity of solids produced per volume of slurry and based on the outlet solution it calculates the pH ouf step 1. 
```python
# Call the calc_step1 method to calculate the necessary values
mfpfr_dat.calc_step1(kps_MgOH, d_mgoh_2)
```

Finally, the second method for step 2 is called. It takes as input the Mg(OH)2 density (_ d_mgoh_2_) and the Ca(OH)2 density (_ d_caoh_2_).  It calculates the flowrates, the concentration of the streams and the requirements of alkaline solution in the 2nd step. 
First, it calculatethe molar flow rate of calcium in the reactor during the 2nd stepin and then based on the concentration, it clculatethe stoichiometric volumetric flow rate of sodium hydroxide for the second step (_QNaOH_2_st_). Then it calculates concentration of the hydroxide ion in mol/L for a ph=13 solution and the added volumetric flow rate of sodium hydroxide needed to reach a pH = 13 (_QNaOH_2_add_). Finally, the the total outlet volumetric flow rate from 2nd step (_Qtot_out_2_), the outlet mass flow rate of calcium and magnesium hydroxide produced during the 2nd step (__M_CaOH2_2, M_MgOH2_2_), the ph of solution during the second step (_ph_2_), and the total outlet volumetric flow rate for 2nd step (_Qout_2_) and its concentration are calculated. 
```python
# Call the calc_step2 method to calculate the necessary values
mfpfr_dat.calc_step2(d_mgoh_2, d_caoh_2 )
```

### 2.3. Assigned the results to output parameters 
After the calculation for the two precipitation steps, you can assigned the results to output parameters: 
```python
ph_2=mfpfr_dat.ph_2

# Calculate the total outlet concentration and the concentration of sulfate ions
Cour_mfpfr = sum([mfpfr_dat.CNa_out_2, mfpfr_dat.CCa_out_2, mfpfr_dat.CCl_out_2, mfpfr_dat.CK_out_2, mfpfr_dat.CMg_out_2, mfpfr_dat.CSO4_out_2]) # mol/l
CSO4_out_2 = mfpfr_dat.CSO4_out_2 # mol/l

# Calculate the concentration of sodium chloride ions
CNa_out_2 = mfpfr_dat.CNa_out_2
Cnacl_out = CNa_out_2 - 2 * CSO4_out_2

# Get the molar masses of the compounds
M_MgOH2_1 = mfpfr_dat.M_MgOH2_1
M_CaOH2 = mfpfr_dat.M_CaOH2_2
M_MgOH2 = mfpfr_dat.M_MgOH2_2

# Create a list of the outlet concentrations in mol/l
Cout_all_m = [mfpfr_dat.CNa_out_2, mfpfr_dat.CCl_out_2, mfpfr_dat.CK_out_2, mfpfr_dat.CMg_out_2, mfpfr_dat.CCa_out_2, mfpfr_dat.CSO4_out_2]

#Outlet flow rate 
Qout_2 = mfpfr_dat.Qout_2

# Calculate the chemical consumption of NaOH
QNAOH = mfpfr_dat.QNaOH_1 + mfpfr_dat.QNaOH_2_add + mfpfr_dat.QNaOH_2_st # convert to kg
```
You can calculate the outlet concentration in g/l 
```python
# Calculate the outlet concentrations in g/l
Cout_mfpfr_g = [Cout_all_m[i] * MW[i] for i in range(len(Cout_all_m))] # g/l
```

### 2.4. Print results 
You can print results from mass calculations 
```python
print("Mg(OH)2 mass flow rate is "+str(round(M_MgOH2_1,2))+"kg/hr")
print("Ca(OH)2 mass flow rate is "+str(round(M_CaOH2,2))+"kg/hr")
```
Mg(OH)2 mass flow rate is 9.8kg/hr        
Ca(OH)2 mass flow rate is 2.06kg/hr     

## 3. Use HClAddition class 
'HClAddition' is a class used to represent the calculation of HCl addition in MF-PFR Unit. It calculates amount of HCl for pH neutralization and the final outlet concentration from MF-PFR unit after pH neutralization. For this calculation, you need to use the outlet flow rate from step 2 (_Qout_2_), the molar concentration of the the outlet flow rate from step 2 (_Cout_all_m), Cl molacular weight (_MW_Cl_), the pH of the solution after the 2nd precipitation step (_ph_2_), and the concentration of the acid solution used for the pH neutralization (_HCl_conc_). The class _returns the flow rate of the required acid solution (_QHCl_) and the new ion concentration in g/l after the pH neutralization (_Cout_mfpfr_g_).   
### 3.1. Oveview
The following attributes are available within theHClAddition class:  
-  `C1, C2 C3, C4, C5, C6 `: (float) Concentration of ions in the solution (mol/L).
-  `z1, z2,z3, z4, z5, z6`: (int) Charge of ions in the solution.

The OsmoticPressure class provides the following method:
```python
calculate_HCl_addition()
```

### 3.2. Create HClAddition objectives and calculate the amount of HCl added and the new outlet concentration 

```python
    # Calculate HClAdditione for  
unit = HClAddition(Qout_2, Cout_all_m, MW_Cl, ph_2, HCl_conc)
QHCl, Cout_mfpfr_g = unit.calculate_HCl_addition(Cout_mfpfr_g)
```
The above line assigns also the results to output parameters.

### 3.4. Print results 
You can print results from HCl addition and the new effluent flow rate.  
```python
# Print the volume of HCl added and the outlet concentration of chloride ions
print(f"HCl flow rate is {round(QHCl,2)} l/hr")

print("Total effluent flow rate is "+str(round(Mout_f,2))+"kg/hr")
print("Total effluent flow rate is "+str(round(Qout_f,2))+"kg/hr")
```
HCl flow rate is 152.66 l/hr  
NaOH flow rate is 531.95 l/hr  

Total effluent flow rate is 1737.65kg/hr   
Total effluent flow rate is 1679.27kg/hr   

## 4. Use energycons class
'energycons' is a class used to represent the calculation of energy consumption and the specific energy consumption for MF-PFR Unit. The takes as input the total volumetric flow rate (_Qtot_), the volumetric flow rate of sodium hydroxide (_QNaOH_), the volumetric flow rate of input (_Qin_), the volumetric flow rate of sodium hydroxide in the first step (_QNaOH_1_), the Added volumetric flow rate of sodium hydroxide needed to reach a pH of 13 (_QNaOH_2_add_), the stoichiometric volumetric flow rate of sodium hydroxide for the second step (_QNaOH_2_st_), the expected pressure drop (_dp_) and the pump efficiency (_npump_). The class 'energycons' _returns the energy consumption for pumping in the two steps (_Epump_1, Epump_2_)_.

### 4.1. Oveview 
The following attributes are available within the 'energycons' class:  
- `Qtot`: Total volumetric flow rate (L/h)
- `QNaOH`: Volumetric flow rate of sodium hydroxide (L/h)
- `Qin`: Volumetric flow rate of input (L/h)
- `QNaOH_1`: Volumetric flow rate of sodium hydroxide in the first step (L/h)
- `QNaOH_2_add`: Added volumetric flow rate of sodium hydroxide needed to reach a pH of 13 (L/h)
- `QNaOH_2_st`: Stoichiometric volumetric flow rate of sodium hydroxide for the second step (L/h)
- `dp`: Pressure drop (bar)
- `npump`: Pump efficiency (-)

The  energycons class provides the following method:
```python
calculate_energy_consumption()
```
This method calculates the  the energy consumption for pumping in the two steps (_Epump_1, Epump_2_).

### 4.2. Create nf_energy objectives and assign the results to output parameters
The following objective is created for energy consumption. Assumptions for pressure drop and pump efficiency need to be made. 

```python
   # Create an instance of the inputpar class with the defined parameters
Epump_1, Epump_2=energycons.energycalc(mfpfr_dat.Qout_2, QNAOH, Qin_mfpfr, mfpfr_dat.QNaOH_1, mfpfr_dat.QNaOH_2_add, mfpfr_dat.QNaOH_2_st, dp, npump)
```
You can calculate total energy consumption for pumping and the Specific energy consumptions. 
```python
    #Electricity consumption for pumping , KWh
E_el_mfpf=(Epump_1+Epump_2+(QHCl*dp_HCl)*1e5/3600/(1000*npump))/1000

    #Specific energy consumption per kg of Mg(OH)2, KWh/kg of Mg(OH)2
SEC_el_prod=(E_el_mfpf)/(M_MgOH2)

    #Specific energy consumption per feed, KWh/m3 of feed
SEC_el_feed=(E_el_mfpf)/(Qin_mfpfr/1000)
```
You can add calculation for the filtration unit, if there are available data and equations. 
### 4.4. Print results 
You can print results from energy calculations. The specific energy consumption is also calculated so you can validate easier the results. 
```python
print("Total electricity energy consumption is "+str(round(E_el_mfpf,2))+ " KW")
print("Specific energy consumption per product is "+str(round(SEC_el_prod,2))+" KWh/kg product ")
print("Specific energy consumption per brine intake is "+str(round(SEC_el_feed,2))+" KWh/m3 of feed ")
```
Total electricity energy consumption is 0.05 KW  
Specific energy consumption per product is 2.88 KWh/kg product  
Specific energy consumption per brine intake is 1.49 KWh/m3 of feed  



## References 


