# Tutotrial: Multi-plug flow reactor (MF-PFR) unit

Nanofiltration (NF) is a membrane liquid-separation technology sharing many characteristics with reverse osmosis (RO). Unlike RO, which has high rejection of virtually all dissolved solutes, 
NF provides high rejection of multivalent ions, such as calcium, and low rejection of monovalent ions, such as chloride [@dupontwebsite]. 



In DesalSim package, the MF-PFR unit is used to model the operation of a chemical precipitation technology. Upon simulation, it will generate the influent/effluent mass flows and their concentrations, the chemical requirements, and the energy requirements.
The MF-PFR function consists of three classes: [NFMass](#use-nfmass-class), [energycons](#use-energycons-class) and [HClAddition](#use-hcladdition-class).  
In this tutorial, we will focus on how to use the three classes. 

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
### 1.3. Set Membrane Characteristics  
You can set membrane characteristics, ion rejection rates and Water recovery. 
```python
    # Ions rejection rates based on membrane characteristics (units: -)
rjr_values = [0.16, 0.29, 0.21, 0.98, 0.95, 0.98]
    # Water recovery based on membrane characteristics (units: -)
Wrec = 0.7 
```
After setting all the required inputs, then you can create the functions' objectives. 

## 2. Use MFPFRCALC class   
MFPFRCALC is a class used to represent Mass Balance for MF-PFR Unit. In particular, it calculates the permeate and concentrate flow rates, and their ion concentrations. 
MFPFRCALC takes as input the names of components (_comp_), the ion concentration in the feed (_C_in_), the rejection rates of the ions (_rjr_values_), the % of water recovery (_Wrec_) and the feed flow rate (_Qf_).  
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
 


TheMFPFRCALC class provides the following method:
```python
calculate_perm()
```
This method calculates the permeate and concentrate flow rates, as well as their corresponding ion concentrations based on the provided attributes. It is automatically called upon initialization of the class instance.

### 2.2. Create MFPFRCALC objects
MFPFRCALC takes as input the names of components (_comp_), the ion concentration in the feed (_C_in_), the rejection rates of the ions (_rjr_values_), the % of water recovery (_Wrec_) and the feed flow rate (_Qf_).
```python

```

### 2.3. Assigned the results to output parameters 
After the calculation of the permeate and concentrate flow rates, as well as their corresponding ion concentrations based on the provided attributes, you can assigned the results to output parameters: 
```python

```

### 2.4. Print results 
You can print results from mass calculations 
```python

```



## 3. Use HClAddition class 
HClAddition is a class used to represent the calculation of osmotic pressure for Nanofiltration Unit. For the calculation of the energy consumption, first the Osmotic pressure for the three streams (feed, concentrate, permeate) need to be calculated. For this calculation, you need to use the ion concentration of the stream (_Ci_in_, _Cperm_, _Cconc_) the ionelectric charge (_z_values_), and the stream temperature (_T_). The class _returns the Osmotic pressure_ of the solution.   
### 3.1. Oveview
The following attributes are available within theHClAddition class:  
-  `C1, C2 C3, C4, C5, C6 `: (float) Concentration of ions in the solution (mol/L).
-  `z1, z2,z3, z4, z5, z6`: (int) Charge of ions in the solution.

The OsmoticPressure class provides the following method:
```python
osmotic_pressure_calculation()
```
This method calculates the osmotic pressure of a solution.

### 3.2. Create HClAddition objectives and calculate the amount of HCl added and the new outlet concentration 

```python
    # Calculate HClAdditione for  
unit = HClAddition(Qout_2, Cout_all_m, MW_Cl, ph_2, HCl_conc)
QHCl, Cout_mfpfr_g = unit.calculate_HCl_addition(Cout_mfpfr_g)
```
The above line assigns also the results to output parameters.

### 3.4. Print results 
You can print results from mass calculations 
```python
# Print the volume of HCl added and the outlet concentration of chloride ions
print(f"HCl flow rate is {round(QHCl,2)} l/hr")
```
## 4. Use energycons class
energycons is a class used to represent the calculation of energy consumption and the specific energy consumption for MF-PFR Unit. For this calculation, takes the Osmotic pressure for the three streams (feed, concentrate, permeate). In addition, the NfEnergy takes as input the expected pressure drop in each stream (_dp_, d_p_, d_in_) and the pump efficiency (_n). The class _returns the Applied pressure, power for applied pressure, the total energy consumption_ and the _specific energy consumption per m<sup>3</sup> permeate_ and _m<sup>3</sup> feed_.
### 4.1. Oveview 
The following attributes are available within the NfEnergy class:  
- `P_osmo_c`: (float) Osmotic pressure of concentrate stream (bar).
- `P_osmo_f`: (float) Osmotic pressure of feed stream (bar).
- `P_osmo_p`: (float) Osmotic pressure of permeate stream (bar).
- `dp`: (float) Pressure drop (bar).
- `d_p`: (float) Permeate stream density (kg/m³).
- `Qperm`: (float) Permeate flow rate (kg/h).
- `Qf`: (float) Concentrate flow rate (kg/h).
- `d_in`: (float) Feed stream density (kg/m³).
- `n`: (float) Pump efficiency (-).

The  energycons class provides the following method:
```python
calculate_energy_consumption()
```
This method calculates the the Applied pressure, power for applied pressure, the total energy consumption and the specific energy consumption per m<sup>3</sup> permeate and _m<sup>3</sup> feed

### 4.2. Create nf_energy objectives and calculate Energy Consumption
The following objective is created for energy consumption. Assumptions for pressure drop and pump efficiency need to be made. 

```python
nf_energy=NfEnergy(P_osmo_c, P_osmo_f, P_osmo_p, dp=2, d_p, Qperm, Qf_nf, d_in,n=0.8) # dp: pressure drop (units: bar) and n: pump efficiency (units: -)
result=nf_energy.calculate_energy_consumption()
```
### 4.3. Assigned the results to output parameters 
```python
E_el_nf = nf_energy.E_el_nf
```
### 4.4. Print results 
You can print results from energy calculations. The specific energy consumption is also calculated so you can validate easier the results. 
```python
for key, value in result.items():
        print(f"{key}: {value}")
```
Applied pressure (Bar): 24.45  
Power for pump (KW): 60.01  
E_el_nf (KW): 75.02  
Specific Energy Consumption (KWh/m3 of permeate): 0.85



## References 


