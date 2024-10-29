# Tutotrial: Nanofiltration (NF) unit

Nanofiltration (NF) is a membrane liquid-separation technology sharing many characteristics with reverse osmosis (RO). Unlike RO, which has high rejection of virtually all dissolved solutes, 
NF provides high rejection of multivalent ions, such as calcium, and low rejection of monovalent ions, such as chloride [@dupontwebsite]. 
![nf](https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/assets/150446818/1d41d6eb-90a7-4b68-ab1b-0ae31f83eb78)


In desalsim package, the Nanofiltration unit is used to model the operation of a nanofiltration (or a membrane-based) technology. Upon simulation, it will generate the influent/effluent mass flows and their concentrations, the applied pressure, and the energy requirements.
The nanofiltration function consists of three classes: [NFMass](#use-nfmass-class), [OsmoticPressure](#use-osmoticpressure-class) , and [NfEnergy](#use-nfenergy-class).  
In this tutorial, we will focus on how to use the three classes. 

**Table 1** gives an overview of the main inputs and outputs for each process unit of Nanofiltration. 
| Process                                   | Input                                       | Output                                                |
|-------------------------------------------|---------------------------------------------|-------------------------------------------------------|
| Nanofiltration                            | Feed flow rate [m³/h]                       | Permeate flow rate and composition [g/L]              |
|                                           | Ion concentration [g/L]                     | Concentrate flow rate and composition [g/L]           |
|                                           | Water recovery [%]                         | Electrical requirements [kWhel]                       |
|                                           |  Ion rejection [-]                         | Osmotic pressure [bar]                                    |


The mathematical description of Nanofiltration technology is given in [Mathematical description](https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/tree/main/paper/Mathematical_description.pdf), see Section A.1. 
## 1. Getting started 
### 1.1. Import classes 
```python
import desalsim
```
Then import the three classes:  
```python
from desalsim.nanofiltration_unit_f import OsmoticPressure
from desalsim.nanofiltration_unit_f import NFMass
from desalsim.nanofiltration_unit_f import NfEnergy
```
Additionally, function for calculating density (`density_calc.py`) or constants (`comparison.py`) where user can add constant values like MW, prices etc, need to be imported. 
```python
from desalsim.density_calc import density_calc
from desalsim import constants
```
### 1.2. Define feed characteristics
You can initialize the feed solution by setting the flow rate, specifying the focus components and their concentration. 
```python
    # Feed concentration
components = ['Na', 'Cl', 'K', 'Mg', 'Ca', 'SO4']
Ci_in = [12.33, 21.67, 0.45, 1.39, 0.45, 3.28]
z_values = [1, -1, 1, 2, 2, -2]
```
> [!NOTE]
> Note that if you want to add more components, you need to update the components list and include the concentration of the new component in the _Ci_in_

You can calculate the density of the feed solution:
```python
mg_in = sum(Ci_in)
T=20+273 #Operating temperature (units: K)

    # Feed flow density 
d_in = density_calc(T-273, mg_in)  # kg/m3
```
```python
    # Feed flowrate
Qf_nf = 3000 / 24 * d_in #m3/d
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

## 2. Use NFMass class   
NFMass is a class used to represent Mass Balance for Nanofiltration Unit. In particular, it calculates the permeate and concentrate flow rates, and their ion concentrations. 
NFMass takes as input the names of components (_comp_), the ion concentration in the feed (_C_in_), the rejection rates of the ions (_rjr_values_), the % of water recovery (_Wrec_) and the feed flow rate (_Qf_).  
### 2.1. Overview 
The following attributes are available within the NFMass class:  
- `comp`: (str) Component name  
- `Cfeedi`: (float) Ion concentration in the feed (g/L)  
- `rjr`: (float) Rejection of the ion by the membrane 
- `Wrec`: (float) Water recovery in the first pass  
- `Qf`: (float) Feed flow rate (kg/h)  
- `Qperm`: (float) Permeate flow rate (kg/h)  
- `Qconc`: (float) Concentrate flow rate (kg/h)  
- `Cpermi`: (float) Ion concentration in the permeate (g/L)  
- `Cconci`: (float) Ion concentration in the concentrate (g/L) 

The NFMass class provides the following method:
```python
calculate_perm()
```
This method calculates the permeate and concentrate flow rates, as well as their corresponding ion concentrations based on the provided attributes. It is automatically called upon initialization of the class instance.

### 2.2. Create NFMass objects
NFMass takes as input the names of components (_comp_), the ion concentration in the feed (_C_in_), the rejection rates of the ions (_rjr_values_), the % of water recovery (_Wrec_) and the feed flow rate (_Qf_).
```python
    # Function to create NFMass objects for different components
def create_nfmass_objects(components, C_in, rjr_values, Wrec, Qf):
    return [NFMass(comp, Ci, rjr, Wrec, Qf) for comp, Ci, rjr in zip(components, C_in, rjr_values)]

    # Create NFMass objects for different components
nfmass_objects = create_nfmass_objects(components, Ci_in, rjr_values, Wrec, Qf_nf)
```

### 2.3. Assigned the results to output parameters 
After the calculation of the permeate and concentrate flow rates, as well as their corresponding ion concentrations based on the provided attributes, you can assigned the results to output parameters: 
```python
    # Components concentrattion in concentrate stream 
Cconc = [nf_mass.Cconci for nf_mass in nfmass_objects]
    # Components concentrattion in permeate stream 
Cperm = [nf_mass.Cpermi for nf_mass in nfmass_objects]
    # Permeate stream mass flow rate
Qperm = nfmass_objects[0].Qperm  # kg/hr
    # Concentrate stream mass flow rate
Qconc = nfmass_objects[0].Qconc  # kg/hr
```

### 2.4. Print results 
You can print results from mass calculations 
```python
print("Permeate stream flow rate is "+str(round(Qperm,2))+"kg/hr")
print("Permeate stream total concentration is "+str(round(sum(Cperm),2))+"g/l")
print("Concentrate stream flow rate is "+str(round(Qconc,2))+"kg/hr")
print("Concentrate stream total concentration is "+str(round(sum(Cconc),2))+"g/l")
```

Permeate stream flow rate is 89974.58kg/hr
Permeate stream total concentration is 26.21g/l
Concentrate stream flow rate is 38560.54kg/hr
Concentrate stream total concentration is 70.73g/l

## 3. Use OsmoticPressure class 
OsmoticPressure is a class used to represent the calculation of osmotic pressure for Nanofiltration Unit. For the calculation of the energy consumption, first the Osmotic pressure for the three streams (feed, concentrate, permeate) need to be calculated. For this calculation, you need to use the ion concentration of the stream (_Ci_in_, _Cperm_, _Cconc_) and the Ions molar mass (_MW_values_). The class _returns the Osmotic pressure_ of the solution.   
### 3.1. Oveview
The following attributes are available within the OsmoticPressure class:  
-  `Ci_in `: (float) Concentration of ions in the solution (mol/L).
-  `MW_values: (float) ons molar mass in g/mol.

The OsmoticPressure class provides the following method:
```python
calculate_osmotic_pressure()
```
This method calculates the osmotic pressure of a solution based on the Gibbs equation: 

$\pi = -\left( \frac{RT}{V} \right) \ln(a_w)$


### 3.2. Create OsmoticPressure objectives and calculate Osmotic Pressure

```python
    # Calculate Osmotic Pressure for the three streams 
P_osmo_f = OsmoticPressure(Ci_in, MW_values).calculate_osmotic_pressure()
P_osmo_p = OsmoticPressure(Cperm, MW_values).calculate_osmotic_pressure()
P_osmo_c = OsmoticPressure(Cconc, MW_values).calculate_osmotic_pressure()
```
## 4. Use NfEnergy class
NfEnergy is a class used to represent the calculation of energy consumption and the specific energy consumption for Nanofiltration Unit. For this calculation, the Osmotic pressure for the three streams (feed, concentrate, permeate) is used. In addition, the NfEnergy takes as input the expected pressure drop in each stream (_dp_, d_p_, d_in_) and the pump efficiency (_n_). The class _returns the Applied pressure, power for applied pressure, the total energy consumption_ and the _specific energy consumption per m<sup>3</sup> permeate_ and _m<sup>3</sup> feed_.
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

The  NfEnergy class provides the following method:
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
Applied pressure (Bar): 21.16  
Power for pump (KW): 51.93 
E_el_nf (KW): 64.92  
Specific Energy Consumption (KWh/m3 of permeate): 0.73

## References 

