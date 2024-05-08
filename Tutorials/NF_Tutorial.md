# Tutotrial: Nanofiltration (NF) unit

Nanofiltration (NF) is a membrane liquid-separation technology sharing many characteristics with reverse osmosis (RO). Unlike RO, which has high rejection of virtually all dissolved solutes, 
NF provides high rejection of multivalent ions, such as calcium, and low rejection of monovalent ions, such as chloride.



The Nanofiltration unit is used to model the operation of a nanofiltration (or a membrane-based) technology. Upon simulation, it will generate the influent/effluent mass flows and their concentrations, the applied pressure, and the energy requirements.
The nanofiltration function consists of three classes: _NFMass_, _OsmoticPressure_, and _NfEnergy_.  
In this tutorial, we will focus on how to use the three classes. 
NFMass is a class used to represent Mass Balance for Nanofiltration Unit. In particular, it calculates the permeate and concentrate flow rates, and their ion concentrations. 
## 2. Installation 
The easiest way is through pip, in command-line interface:   
```
pip install DesalSim==1.0.1
```

If you want to install the latest GitHub verstion:
1. Download the repository to your local machine:
```
https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-
```
2. Install the required dependencies:
 ```
pip install -r requirements.txt
 ```

```python
import Desalsim
```
Then:  
```python
from Desalsim.nanofiltration_unit_f import OsmoticPressure
from Desalsim.nanofiltration_unit_f import NFMass
from Desalsim.nanofiltration_unit_f import NfEnergy
```
Similarly for the other process units. Additionally, function for calculating density (`density_calc.py`) or constants (`comparison.py`) where user can add constant values like MW, prices etc, need to be imported. 
```python
from Desalsim.density_calc import density_calc
import Desalsim.constants
import Desalsim.scaleup
```

### 3.1.1. Define feed characteristics
You can initialize the feed solution by setting the flow rate, specifying the focus components and their concentration. 
```python
    # Feed concentration
components = ['Na', 'Cl', 'K', 'Mg', 'Ca', 'SO4']
Ci_in = [12.33, 21.67, 0.45, 1.39, 0.45, 3.28]
z_values = [1, -1, 1, 2, 2, -2]

    # Feed flowrate
Qsw = 3000 / 24 * d_in #m3/d
```
Note that if you want to add more components, you need to update the components list and include the concentration of the new component in the _Ci_in_

You can calculate the density of the feed solution:
```python
mg_in = sum(Ci_in)
T=20+273 #Operating temperature (units: K)

    # Feed flow density 
d_in = density_calc(T-273, mg_in)  # kg/m3
```

## 3.2. Use process unit model

### 3.2.1. Nanofiltration 

To run simulation model for Nanofiltration unit, you need to implement the following steps. 

**Table 1** gives an overview of the main inputs and outputs for each process unit of Nanofiltration. 
| Process                                   | Input                                       | Output                                                |
|-------------------------------------------|---------------------------------------------|-------------------------------------------------------|
| Nanofiltration                            | Feed flow rate [m³/h]                       | Permeate flow rate and composition [g/L]              |
|                                           | Ion concentration [g/L]                     | Concentrate flow rate and composition [g/L]           |
|                                           | Water recovery [%]                         | Electrical requirements [kWhel]                       |
|                                           |  Ion rejection [-]                         | Osmotic pressure [bar]                                    |


##### Setting Membrane Characteristics
You can set membrane characteristics, ion rejection rates and Water recovery. 
```python
    # Ions rejection rates based on membrane characteristics (units: -)
rjr_values = [0.16, 0.29, 0.21, 0.98, 0.95, 0.98]
    # Water recovery based on membrane characteristics (units: -)
Wrec = 0.7 
```
##### Create NFMass Objects
After setting all the required inputs, then you can create the functions' objectives. 
NFMass is a class used to represent Mass Balance for Nanofiltration Unit. In particular, it calculates the permeate and concentrate flow rates, and their ion concentrations. 
NFMass takes as input the names of components (_comp_), the ion concentration in the feed (_C_in_), the rejection rates of the ions (_rjr_values_), the % of water recovery (_Wrec_) and the feed flow rate (_Qf_).  

```python
    # Function to create NFMass objects for different components
def create_nfmass_objects(components, C_in, rjr_values, Wrec, Qf):
    return [NFMass(comp, Ci, rjr, Wrec, Qf) for comp, Ci, rjr in zip(components, C_in, rjr_values)]

    # Create NFMass objects for different components
nfmass_objects = create_nfmass_objects(components, Ci_in, rjr_values, Wrec, Qf_nf)
```
Assigned the results to output parameters 

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

##### Calculate Osmotic Pressure
For the calculation of the energy consumption, first the Osmotic pressure for the three streams (feed, concentrate, permeate) need to be calculated. For this calculation, you need to use the ion concentration of the stream (_Ci_in_, _Cperm_, _Cconc_) the ionelectric charge (_z_values_), and the stream temperature (_T_).
```python
    # Calculate Osmotic Pressure
P_osmo_f = OsmoticPressure(Ci_in, z_values, T).osmotic_pressure_calculation()
P_osmo_p = OsmoticPressure(Cperm, z_values, T).osmotic_pressure_calculation()
P_osmo_c = OsmoticPressure(Cconc, z_values, T).osmotic_pressure_calculation()
```

##### Calculate Energy Consumption
The following objective is created for energy consumption. Assumptions for pressure drop and pump efficiency need to be made. 
```python
nf_energy=NfEnergy(P_osmo_c, P_osmo_f, P_osmo_p, dp=2, d_p, Qperm, Qf_nf, d_in,n=0.8) # dp: pressure drop (units: bar) and n: pump efficiency (units: -)
result=nf_energy.calculate_energy_consumption()
E_el_nf = nf_energy.E_el_nf
```
You can print results from energy calculations. The specific energy consumption is also calculated so you can validate easier the results. 
```python
for key, value in result.items():
        print(f"{key}: {value}")
```

Applied pressure (Bar): 24.45  
Power for pump (KW): 60.01  
E_el_nf (KW): 75.02  
Specific Energy Consumption (KWh/m3 of permeate): 0.85
