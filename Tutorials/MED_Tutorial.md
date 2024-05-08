# Tutotrial: Multi-effect distillation (MED) unit

MED is a thermal based process that is used to desalinate water. In this work, MED aims to recover high quality water and concentrate further the brine stream.
![med](https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/assets/150446818/54c49681-d70f-4cf4-ab5f-4a5d7d8791b7)


In **DesalSim** package, the MED unit is used to model the operation of a Multi-effect distillation technology. Upon simulation, it will generate the influent/effluent mass flows and their concentrations, the cooling water, and the energy requirements.
The MED function consists of one classes: [MEDCalculator](#use-medalculator-class) that constsis.  
In this tutorial, we will focus on how to use the the class and their methods. 

**Table 1** gives an overview of the main inputs and outputs for each process unit of Multi-effect distillation. 
| Process                                   | Input                                       | Output                                                |
|-------------------------------------------|---------------------------------------------|-------------------------------------------------------|
| Multi-effect distillation             | Feed flow rate [m3/h]                        | Flow rate of water [m3/h]                                  |
|                                        | Ion concentration [g/L]                      | Effluent flow rate and composition [g/L]                   |
|                                        | Feed temperature [°C]                        | Electrical [kWhel] and thermal [kWhth] requirements        |
|                                        | Steam temperature [°C]                       | Cooling water flow rate [m3/h]                             |

The mathematical description of Nanofiltration technology is given in [Mathematical description](https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/tree/main/paper/Mathematical_description.pdf), see Section A.2. 

## 1. Getting started 
### 1.1. Import class
```python
import Desalsim
```
Then import the class:  
```python
from Desalsim.med_unit_f import MEDCalculator
```
Additionally, function for calculating density (`density_calc.py`) need to be imported. 
```python
from Desalsim.density_calc import density_calc
```
### 1.2. Define feed characteristics
You can initialize the feed solution by setting the flow rate, specifying the focus components and their concentration. 
```python
  #Feed concentration
components = ['Na', 'Cl', 'K', 'Mg', 'Ca', 'SO4']
Cin_med = [10.36, 15.39, 0.36, 0.028, 0.02, 0.07]

#Feed flow rate 
Qf_med =1000 #l/hr

#input conditions
T=20
```
Note that if you want to add more components, you need to update the components list and include the concentration of the new component in the _Ciin_med_. 

You can calculate the density of the feed solution and the mass flow rate:
```python
d=density_calc(T, sum(Cin_med))

# Mass flow rate (units: kg/hr)
Mf_med=Qf_med*d/1000 
```
### 1.3. Set operating assumptions  
You need to set operating assumptions related to temperatures such as the temperature in the last effect, the intake/outake cooling water temperature etc.  
```python
# Assumptions:
T_in=40 #(oC)
T_N=45 #Temperature in the last effect (oC)
T_cw_in=25 #intake cooling water temperature (oC)
T_cw_out=35 #out cooling water temperature (oC)
T_s=70 #steam temperature oC
DT_loss=1 #temperature difference (oC)
T3=69
```
#### Calculate latent heat of motive steam
The latent heat of motive steam is calculated based on the set steam temperature and steam tables. 
```python
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
```

Additionally, you need to define the aimed brine vonventration leaving effect n (_Cb_out_) and the brine circulation flow rate (_Xr_), and the number of effects (_N_). 
```python
Cb_out=200 #Brine Concentration leaving effect n (unit: g/l)
Xr=5.5 # brine circulation flow rate (units: -)
N=2 #Number of effects (-)
```
Finally, you need to set assumptions related to pumping like pressure drop (_dp_) and pump efficiency (_npump_. 
```python
dp=0.1  # pressure drop (units: bar)
dp_slurry=1 # pressure drop (units: bar)
npump=0.8 #pump efficiency (units: -)
```
### 1.4. Set constants 
You need to set constant parameters like the specific heat capacity of water
```python
Cp_w=4182 # specific heat capacity of water (j/kgC)
cp_sol=4184 # specific heat capacity of solution (j/kgC)
```
After setting all the required inputs, then you can create the functions' objectives. 

## 2. Use MEDCalculator class   
MEDCalculator is a class used to represent Mass Balance for Nanofiltration Unit. In particular, it calculates the permeate and concentrate flow rates, and their ion concentrations. 
MEDCalculator takes as input the names of components (_comp_), the ion concentration in the feed (_C_in_), the rejection rates of the ions (_rjr_values_), the % of water recovery (_Wrec_) and the feed flow rate (_Qf_).  
### 2.1. Overview 
The following attributes are available within the NFMass class:  
- `Qf`: (float) Flow rate (m^3/s).
- `CNa_in, CCl_in, CK_in, CMg_in, CCa_in, CSO4_in` : float
        Initial concentrations of various ions (g/l).

The MEDCalculator class provides the following method:
```python
calculate_perm()
```
This method calculates the permeate and concentrate flow rates, as well as their corresponding ion concentrations based on the provided attributes. It is automatically called upon initialization of the class instance.

### 2.2. Create MEDCalculator objects
MEDCalculator takes as input the feed volumetric flow rate (_Qf_med_) and mass flow rate (_Mf_med_), the concentration in the feed for the components (_CNa_in, CCl_in, CK_in, CMg_in, CCa_in, CSO4_in_) and the 
 Cin_med[0], Cin_med[1], Cin_med[2], Cin_med[3], Cin_med[4], Cin_med[5], and the inlet temperature (_T_in_). 
 
```python
# Create an instance of the MEDCalculator class
med_dat = MEDCalculator(Qf_med, Mf_med, Cin_med[0], Cin_med[1], Cin_med[2], Cin_med[3], Cin_med[4], Cin_med[5], T_in)

```

### 2.3. Use 'salinity_calc' method 
```python

```
```python

```
```python

```
### 2.3.1. Assigned the results to output parameters 
```python

```
### 2.4. Use 'mass_balance_med' method 

```python

```

### 2.4.1. Assigned the results to output parameters 
```python

```

### 2.5. Use 'temperature_calc' method 
```python

```

### 2.5.1. Assigned the results to output parameters 
```python

```

### 2.6. Use 'performance_parameters' method 
```python

```

### 2.6.1. Assigned the results to output parameters 
```python

```

### 2.7. Use 'output_concentration' method 
```python

```
### 2.7.1. Assigned the results to output parameters 

```python

```
After the calculation of , you can assigned the results to output parameters: 
```python

```

### 2.8. Print results 
You can print results from mass calculations 
```python

```


