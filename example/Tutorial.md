# Tutorial 
## 1. Introduction 
This tutorial provides a comprehensive overview and guide to utilizing a simulation package tailored for analyzing desalination and brine treatment technologies. Here's what you'll find:  
**2. Usage:** Instructions on how to use the simulation models, including input parameters and result interpretation.  
**3. Technical Process Models:** Detailed descriptions of each technology model, including input-output relationships and simulation steps.  
**4. Economic Models:** Explanation of economic models for evaluating operating and investment costs.  
**5. Treatment Chains Comparison:** Guidance on comparing different treatment chains using provided tools.  

## 2. Installation 
The easiest way is through pip, in command-line interface:   
```
pip install DesalSim==0.5
```

If you want to install the latest GitHub verstion:
1. Clone the repository to your local machine:
```
git clone https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-
```
2. Install the required dependencies:
 ```
pip install -r requirements.txt
 ```
## 3. Usage 
Each simulation model serves as a standalone tool for analyzing the performance of a specific desalination or brine treatment technology. Before running the simulation, ensure that you have provided the necessary input parameters, such as feed flow rates, salinity levels, membrane properties, heat sources, and operating conditions.

The simulation results, including salt concentration profiles, ion fluxes, energy consumption, chemical consumption, and operational costs, will be generated based on the specified inputs and displayed in the console output or saved to output files for further analysis.

However, simulation models of more than one technology can be combined to simulate and evaluate the performance of a treatment chain (desalination and brine treatment system). In this case, the output flow rates and stream concentrate are the input data for the next technology. 

Additionally, two example files are provided to demonstrate the usage of the simulation suite (see [Example 1](example_1.py) and [Example 2](example_2.py)). These examples simulate and evaluate two different treatment chains, showcasing the integration of multiple technologies. The economic evaluation of the treatment chain is given in [Example 1](example_1.py) and in the [Economic Tutorial](Economic_Tutorial.md). 
Furthermore, a [comparison file](comparison.py) is included, where the results of the two examples are compared in terms of various parameters. Users can extend this comparison by adding more indicators as needed.

**<u>Followed steps:</u>**<br>
Step 1: Import required fucntions for process units in the treatment chain.<br>
Step 2: Set input data like feed flow rate, ion concentration, relevant ions for the feed solution.<br>
Step 3: Set input parameters for each process unit as shown in **Table 1** and for economic model as shown in **Table 2** and **Table 3**.<br>
Step 4: Call function of each process unit, creat objects for each calculation.<br>
Step 5: Results interpretation. <br>
<br>

### 4. Technical process models 
For more detailed steps and instructions see [Tutorial for Example 1](Example_1_Tutorial.md).  
The mathematical description of each technology is given in  [Mathematical description](../paper/Mathematical_description.pdf).  

**Table 1** gives an overview of the main inputs and outputs for each process unit of example 1. 
| Process                                | Input                                        | Output                                                     |
|----------------------------------------|----------------------------------------------|------------------------------------------------------------|
| Nanofiltration                         | Feed flow rate [m3/h]                        | Permeate flow rate and composition [g/L]                   |
|                                        | Ion concentration [g/L]                      | Concentrate flow rate and composition [g/L]                |
|                                        | Osmotic pressure [bar]                       | Electrical requirements [kWhel]                            |
|                                        | Water recovery [%]                          |                                        |
|                                        | Ion rejection [-]                          |                                       |
| Multi-effect distillation             | Feed flow rate [m3/h]                        | Flow rate of water [m3/h]                                  |
|                                        | Ion concentration [g/L]                      | Effluent flow rate and composition [g/L]                   |
|                                        | Feed temperature [°C]                        | Electrical [kWhel] and thermal [kWhth] requirements        |
|                                        | Steam temperature [°C]                       | Cooling water flow rate [m3/h]                             |
| Thermal crystallizer                  | Feed flow rate [m3/h]                        | Flow rate of water [kg/h]                                  |
|                                        | Ion concentration [g/L]                      | Flow rate of NaCl [kg/h]                                   |
|                                        | Feed temperature [°C]                        | Cooling water flow rate [m3/h]                             |
|                                        | Steam temperature [°C]                       | Electrical [kWhel] and thermal [kWhth] requirements        |
| Multi-plug flow reactor               | Feed flow rate [m3/h]                        | Alkaline solution flow rate [L/h]                          |
|                                        | Ion concentration [g/L]                      | Flow rate of Mg(OH)₂ [kg/h]                               |
|                                        | Concentration of the alkaline solution (NaOH) [M] | Flow rate of Ca(OH)₂ [kg/h]                          |
|                                        | Concentration of the acid solution (HCl) [M] | Acid solution flow rate [L/h]                             |
|                                        |                                            | Effluent flow rate [m3/h] and composition [g/L]            |
|                                        |                                            | Electricity requirements [kWhel]                          |
| Eutectic freeze crystallizer          | Feed flow rate [m3/h]                        | Flow rate of Na2SO4 [kg/h]                                |
|                                        | Ion concentration [g/L]                      | Flow rate of ice [kg/h]                                    |
|                                        | Feed temperature [°C]                        | Effluent flow rate [m3/h] and composition [g/L]            |
|                                        |                                            | Electricity requirements [kWhel]                          |
| Electrodialysis with bipolar membranes| Feed flow rate [m3/h]                        | Flow rate of acid [m3/h] and composition [g/L]             |
|                                        | Ion concentration [g/L]                      | Flow rate of base [m3/h] and composition [g/L]             |
|                                        | Electric density                            | Flow rate of salt [m3/h] and composition [g/L]             |
|                                        |                                            | Electricity requirements [kWhel]                          |
| Electrodialysis                        | Feed flow rate [m3/h]                        | Flow rate of diluted stream [m3/h] and composition [g/L]   |
|                                        | Ion concentration [g/L]                      | Flow rate of concentrate stream [m3/h] and composition [g/L]|
|                                        | Electric density                            | Electricity requirements [kWhel]                          |


### 5. Economic models 

For more detailed steps and instructions see [Economic Tutorial](Economic_Tutorial.md).  
The mathematical description of economic model is given also in [Mathematical description](paper/Mathematical_description.pdf).  

**Table 2** gives an overview of the main inputs and outputs of economic model (`economic_f.py`). 

|  Input                                     | Output                                    |
|-------------------------------------------|-------------------------------------------|
| Selling price for products [€/ton] or [€/m<sup>3</sup>] | Operating cost (OPEX) [€/year]          |
| Prices for energy [€/KWh], input chemicals [€/m<sup>3</sup>], cooling water [€/m<sup>3</sup>] | Investment cost (CAPEX) [€]               |
| Operating hours, lifetime                 | Revenues from selling products [€/year] |
| Interest rate, Inflation rate             |                                         |
|Equipment cost [€]  |                                          |
| Assumptions on CAPEX and OPEX calculations |                                          |


For the economic analysis of a full-scale desalination plant, the equipment costs of pilot-scale units are scaled-up to a capacity of 30000 m<sup>3</sup>/d. The equipment (material) costs of the full-scale plant are derived from the cost of the same equipment in the pilot plant with known capacity using function `scaleup.py`. 


**Table 3** gives an overview of the main assumptions made to calculate the CAPEX and OPEX. 
| CAPEX                             | Annual OPEX                                    |
|-----------------------------------|------------------------------------------------|
| Installation: 25% of purchased equipment cost| Maintenance: 3% of the fixed-capital investment            |
| Buildings, process, and auxiliary: 20% of purchased equipment cost| Operating Supplies: 5% of maintenance |
| Land: 6% of purchased equipment cost  | Operating Labor: 15% of annual OPEX                             |
| Indirect costs: 15% of direct cost                   | Direct supervisory and clerical labor:15% of operating labor                         |
| Working capital: 20% of total investment cost  | Laboratory charges: 15% of operating labor                         |
|                                   | Patents and royalties: 3% of annual OPEX                          |
|                                   | Fixed charges: 5% of annual OPEX                                  |
|                                   | Plant overhead costs: 5% of annual OPEX                           |


### 6. Treatment chains comparison 
In comparison file, results from different treatment chains are summarised. Indicators are formulated to compare the treatment chains. 

#### Import results
First, import the results from the two examples.

```python
#Import results 
import example_1 as sc1
import example_2 as sc2
```
Import required functions: 
```python
import numpy as np
import pandas as pd
import numpy as np 
```
#### Create lists with results
```python
X = ['Example 1', 'Example 2']
X_axis = np.arange(len(X))
```
#### Electrical consumption Vs Thermal consumption  
For instance, the two examples are compared based on their electrical and thernal energy requirements. 

```python
  # Create lists for OPEX and asigned reuslts for Electrical consumption and thermal consumption 
Eel = [ sc1.sum_el_en, sc2.sum_el_en]
Eth = [sc1.sum_th_en,   sc2.sum_th_en]
  # Yearly calculation 
Eel = [i * hr/1e6 for i in Eel] # Total electrical energy consumption 
Eth = [i * hr/1e6 for i in Eth] # Total thermal energy consumption 
```
#### Visualization 
For the visualization, a bar figure is created and saved.  

```python
 # Create Figure 1: Electrical consumption Vs thermal consumption 
plt.bar(X_axis - 0.2, Eel, 0.4, color="#00516a", label = 'Electrical (GWel)')
plt.bar(X_axis + 0.2, Eth, 0.4, color="sandybrown", label = 'Thermal (GWth)')
plt.xticks(X_axis, X)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel("Scenarios")
plt.ylabel("Eenergy consumption (GW)")
plt.legend()
plt.savefig('electricVSthermal.png')
plt.show()
```
<figure>
  <img src="https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/assets/150446818/640b2dde-d5c6-439d-a44d-fc0e73dcf342" alt="Image" style="width:600px;">
</figure>

#### Operating costs (OPEX)
For instance, the two examples are compared based on the operating costs (OPEX). 

```python
  # Create lists for OPEX and asigned reuslts 
OPEX = [sc1.OPEX, sc2.OPEX]
# Yearly calculation 
OPEX = [i/1e6 for i in OPEX]
```
#### Visualization 
For the visualization, a bar figure is created and saved.  

```python
# Create Figure 2: OPEX 
plt.bar(X_axis - 0.4, OPEX, 0.4, color="#00516a")    
plt.xticks(X_axis, X)
plt.xlabel("Scenarios")
plt.ylabel("OPEX (M€/year)")
plt.savefig('OPEX.png')
plt.show()

```

<figure>
  <img src="https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/assets/150446818/71c0fa99-49f8-4b69-916a-48cdaa6d303e" alt="Image" style="width:600px;">
</figure>



