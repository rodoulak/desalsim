# Tutorial 
## Introduction 

## Usage 
Each simulation model serves as a standalone tool for analyzing the performance of a specific desalination or brine treatment technology. Before running the simulation, ensure that you have provided the necessary input parameters, such as feed flow rates, salinity levels, membrane properties, heat sources, and operating conditions.

The simulation results, including salt concentration profiles, ion fluxes, energy consumption, chemical consumption, and operational costs, will be generated based on the specified inputs and displayed in the console output or saved to output files for further analysis.

However, simulation models of more than one technology can be combined to simulate and evaluate the performance of a treatment chain (desalination and brine treatment system). In this case, the output flow rates and stream concentrate are the input data for the next technology. 

Additionally, two example files are provided to demonstrate the usage of the simulation suite (see [Example Folder](example), [Example 1](example/example_1.py) and [Example 2](example/example_2.py)). These examples simulate and evaluate two different treatment chains, showcasing the integration of multiple technologies. Furthermore, a [comparison file](example/comparison.py) is included, where the results of the two examples are compared in terms of various parameters. Users can extend this comparison by adding more indicators as needed.

**<u>Followed steps:</u>**<br>
Step 1: Import required fucntions for process units in the treatment chain.<br>
Step 2: Set input data like feed flow rate, ion concentration, relevant ions for the feed solution.<br>
Step 3: Set input parameters for each process unit as shown in **Table 1** and for economic model as shown in **Table 2** and **Table 3**.<br>
Step 4: Call function of each process unit, creat objects for each calculation.<br>
Step 5: Results interpretation. <br>
<br>

### Technical process models 
For more detailed steps and instructions see [Tutorial for Example 1](example/Example_1_Tutorial.md).  
The mathematical description of each technology is given in  [Mathematical description](paper/Mathematical_description.pdf).  

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


### Economic models 

For more detailed steps and instructions see [Tutorial for Economic model](example/Economic_Tutorial.md).  
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


### Treatment chains comparison 
In comparison file, results from different treatment chains are summarised. Indicators are formulated to compare the treatment chains. 
For example, the two examples are compared based on the operating costs (OPEX). 
<figure>
  <img src="https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/assets/150446818/71c0fa99-49f8-4b69-916a-48cdaa6d303e" alt="Image" style="width:600px;">
</figure>



