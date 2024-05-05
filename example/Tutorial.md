# Tutorial 

## Usage 
Each simulation model serves as a standalone tool for analyzing the performance of a specific desalination or brine treatment technology. Before running the simulation, ensure that you have provided the necessary input parameters, such as feed flow rates, salinity levels, membrane properties, heat sources, and operating conditions.

The simulation results, including salt concentration profiles, ion fluxes, energy consumption, chemical consumption, and operational costs, will be generated based on the specified inputs and displayed in the console output or saved to output files for further analysis.

However, simulation models of more than one technology can be combined to simulate and evaluate the performance of a treatment chain (desalination and brine treatment system). In this case, the output flow rates and stream concentrate are the input data for the next technology. 

Additionally, two example files are provided to demonstrate the usage of the simulation suite. These examples simulate and evaluate two different treatment chains, showcasing the integration of multiple technologies. Furthermore, a comparison file is included, where the results of the two examples are compared in terms of various parameters. Users can extend this comparison by adding more indicators as needed.

### Example 1
**Figure 1** presents the process flow diagram of example 1 which consists of four technologies: Nanofiltration (NF), Multiple Feed Plug Flow Reactor (MF-PFR), Electrodialysis (ED), Electrodialysis With Bipolar Membranes (EDBM). The treatment chain represents an MLD system aiming to maximize valuable resources recovery from brine, such as Mg(OH)<sub>2</sub>, Ca(OH)<sub>2</sub>, HCl, and NaOH. The seawater stream or concentrate stream from a Reverse Osmosis plant (RO) first goes to the NF unit. The NF unit is separated into two different streams: one that is high in monovalent ions and one that is high in multi-valent ions. The latter stream from NF, high in monovalent ions, is directed to ED, in which the NaCl stream is concentrated further, and a dilute stream is also recovered. The former is directed to a treatment line comprising selective MF-PFR and EDBM units. In particular, the retentate is sent to the MF-PFR, in which magnesium and calcium are recovered in the form of hydroxide precipitates via a chemical reaction between the NF retentate and an alkaline reactant. Then, the brine stream is free from Mg<sup>2+</sup> and Ca<sup>2+</sup> mixed with the ED concentrate stream. The mixed solution (NaCl rich) is fed to EDBM. EDBM unit recovers, and the saline solution (low concentration) can be recycled back into the treatment chain. 
<figure>
  <img src="https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/assets/150446818/55cc6b6f-dde8-4b12-ae61-fa23665c288e" alt="Image" style="width:600px;">
</figure>

**Figure 1**. Process flow diagram of example 1.
<br>

**<u>Followed steps:</u>**<br>
Step 1: Import required fucntions for process units in the treatment chain.<br>
Step 2: Set input data like feed flow rate, ion concentration, relevant ions for the feed solution.<br>
Step 3: Set input parameters for each process unit as shown in **Table 1** and for economic model as shown in **Table 2** and **Table 3**.<br>
Step 4: Call function of each process unit, creat objects for each calculation.<br>
Step 5: Results interpretation. <br>
<br>

For more details see 
**Table 1** gives an overview of the main inputs and outputs for each process unit of example 1. 
| Process                                   | Input                                       | Output                                                |
|-------------------------------------------|---------------------------------------------|-------------------------------------------------------|
| Nanofiltration                            | Feed flow rate [m³/h]                       | Permeate flow rate and composition [g/L]              |
|                                           | Ion concentration [g/L]                     | Concentrate flow rate and composition [g/L]           |
|                                           | Osmotic pressure [bar]                      | Electrical requirements [kWhel]                       |
|                                           | Water recovery [%]                          | Ion rejection [-]                                     |
| Multi-plug flow reactor                   | Feed flow rate [m³/h]                       | Alkaline solution flow rate [L/h]                    |
|                                           | Ion concentration [g/L]                     | Flow rate of Mg(OH)₂ [kg/h]                          |
|                                           | Concentration of the alkaline solution [M] | Flow rate of Ca(OH)₂ [kg/h]                          |
|                                           | Concentration of the acid solution [M]     | Acid solution flow rate [L/h]                        |
|                                           | Products characteristics e.g. solubility...| Effluent flow rate [m³/h] and composition [g/L]      |
|                                           |                                             | Electricity requirements [kWhel]                     |
| Electrodialysis with bipolar membranes   | Feed flow rate [m³/h]                       | Flow rate of acid [m³/h] and composition [g/L]       |
|                                           | Ion concentration [g/L]                     | Flow rate of base [m³/h] and composition [g/L]       |
|                                           | Current density [A/m²]                      | Flow rate of salt [m³/h] and composition [g/L]       |
|                                           | Number of triplets and Membrane area and other characteristics      | Electricity requirements [kWhel]                     |
| Electrodialysis                          | Feed flow rate [m³/h]                       | Flow rate of diluted stream [m³/h] and composition [g/L]|
|                                           | Ion concentration [g/L]                     | Flow rate of concentrate stream [m³/h] and composition [g/L]        |
|                                           | Current density [A/m²]                      | Electricity requirements [kWhel]                     |


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

### Comparison file 
In comparison file, results from different treatment chains are summarised. Indicators are formulated to compare the treatment chains. 
