# Desalination and Brine Treatment Simulation Suite

 
## Overview 
This project comprises a suite of simulation models for desalination and brine treatment technologies, including reverse osmosis, nanofiltration, multi-effect distillation, chemical precipitation, eutectic freeze crystallization, electrodialysis, and thermal crystallization. These technologies are utilized for desalination, ion separation, and salt recovery processes in various industrial and environmental applications.

The simulation models implemented here calculate various parameters such as salt concentration profiles, ion fluxes, energy consumption, chemical consumption, and operational costs. They provide insights into the performance and economics of the technologies under different operating conditions and input parameters. Additionally, the models allow for the integration of different technologies in various configurations to optimize process efficiency and resource utilization.

## Purpose 
The purpose of this software suite is to provide researchers, engineers, and policymakers with a powerful tool for evaluating the performance and economics of desalination and brine treatment systems. By integrating technical process models with economic and environmental analyses, the suite enables users to make informed decisions about technology integration, process optimization, and resource management.

## Installation 
To run the simulation models, follow these steps:
1. Clone the repository to your local machine:
```
git clone https://github.com/your-username/desalination-brine-treatment-simulation.git
```
2. Install the required dependencies:
 ```
pip install -r requirements.txt
 ```
3.  Choose the simulation model for the desired technology (e.g., electrodialysis, reverse osmosis, thermal desalination) and run the corresponding script.

## Usage 
Each simulation model serves as a standalone tool for analyzing the performance of a specific desalination or brine treatment technology. Before running the simulation, ensure that you have provided the necessary input parameters, such as feed flow rates, salinity levels, membrane properties, heat sources, and operating conditions.

The simulation results, including salt concentration profiles, ion fluxes, energy consumption, chemical consumption, and operational costs, will be generated based on the specified inputs and displayed in the console output or saved to output files for further analysis.

However, simulation models of more than one technology can be combined to simulate and evaluate the performance of a treatment chain (desalination and brine treatment system). In this case, the output flow rates and stream concentrate are the input data for the next technology. 

Additionally, two example files are provided to demonstrate the usage of the simulation suite. These examples simulate and evaluate two different treatment chains, showcasing the integration of multiple technologies. Furthermore, a comparison file is included, where the results of the two examples are compared in terms of various parameters. Users can extend this comparison by adding more indicators as needed.

### Example 1
**Figure 1** presents the process flow diagram of example 1 which consists of four technologies: Nanofiltration (NF), Multiple Feed Plug Flow Reactor (MF-PFR), Electrodialysis (ED), Electrodialysis With Bipolar Membranes (EDBM). The treatment chain represents an MLD system aiming to maximize valuable resources recovery from brine, such as Mg(OH)2, Ca(OH)2, HCl, and NaOH. The seawater stream or concentrate stream from a Reverse Osmosis plant (RO) first goes to the NF unit. The NF unit is separated into two different streams: one that is high in monovalent ions and one that is high in multi-valent ions. The latter stream from NF, high in monovalent ions, is directed to ED, in which the NaCl stream is concentrated further, and a dilute stream is also recovered. The former is directed to a treatment line comprising selective MF-PFR and EDBM units. In particular, the retentate is sent to the MF-PFR, in which magnesium and calcium are recovered in the form of hydroxide precipitates via a chemical reaction between the NF retentate and an alkaline reactant. Then, the brine stream is free from Mg2+ and Ca2+ mixed with the ED concentrate stream. The mixed solution (NaCl rich) is fed to EDBM. EDBM unit recovers, and the saline solution (low concentration) can be recycled back into the treatment chain. 
<figure>
  <img src="https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/assets/150446818/55cc6b6f-dde8-4b12-ae61-fa23665c288e" alt="Image" style="width:600px;">
  <figcaption style="text-align:left;">Figure 1. Process flow diagram of example 1.</figcaption>
</figure>

**Table 1** gives an overview of the main inputs and outputs for each process unit of example 1. 

![image](https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/assets/150446818/0f1c7b73-1a72-459a-9b72-08e6ef7cc8b8)

**Table 2** gives an overview of the main inputs and outputs of economic model. 
![image](https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/assets/150446818/459ed1f1-65d7-4ca1-b712-bf2d8eb9a7db)

 Header 1 | Header 2 | Header 3 
----------|----------|----------
 Row 1, Col 1 | Row 1, Col 2 | Row 1, Col 3 
 Row 2, Col 1 | Row 2, Col 2 | Row 2, Col 3 

## Key features 
- Simulation models for various desalination and brine treatment technologies
- Analysis of salt concentration profiles, ion fluxes, energy consumption, and chemical usage
- Integration of different technologies to optimize process efficiency and resource utilization
- Economic models for technologies and integrated systems 
- Output visualization and data export for further analysis

## Contributing
Contributions to the project are welcome! If you'd like to contribute, please follow the standard GitHub workflow:
1. Fork the repository.
2. Create a new branch (git checkout -b feature/improvement).
3. Make your changes and commit them (git commit -am 'Add new feature').
4. Push to the branch (git push origin feature/improvement).
5. Create a new Pull Request.

## Authors and Acknowledgements
The software was developed by Rodoula Ktori, with theoretical support by all co-authors. The validation of each simulation model for the respective technology was conducted by technological experts: 
- Nanofiltration: Dionysia Diamantidou and Niels van Linden,
- Multi-effect Distillation: Alessandro Trezzi,
- MF-PFR: Fabrizio Vassallo, Carmelo Morgante, Andrea Cipollina,
- EDBM: Calogero Cassaro, Andrea Culcasi,
- EFC: Marcos Rodriguez Pascual.

This work was supported by the EU within the WATER MINING (Next generation water-smart management systems: large scale demonstrations for a circular economy and society) - Horizon 2020 research and innovation programme under grant agreement No 869474.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

