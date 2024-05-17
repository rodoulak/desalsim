# DesalSim: Desalination and Brine Treatment simulation model

 
## Overview 
This project comprises a package of simulation models for desalination and brine treatment technologies, including reverse osmosis, nanofiltration, multi-effect distillation, chemical precipitation, eutectic freeze crystallization, electrodialysis, and thermal crystallization. These technologies are utilized for desalination, ion separation, and salt recovery processes in various industrial and environmental applications.

The simulation models implemented here calculate various parameters such as salt concentration profiles, ion fluxes, energy consumption, chemical consumption, and operational costs. They provide insights into the performance and economics of the technologies under different operating conditions and input parameters. Additionally, the models allow for the integration of different technologies in various configurations to optimize process efficiency and resource utilization.

![system description](https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/assets/150446818/bb10e07d-b878-45c8-878a-0c56222546cf)

## Table of contents
* [Purpose](#purpose)
* [Key features](#key-features)
* [Installation](#installation)
* [Usage](#usage)
* [Contributing](#contributing)
* [Authors and Acknowledgements](#authors-and-acknowledgements)
* [License](#license)

## Purpose 
The purpose of this software suite is to provide researchers, engineers, and policymakers with a powerful tool for evaluating the performance and economics of desalination and brine treatment systems. By integrating technical process models with economic and environmental analyses, the suite enables users to make informed decisions about technology integration, process optimization, and resource management.


## Installation  
The easiest way is through pip, in command-line interface:   
```
pip install desalsim==1.0.2
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
## Key features 
- Simulation models for various desalination and brine treatment technologies
- Analysis of salt concentration profiles, ion fluxes, energy consumption, and chemical usage
- Integration of different technologies to optimize process efficiency and resource utilization
- Economic models for technologies and integrated systems 
- Output visualization and data export for further analysis 

### File Descriptions

A brief overview of each file in the project can be found in (https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/tree/main/file-description). 

## Usage 
Each simulation model serves as a standalone tool for analyzing the performance of a specific desalination or brine treatment technology. Before running the simulation, ensure that you have provided the necessary input parameters, such as feed flow rates, salinity levels, membrane properties, heat sources, and operating conditions.

The simulation results, including salt concentration profiles, ion fluxes, energy consumption, chemical consumption, and operational costs, will be generated based on the specified inputs and displayed in the console output or saved to output files for further analysis.

However, simulation models of more than one technology can be combined to simulate and evaluate the performance of a treatment chain (desalination and brine treatment system). In this case, the output flow rates and stream concentrate are the input data for the next technology. 

Additionally, two example files are provided to demonstrate the usage of the simulation suite (see [Example Folder](https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/tree/main/example), [Example 1](https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/tree/main/example/example_1.py) and [Example 2](https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/tree/main/example/example_2.py)). These examples simulate and evaluate two different treatment chains, showcasing the integration of multiple technologies. The economic evaluation of the treatment chain is given in [Example 1](https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/tree/main/example/example_1.py) and in the [Economic Tutorial](https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/tree/main/example/Economic_Tutorial.md). 

Furthermore, a [comparison file](https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/tree/main/example/comparison.py) is included, where the results of the two examples are compared in terms of various parameters. Users can extend this comparison by adding more indicators as needed.

For more details on input/output parameters and assumption see [Link to Tutorial File](https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/tree/main/example/Tutorial.md).

### Documentation 
You can find Tutorials and documents at [Tutorial File](https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/tree/main/example/Tutorial.md). Additionally, you can find tests for every process unit and the economic model in the [tests folder](https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/tree/main/tests) that verify that the code is running properly. 


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
- EDBM: Calogero Cassaro, Andrea Culcasi, Giovanni Virruso
- EFC: Marcos Rodriguez Pascual.

This work was supported by the EU within the WATER MINING (Next generation water-smart management systems: large scale demonstrations for a circular economy and society) - Horizon 2020 research and innovation programme under grant agreement No 869474.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

