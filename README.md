# Desalination and Brine Treatment simulation model

 
## Overview 
This project comprises a package of simulation models for desalination and brine treatment technologies, including reverse osmosis, nanofiltration, multi-effect distillation, chemical precipitation, eutectic freeze crystallization, electrodialysis, and thermal crystallization. These technologies are utilized for desalination, ion separation, and salt recovery processes in various industrial and environmental applications.

The simulation models implemented here calculate various parameters such as salt concentration profiles, ion fluxes, energy consumption, chemical consumption, and operational costs. They provide insights into the performance and economics of the technologies under different operating conditions and input parameters. Additionally, the models allow for the integration of different technologies in various configurations to optimize process efficiency and resource utilization.

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


## Key features 
- Simulation models for various desalination and brine treatment technologies
- Analysis of salt concentration profiles, ion fluxes, energy consumption, and chemical usage
- Integration of different technologies to optimize process efficiency and resource utilization
- Economic models for technologies and integrated systems 
- Output visualization and data export for further analysis 

### File Descriptions

Here is a brief overview of each file in the project:

- `README.md`: This README file providing an overview of the project, usage instructions, and other details.
- `nanofiltration_unit_f.py`: This file contains the main script for running the simulation model of Nanofiltration (NF) unit.
- `ro_unit.py`: This file contains the main script for running the simulation model of Reverse Osmosis (RO) unit.
- `med_unit_f.py`: This file contains the main script for running the simulation model of Multi-effect distillation (MED) unit.
- `mfpfr_unit.py`: This file contains the main script for running the simulation model of Multiple Feed Plug Flow Reactor (MF-PFR) unit.
- `economic_f.py`: This file contains the main script for running the economic model of a technology or a system of technologies.
- `ed_unit_f.py`: This file contains the main script for running the simulation model of Electrodialysis (ED) unit.
- `edbm_unit_f.py`: This file contains the main script for running the simulation model of Electrodialysis With Bipolar Membranes (EDBM) unit.
- `efc_unit_f.py`: This file contains the main script for running the simulation model of Eutectic freeze crystallization (EFC) unit.
- `thermal_cryst_unit_f.py`: This file contains the main script for running the simulation model of Thermal crystallizer (Tcr) unit.
- `example1.py`: Example script demonstrating the usage of the simulation suite for a specific treatment chain.
- `example2.py`: Another example script showcasing a different treatment chain simulation.
- `comparison.py`: Script for results intrepetation. It compares the results of example1 and example2 simulations.
- `constants.py`: This file contains the main script for input contant parameters.
- `density_calc.py`: This file contains the main script for calculating density of a solution based on the tempretaure and the concentration of the solution.
- `scaleup.py`: This file contains the main script for scaling-up a technology with known capacity using the six-tenths factor rule (m=0.6).
- `LICENSE`: License file containing the MIT License terms.

## Installation  
The easiest way is through pip, in command-line interface:   
```
pip install desalination-and-brine-treatment==0.2
```

If you want to install the latest GitHub verstion:
1. Clone the repository to your local machine:
```
git clone https://github.com/your-username/desalination-brine-treatment-simulation.git
```
2. Install the required dependencies:
 ```
pip install -r requirements.txt
 ```

## Usage 
Each simulation model serves as a standalone tool for analyzing the performance of a specific desalination or brine treatment technology. Before running the simulation, ensure that you have provided the necessary input parameters, such as feed flow rates, salinity levels, membrane properties, heat sources, and operating conditions.

The simulation results, including salt concentration profiles, ion fluxes, energy consumption, chemical consumption, and operational costs, will be generated based on the specified inputs and displayed in the console output or saved to output files for further analysis.

However, simulation models of more than one technology can be combined to simulate and evaluate the performance of a treatment chain (desalination and brine treatment system). In this case, the output flow rates and stream concentrate are the input data for the next technology. 

Additionally, two example files are provided to demonstrate the usage of the simulation suite (see [Example Folder](example), [Example 1](example/example_1.py) and [Example 2](example/example_2.py)). These examples simulate and evaluate two different treatment chains, showcasing the integration of multiple technologies. The economic evaluation of the treatment chain is given in [Example 1](example/example_1.py) and in the [Economic Tutorial](example/Economic_Tutorial.md). 

Furthermore, a [comparison file](example/comparison.py) is included, where the results of the two examples are compared in terms of various parameters. Users can extend this comparison by adding more indicators as needed.

For more details on input/output parameters and assumption see [Link to Tutorial File](example/Tutorial.md).

### Documentation 
You can find Tutorials and documents at: 
1. [Tutorial File](example/Tutorial.md)
2. [Tutorial for Example 1](example/Example_1_Tutorial.md)
3. [Economic Tutorial](example/Economic_Tutorial.md)
4. The mathematical description of each technology is given in [Mathematical description](paper/Mathematical_description.pdf)
5. [Example 1](example/example_1.py)
6. [Example 2](example/example_2.py)


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

