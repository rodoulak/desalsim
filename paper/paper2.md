---
title: 'Desalination and brine treatment systems integrated modelling framework: simulation and evaluation of water and resource recovery'
tags:
  - Python
  - Desalination
  - Brine treatment
  - resource recovery
  - Techno-economic assessment 
authors:
  - name: Rodoula Ktori
    affiliation: 1 
  - name: Fabrizio Vassallo
    affiliation: 2
  - name: Giovanni Virruso
    affiliation: 2
  - name: Carmelo Morgante
    affiliation: 2
  - name: Andrea Culcasi
    affiliation: 2
  - name: Dionysia Diamantidou
    affiliation: 3
  - name: Niels Van Linden
    affiliation: 3
  - name: Alessandro Trezzi
    affiliation: 4
  - name: Adithya Krishnan
    affiliation: 5
  - name: Andrea Cipollina
    affiliation: 2
  - name: Giorgio Micale
    affiliation: 2
  - name: Mark C.M. van Loosdrecht
    affiliation: 1
  - name: Dimitrios Xevgenos
    affiliation: 6
affiliations:
  - name: Department of Biotechnology, Delft University of Technology, Van der Maasweg 9, 2629 HZ, Delft, The Netherlands
    index: 1
  - name: Dipartimento di Ingegneria, Università degli Studi di Palermo - viale delle Scienze Ed.6, 90128 Palermo, Italy
    index: 2
  - name: Lenntech BV, Distributieweg 3, 2645 EG Delfgauw, The Netherlands
    index: 3
  - name: Sofinter S.p.A, Piazza Francesco Buffoni, 3, 21013 Gallarate VA, Italy
    index: 4
  - name: Water & Energy Intelligence BV, the Netherlands
    index: 5
  - name: Department of Engineering Systems and Services, Delft University of Technology, Jaffalaan 5, 2628 BX Delft, The Netherlands
    index: 6
    
bibliography: paper.bib
---

# Summary 

Desalination plays a crucial role in addressing the growing challenges of water scarcity. In recent years, the integration of desalination and brine treatment technologies has been increasingly studied, aiming to develop sustainable solutions for resource recovery from seawater. However, designing treatment chains and optimizing these processes for maximum efficiency, sustainability, and cost-effectiveness are complex tasks that require data, sophisticated analysis and decision-making strategies.
Our software offers a comprehensive modelling framework for simulating desalination and minerals recovery systems. Integrating technical process models with economic and environmental analysis provides valuable insights into the integration of these technologies and their impact on production efficiency, energy consumption, and environmental performance.
Through our software's simulations, researchers, engineers, and policymakers gain the power to evaluate the resource recovery potential, economics, and greenhouse gas emissions associated with different configuration combinations. This empowerment with crucial information for early-stage design and strategic planning is a significant step toward fostering more sustainable water management practices.

# Statement of need
Traditionally, simulation models were developed to evaluate the influence of certain parameters on the characteristics of the recovered products and the performance of the technology in terms of energy, chemicals, and water consumption. In the desalination field, there is a lack of open-access simulation models. While commercial software programs, like WAVE [@dupontwebsite], exist for membranes, and numerous publications discuss techno-economic models for desalination [@el2002fundamentals] and brine treatment technologies [@xevgenos2015design, @micari2020techno, @chen2021zero, @panagopoulos2021beneficiation, @poirier2022techno, @morgante2022valorisation], there is a noticeable absence of open-access simulation tools in the literature. With the shift towards circular systems and integrated desalination and brine treatment technologies for resource recovery, there is a need for a unified tool that integrates simulation models for each technology. 
Simulation and evaluation are essential for informing future design and decision-making processes. An open-access simulation tool is essential in the field of desalination to improve the credibility of evaluations, facilitate study repeatability, support validation efforts, and enable comparisons.
Therefore, while the modelling of desalination and brine treatment processes has been extensively studied in the literature, our software offers an integrated platform for integrating these models and facilitating seamless data exchange. Based on this need, we provide, for the first time, open-source software that brings together a range of technical processes and economic models developed over the last decades. By making these models accessible and transparent, our software aims to advance research, engineering, and decision-making in the field of desalination and resource recovery.
The proposed software provides a variety of examples to help modellers design and evaluate different technical configurations. 

# Limitations 
The proposed software is not intended to be a substitute for sophisticated physical models or a system dynamics approach. In cases requiring detailed process representations, the process model may need to be enhanced to provide more detailed results and optimization opportunities. This work highlights that the proposed software is particularly valuable for evaluating the integration of different processes and preliminary designs, capturing the technical, economic, and environmental impacts of technology integration. 

# Acknowledgements 
The software was developed by Rodoula Ktori, with theoretical support from all co-authors. Technological experts conducted the validation of each simulation model for the respective technology: Nanofiltration: Dionysia Diamantidou, Niels van Linden; Multi-effect Distillation: Alessandro Trezzi; MF-PFR: Fabrizio Vassallo, Carmelo Morgante, Andrea Cipollina; EDBM: Giovanni Virruso, Andrea Culcasi; EFC: Marcos Rodriguez Pascual.

The technical process models were developed based on the report from [@xevgenos2023deliverable] and the following literature [@nayar2019cost, @morgante2022valorisation, @cassaro2023electrodialysis]. Then they were valisated with experimental results from [@morgante2024pioneering]. The development of economic models were influenced by [@peters1968plant, @abraham2011socio, @bilton2011feasibility, @choi2015economic, @kesieme2013economic]. 

This work was supported by the EU within the WATER MINING (Next generation water-smart management systems: large scale demonstrations for a circular economy and society) - Horizon 2020 research and innovation programme under grant agreement No 869474.

# References 

