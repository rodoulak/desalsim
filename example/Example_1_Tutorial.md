# Example 1: Tutorial 

## 1. Introduction
Welcome to our comprehensive tutorial on running simulation models and evaluating the 
performance of treatment chains for water purification processes. In this tutorial, 
we provide step-by-step instructions on how to [create treatment chains](#create-treatment-chain), 
[define feed characteristics](#define-feed-characteristics), [use process unit models](#use-process-unit-model), and [analyze the results obtained from the simulation](#results-evaluation ), using **Example 1** as case study.
Each section guides you through setting up the simulation environment, running the models, and interpreting the results. Additionally, we discuss technical, economic, and environmental indicators to evaluate the performance of the treatment chain. 

### Example description 
**Figure 1** presents the process flow diagram of example 1 which consists of four technologies: Nanofiltration (NF), Multiple Feed Plug Flow Reactor (MF-PFR), Electrodialysis (ED), Electrodialysis With Bipolar Membranes (EDBM). The treatment chain represents an MLD system aiming to maximize valuable resources recovery from brine, such as Mg(OH)<sub>2</sub>, Ca(OH)<sub>2</sub>, HCl, and NaOH. The seawater stream or concentrate stream from a Reverse Osmosis plant (RO) first goes to the NF unit. The NF unit is separated into two different streams: one that is high in monovalent ions and one that is high in multi-valent ions. The latter stream from NF, high in monovalent ions, is directed to ED, in which the NaCl stream is concentrated further, and a dilute stream is also recovered. The former is directed to a treatment line comprising selective MF-PFR and EDBM units. In particular, the retentate is sent to the MF-PFR, in which magnesium and calcium are recovered in the form of hydroxide precipitates via a chemical reaction between the NF retentate and an alkaline reactant. Then, the brine stream is free from Mg<sup>2+</sup> and Ca<sup>2+</sup> mixed with the ED concentrate stream. The mixed solution (NaCl rich) is fed to EDBM. EDBM unit recovers, and the saline solution (low concentration) can be recycled back into the treatment chain. 

<figure>
  <img src="https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/assets/150446818/55cc6b6f-dde8-4b12-ae61-fa23665c288e" alt="Image" style="width:600px;">
</figure>

**Figure 1**. Process flow diagram of example 1.
<br>

## 2. Installation
Instructions on how to install the required dependencies or software.

## 3. Running the Simulation models

### 3.1. Create treatment chain 
To create the treatment chain, the required units have to be imported. 
For **Example 1** which consists of four technologies:
- Nanofiltration (NF),
- Multiple Feed Plug Flow Reactor (MF-PFR),
- Electrodialysis (ED),
- Electrodialysis With Bipolar Membranes (EDBM)
their functions are imported:
```python
import nanofiltration_unit_f
```
or
```python
from nanofiltration_unit_f import OsmoticPressure
from nanofiltration_unit_f import NFMass
from nanofiltration_unit_f import NfEnergy
```
Similarly for the other process units. Additionally, function for calculating density (`density_calc.py`) or constants (`comparison.py`) where user can add constant values like MW, prices etc, need to be imported. 
```python
from density_calc import density_calc
import constants
import scaleup
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


### 3.2.2.  Multi-plug flow reactor

To run simulation model for Multi-plug flow reactor (MFPFR) unit, you need to implement the following steps. 

**Table 2** gives an overview of the main inputs and outputs for each process unit of Multi-plug flow reactor. 
| Process                                   | Input                                       | Output                                                |
|-------------------------------------------|---------------------------------------------|-------------------------------------------------------|
| Multi-plug flow reactor                   | Feed flow rate [m³/h]                       | Alkaline solution flow rate [L/h]                    |
|                                           | Ion concentration [g/L]                     | Flow rate of Mg(OH)₂ [kg/h]                          |
|                                           | Concentration of the alkaline solution [M] | Flow rate of Ca(OH)₂ [kg/h]                          |
|                                           | Concentration of the acid solution [M]     | Acid solution flow rate [L/h]                        |
|                                           | Products characteristics e.g. solubility...| Effluent flow rate [m³/h] and composition [g/L]      |
|                                           |                                             | Electricity requirements [kWhel]                     |

##### Setting input flow rate and ion concentration 
First, the results from the previous process unit (in this case, Nanofiltration) need to be assigned as input parameters for the next process (MFPFR).

```python
    # The feed ion concentration is concentration of NFconcentrate stream.
Cin_mfpfr = Cconc  # Concentrations of [Na, Cl, K, Mg, Ca, SO4]
    # Flow rate in l/hr
Qin_mfpfr = Qconc
    # Calculate the density of the input
d_in = density_calc(25, sum(Cin_mfpfr)) / 1000
```
##### Setting other input parameters
Then the required input for MFPFR unit need to be added from unser. 

First, the concentration of the alkaline solution (NaOH) and acid solution (HCl) are import. Note that different chemicals and concentrations can be used for the percicipation and the pH neutralization.

```python
    # Concentration of NaOH solution for step 1 and step 2 in MOL/L
C_NaOH = [1, 1]
    # Concentration of HCl in MOL/L
HCl_conc=1 
```
Then the conversion rate of Mg in every step has to be assumed. The assumption relies on experimental data.
```python
    # Conversion rate for step 1 and step 2
conv = [95, 93]  
```

Finally, the products characteristics need to be set. 

```python
    # Product solublity of Mg(OH)2
kps_MgOH=5.61*0.000000000001
    # Product solublity of Ca(OH)2
kps_CaOH=5.5*0.000001
    # Mg(OH)2 density (units: kg/l)
d_mgoh_2=2.34
    # Ca(OH)2 density (units: kg/l)
d_caoh_2=2.211 
```
##### Create MFPFRCALC Objects
```python
    # Create an instance of the inputpar class with the defined parameters
mfpfr_dat = MFPFRCALC(Qin_mfpfr, Cin_mfpfr, *C_NaOH, *conv)
```

##### Calculations for precipitation steps
```python
    # Call the calc_step1 and calc_step2 methods to calculate the necessary values
mfpfr_dat.calc_step1(kps_MgOH, d_mgoh_2)
mfpfr_dat.calc_step2(d_mgoh_2, d_caoh_2 )
ph_2=mfpfr_dat.ph_2
```

Get total outlet flowtate after the second step of precipitation
```python
    # Outlet flow rate 
Qout_2 = mfpfr_dat.Qout_2
```

Get the total mass of the recovered products from each precipitation step
```python
    # Get the masses of Mg(OH)2 in the first step 
M_MgOH2_1 = mfpfr_dat.M_MgOH2_1

    # Get the masses of Mg(OH)2 and Ca(OH)2 in the second step 
M_CaOH2 = mfpfr_dat.M_CaOH2_2
M_MgOH2 = mfpfr_dat.M_MgOH2_2

print("Mg(OH)2 mass flow rate is "+str(round(M_MgOH2_1,2))+"kg/hr")
print("Ca(OH)2 mass flow rate is "+str(round(M_CaOH2,2))+"kg/hr")
```
Mg(OH)2 mass flow rate is 401.57kg/hr  
Ca(OH)2 mass flow rate is 95.96kg/hr  

##### Calculate outlet concentration 
```python
    # Create a list of the outlet concentrations in mol/l
Cout_all_m = [mfpfr_dat.CNa_out_2, mfpfr_dat.CCl_out_2, mfpfr_dat.CK_out_2, mfpfr_dat.CMg_out_2, mfpfr_dat.CCa_out_2, mfpfr_dat.CSO4_out_2]

    # Calculate the outlet concentrations in g/l
MW = [MW_Na, MW_Ca, MW_Cl, MW_K, MW_Mg, MW_SO4]
Cout_mfpfr_g = [Cout_all_m[i] * MW[i] for i in range(len(Cout_all_m))] # g/l
```
##### Calculate NaOH consumption 
```python
    # Calculate the chemical consumption of NaOH
QNAOH = mfpfr_dat.QNaOH_1 + mfpfr_dat.QNaOH_2_add + mfpfr_dat.QNaOH_2_st # convert to kg
print(f"NaOH flow rate is {round(QNAOH,2)} l/hr")
```
NaOH flow rate is 21918.92 l/hr  

##### Calculate amount of HCl for pH neutralization and the final outlet concentration from MF-PFR unit after pH neutralization
```python
    # Create an instance of the HCladdition class
unit = HClAddition(Qout_2, Cout_all_m, MW_Cl, ph_2, HCl_conc)

    # Call the calculate_HCladdition method
QHCl, Cout_mfpfr_g = unit.calculate_HCl_addition(Cout_mfpfr_g)

    # Print the volume of HCl added 
print(f"HCl flow rate is {round(QHCl,2)} l/hr")
```
HCl flow rate is 6025.53 l/hr  

##### Calculate final outlet flow rate from MF-PFR unit 
```python
    # Volumetric flow rate 
Qout_f=Qout_2+QHCl #l/h
# Density calculation 
d_out_s=density_calc(25,round(sum(Cout_mfpfr_g),2))/1000 #kg/m3

    # Mass flow rate 
M_mfpfr_out=Qout_f*d_out_s #kg/h
print("Total effluent flow rate is "+str(round(M_mfpfr_out,2))+"kg/hr")
print("Total effluent flow rate is "+str(round(Qout_f,2))+"kg/hr")
```
Total effluent flow rate is 68505.96kg/hr  
Total effluent flow rate is 66280.95kg/hr

##### Calculate Energy consumption 
First, create an instance of the inputpar class with the defined parameters. 
```python
    # Create an instance of the inputpar class with the defined parameters
Epump_1, Epump_2=energycons.energycalc(mfpfr_dat.Qout_2, QNAOH, Qin_mfpfr, mfpfr_dat.QNaOH_1, mfpfr_dat.QNaOH_2_add, mfpfr_dat.QNaOH_2_st, dp, npump)
```
Calculate the total pumping energy including the HCl stream
```python
    # Electricity consumption for pumping , KWh
E_el_mfpf=(Epump_1+Epump_2+(QHCl*dp_HCl)*1e5/3600/(1000*npump))/1000
print("Total electricity energy consumption is "+str(round(E_el_mfpf,2))+ " KW")
```
Note that you can add a calculation for filtration unit and then sum the energy requirements. 
Specific energy consumption can also be calculated: 
```python
    # Specific energy consumption per kg of Mg(OH)2, KWh/kg of Mg(OH)2
SEC_el_prod=(E_el_mfpf)/(M_MgOH2)
print("Specific energy consumption per product is "+str(round(SEC_el_prod,2))+" KWh/kg product ")

    # Specific energy consumption per feed, KWh/m3 of feed
SEC_el_feed=(E_el_mfpf)/(Qin_mfpfr/1000)
print("Specific energy consumption per brine intake is "+str(round(SEC_el_feed,2))+" KWh/m3 of feed ")
```
Specific energy consumption per product is 2.88 KWh/kg product  
Specific energy consumption per brine intake is 1.49 KWh/m3 of feed 

### 3.2.2.  Other units
You need to follow similar steps for the other two processes. 
**Table 3** gives an overview of the main inputs and outputs for each process unit of Electrodialysis with bipolar membranes and Electrodialysis. 
| Process                                   | Input                                       | Output                                                |
|-------------------------------------------|---------------------------------------------|-------------------------------------------------------|
| Electrodialysis with bipolar membrane (EDBM)  | Feed flow rate [m³/h]                       | Flow rate of acid [m³/h] and composition [g/L]       |
|                                           | Ion concentration [g/L]                     | Flow rate of base [m³/h] and composition [g/L]       |
|                                           | Current density [A/m²]                      | Flow rate of salt [m³/h] and composition [g/L]       |
|                                           | Number of triplets and Membrane area and other characteristics      | Electricity requirements [kWhel]                     |
| Electrodialysis (ED)                         | Feed flow rate [m³/h]                       | Flow rate of diluted stream [m³/h] and composition [g/L]|
|                                           | Ion concentration [g/L]                     | Flow rate of concentrate stream [m³/h] and composition [g/L]        |
|                                           | Current density [A/m²]                      | Electricity requirements [kWhel]                     |

_Note that the feed flow rate and concentration of the units are the effluent flow rate and ions concentration of the unit before in the treatment chain._ 
In this treatment chain, Electrodialysis with bipolar membrane has two streams as feed for the salt channel. The two streams are mixed. For this the following calculations are required to calculate the new flow rate and concentration after the mixing. 
```python
    # Feed flow rate L/h
Q_in_edbm=M_mfpfr_out+Mc # Where M_mfpfr_out is the effluent from MF-PFR and Mc the effluent from ED

    # Feed concentration g/L
Cin_edbm=sum(Cout_mfpfr_g)*M_mfpfr_out/Q_in_edbm+Sc_o*Mc/Q_in_edbm # Where Cout_mfpfr_g is the effluent from MF-PFR and Sc_o the effluent from ED
C_in_mix=[]
for i in range(len(Cconc)):
    C_in_mix.append(Cout_mfpfr_g[i]*M_mfpfr_out/Q_in_edbm+Sc_out[i]*Mc/Q_in_edbm)
```

## 4. Results evaluation 
After the simulation of the treatment chain, the performance of the process units needs to be evaluated individually and overall as a system. 
### 4.1. Summarise results
The following code can be used to summarise the most important results from each process unit. Note that more results can be added. 
```python
class indic:
    
    def __init__(self, tech, Qout, Qin, Qprod_1, prod_1_name, Qprod_2, prod_2_name, E_el, E_th, Cin, Cout, chem1,chem2):
        self.tech=tech
        self.Qout=Qout
        self.Qin=Qin
        self.E_el=E_el
        self.E_th=E_th
        self.Cin=Cin
        self.Cout=Cout
        self.Qprod_1=Qprod_1
        self.Qprod_2=Qprod_2
        self.prod_1_name=prod_1_name
        self.prod_2_name=prod_2_name
        self.chem=chem1+chem2
        
        
    def techn_indi(self):        
        if self.prod_1_name=="water": 
            self.Water_rr=self.Qprod_1/self.Qin*100 # %
            self.Qwater=self.Qprod_1
        elif self.prod_2_name=="water":
            self.Water_rr=self.Qprod_2/self.Qin*100 # %
            self.Qwater=self.Qprod_2
        else:
            self.Qwater=0
            self.Water_rr=0
            print("no water production by " + self.tech)
```
Let's use Nanofiltration unit as example, here is how it can be used: 
```python
tec1=indic("NF", Qconc, Qf_nf, Qperm, "none", 0, "none", E_el_nf, 0, Ci_in, Cconc, QHCl_nf, Qantsc_nf)   
tec1.techn_indi()
```
##### Create lists with important results
After collecting all the important results, they can summarise in lists: 
```python
    # List results 
tec_names=[tec1.tech,tec2.tech, tec3.tech,  tec4.tech]
E_el_all=[tec1.E_el,tec2.E_el, tec3.E_el, tec4.E_el]
E_th_all=[tec1.E_th,tec2.E_th, tec3.E_th,  tec4.E_th]
Cout_all=[tec1.Cout,tec2.Cout, tec3.Cout,  tec4.Cout]
Qout_all=[tec1.Qout, tec2.Qout, tec3.Qout,  tec4.Qout]
Qchem_all=[tec1.chem,tec2.chem,tec3.chem, tec4.chem]
```
### 4.2. Formulate performance indicators
In Example 1, technical, economic and environmental indicators are used to evaluate the treatment chain. 
### 4.2.1. Technical indicators 
For instance, a technical indicator is the efficiency of the system in terms of **Water recovery (%)**. 

$$
\text{Water recovery} = \frac{\sum\limits_{i=1}^{n} Qw_i}{Qsw} \times 100
$$  


Where i is the number of technologies. 
```python
    # Calculate the toal quantity of water production 
Qw_tot=tec1.Qwater+tec2.Qwater+ tec3.Qwater+ tec4.Qwater

    # Calculate water recovery (%)
rec=Qw_tot/Qsw*100 #%
```
Another technical indicator is the **Specific energy consumption**.  

$$
\text{SEC} = \frac{\sum\limits_{i=1}^{n} E_{\text{el}}}{Qw_{\text{tot}}}
$$

```python
    # Energy performance: Calculate specific electrical energy consumption #kwh/kg of desalinated water 
SEC=sum(E_el_all)/Qw_tot
print("specific electrical consumption is " +str(SEC)+ " Kwh/kg of desalinated water")
```

### 4.2.2. Economic indicators 
For the economic analysis, the total amount of **Revenues** from selling products is used to evaluate the economic performance of the treatment chain. 

$$
\text{revenues} = \sum_{i=1}^{n} Q_{\text{product}_i} \times \text{Selling price of product}_i
$$

#### Set input parameters
First, the updated market prices of the recovered products need to be set. 

```python
    # Quantity of recovered products
prd=[Qw_tot, M_MgOH2_1, Q_b_out, Q_a_out]

    # Specify products 
prd_name= ["Water",  "Mg(OH)2", "NaOH", "HCl"]

    # Market prices
hcl_pr=5.78 #euro/l 1M solution 
w_pr=0.001 #euro/kg
mgoh2_pr=1.0 #euro/kg
naoh_pr=7.2 #euro/L 1M NaOH solution  
```

#### Revenue calculation
After the set of input parameters, the **Revenues** of the treatment chain are calculated: 
```python
    # Initialize lists
reve_t=0
reve_list=[]
    # Revenue calculation
for i in range(len(prd)):
    rev_calc=revenue(prd[i], prd_name[i])    
    rev_calc.rev(hr, w_pr, nacl_pr, mgoh2_pr,na2so4_pr, naoh_pr, hcl_pr)
    print("Revenues from selling product " + prd_name[i]+" are " + str(round(rev_calc.rev_prd,2))+" Euro/year")
    reve_t = reve_t+rev_calc.rev_prd
    reve_list.append(rev_calc.rev_prd)
```
**_Note that a detailed description of the economic model and more economic indicators can be found in ._** 

### 4.2.3. Environmental indicators 
A simple indicator to evaluate the environmental performance of the treatment chain is the **Carbon dioxide emissions**.  

$$
\text{Emissions} = \sum_{i=1}^{n} E_i \times \text{CO2}\text{el}  + \sum_{i=1}^{n} E_{\text{th}} \times \text{CO2}\text{st}
$$  

#### Set input parameters
First, the emissions factor from electrical and thermal energy consumption must be set based on the case study's location and fuel. 
```python
  # Emission factor for electricity consumption 
CO2_el=0.275 # kg/kwh
  # Emission factor for thermal energy consumption 
CO2_st=0.255 # kg/kwh
```
#### Carbon dioxide emissions
After the set of input parameters, the **Carbon dioxide emissions** of the treatment chain are calculated: 

```python
# Calculate Carbon dioxide emission 
emis=(sum(E_el_all)*CO2_el)+(sum(E_th_all)*CO2_st)
```
### 4.2.4. Export results to excel file 
After running the process and economic models, the calculation of indicators and results can be exported to Excel. 

Here is an example: 
```python
#%% Results to excel
  # Create dataframes 
dfel=pd.DataFrame(E_el_all ,tec_names)
dfeth=pd.DataFrame(E_th_all, tec_names)
dfind=(Qw_tot, abs_Qw, rec, Tot_el, Tot_th, emis_t, Q_w_in,  OPEX, CAPEX)
dfprod=(Qw_tot, M_MgOH2_1, M_CaOH2,  M_b_out, M_a_out)

dfprodn=("Water (kg/hr)",  "Mg(OH)2(kg/hr)", "Ca(OH)2(kg/hr)",  "1M NaOH (kg/hr)","1M HCl(kg/hr)")
ind=np.array(["Water production", "absolute water production", "water recovery","Total electrical consumption", "Total thermal energy consumption", "Carbon dioxide emission kg co2/year",
              "Water footprint","OPEX", "CAPEX"])
units=pd.DataFrame(["kg/hr", "kg/hr", "%","KWh", "KWh"," kg co2/year","kg/hr","€/year", "€"])
dfind_t=pd.DataFrame(data=dfind, index=ind)
dfprodn=pd.DataFrame(data=dfprod, index=dfprodn)

  # Write results in excel document
with pd.ExcelWriter('results_example.xlsx') as writer:
    dfel.to_excel(writer,sheet_name="example", startcol=0, startrow=14, header=True)
    dfeth.to_excel(writer,sheet_name="example", startcol=2, startrow=14, header=True)
    dfprodn.to_excel(writer,sheet_name="example", startcol=0, startrow=23, header=True)
    dfind_t.to_excel(writer,sheet_name="indicators")
    units.to_excel(writer,sheet_name="indicators",startcol=2, startrow=0, header=True)
```

Additionally, table with flowrates and concentrations can be developed (see 'exmple_1.py) and printed in Excel. 

### 4.3. Visualization 
For the visualization of the treatment chain, a **Sankey** diagram is created. For this diagram, the mass flow rates are needed. 

First, you need to import: 
```python
#sankey diagram 
import plotly.graph_objects as go
from plotly.offline import plot
```
Then you can create the figure as below: 
```python
#Figure 1: Sankey diagram for Example: Mass flow rates
fig = go.Figure(data=[go.Sankey(
    node = dict(
      pad = 15,
      thickness = 20,
      line = dict(color = "black", width = 0.5),
      label = ["Seawater", "Water", "NaCl",  "NF", "ED", "MF-PFR", "EDBM", "HCl", "NaOH", "Mg(OH)\u2082", "Ca(OH)\u2082", "Utilities", "Saline solution"],
      color = "blue"
    ),
    link = dict(
      source = [0,3,3,4,4,5,5,5,6,6,6,11,11], 
      target = [3,4,5,6,12,6,9,10,7,8,12,5,6],
      value = [Qsw, Qperm, Qconc,Mc, Md,M_mfpfr_out, M_MgOH2_1,M_CaOH2,Q_a_out, Q_b_out,Q_s_out, QNAOH+QHCl, Q_w_in ]
  ))])
color_for_nodes = ["lightsteelblue","darkcyan","maroon", "midnightblue", "midnightblue", "midnightblue", "maroon"]
fig.update_traces(node_color = color_for_nodes)
fig.update_layout(title_text="Sankey diagram for Example: Mass flow rates ", font_size=15)
fig.show()
plot(fig)
```


<figure>
  <img src="https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/assets/150446818/09d80c5b-3398-4618-8d5c-3b0d2db7ad07" alt="Image" style="width:600px;">
</figure>

**Figure 2**. Sankey diagram of example 1.
<br>

