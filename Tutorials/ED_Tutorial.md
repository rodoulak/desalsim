# Tutotrial: Electrodialysis (ED) unit

ED  is a technique that uses ion exchange membranes and electricity to extract ionic compounds from aqueous solutions, with an electrical field serving as the driving force for solute separation. In this work, ED aims to purify NF permeate. 
![ed](https://github.com/user-attachments/assets/386c055d-9935-41fa-b15c-e790fc169f9a)


In **desalsim** package, the ED unit is used to model the operation of a Electrodialysis technology. Upon simulation, it calculates the flow rate of the concentrate and dilute streams, their ion concentration, and the electricity requirements of the unit.
The ED function consists of one class: [ElectrodialysisCalc](#use-electrodialysiscalc-class ) that constsis of six methods: `Ts_cp`, `w_cp`, `Ls_cpv`, `Lw_cp`, `p_osmo`, `dC` [see Section 2.3](#23-use-ts_cp-w_cp-ls_cpv-lw_cp-p_osmo-dc-methods).  
In this tutorial, we will focus on how to use the the class and their methods. 

**Table 1** gives an overview of the main inputs and outputs for each process unit of Electrodialysis. 
| Process                                   | Input                                       | Output                                                |
|-------------------------------------------|---------------------------------------------|-------------------------------------------------------|
| Electrodialysis (ED)                         | Feed flow rate [m³/h]                       | Flow rate of diluted stream [m³/h] and composition [g/L]|
|                                           | Ion concentration [g/L]                     | Flow rate of concentrate stream [m³/h] and composition [g/L]        |
|                                           | Current density [A/m²]                      | Electricity requirements [kWhel]                     |

The mathematical description of Electrodialysis  is given in [Mathematical description](https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/tree/main/paper/Mathematical_description.pdf), see Section A.8. 

## 1. Getting started 
### 1.1. Import class
```python
import desalsim
```
Then import the class:  
```python
from desalsim.ed_unit_f import ElectrodialysisCalc
```
Additionally, function for calculating density (`density_calc.py`) or constants ('comparison.py') where user can add constant values like MW, prices etc, need to be imported.
```python
from desalsim.density_calc import density_calc 
from desalsim import constants
from desalsim import scaleup
```
### 1.2. Define feed characteristics
You can initialize the feed solution by setting the flow rate, specifying the focus components and their concentration. 
```python
# Input conditions
Sc_i = 43.39  # Salinity at concentrate inlet (g/kg)
Sd_i=Sc_i # Salinity at diluate inlet (g/kg)

  # Feed concentration
components = ['Na', 'Cl', 'K', 'Mg', 'Ca', 'SO4', 'HCO3', 'H', 'OH']
Csw= [17.17, 25.47, 0.57, 0.04, 0.03, 0.10]

#Feed flow rate L/h
Qed_in=1000

# Temperature 
T=20+273.15 #K
```
> [!NOTE]
> Note that if you want to add more components, you need to update the components list and include the concentration of the new component in the _Csw_. 

You can calculate the density of the feed solution and the water quantity in inflow:
```python
d_in_ed=(density_calc(T-273, Sc_i)/1000)
Qed_in=1000/d_in_ed

# Calculate the flowrates for dilute and concentrate streams 
Qed_in_c=Qed_in/17/Ncp
Qed_in_d=Qed_in*16/17/Ncp
```
> [!NOTE]
> In this case we assumed that the concentrate stream is the 1/17 of the total feed flow rate and the dilute stream is the 16/17 of the total feed flow rate.

### 1.3. Set operating assumptions  
You need to set operating assumptions such as the electrical current density.   
```python
# Assumptions:
# The electrical current desnity
Ij=400 #Am2
# Set Number of computational cells per cell-pair
N = 50
# Number of identical parallel cell-pairs
Ncp=1
# Set active area of cell-pair (m^2)
A = 1.1
# Set Membrane efficiency
Mem_eff=0.64
# Set applied voltage (V)
Vcp = 8
# Set Voltage across the electrodes (V)
Vel=2.1

# Effective cell-pair area (m^2)
Acp =scaleup.scaleup(24, 1000, Qed_in) 
Acp_tot=Acp
```

You need to set the salinity at dilute outlet concentrate outlet in g/kg. 
```python
Sc_o = 200  # Salinity at concentrate outlet (g/kg)
Sd_o=20 # Salinity at doluate outlet (g/kg)
```

Finally, you need to set assumptions related to pumping like pressure drop (_dp_) and pump efficiency (_npump_). 
```python
npump=0.8 #pump efficiency (units: -)
dp=1 #pressure drop (units: bar)
```

### 1.4. Set constants 
You need to set constant parameters: 
```python
R = 0.002  # Resistance of rinse stream (ohm)
Rp = 0.015  # Resistance of polarization (ohm)
A = 1.1  # Active area of cell-pair (m^2)
F = 96485.3329  # Faraday constant (C/mol)
rho_w=1000 # water desnity kg/m3

D=1.61e-9 #Diffusion coefficient (m^2/s)
tcu=0.5
veloc=8.9e-7 #m2/s
h=0.5 #mm
Sh=18
```

### 1.5. Initializations
First, you need to initialize the parameters:
```python
# Initializations
Sc = np.zeros(N)
Sd = np.zeros(N)
Ns_c = np.zeros(N)
Ns_d = np.zeros(N)
Nw_c = np.zeros(N)
Nw_d = np.zeros(N)
Js = np.zeros(N)
Jw = np.zeros(N)
Mw_in_d_l=np.zeros(N)
Ms_d=np.zeros(N)
Mw_d=np.zeros(N)
M_d=np.zeros(N)
M_c=np.zeros(N)
Q_c=np.zeros(N)
Q_d=np.zeros(N)
```
Then set the initial values for the concentrate stream:
```python
# Set initial values salt stream 
Sc[0] = Sc_i
Ns_c[0]=Qed_in_c*d_in_ed*Sc[0]/MWs #mol/hr
Ms_in_c=Qed_in_c*d_in_ed*Sc_i/1000 #kg salr/hr
Mw_in_c=Qed_in_c*d_in_ed-Ms_in_c #kg water/hr
Nw_c[0]=Mw_in_c*1000/MWw #mol/hr
M_c[0]=Mw_in_c+Ms_in_c
Q_c[0]=Qed_in_c
```
```python
# Set initial values diluate stream 
Csw= [17.17, 25.47, 0.57, 0.04, 0.03, 0.10]
Sd[0] = Sc_i #g/kg 
Ns_d[0]=Qed_in_d*d_in_ed*Sd[0]/(MWs) #mol/s
Ms_in_d=Qed_in_d*d_in_ed*Sd[0]/1000 #kg salr/hr
Mw_in_d=Qed_in_d*d_in_ed-Ms_in_d #kg water/hr
Nw_d[0]=Mw_in_d*1000/MWw #mol/hr
Mw_in_d_l[0]=Mw_in_d
Ms_d[0]=Ms_in_d
Mw_d[0]=Mw_in_d
M_d[0]=Mw_in_d+Ms_in_d
Q_d[0]=Qed_in_d
```

Finally, initialize the total cell-pair area
```python
# Initialize the Acp_tot array
Acp_tot_j=Acp_tot/N
```
After setting all the required inputs and initialize the values, then you can create the functions' objectives. 

## 2. Use ElectrodialysisCalc class   
`ElectrodialysisCalc` is a class used to represent mass and energy balance for ED Unit. In particular, it calculateshe flowrate in each channel, the outlet concentration in each channel, the external Voltage and power needed. 

### 2.1. Overview 
The following attributes are available within the `ElectrodialysisCalc` class:  
- MWs :Molecular weight of NaCl (g/mol)
- MWw : Molecular weight of water (g/mol)
- R :Resistance of rinse stream (ohm)
- Rp : Resistance of polarization (ohm)
- A : Active area of cell-pair (m^2)
- F : Faraday constant (C/mol)
- T : Temperature in Kelvin
- dp : Parameter dp
- npump : Pump efficiency
- rho_w : Density of water (kg/m^3)
- D : Diffusion coefficient (m^2/s)
- tcu: ntcu parameter
- h : Height (mm)
- Sh : Sh parameter
- Mem_eff : Membrane efficiency
- Ncp : Number of cell-pairs
- Qed_in :  Inlet flow rate (l/h)
- Qed_in_c : Concentrate inlet flow rate (l/h)
- Qed_in_d : Dilute inlet flow rate (l/h)


The ElectrodialysisCalc class provides the following methods:
```python
 # Calculates the transport number for salt in concentrate compartment"
Ts_cp(S)
# Calculates the transport number for water in concentrate compartment"
w_cp(Sc,Sd)
# Permeability for salt in concentrate compartment"
Ls_cp(Sc, Sd)
# Permeability for water in concentrate compartment
Lw_cp(S)
# calculation for Osmotic pressure in concentrate compartment
p_osmo(S, T, MWs)
# Calculate the change in concentration
dC(Ts_cp, tcu, D, Ij,  h, Sh )
```

### 2.2. Create ElectrodialysisCalc objects
ElectrodialysisCalc has no inputs.  
 
```python
#Create an instance of the ElectrodialysisCalc class 
ed_em = ElectrodialysisCalc()
```

### 2.3. Use `Ts_cp`, `w_cp`, `Ls_cpv`, `Lw_cp`, `p_osmo`, `dC`  methods 
The ED system is modeled by adapting a model developed by [Nayar et al](https://www.sciencedirect.com/science/article/pii/S0011916418312761), keeping both the con- centrate and diluate channels fully continuous, with the salinities of both channels varying along the length of the ED stack. The following code is used to simulate the ED unit using the `Ts_cp`, `w_cp`, `Ls_cpv`, `Lw_cp`, `p_osmo`, `dC`   methods. 

```python
# Iterate over cells
for j in range(1, N):
    #Calculate salinity change 
    concentration_diff = Sc[j - 1] - Sd[j - 1]
    Sc[j] = Sc[j - 1] + (Sc_o - Sc_i) / (N - 1)
    Sd[j] = Sd[j - 1] + (Sd_o - Sd_i) / (N - 1)
    
    #Calculate net salt flux 
    Js[j] = (ElectrodialysisCalc.Ts_cp(Sd[j - 1]) * Ij / F - 
             (ElectrodialysisCalc.Ls_cp(Sc[j - 1], Sd[j - 1])) * concentration_diff)
    #calculate net water flux 
    Jw[j] = (ElectrodialysisCalc.Tw_cp(Sc[j - 1], Sd[j - 1]) * Ij / F +
             ElectrodialysisCalc.Lw_cp(Sc[j - 1]) * (ElectrodialysisCalc.p_osmo(Sc[j - 1], T, MWs) - 
                                                      ElectrodialysisCalc.p_osmo(Sd[j - 1], T, MWs)))
    
    #Calculate he total concentrate and dilute molar flow rates
    Ns_c[j] = Ns_c[j - 1] + Acp_tot_j * Js[j]
    Ns_d[j] = Ns_d[j - 1] - Acp_tot_j * Js[j]

    Nw_c[j] = Nw_c[j - 1] + Acp_tot_j * Jw[j]
    Nw_d[j] = Nw_d[j - 1] - Acp_tot_j * Jw[j]

    # Update the flow rates of the concentrate and dilute streams
    Q_c[j] = Nw_c[j] * MWw / (rho_w * (1 - Sc[j] / 1000))
    Q_d[j] = Nw_d[j] * MWw / (rho_w * (1 - Sd[j] / 1000))

```

 
### 2.3.1. Assigned the results to output parameters 
You can assigned the results to output parameters: 
```python
Cc_na_f=Sc[N-1]/MWs*constants.MW_Na
Cc_cl_f=Sc[N-1]/MWs*constants.MW_cl
Sc_out=[Cc_na_f, Cc_cl_f]

```

### 2.4. Calculate the concentrate stream flow ate 

```python
#Calculate the concnetrate stream flow rate 
Mc=(Ns_c[N-1]*MWs/1000+Nw_c[N-1]*MWw/1000) #(kg/hr)

dc_out=density_calc(T-273, Sc[N-1])/1000 #(kg/l)
Qc=Mc/dc_out #concnetrate stream volume flow rate (l/hr)

i=2
for i in range(2,len(Csw)):
    Sc_out.append(Csw[i]*Qed_in_c/Qc) #The total effluent concentration concentrate stream
```
### 2.5. Calculate the dilute stream flow ate 

```python
#Calculations for diluate stream 
Md=(Ns_d[N-1]*MWs/1000+Nw_d[N-1]*MWw/1000) #mass flow rate (kg/hr)


Sd_f=Sd[N-1]
Cd_na_f=Sd_f/MWs*constants.MW_Na
Cd_cl_f=Sd_f/MWs*constants.MW_cl
dd_out=density_calc(T-273, Sd[N-1])/1000 #density of diluate stream
Qd=Md/dd_out #diluate stream volume flow rate (l/hr)


Sd_out=[Cd_na_f, Cd_cl_f]
# The total effluent concentration dilute
for i in range(2,len(Csw)):
    Sd_out.append(Csw[i]*Qed_in/Qd)
```


### 2.6. Calculate energy consumption 
You can calculate the total energy requirements for the ED unit. For this, you can use the voltage applied across an ED cell-pair (_Vcp_), voltage across the electrodes (_Vel_) and the energy for pumping (_Ppump_ed_). 

```python
# Energy consumption 
#power required 
Ws = 0
for j in range(N):
     Ws += Ij * Acp_tot_j * (Ncp * Vcp + Vel)
print("Power required is "+str(round(Ws/1000,2))+"KW")

#Calculate energy consumption for pumping 
Ppump_ed=(Qed_in_d*1+Qed_in_c*1+Qc*2+Qd*1)/1000/3600*1e5/npump
Eel_t_ed=Ws/1000+Ppump_ed/1000

# Specific enrgy consumption
sec_ed=Eel_t_ed/(Qed_in/1000)
```

### 2.7. Print results 
You can print results from the calculations 
```python
print("Mass flowrate concentrate stream is "+str(round(Mc,2))+ " kg/hr")
print("Volume flowrate concentrate stream is "+str(round(Qc,2))+" l/hr")
print("The total effluent concentration concentrate stream  is " + str(round(Sc[N-1],2))+"g/kg")
print("-----------------------------------------")
```
Mass flowrate concentrate stream is 78.5 kg/hr
Volume flowrate concentrate stream is 67.74 l/hr
The total effluent concentration concentrate stream  is 200.0g/kg

```python
print("Mass flowrate of diluate stream is "+str(round(Md,2))+" kg/hr")
print("volume flowrate diluate stream is "+str(round(Qd,2))+" l/hr")
print("The total effluent concentration dilute is " + str(round(Sd[N-1],2))+"g/kg")
print("-----------------------------------------")
```
Mass flowrate of diluate stream is 921.5 kg/hr
volume flowrate diluate stream is 909.36 l/hr
The total effluent concentration dilute is 20.0g/kg

```python
#solid mass balance
bal=Qed_in-Md-Mc
bal=(Qed_in*sum(Csw) -Md*(sum(Sd_out))-Mc*Sc_o)/1000
print("Mass balance difference is "+str(round(bal,2)))
error_perc=abs(bal)/(Qed_in*sum(Csw))*100
print("Balance error percentage is "+str(round(error_perc,2))+"%")
print("-----------------------------------------")
```
Mass balance difference is 7.21
Balance error percentage is 0.02%

```python
# Energy consumption 
print("Power required is "+str(round(Ws/1000,2))+"KW")
print("Total energy consumption is "+str(round(Eel_t_ed,2))+"KW")
print("Specific energy consumption of Electrodialysis (ED) is "+str(round(sec_ed,2))+"KW/m3 feed")

```
Power required is 95.19KW
Total energy consumption is 95.26KW
Specific energy consumption of Electrodialysis (ED) is 98.23KW/m3 feed
