# Tutorial for Economic model 

## 1. Introduction 
Welcome to our comprehensive tutorial on running economic models. The economic model (`economic_f.py`) can be used to perform economic assessments of the treatment chain or of individual units to identify economic hotspots. The economic assessment will help researchers, engineers, and decision-makers improve the system's performance and select the most economically feasible integration of technologies. 

In this tutorial, we provide step-by-step instructions on how to [use economic models](#running-economic-model) and [analyse the results](#results-evaluation) obtained using **Example 1** as case study (see **Figure 1**).  

<figure>
  <img src="https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-/assets/150446818/55cc6b6f-dde8-4b12-ae61-fa23665c288e" alt="Image" style="width:600px;">
</figure>

**Figure 1**. Process flow diagram of example 1.
<br>

## 2. Running Economic model 

**Table 1** gives an overview of the main inputs and outputs of economic model (`economic_f.py`). 

|  Input                                     | Output                                    |
|-------------------------------------------|-------------------------------------------|
| Selling price for products [€/ton] or [€/m<sup>3</sup>] | Operating cost (OPEX) [€/year]          |
| Prices for energy [€/KWh], input chemicals [€/m<sup>3</sup>], cooling water [€/m<sup>3</sup>] | Investment cost (CAPEX) [€]               |
| Operating hours, lifetime                 | Revenues from selling products [€/year] |
| Interest rate, Inflation rate             |                                         |
|Equipment cost [€]  |                                          |
| Assumptions on CAPEX and OPEX calculations |                                          |

### 2.1. Import economic model functions 
First, import the required functions from `economic_f.py`.  
```python
from desalsim.economic_f import econom
from desalsim.economic_f import revenue
```

### 2.2. Set constant values 

```python
    # Operating hours per year
hr=24*300 #hours/year
    # Plant's lifetime 
lf= 20 #years
    # Interest rate 
r=0.06 #(-)
    # Inflation rate 
inf=0.02 #(-)
```
Import 'constants.py' for more constant parameters like molecular weights (MW) etc.  

```python
import desalsim.constants
```

### 2.3. Set input data 
The market prices for utilities (electricity, chemicals, steam, cooling water etc) need to be set. Use updated market values. 

```python
    # Market prices
      # Electricity price 
el_pr=0.253 #euro/kwh
      # Steam price 
s_pr=0  #euro/kwh
      # Chemicals prices 
hcl_pr=5.78  #euro/l 1M solution 
naoh_pr=5.2  #euro/L 1M NaOH solution  
antisc_pr=0.002 #euro/ml
      # Cooling water price 
cw_pr=0.000118  #euro/kg
```
The market prices for focused products need to be set. Use updated market values. 
```python
    # Market prices
      # Chemicals prices 
hcl_pr=5.78  #euro/l 1M solution 
naoh_pr=5.2  #euro/L 1M NaOH solution
      # Industrial water 
w_pr=0.001 #euro/kg
      # Salts prices 
mgoh2_pr=1 # euro/kg
caoh2_pr=0.125 #euro/kg

```
### 2.4. Import results from technical model

Before importing the results from `example_1.py`, ensure that the script has been executed successfully. This can be done either by running the script directly or by making sure the necessary variables are defined and ready for import.

To import the results from `example_1.py`: 
```python
from example_1 import E_el_all, E_th_all, Qchem_all, Mhcl_need, Qnaoh_need, Q_w_in

```
> [!IMPORTANT]
> Ensure that `example_1.py` has been executed so that these variables are properly defined. If these variables do not exist or the script has not been run, you will encounter an error.

After importing, you can create lists with the results as shown below:  
```python
el_conc=E_el_all
s_conc=E_th_all
chem1_conc=Qchem_all
chem1_pr=[0.002,0,constants.naoh_pr_s,0]
chem2_conc=[0,0,Mhcl_need,0]
chem2_pr=[0.002,0,constants.hcl_pr_s,0]
wat_conc=[0,0, Qnaoh_need, Q_w_in]
cw_conc=[0,0, 0, 0]
```

### 2.5 Equipment cost 
You can set equipment cost if it is known. In case of preliminary study, where equipment cost is still unkown scale-up function (`scaleup.py`) can be used. 

$$
\frac{{\text{Cost of purchased equipment (Plant A)}}}{{\text{Cost of purchased equipment (Plant B)}}} = \left( \frac{{\text{Capacity of Plant A}}}{{\text{Capacity of Plant B}}} \right)^m
$$

Where _m_ is six-tenths factor rule (m=0.6 or m=0.8) depends on the nature of technology. 

In this example, equipment costs from pilot units are used as reference. The equipment costs are included in the 'constantc.py' so the user can call the function in case of scaleup and unkown equipment costs. 

```python
#Units: nf,  mfpfr,  edbm, ed
Mf_basic_sc=[2290.86, 595.03, 229.01, 253]
```
After setting the reference scenario, the equipment cost of the scaled-up unit can be calculated as below: 

```python
    # Import equipment cost for reference scenario from constants function
eq_c=[constants.eq_c_nf, constants.eq_c_ed, constants.eq_c_mfpfr, constants.eq_c_edbm]
    # Capacity of reference scenario
Mf_basic_sc=[constants.Mf_basic_sc[0],constants.Mf_basic_sc[8], constants.Mf_basic_sc[3], constants.Mf_basic_sc[6]]

    # Capacity of evaluated treatment chain 
Mf_sce=[Qsw,   Qed_in, Qin_mfpfr, Q_in_edbm]

    # Calculation of the new equipment cost
for i in range(len(eq_c)):
    if Mf_basic_sc[i]!=Mf_sce[i]:
        eq_c[i]=scaleup.scaleup_eq(eq_c[i],Mf_basic_sc[i],Mf_sce[i],tec_names[i])
```

### 2.6. Calculate Capital costs (CAPEX)
The CAPEX consists of fixed-capital investment and working capital, and the former one includes hardware costs, costs of buildings, process, and auxiliary, land, working capital and other indirect costs. 

**Table 2** gives an overview of the main assumptions made to calculate the CAPEX.
| CAPEX                             | 
|-----------------------------------|
| Installation: 25% of purchased equipment cost| 
| Buildings, process, and auxiliary: 20% of purchased equipment cost| 
| Land: 6% of purchased equipment cost  | 
| Indirect costs: 15% of direct cost                   | 
| Working capital: 20% of total investment cost  | 

The assumptions in **Table 2** are set. Note that they can change for different case studies. 
```python
    # Set assumptions for Capital investment cost 
inst_percent=0.25 # % of purchased equipment cost
buildings_percent=0.2 # of purchased equipment cost
land_percent=0.06 # of purchased equipment cost
indirect_c_percent=0.15 # of direct cost   
workinf_c_percent=0.2 #of total investment cost
```
Based on those assumptions the normalized factor is calculated and all the assumptions are assigned to the _capex_assumptions_ list.  
```python
capex_assumptions=[inst_percent,buildings_percent, land_percent, indirect_c_percent, workinf_c_percent]
```

Then the **CAPEX** is calculated: 
```python
     # Create list with capital cost for each unit. This will help us to identify the hotspots in capital investment.  
capex_list=[]
      # Initialize value
CAPEX=0
      # Calculate investment cost (CAPEX) 
for i in range(len(eq_c)):
    total_econom=econom(eq_c[i], el_conc[i], s_conc[i], chem1_conc[i], chem1_pr[i],chem2_conc[i], chem2_pr[i], cw_conc[i], wat_conc[i])
    total_econom.capex_calc(capex_assumptions)
    CAPEX=total_econom.t_capital_inv
    capex_list.append(total_econom.t_capital_inv)    
print("Total investment cost (CAPEX) of system is " + str(round(CAPEX))+ " Euro")    
```

#### 2.6.1 Calculate annual capital costs 
For the calculation of the annualized CAPEX, the amortization factor (α) is used: 

$$
\text{Annualized CAPEX} = \text{CAPEX} \times a
$$  

where a is the amortisation factor. 

Amortisation factor is calculated as: 

$$
a = \frac{r \cdot (1+r)^{lf}}{(1+r)^{lf} - 1}
$$  

where r is discount rate, lf is plant lifetime (year).

```python
    # Calculate amortisation factor
lf= 20 #years
r=0.06 # interest rate
a=(r*(1+r)**lf)/((1+r)**lf-1)
```
### 2.7. Calculate Operating costs (OPEX)
The OPEX refers to expenditure directly generated by manufacturing operation or connected to the equipment of a technical unit. Table 3 gives an overview of the costs that constitute OPEX (Peters, Timmerhaus and West, 2003). In this study, the utilities in this system are mainly energy, chemicals, and water costs. The calculation of yearly electrical (Cel) and thermal (Cth) energy costs follows equations: 

$$
C_{el} = Etot_{\text{el}} \cdot t_{\text{operation}} \cdot P_{\text{el}}
$$

$$
C_{th} = Etot_{\text{th}} \cdot t_{\text{operation}} \cdot P_{\text{st}}
$$

$$
C_{e} = C_{el} + C_{th}
$$  

Where:  
E<sub>el</sub> and E<sub>th</sub> are the total energy consumption per operating hour (in kWh/hr),   
t<sub>operation</sub> is the total operation time in one year (in hr),  
P<sub>el</sub> and P<sub>steam</sub> are the prices of electricity and steam, respectively (in €/kWh).  

The calculation of chemicals and water costs is similar to the energy cost, multiplying the amount of consumption every year by their price. The other costs are calculated in _opex_calc_ class (see 'economic.py').   

**Table 3** gives an overview of the main assumptions made to calculate the OPEX. 
| Annual OPEX                                    |
|------------------------------------------------|
| Maintenance: 3% of the fixed-capital investment            |
| Operating Supplies: 5% of maintenance |
| Operating Labor: 15% of annual OPEX                             |
| Direct supervisory and clerical labor: 15% of operating labor                         |
| Laboratory charges: 15% of operating labor                         |
| Patents and royalties: 3% of annual OPEX                          |
| Fixed charges: 5% of annual OPEX                                  |
| Plant overhead costs: 5% of annual OPEX                           |

```python
   # Set assumptions for Operating cost (OPEX)
main_c_percent=0.03  # % of the fixed-capital investment
oper_sup_c_percent=0.05  # % of maintenance 
oper_lab_c_percent=0.15    # % of annual OPEX
super_c_percent=0.15 # % of operating labor 
lab_c_percent=0.15 # % of operating labor 
pat_c_percent=0.03 # % of annual OPEX
fix_char_percent=0.05 # % of annual OPEX
over_c_percent=0.05 # %of annual OPEX
```
Based on those assumptions the normalized factor is calculated and all the assumptions are assigned to the _economic assumptions_ list.  
```python
    # Calculate normalized factor 
norm_factor =(1 -oper_lab_c_percent-oper_lab_c_percent*super_c_percent- oper_lab_c_percent*lab_c_percent-pat_c_percent-fix_char_percent-over_c_percent)

    # Create a list for economic assumptions
economic_assumptions=[main_c_percent,oper_sup_c_percent, oper_lab_c_percent, super_c_percent, lab_c_percent, pat_c_percent, 
                      fix_char_percent, over_c_percent, norm_factor ]
```
Then the **OPEX** is calculated: 
```python
  # Initialize values
OPEX=0

   # Create list with operating for each unit. This will help us to identify the hotspots in operating costs.  
opex_list=[]

  # Calculate OPEX 
for i in range(len(eq_c)):
    total_econom=econom(eq_c[i], el_conc[i], s_conc[i], chem1_conc[i], chem1_pr[i],chem2_conc[i], chem2_pr[i], cw_conc[i], wat_conc[i])
    total_econom.opex_calc(hr, el_pr, s_pr, cw_pr, w_pr, economic_assumptions)
    OPEX=total_econom.opex
    capex_list.append(total_econom.t_capital_inv)    
print("Total operating cost (OPEX) is "+str(OPEX)+ " Euro/year")   
```
### 2.8. Calculate Revenues 
The total amount of **Revenues** from selling products is used to evaluate the economic performance of the treatment chain. 

$$
\text{revenues} = \sum_{i=1}^{n} Q_{\text{product}_i} \times \text{Selling price of product}_i
$$

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

## 3. Results evaluation 
Running the economic model for **Example 1**, the following results are obtained. 

Total operating cost (OPEX) is 3714921 Euro/year  
Total investment cost (CAPEX) of the system is 6017026 Euro    
Annual capital investment cost of the system is 524592 Euro/year   

Revenues from selling product Water are 0.0 Euro/year  
Revenues from selling product Mg(OH)2 are 2891274 Euro/year  
Revenues from selling product NaOH are 3943229426 Euro/year  
Revenues from selling product HCl are 3274716868 Euro/year  
