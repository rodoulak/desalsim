#scenarios comparison 

#Import results 
import example_1 as sc1
import example_1 as sc2

#import 
import numpy as np
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt 
from pandas.plotting import register_matplotlib_converters
import matplotlib.ticker as ticker
register_matplotlib_converters()

hr=300*24 #hours   

#%%figures 

X = ['Sce 1', 'Sce 2']
X_axis = np.arange(len(X))

#Figure 1: Electrical consumption Vs thermal consumption 
Eel = [ sc1.sum_el_en, sc2.sum_el_en]
Eth = [sc1.sum_th_en,   sc2.sum_th_en]

Eel = [i * hr/1e6 for i in Eel]
Eth = [i * hr/1e6 for i in Eth]
plt.bar(X_axis - 0.2, Eel, 0.4, color="#00516a", label = 'Electrical (GWel)')
plt.bar(X_axis + 0.2, Eth, 0.4, color="sandybrown", label = 'Thermal (GWth)')
plt.xticks(X_axis, X)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel("Scenarios")
plt.ylabel("Eenergy consumption (GW)")
# plt.title("electrical energy vs thermal energy")
plt.legend()
plt.savefig('electricVSthermal.png')
plt.show()

#Figure 2: water production in kg/hr
Qw = [sc1.abs_Qw, sc2.abs_Qw]

Qw_y = [i * hr/1E6 for i in Qw]
# Qw_y=Qw*constants.hr/1000
plt.bar(X_axis - 0.4, Qw_y, 0.4, color="#00516a")  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xticks(X_axis, X)
plt.xlabel("Scenarios")
plt.ylabel("Water production in 1000 m\u00B3/year")
plt.savefig('water_production.png')
plt.show()

#Figure 3: Capex 
Cap = [sc1.CAPEX, sc2.CAPEX]
Cap = [i/1e6 for i in Cap]
plt.bar(X_axis - 0.4, Cap, 0.4, color="#00516a")  
plt.xticks(X_axis, X)
plt.xlabel("Scenarios")
plt.ylabel("Capex (M€)")
plt.savefig('capex.png')
plt.show()

#Figure 4: OPEX 
OPEX = [sc1.OPEX, sc2.OPEX]
OPEX = [i/1e6 for i in OPEX]
plt.bar(X_axis - 0.4, OPEX, 0.4, color="#00516a")    
plt.xticks(X_axis, X)
plt.xlabel("Scenarios")
plt.ylabel("OPEX (M€/year)")
plt.savefig('OPEX.png')
plt.show()

#Figure 5: Revenues 
Rev= [sc1.Rev,  sc2.Rev]
Rev = [i/1e6 for i in Rev]
plt.bar(X_axis - 0.4, Rev, 0.4, color="#00516a")  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xticks(X_axis, X)
plt.xlabel("Scenarios")
plt.ylabel("Revenues (M€/year)")

plt.savefig('revenues.png')
plt.show()

#Figure 12: CO2 emissions 
co2_em=[sc1.emis_t, sc2.emis_t]
co2_em = [i * hr/1e6 for i in co2_em]
plt.bar(X_axis - 0.4, co2_em, 0.4, color="#00516a")  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xticks(X_axis, X)
plt.xlabel("Scenarios")
plt.ylabel("CO\u2082 emissions (KTon/year)")
plt.savefig('CO2emissions.png')
plt.show()


#%%
#create dataframes 
ind=np.array(["Water production", "Total electrical consumption (GWh)", "Total thermal energy consumption (GWh)","Carbon dioxide emission (Kton co2/year) ",
              "OPEX (M€/year)", "CAPEX (M€)", "Resource efficiency (%)", "Production efficiency (%)", "Water footprint (m3/year)", 
              "Brine production (ton/year)", "levelized cost of water €/m\u00B3"])
sum_res_1=np.array([Qw_y[1],Eel[1], Eth[1],co2_em[1], OPEX[1], Cap[1]])
df1=pd.DataFrame(sum_res_1, ind)

sum_res_2=np.array([Qw_y[2],Eel[2], Eth[2],co2_em[2], OPEX[2], Cap[2]])
df2=pd.DataFrame(sum_res_2, ind)

#all
df=pd.DataFrame({'sce 1': [Qw_y[0],Eel[0], Eth[0],co2_em[0], OPEX[0], Cap[0]], 
                 'sce 2': [Qw_y[1],Eel[1], Eth[1],co2_em[1], OPEX[1], Cap[1]],
                 })

sc1_dfel=pd.DataFrame(sc1.E_el_all , sc1.tec_names)
sc1_dfeth=pd.DataFrame(sc1.E_th_all, sc1.tec_names)
sc2_dfel=pd.DataFrame(sc2.E_el_all , sc2.tec_names)
sc2_dfeth=pd.DataFrame(sc2.E_th_all, sc2.tec_names)

#%%export to excel
with pd.ExcelWriter('results_scenario.xlsx') as writer:
    df1.to_excel(writer,sheet_name="example 1")
    df2.to_excel(writer,sheet_name="example 2")
    df.to_excel(writer,sheet_name="all scenarios")
    
    sc1.sum_table_C.to_excel(writer,sheet_name="example 1", startcol=0, startrow=20, header=False)
    sc1.sum_table_f.to_excel(writer,sheet_name="example 1", startcol=0, startrow=18, header=True, index=True )
    sc1.sum_table_d.to_excel(writer,sheet_name="example 1", startcol=0, startrow=29, header=True)

    sc2.sum_table_C.to_excel(writer,sheet_name="example 2", startcol=0, startrow=20, header=False)
    sc2.sum_table_f.to_excel(writer,sheet_name="example 2", startcol=0, startrow=18, header=True)
    sc2.sum_table_d.to_excel(writer,sheet_name="example 2", startcol=0, startrow=30, header=False)
    
    sc2_dfel.to_excel(writer,sheet_name="example 2", startcol=0, startrow=33, header=True)
    sc2_dfeth.to_excel(writer,sheet_name="example 2", startcol=2, startrow=33, header=True)
    
    sc1_dfel.to_excel(writer,sheet_name="example 1", startcol=0, startrow=33, header=True)
    sc1_dfeth.to_excel(writer,sheet_name="example 1", startcol=2, startrow=33, header=True)


    