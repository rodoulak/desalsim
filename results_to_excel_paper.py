#scenarios to excel 
import numpy as np
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt 
import sc2_water_mining as sc2
import sc1_swro_as_feed as sc1
import constants
import sc3_water_recovery as sc3
import sc4_swro_ED as sc4
#import scenario_4_nf_med_mfpfr_edbm_evaporaponds
#import scenario_3_nf_med_thcryst_mfpfr_edbm
#import sc1_swro_as_feed_nf_ro_med
from pandas.plotting import register_matplotlib_converters
import matplotlib.ticker as ticker
register_matplotlib_converters()
# %matplotlib inline
hr=300*24 #hours   
#%%

# df=pd.DataFrame({'sce 1': [sc2_water_mining.Qw_tot,sc2_water_mining.SEC,sc2_water_mining.Q_br_prod,sc2_water_mining.emis_t,sc2_water_mining.Chem_cons, sc2_water_mining.Min_waste, sc2_water_mining.OPEX_sc1, sc2_water_mining.CAPEX_sc1, sc2_water_mining.Pr_c_ext,sc2_water_mining.Pr_c,sc2_water_mining.Ec_m, sc2_water_mining.sum_el_en, sc2_water_mining.sum_th_en],
#                  'sce 2a':[sc1_swro_as_feed.Qw_tot,sc1_swro_as_feed.SEC,sc1_swro_as_feed.Q_br_prod,sc1_swro_as_feed.emis_t,sc1_swro_as_feed.Chem_cons, sc1_swro_as_feed.Min_waste, sc1_swro_as_feed.OPEX_sc2, sc1_swro_as_feed.CAPEX_sc2, sc1_swro_as_feed.Pr_c_ext,sc1_swro_as_feed.Pr_c,sc1_swro_as_feed.Ec_m, sc1_swro_as_feed.sum_el_en, sc1_swro_as_feed.sum_th_en], 
#                  'sce 2b': [sc1_swro_as_feed_nf_ro_med.Qw_tot,sc1_swro_as_feed_nf_ro_med.SEC,sc1_swro_as_feed_nf_ro_med.Q_br_prod,sc1_swro_as_feed_nf_ro_med.emis_t,sc1_swro_as_feed_nf_ro_med.Chem_cons, sc1_swro_as_feed_nf_ro_med.Min_waste, sc1_swro_as_feed_nf_ro_med.OPEX_sc2, sc1_swro_as_feed_nf_ro_med.CAPEX_sc2, sc1_swro_as_feed_nf_ro_med.Pr_c_ext,sc1_swro_as_feed_nf_ro_med.Pr_c,sc1_swro_as_feed_nf_ro_med.Ec_m, sc1_swro_as_feed_nf_ro_med.sum_el_en, sc1_swro_as_feed_nf_ro_med.sum_th_en], 
#                  'sce 3': [sc3_water_recovery.Qw_tot,sc3_water_recovery.SEC,sc3_water_recovery.Q_br_prod,sc3_water_recovery.emis_t,sc3_water_recovery.Chem_cons, sc3_water_recovery.Min_waste, sc3_water_recovery.OPEX_sc3, sc3_water_recovery.CAPEX_sc3, sc3_water_recovery.Pr_c_ext,sc3_water_recovery.Pr_c,sc3_water_recovery.Ec_m, sc3_water_recovery.sum_el_en, sc3_water_recovery.sum_th_en], 
#                  'sce 4': [scenario_3_nf_med_thcryst_mfpfr_edbm.Qw_tot,scenario_3_nf_med_thcryst_mfpfr_edbm.SEC,scenario_3_nf_med_thcryst_mfpfr_edbm.Q_br_prod,scenario_3_nf_med_thcryst_mfpfr_edbm.emis_t,scenario_3_nf_med_thcryst_mfpfr_edbm.Chem_cons, scenario_3_nf_med_thcryst_mfpfr_edbm.Min_waste, scenario_3_nf_med_thcryst_mfpfr_edbm.OPEX_sc4, scenario_3_nf_med_thcryst_mfpfr_edbm.CAPEX_sc4,scenario_3_nf_med_thcryst_mfpfr_edbm.Pr_c_ext,scenario_3_nf_med_thcryst_mfpfr_edbm.Pr_c,scenario_3_nf_med_thcryst_mfpfr_edbm.Ec_m, scenario_3_nf_med_thcryst_mfpfr_edbm.sum_el_en, scenario_3_nf_med_thcryst_mfpfr_edbm.sum_th_en ],
#                  'sce 5': [scenario_4_nf_med_mfpfr_edbm_evaporaponds.Qw_tot,scenario_4_nf_med_mfpfr_edbm_evaporaponds.SEC,scenario_4_nf_med_mfpfr_edbm_evaporaponds.Q_br_prod,scenario_4_nf_med_mfpfr_edbm_evaporaponds.emis_t,scenario_4_nf_med_mfpfr_edbm_evaporaponds.Chem_cons, scenario_4_nf_med_mfpfr_edbm_evaporaponds.Min_waste, scenario_4_nf_med_mfpfr_edbm_evaporaponds.OPEX_sc5, scenario_4_nf_med_mfpfr_edbm_evaporaponds.CAPEX_sc5,scenario_4_nf_med_mfpfr_edbm_evaporaponds.Pr_c_ext,scenario_4_nf_med_mfpfr_edbm_evaporaponds.Pr_c,scenario_4_nf_med_mfpfr_edbm_evaporaponds.Ec_m, scenario_4_nf_med_mfpfr_edbm_evaporaponds.sum_el_en, scenario_4_nf_med_mfpfr_edbm_evaporaponds.sum_th_en ]
# })
# print(df)

#%%figures 
#
X = ['Sce 1', 'Sce 2', 'Sce 3', 'Sce 4']#,'Sce 4', 'Sce 5']
# X = ['sce 1a','sce 1b', 'sce 2', 'sce 3','sce 4', 'sce 5']
X_axis = np.arange(len(X))

#Figure 1: Electrical consumption Vs thermal consumption 
Eel = [ sc3.sum_el_en, sc2.sum_el_en, sc1.sum_el_en,  sc4.sum_el_en]#, scenario_3_nf_med_thcryst_mfpfr_edbm.sum_el_en, scenario_4_nf_med_mfpfr_edbm_evaporaponds.sum_el_en ]
Eth = [sc3.sum_th_en,   sc2.sum_th_en, sc1.sum_th_en, sc4.sum_th_en]#, scenario_3_nf_med_thcryst_mfpfr_edbm.sum_th_en, scenario_4_nf_med_mfpfr_edbm_evaporaponds.sum_th_en]

# Eel = [sc1_swro_as_feed.sum_el_en, sc1_swro_as_feed_nf_ro_med.sum_el_en, sc2_water_mining.sum_el_en, sc3_water_recovery.sum_el_en, scenario_3_nf_med_thcryst_mfpfr_edbm.sum_el_en, scenario_4_nf_med_mfpfr_edbm_evaporaponds.sum_el_en ]
# Eth = [sc1_swro_as_feed.sum_th_en, sc1_swro_as_feed_nf_ro_med.sum_th_en, sc2_water_mining.sum_th_en,sc3_water_recovery.sum_th_en, scenario_3_nf_med_thcryst_mfpfr_edbm.sum_th_en, scenario_4_nf_med_mfpfr_edbm_evaporaponds.sum_th_en]

Eel = [i * constants.hr/1e6 for i in Eel]
Eth = [i * constants.hr/1e6 for i in Eth]
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
Qw = [sc3.abs_Qw, sc2.abs_Qw, sc1.Qw_tot, sc4.Qw_tot]#, scenario_3_nf_med_thcryst_mfpfr_edbm.Qw_tot, scenario_4_nf_med_mfpfr_edbm_evaporaponds.Qw_tot]
#Qw = [sc1_swro_as_feed.abs_Qw, sc1_swro_as_feed_nf_ro_med.abs_Qw, sc2_water_mining.abs_Qw, sc3_water_recovery.Qw_tot, scenario_3_nf_med_thcryst_mfpfr_edbm.Qw_tot, scenario_4_nf_med_mfpfr_edbm_evaporaponds.Qw_tot]

Qw_y = [i * constants.hr/1E6 for i in Qw]
# Qw_y=Qw*constants.hr/1000
plt.bar(X_axis - 0.4, Qw_y, 0.4, color="#00516a")  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xticks(X_axis, X)
plt.xlabel("Scenarios")
plt.ylabel("Water production in 1000 m\u00B3/year")
# plt.title("water production")
# plt.legend()
plt.savefig('water_production.png')
plt.show()

#Figure 3: Capex 
Cap = [sc3.CAPEX, sc2.CAPEX, sc1.CAPEX, sc4.CAPEX]#, scenario_3_nf_med_thcryst_mfpfr_edbm.CAPEX_sc4, scenario_4_nf_med_mfpfr_edbm_evaporaponds.CAPEX_sc5]
#Cap = [sc1_swro_as_feed.CAPEX_sc2, sc1_swro_as_feed_nf_ro_med.CAPEX_sc2, sc2_water_mining.CAPEX_sc1, sc3_water_recovery.CAPEX_sc3, scenario_3_nf_med_thcryst_mfpfr_edbm.CAPEX_sc4, scenario_4_nf_med_mfpfr_edbm_evaporaponds.CAPEX_sc5]
Cap = [i/1e6 for i in Cap]
plt.bar(X_axis - 0.4, Cap, 0.4, color="#00516a")  
plt.xticks(X_axis, X)
plt.xlabel("Scenarios")
plt.ylabel("Capex (M€)")
# plt.title("capex")
# plt.legend()
plt.savefig('capex.png')
plt.show()

#Figure 4: OPEX 
OPEX = [sc3.OPEX, sc2.OPEX, sc1.OPEX, sc4.OPEX]#, scenario_3_nf_med_thcryst_mfpfr_edbm.OPEX_sc4, scenario_4_nf_med_mfpfr_edbm_evaporaponds.OPEX_sc5]
#OPEX = [sc1_swro_as_feed.OPEX_sc2, sc1_swro_as_feed_nf_ro_med.OPEX_sc2, sc2_water_mining.OPEX_sc1, sc3_water_recovery.OPEX_sc3, scenario_3_nf_med_thcryst_mfpfr_edbm.OPEX_sc4, scenario_4_nf_med_mfpfr_edbm_evaporaponds.OPEX_sc5]
OPEX = [i/1e6 for i in OPEX]
plt.bar(X_axis - 0.4, OPEX, 0.4, color="#00516a")    
plt.xticks(X_axis, X)
plt.xlabel("Scenarios")

plt.ylabel("OPEX (M€/year)")
# plt.title("OPEX")
# plt.legend()
plt.savefig('OPEX.png')
plt.show()

#Figure 5: Revenues 
Rev= [sc3.Rev,  sc2.Rev, sc1.Rev, sc4.Rev]#, scenario_3_nf_med_thcryst_mfpfr_edbm.Rev, scenario_4_nf_med_mfpfr_edbm_evaporaponds.Rev]
#Rev= [sc1_swro_as_feed.Rev, sc1_swro_as_feed_nf_ro_med.Rev, sc2_water_mining.Rev, sc3_water_recovery.Rev, scenario_3_nf_med_thcryst_mfpfr_edbm.Rev, scenario_4_nf_med_mfpfr_edbm_evaporaponds.Rev]

Rev = [i/1e6 for i in Rev]
plt.bar(X_axis - 0.4, Rev, 0.4, color="#00516a")  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xticks(X_axis, X)
plt.xlabel("Scenarios")
plt.ylabel("Revenues (M€/year)")
# plt.title("REVENUES")
# plt.legend()
plt.savefig('revenues.png')
plt.show()

#Figure 6: Annual cost (OPEX +annualize capex) VS Revenue
operc_c_t=[ sc3.oper_c_t,  sc2.oper_c_t, sc1.oper_c_t, sc4.oper_c_t]#, scenario_3_nf_med_thcryst_mfpfr_edbm.oper_c_t, scenario_4_nf_med_mfpfr_edbm_evaporaponds.oper_c_t]
#operc_c_t=[ sc1_swro_as_feed.oper_c_t, sc1_swro_as_feed_nf_ro_med.oper_c_t,  sc2_water_mining.oper_c_t, sc3_water_recovery.oper_c_t, scenario_3_nf_med_thcryst_mfpfr_edbm.oper_c_t, scenario_4_nf_med_mfpfr_edbm_evaporaponds.oper_c_t]
  
plt.bar(X_axis - 0.2, Rev, 0.4, color="#00516a", label = 'Revenue')
plt.bar(X_axis + 0.2, OPEX, 0.4, color="sandybrown", label = 'Production cost')
plt.xticks(X_axis, X)
plt.xlabel("Scenarios")
plt.ylabel("M€/year")
# plt.title("Revenue vs Annual cost")
plt.legend(loc='upper right')
plt.savefig('revenueVScost.png')
plt.show()

#for products 
#Figure 7: mgoh2 production 
M_magoh2=[ 0,sc2.M_MgOH2, sc1.M_MgOH2_1, sc4.M_MgOH2_1]#, scenario_3_nf_med_thcryst_mfpfr_edbm.M_MgOH2_1, scenario_4_nf_med_mfpfr_edbm_evaporaponds.M_MgOH2_1]
#M_magoh2=[sc1_swro_as_feed.M_MgOH2_1, sc1_swro_as_feed_nf_ro_med.M_MgOH2_1, sc2_water_mining.M_MgOH2, 0, scenario_3_nf_med_thcryst_mfpfr_edbm.M_MgOH2_1, scenario_4_nf_med_mfpfr_edbm_evaporaponds.M_MgOH2_1]
plt.bar(X_axis - 0.4, M_magoh2, 0.4, color="#00516a")  
  
plt.xticks(X_axis, X)
plt.xlabel("Scenarios")
plt.ylabel("Mg(OH)\u2082 Production (kg/year)")
# plt.title("M_magoh2 Production")
# plt.legend()
plt.savefig('M_magoh2Production.png')
plt.show()

#Figure 8: nacl production 
M_nacl=[sc3.M_Nacl, sc2.M_Nacl, sc1.M_Nacl,0]#, scenario_3_nf_med_thcryst_mfpfr_edbm.M_Nacl, scenario_4_nf_med_mfpfr_edbm_evaporaponds.Mnacl]
#M_nacl=[sc1_swro_as_feed.M_Nacl, sc1_swro_as_feed_nf_ro_med.M_Nacl, sc2_water_mining.M_Nacl, sc3_water_recovery.M_Nacl, scenario_3_nf_med_thcryst_mfpfr_edbm.M_Nacl, scenario_4_nf_med_mfpfr_edbm_evaporaponds.Mnacl]

plt.bar(X_axis - 0.4, M_nacl, 0.4, color="#00516a")  
plt.xticks(X_axis, X)
plt.xlabel("Scenarios")
plt.ylabel("NaCl Production (kg/year)")
# plt.title("NaCl Production")
# plt.legend()
plt.savefig('NaClProduction.png')
plt.show()

#Figure 9: na2so4 production 

M_na2so4=[0, sc2.M_Na2so4, sc1.M_Na2so4,0]
#M_na2so4=[sc1_swro_as_feed.M_Na2so4, sc1_swro_as_feed_nf_ro_med.M_Na2so4, sc2_water_mining.M_Na2so4, 0, 0, 0]

plt.bar(X_axis - 0.4, M_na2so4, 0.4, color="#00516a")  
plt.xticks(X_axis, X)
plt.xlabel("Scenarios")
plt.ylabel("Na\u2082SO\u2084 Production (kg/hr)")
# plt.title("na2so4 Production")
# plt.legend()
plt.savefig('na2so4Production.png')
plt.show()

#Figure 10: naoh production 
M_naoh=[0,  sc2.Qnaoh_s, sc1.Qnaoh_s, sc4.Qnaoh_s]#, scenario_3_nf_med_thcryst_mfpfr_edbm.Qnaoh_s, scenario_4_nf_med_mfpfr_edbm_evaporaponds.Qnaoh_s]
#M_naoh=[sc1_swro_as_feed.Qnaoh_s, sc1_swro_as_feed_nf_ro_med.Qnaoh_s, sc2_water_mining.Qnaoh_s, 0, scenario_3_nf_med_thcryst_mfpfr_edbm.Qnaoh_s, scenario_4_nf_med_mfpfr_edbm_evaporaponds.Qnaoh_s]

plt.bar(X_axis - 0.4, M_naoh, 0.4, color="#00516a")  
plt.xticks(X_axis, X)
plt.xlabel("Scenarios")
plt.ylabel("NaOH Production (kg/hr)")
# plt.title("NaOH Production")
# plt.legend()
plt.savefig('NaOHProduction.png')
plt.show()

#Figure 11: hcl production 
M_hcl=[0,  sc2.Qhcl_s, sc1.Qhcl_s,  sc4.Qhcl_s]#, scenario_3_nf_med_thcryst_mfpfr_edbm.Qhcl_s, scenario_4_nf_med_mfpfr_edbm_evaporaponds.Qhcl_s]
#M_hcl=[sc1_swro_as_feed.Qhcl_s, sc1_swro_as_feed_nf_ro_med.Qhcl_s, sc2_water_mining.Qhcl_s, 0, scenario_3_nf_med_thcryst_mfpfr_edbm.Qhcl_s, scenario_4_nf_med_mfpfr_edbm_evaporaponds.Qhcl_s]

plt.bar(X_axis - 0.4, M_hcl, 0.4, color="#00516a")  
plt.xticks(X_axis, X)
plt.xlabel("Scenarios")
plt.ylabel("HCl Production (kg/hr)")
# plt.title("hCl Production")
# plt.legend()
plt.savefig('HClProduction.png')
plt.show()

#Figure 12: CO2 emissions 
co2_em=[sc3.emis_t, sc2.emis_t, sc1.emis_t, sc4.emis_t]#, scenario_3_nf_med_thcryst_mfpfr_edbm.emis_t, scenario_4_nf_med_mfpfr_edbm_evaporaponds.emis_t]
#co2_em=[sc1_swro_as_feed.emis_t, sc1_swro_as_feed_nf_ro_med.emis_t, sc2_water_mining.emis_t, sc3_water_recovery.emis_t, scenario_3_nf_med_thcryst_mfpfr_edbm.emis_t, scenario_4_nf_med_mfpfr_edbm_evaporaponds.emis_t]
co2_em = [i * constants.hr/1e6 for i in co2_em]
plt.bar(X_axis - 0.4, co2_em, 0.4, color="#00516a")  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xticks(X_axis, X)
plt.xlabel("Scenarios")
plt.ylabel("CO\u2082 emissions (KTon/year)")
# plt.title("co2 emissions")
# plt.legend()
plt.savefig('CO2emissions.png')
plt.show()

#Figure 13: emissions vs thermal consumption 
fig,ax1=plt.subplots()
ax1.bar(X_axis- 0.4, co2_em , 0.4, color="#00516a", label='Emissions')  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_ylabel("CO\u2082 emissions (KTon/year)")
ax1.set_ylim(0)
plt.xticks(X_axis, X)
plt.legend(loc='upper left')
ax2=ax1.twinx()
ax2.bar(X_axis, Eth, width=0.4, color="sandybrown", label = 'Thermal energy consumption')
ax2.set_ylabel("Thermal energy consumption (GWh/year)")
# ax2.set_ylim(0,5.0*1e7)
plt.legend(loc='upper right')
fig.set_size_inches(8, 6)
plt.tight_layout()
plt.savefig('CO2vsthermal.png')
# plt.title("thermal cons vs co2 emissions")
plt.show()



#Figure 14: Electrical consumption Vs thermal consumption 
fig,ax1=plt.subplots()
ax1.bar(X_axis- 0.4,  Eel , 0.4, color="#00516a", label='Electrical (GWel)')  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_ylabel("Electrical energy consumption (GW)")
ax1.set_ylim(0)
plt.xticks(X_axis, X)
ax2=ax1.twinx()
ax2.bar(X_axis, Eth, width=0.4,  color="#008B8B", label = 'Thermal (GWth)')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_ylabel("Thermal energy consumption (GWh/year)")
# plt.title("electrical energy vs thermal energy")

lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

# Position the legend outside the figure
ax2.legend(lines + lines2, labels + labels2, loc='upper right', bbox_to_anchor=(1.5, 0.6))

plt.savefig('electricVSthermal.png', bbox_inches='tight')
plt.show()


#Figure 15: OPEX Vs CAPEX
fig,ax1=plt.subplots()
ax1.bar(X_axis- 0.4,  OPEX , 0.4, color="#00516a", label='OPEX')  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_ylabel("OPEX (M€/year)")
ax1.set_ylim(0)
plt.xticks(X_axis, X)
plt.legend(loc='upper left')
ax2=ax1.twinx()
ax2.bar(X_axis, Cap, width=0.4,  color="sandybrown", label = 'Capex ')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_ylabel("CAPEX (M€)")
plt.legend(loc='upper right')


plt.savefig('CAPEXVSOPEX.png')
plt.show()

#Figure 16: OPEX Vs Eel

fig,ax1=plt.subplots()
ax1.bar(X_axis- 0.4,  OPEX , 0.4, color="#00516a", label='OPEX')  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_ylabel("OPEX (M€/year)")
ax1.set_ylim(0)
plt.xticks(X_axis, X)
plt.legend(loc='upper left')
ax2=ax1.twinx()
ax2.bar(X_axis, Eel, width=0.4,  color="#008B8B", label = 'Electr ')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_ylabel("Electrical energy consumption (GW)")
plt.legend(loc='upper right')


plt.savefig('EnergyVSOPEX.png')
plt.show()

#Figure 17: OPEX Vs Eth
fig,ax1=plt.subplots()
ax1.bar(X_axis- 0.4,  OPEX , 0.4, color="#00516a", label='OPEX')  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_ylabel("OPEX (M€/year)")
ax1.set_ylim(0)
plt.xticks(X_axis, X)
plt.legend(loc='upper left')
ax2=ax1.twinx()
ax2.bar(X_axis, Eth, width=0.4,  color="sandybrown", label = 'Eth ')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_ylabel("Thermal energy consumption (GW)")
plt.legend(loc='upper right')


plt.savefig('thermalVSOPEX.png')
plt.show()

#Figure 18: Resource recovery Vs Eel
Rr=[sc3.Rr, sc2.Rr, sc1.Rr, sc4.Rr]
Rr = [i * 100 for i in Rr]
fig,ax1=plt.subplots()
ax1.bar(X_axis-0.4,  Rr , 0.4, color="#00516a", label='Rr')  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_ylabel("Resource efficiency")
ax1.set_ylim(0)
plt.xticks(X_axis, X)
plt.legend(loc='upper left')
ax2=ax1.twinx()
ax2.bar(X_axis, Eel, width=0.4,  color="#008B8B", label = 'Electr ')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_ylabel("Electrical energy consumption (GW)")
plt.legend(loc='upper right')


plt.savefig('ResourceVSEel.png')
plt.show()

#Figure 19: Resource recovery Vs Eel
Rr2=[sc3.Rr2, sc2.Rr2, sc1.Rr2, sc4.Rr2]
Rr2 = [i * 100 for i in Rr2]
fig,ax1=plt.subplots()
ax1.bar(X_axis-0.4,  Rr2 , 0.4, color="#00516a", label='Resource recovery')  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_ylabel("Resource efficiency")
ax1.set_ylim(0)
plt.xticks(X_axis, X)

ax2=ax1.twinx()
ax2.bar(X_axis, Eel, width=0.4,  color="#008B8B", label = 'Electrical cons. ')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_ylabel("Electrical energy consumption (GW)")


lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

# Position the legend outside the figure
ax2.legend(lines + lines2, labels + labels2, loc='upper right', bbox_to_anchor=(0.8, -0.12))


plt.savefig('ResourceVSEe2l.png', bbox_inches='tight')
plt.show()

#Figure 19b: Resource recovery Vs Eth
Rr2=[sc3.Rr2, sc2.Rr2, sc1.Rr2, sc4.Rr2]
Rr2 = [i * 100 for i in Rr2]
fig,ax1=plt.subplots()
ax1.bar(X_axis-0.4,  Rr2 , 0.4, color="#00516a", label='Resource recovery')  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_ylabel("Resource efficiency")
ax1.set_ylim(0)
plt.xticks(X_axis, X)

ax2=ax1.twinx()
ax2.bar(X_axis, Eth, width=0.4,  color="#008B8B", label = 'Thermal Energy cons. ')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_ylabel("Thermal energy consumption (GW)")


lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

# Position the legend outside the figure
ax2.legend(lines + lines2, labels + labels2, loc='upper right', bbox_to_anchor=(0.8, -0.12))


plt.savefig('ResourceVSEthermal.png', bbox_inches='tight')
plt.show()


#Figure 20: Resource recovery 
plt.bar(X_axis - 0.4, Rr, 0.4, color="#00516a")  
plt.xticks(X_axis, X)
plt.xlabel("Scenarios")
plt.ylabel("Resource efficiency (%)")
plt.savefig('Resourceeff.png')
plt.show()

#Figure 21: Salt recovery 
Salt_r = [sc3.Salt_r, sc2.Salt_r, sc1.Salt_r, sc4.Salt_r]#
Salt_r_y = [i * constants.hr/1000 for i in Salt_r]
plt.bar(X_axis - 0.4, Salt_r_y, 0.4, color="#00516a")  
plt.xticks(X_axis, X)
plt.xlabel("Scenarios")
plt.ylabel("Salt recovery")
plt.savefig('saltrecov.png')
plt.show()

# add the third bar
#ax1.bar(X_axis + 0.4, Qw_y, 0.4, color="#a6cee3", label='water production')

#Figure 22: Salt recovery Vs CO2 
fig,ax1=plt.subplots()
ax1.bar(X_axis-0.4,  Salt_r_y , 0.4, color="#00516a", label='Salt recovery')  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_ylabel("Salt recovery (ton/year)")
ax1.set_ylim(0)
plt.xticks(X_axis, X)
plt.legend(loc='upper left')
ax2=ax1.twinx()
ax2.bar(X_axis, co2_em, width=0.4,  color="sandybrown", label = 'CO\u2082 emissions ')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_ylabel("CO\u2082 emissions (KTon/year)")
plt.legend(loc='upper right')


plt.savefig('saltVSco2_b.png')
plt.show()

#Figure 23: water recovery Vs CO2 
fig,ax1=plt.subplots()
ax1.bar(X_axis-0.4,  Qw_y, 0.4, color="#a6cee3", label='Water production')  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_ylabel("Water production in 1000 m\u00B3/year")
ax1.set_ylim(0)
plt.xticks(X_axis, X)
plt.legend(loc='upper left')
ax2=ax1.twinx()
ax2.bar(X_axis, co2_em, width=0.4,  color="sandybrown", label = 'CO\u2082 emissions ')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_ylabel("CO\u2082 emissions (KTon/year)")
plt.legend(loc='upper right')


plt.savefig('waterVSco2.png')
plt.show()

#Figure 24: salt recovery +water Vs CO2 
# define the width of each bar
fig, ax1 = plt.subplots()
bar_width = 0.25
space_width = 0.05

X_axis = np.arange(len(X))

ax1.bar(X_axis - bar_width - space_width, Qw_y, bar_width, color="#a6cee3", label='Water production')
ax1.bar(X_axis, Salt_r_y, bar_width, color="#00516a", label='Salt recovery')
ax1.set_ylabel("Production (ton/year)")
ax1.set_ylim(0)
ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_xticks(X_axis)
ax1.set_xticklabels(X)
ax1.tick_params(axis='x', which='major', pad=10)

ax2 = ax1.twinx()
ax2.bar(X_axis + bar_width + space_width, co2_em, width=bar_width, color="sandybrown", label='CO\u2082 emissions')
ax2.set_ylabel("CO\u2082 emissions KTon/year")
ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

handles1, labels1 = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()
handles = handles1 + handles2
labels = labels1 + labels2
ax1.legend(handles, labels, loc='upper left', fontsize=10)

# adjust x-axis limits to add space between the bars
ax1.set_xlim(X_axis[0] - bar_width - space_width, X_axis[-1] + 2*bar_width + space_width)

plt.tight_layout()
plt.savefig('salt_waterVSco2.png')
plt.show()



#Figure 25: Salt recovery Vs Cost (production cost or OPEX) 
fig,ax1=plt.subplots()
ax1.bar(X_axis-0.4,  Salt_r_y , 0.4, color="#00516a", label='Salt recovery')  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_ylabel("Salt recovery (ton/year)")
ax1.set_ylim(0)
plt.xticks(X_axis, X)
plt.legend(loc='upper left')
ax2=ax1.twinx()
ax2.bar(X_axis, OPEX, width=0.4,  color="sandybrown", label = 'OPEX ')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_ylabel("OPEX (M€/year)")
plt.legend(loc='upper right')


plt.savefig('saltVScost.png')
plt.show()

#Figure 26: Salt recovery Vs Eel consumption 
fig,ax1=plt.subplots()
ax1.bar(X_axis-0.4,  Salt_r_y , 0.4, color="#00516a", label='Salt recovery')  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_ylabel("Salt recovery (ton/year)")
ax1.set_ylim(0)
plt.xticks(X_axis, X)

ax2=ax1.twinx()
ax2.bar(X_axis, Eel, width=0.4,  color="#008B8B", label = 'Electrical cons. ')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_ylabel("Electrical energy consumption (GW)")
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

# Position the legend outside the figure
ax2.legend(lines + lines2, labels + labels2, loc='upper right', bbox_to_anchor=(0.8, -0.12))
plt.savefig('saltVSEel.png', bbox_inches='tight')
plt.show()

#Figure 27: resource efficiency Vs feasibility  
Ec_m = [sc1.Ec_m, sc2.Ec_m, sc3.Ec_m, sc4.Ec_m]#economic margin 
fig,ax1=plt.subplots()
ax1.bar(X_axis-0.4,  Rr , 0.4, color="#00516a", label='Rr')  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_ylabel("Resource efficiency")
ax1.set_ylim(0)
plt.xticks(X_axis, X)
plt.legend(loc='upper left')
ax2=ax1.twinx()
ax2.bar(X_axis, Ec_m, width=0.4,  color="sandybrown", label = 'Economic Margin ')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_ylabel("Economic Margin euro/m\u00B3 of water")
plt.legend(loc='upper right')


plt.savefig('ResourceVSEconomic.png')
plt.show()

#Figure 28: resource efficiency Vs revenues 
Ec_m = [sc3.Ec_m, sc2.Ec_m, sc1.Ec_m, sc4.Ec_m]#economic margin 
fig,ax1=plt.subplots()
ax1.bar(X_axis-0.4,  Rr , 0.4, color="#00516a", label='Rr')  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_ylabel("Resource efficiency")
ax1.set_ylim(0)
plt.xticks(X_axis, X)
plt.legend(loc='upper left')
ax2=ax1.twinx()
ax2.bar(X_axis, Rev, width=0.4,  color="sandybrown", label = 'Revenues ')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_ylabel("Revenues (M€/year)")
plt.legend(loc='upper right')


plt.savefig('ResourceVSEcrevenues.png')
plt.show()

#Figure 29: resource efficiency Vs OPEX
fig,ax1=plt.subplots()
ax1.bar(X_axis-0.4, Rr, 0.4, color="#00516a", label='Resource efficiency (%)')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
ax1.set_ylabel("Resource efficiency")
ax1.set_ylim(0)
plt.xticks(X_axis, X)

ax2 = ax1.twinx()
ax2.bar(X_axis, OPEX, width=0.4, color="sandybrown", label='OPEX')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
ax2.set_ylabel("OPEX (M€/year)")

lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

# Position the legend outside the figure
ax2.legend(lines + lines2, labels + labels2, loc='upper right', bbox_to_anchor=(0.8, -0.12))

# Save the figure with the legend outside
plt.savefig('Resource_effVSOPEX.png', bbox_inches='tight')
plt.show()

#Figure 30: resource efficiency Vs CAPEX
fig,ax1=plt.subplots()
ax1.bar(X_axis-0.4, Rr, 0.4, color="#00516a", label='Resource efficiency (%)')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
ax1.set_ylabel("Resource efficiency")
ax1.set_ylim(0)
plt.xticks(X_axis, X)

ax2 = ax1.twinx()
ax2.bar(X_axis, Cap, width=0.4, color="sandybrown", label='CAPEX')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
ax2.set_ylabel("CAPEX (M€)")

lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

# Position the legend outside the figure
ax2.legend(lines + lines2, labels + labels2, loc='upper right', bbox_to_anchor=(0.8, -0.12))

# Save the figure with the legend outside
plt.savefig('Resource_effVSCAPEX.png', bbox_inches='tight')
plt.show()

#figure 31 production efficiency 
prd_eff= [sc3.prd_eff, sc2.prd_eff, sc1.prd_eff, sc4.prd_eff]
plt.bar(X_axis - 0.4, prd_eff, 0.4, color="#00516a")  
plt.xticks(X_axis, X)
plt.xlabel("Scenarios")
plt.ylabel("Production efficiency (%)")
plt.savefig('productioneff.png')
plt.show()


#figure 32 water consumption -> water footprint 

Water_cons=[sc3.Q_w_in,sc2.Q_w_in, sc1.Q_w_in, sc4.Q_w_in]
Water_cons = [i * constants.hr/1000 for i in Water_cons]
plt.bar(X_axis - 0.4, Water_cons, 0.4, color="#00516a")  
plt.xticks(X_axis, X)
plt.xlabel("Scenarios")
plt.ylabel("Water footprint (m\u00B3/year)")
plt.savefig('waterfootprint.png')
plt.show()


#figure 33 water consumption -> water footprint 

Q_br_prod=[sc3.Q_br_prod,sc2.Q_br_prod, sc1.Q_br_prod, sc4.Q_br_prod]
Q_br_prod = [i * constants.hr/1000 for i in Q_br_prod]
plt.bar(X_axis - 0.4, Q_br_prod, 0.4, color="#00516a")  
plt.xticks(X_axis, X)
plt.xlabel("Scenarios")
plt.ylabel("Brine production (ton/year)")
plt.savefig('brineproduction.png')
plt.show()

#Figure 34: production efficiency (revenue/cost) vs resource efficiency 

fig,ax1=plt.subplots()
ax1.bar(X_axis-0.4,  Rr2 , 0.4, color="#00516a", label='Resource recovery')  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_ylabel("Resource recovery efficiency")
ax1.set_ylim(0)
plt.xticks(X_axis, X)
ax2=ax1.twinx()
ax2.bar(X_axis, prd_eff, width=0.4,  color="sandybrown", label = 'Production efficiency ')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_ylabel("Production efficiency (€/€)")
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

# Position the legend outside the figure
ax2.legend(lines + lines2, labels + labels2, loc='upper right', bbox_to_anchor=(0.8, -0.12))
plt.savefig('ResourceVSproduction_Eff.png', bbox_inches='tight')
plt.show()

#Figure 35: production efficiency (revenue/cost) vs resource efficiency 

fig,ax1=plt.subplots()
ax1.bar(X_axis-0.4,  Rr2 , 0.4, color="#00516a", label='Rr')  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_ylabel("Resource recovery efficiency")
ax1.set_ylim(0)
plt.xticks(X_axis, X)
plt.legend(loc='upper left')
ax2=ax1.twinx()
ax2.bar(X_axis, Rev, width=0.4,  color="sandybrown", label = 'Revenue ')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_ylabel("Revenue (M€/year)")
plt.legend(loc='upper right')


plt.savefig('ResourceVSrevenue.png')
plt.show()

#Figure 36: Salt recovery Vs Eth consumption 
fig,ax1=plt.subplots()
ax1.bar(X_axis-0.4,  Salt_r_y , 0.4, color="#00516a", label='Salt recovery')  
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_ylabel("Salt recovery (ton/year)")
ax1.set_ylim(0)
plt.xticks(X_axis, X)

ax2=ax1.twinx()
ax2.bar(X_axis, Eth, width=0.4,  color="#008B8B", label = 'Thermal cons. ')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_ylabel("Thermal energy consumption (GW)")
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

# Position the legend outside the figure
ax2.legend(lines + lines2, labels + labels2, loc='upper right', bbox_to_anchor=(0.8, -0.12))
plt.savefig('saltVSEth.png', bbox_inches='tight')
plt.show()

#Figure 37: Water recovery Vs Eth consumption 
fig,ax1=plt.subplots()
ax1.bar(X_axis-0.4,  Qw_y, 0.4, color="#a6cee3", label='Water production') 
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_ylabel("Water recovery (1000 m\u00B3/year)")
ax1.set_ylim(0)
plt.xticks(X_axis, X)

ax2=ax1.twinx()
ax2.bar(X_axis, Eth, width=0.4,  color="#008B8B", label = 'Thermal cons. ')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_ylabel("Thermal energy consumption (GW)")
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

# Position the legend outside the figure
ax2.legend(lines + lines2, labels + labels2, loc='upper right', bbox_to_anchor=(0.8, -0.12))
plt.savefig('waterVSEth.png', bbox_inches='tight')
plt.show()

#Figure 38: Water recovery Vs Eel consumption 
fig,ax1=plt.subplots()
ax1.bar(X_axis-0.4,  Qw_y, 0.4, color="#a6cee3", label='Water production') 
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_ylabel("Water recovery (1000 m\u00B3/year)")
ax1.set_ylim(0)
plt.xticks(X_axis, X)

ax2=ax1.twinx()
ax2.bar(X_axis, Eel, width=0.4,  color="#008B8B", label = 'Electrical cons. ')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_ylabel("Electrical energy consumption (GW)")
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

# Position the legend outside the figure
ax2.legend(lines + lines2, labels + labels2, loc='upper right', bbox_to_anchor=(0.8, -0.12))
plt.savefig('waterVSEel.png', bbox_inches='tight')
plt.show()
#export data to excel and make tables 

#%%levelized cost for each product #%%levelized cost for each product 
LPCw_1=[sc3.LPCw, sc2.LPCw, sc1.LPCw, sc4.LPCw_1]
LPCnacl_1=[sc3.LPCnacl_1, sc2.LPCnacl_1, sc1.LPCnacl_1, sc4.LPCnacl_1]
LPCmg_1=[sc3.LPCmg_1,sc2.LPCmg_1,sc1.LPCmg_1,sc4.LPCmg_1]
LPCna2so4_1=[sc3.LPCna2so4_1, sc2.LPCna2so4_1, sc1.LPCna2so4_1, sc4.LPCna2so4_1]
LPChcl_1=[sc3.LPChcl_1, sc2.LPChcl_1, sc1.LPChcl_1, sc4.LPChcl_1]
LPCnaoh_1=[sc3.LPCnaoh_1, sc2.LPCnaoh_1, sc1.LPCnaoh_1, sc4.LPCnaoh_1]

dfLPC1=pd.DataFrame({'sce 1': [LPCw_1[0],LPCnacl_1[0], LPCmg_1[0],LPCna2so4_1[0],LPChcl_1[0] , LPCnaoh_1[0]], 
                 'sce 2': [LPCw_1[1],LPCnacl_1[1], LPCmg_1[1],LPCna2so4_1[1], LPChcl_1[1], LPCnaoh_1[1]],
                 'sce 3': [LPCw_1[2],LPCnacl_1[2], LPCmg_1[2],LPCna2so4_1[2], LPChcl_1[2], LPCnaoh_1[2]], 
                 'sce 4': [LPCw_1[3],LPCnacl_1[3], LPCmg_1[3],LPCna2so4_1[3], LPChcl_1[3], LPCnaoh_1[3]], })

LPCw_2=[sc3.LPCw, sc2.LPCw, sc1.LPCw, sc4.LPCw_2]
LPCnacl_2=[sc3.LPCnacl_2, sc2.LPCnacl_2, sc1.LPCnacl_2, sc4.LPCnacl_2]
LPCmg_2=[sc3.LPCmg_2,sc2.LPCmg_2,sc1.LPCmg_2,sc4.LPCmg_2]
LPCna2so4_2=[sc3.LPCna2so4_2, sc2.LPCna2so4_2, sc1.LPCna2so4_2, sc4.LPCna2so4_2]
LPChcl_2=[sc3.LPChcl_2, sc2.LPChcl_2, sc1.LPChcl_2, sc4.LPChcl_2]
LPCnaoh_2=[sc3.LPCnaoh_2, sc2.LPCnaoh_2, sc1.LPCnaoh_2, sc4.LPCnaoh_2]

dfLPC2=pd.DataFrame({'sce 1': [LPCw_2[0],LPCnacl_2[0], LPCmg_2[0],LPCna2so4_2[0],LPChcl_2[0] , LPCnaoh_2[0]], 
                 'sce 2': [LPCw_2[1],LPCnacl_2[1], LPCmg_2[1],LPCna2so4_2[1], LPChcl_2[1], LPCnaoh_2[1]],
                 'sce 3': [LPCw_2[2],LPCnacl_2[2], LPCmg_2[2],LPCna2so4_2[2], LPChcl_2[2], LPCnaoh_2[2]], 
                 'sce 4': [LPCw_2[3],LPCnacl_2[3], LPCmg_2[3],LPCna2so4_2[3], LPChcl_2[3], LPCnaoh_2[3]], })

#%%
#create dataframes 

sum_res_2=np.array([Qw_y[2],Eel[2], Eth[2],co2_em[2], OPEX[2], Cap[2],  Rr[2], prd_eff[2], Water_cons[2],Q_br_prod[2], LPCw_1[2]])
ind=np.array(["Water production", "Total electrical consumption (GWh)", "Total thermal energy consumption (GWh)","Carbon dioxide emission (Kton co2/year) ",
              "OPEX (M€/year)", "CAPEX (M€)", "Resource efficiency (%)", "Production efficiency (%)", "Water footprint (m3/year)", 
              "Brine production (ton/year)", "levelized cost of water €/m\u00B3"])
df2=pd.DataFrame(sum_res_2, ind)
sum_res_1=np.array([Qw_y[1],Eel[1], Eth[1],co2_em[1], OPEX[1], Cap[1], Rr[1], prd_eff[1], Water_cons[1],Q_br_prod[1], LPCw_1[1]])
df1=pd.DataFrame(sum_res_1, ind)
sum_res_3=np.array([Qw_y[0],Eel[0], Eth[0],co2_em[0], OPEX[0], Cap[0], Rr[0], prd_eff[0], Water_cons[0],Q_br_prod[0], LPCw_1[0]])
df3=pd.DataFrame(sum_res_3, ind)
sum_res_4=np.array([ Qw_y[3],Eel[3], Eth[3],co2_em[3], OPEX[3], Cap[3], Rr[3], prd_eff[3], Water_cons[3],Q_br_prod[3], LPCw_1[3]])
df4=pd.DataFrame(sum_res_4, ind)

#all
df=pd.DataFrame({'sce 1': [Qw_y[0],Eel[0], Eth[0],co2_em[0], OPEX[0], Cap[0], Rr[0], prd_eff[0], Water_cons[0],Q_br_prod[0], LPCw_1[0]], 
                 'sce 2': [Qw_y[1],Eel[1], Eth[1],co2_em[1], OPEX[1], Cap[1],Rr[1], prd_eff[1], Water_cons[1],Q_br_prod[1], LPCw_1[1]],
                 'sce 3': [Qw_y[2],Eel[2], Eth[2],co2_em[2], OPEX[2], Cap[2], Rr[2],prd_eff[2],Water_cons[2],Q_br_prod[2], LPCw_1[2]], 
                 'sce 4': [Qw_y[3],Eel[3], Eth[3],co2_em[3], OPEX[3], Cap[3], Rr[3], prd_eff[3],Water_cons[3],Q_br_prod[3], LPCw_1[3]], })

sc1_dfel=pd.DataFrame(sc1.E_el_all , sc1.tec_names)
sc1_dfeth=pd.DataFrame(sc1.E_th_all, sc1.tec_names)
sc2_dfel=pd.DataFrame(sc2.E_el_all , sc2.tec_names)
sc2_dfeth=pd.DataFrame(sc2.E_th_all, sc2.tec_names)
sc3_dfel=pd.DataFrame(sc3.E_el_all , sc3.tec_names)
sc3_dfeth=pd.DataFrame(sc3.E_th_all, sc3.tec_names)
sc4_dfel=pd.DataFrame(sc4.E_el_all , sc4.tec_names)
sc4_dfeth=pd.DataFrame(sc4.E_th_all, sc4.tec_names)
#%%export to excel
with pd.ExcelWriter('results_scenario.xlsx') as writer:
    df1.to_excel(writer,sheet_name="sc2_water_mining")
    df2.to_excel(writer,sheet_name="sc1_water_recovery_scenario")
    df3.to_excel(writer,sheet_name="sc3_swro_as_feed")
    df4.to_excel(writer,sheet_name="sc4_swro_ED")
    df.to_excel(writer,sheet_name="all scenarios")
    dfLPC1.to_excel(writer,sheet_name="levelized costs")
    dfLPC2.to_excel(writer,sheet_name="levelized costs", startcol=0, startrow=10, header=True)
    sc2.sum_table_C.to_excel(writer,sheet_name="sc2_water_mining", startcol=0, startrow=20, header=False)
    sc2.sum_table_f.to_excel(writer,sheet_name="sc2_water_mining", startcol=0, startrow=18, header=True, index=True )
    sc2.sum_table_d.to_excel(writer,sheet_name="sc2_water_mining", startcol=0, startrow=29, header=True)
    sc3.sum_table_C.to_excel(writer,sheet_name="sc1_water_recovery_scenario", startcol=0, startrow=20, header=False)
    sc3.sum_table_f.to_excel(writer,sheet_name="sc1_water_recovery_scenario", startcol=0, startrow=18, header=True)
    sc3.sum_table_d.to_excel(writer,sheet_name="sc1_water_recovery_scenario", startcol=0, startrow=29, header=True)
    sc1.sum_table_C.to_excel(writer,sheet_name="sc3_swro_as_feed", startcol=0, startrow=20, header=False)
    sc1.sum_table_f.to_excel(writer,sheet_name="sc3_swro_as_feed", startcol=0, startrow=18, header=True)
    sc1.sum_table_d.to_excel(writer,sheet_name="sc3_swro_as_feed", startcol=0, startrow=30, header=False)
    sc4.sum_table_C.to_excel(writer,sheet_name="sc4_swro_ED", startcol=0, startrow=20, header=False)
    sc4.sum_table_f.to_excel(writer,sheet_name="sc4_swro_ED", startcol=0, startrow=18, header=True)
    sc4.sum_table_d.to_excel(writer,sheet_name="sc4_swro_ED", startcol=0, startrow=30, header=False)
    sc1_dfel.to_excel(writer,sheet_name="sc3_swro_as_feed", startcol=0, startrow=33, header=True)
    sc1_dfeth.to_excel(writer,sheet_name="sc3_swro_as_feed", startcol=2, startrow=33, header=True)
    sc2_dfel.to_excel(writer,sheet_name="sc2_water_mining", startcol=0, startrow=33, header=True)
    sc2_dfeth.to_excel(writer,sheet_name="sc2_water_mining", startcol=2, startrow=33, header=True)
    sc3_dfel.to_excel(writer,sheet_name="sc1_water_recovery_scenario", startcol=0, startrow=33, header=True)
    sc3_dfeth.to_excel(writer,sheet_name="sc1_water_recovery_scenario", startcol=2, startrow=33, header=True)
    sc4_dfel.to_excel(writer,sheet_name="sc4_swro_ED", startcol=0, startrow=33, header=True)
    sc4_dfeth.to_excel(writer,sheet_name="sc4_swro_ED", startcol=2, startrow=33, header=True)
    