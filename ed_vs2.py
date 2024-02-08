import numpy as np
import nanofiltration_unit_swro as nf
import constants
import math
import density_calc
# Constants
MWs = 58.44  # Molecular weight of NaCl (g/mol)
MWw=18.01528
Acp = 0.395  # Effective cell-pair area (m^2)
Sc_i = sum(nf.Cconc_out)  # Salinity at concentrate inlet (g/kg)
Sc_o = 200  # Salinity at concentrate outlet (g/kg)
N = 50  # Number of computational cells per cell-pair
Vcp = 8  # Applied voltage (V)
Vel=2.1
Ij=300 #Am2
Mem_eff=0.64

R = 0.002  # Resistance of rinse stream (ohm)
Rp = 0.015  # Resistance of polarization (ohm)
A = 1.1  # Active area of cell-pair (m^2)
F = 96485.3329  # Faraday constant (C/mol)
T=20+273.15 #K
dp=1
npump=0.8
# Transport and permeability numbers
def Ts_cp(S):
    return -4e-6*S**2 + 4e-5*S + 0.96

def Tw_cp(Sc,Sd):
    return -4e-5*Sc**2 - 1.9e-2*Sd + 11.2

def Ls_cp(Sc, Sd):
    return min(2e-12*Sd**2 - 3e-10*Sd + 6e-8, 2e-12*Sc**2 - 3e-10*Sc + 6e-8)

def Lw_cp(S):
    return 5*S**(-0.416)

def p_osmo(S):
    
    C1=S/MWs*constants.MW_Na
    z1=1
    C2=S/MWs*constants.MW_cl
    z2=1
    sum_Ci=sum([C1, C2])
    zi_2=[z1**2, z2**2]
    conc=[C1/constants.MW_Na, C2/constants.MW_cl] #mol/l
    mi=[C1*1000/(constants.MW_Na*1000*((1e+6-sum_Ci*1000)/1e+6)), C2*1000/(constants.MW_cl*1000*((1e+6-sum_Ci*1000)/1e+6))]
    mizi_2=[]
    for i in range(2) :
        mizi_2.append(mi[i]*zi_2[i])
        SI=sum(mizi_2)/2
        B=-348.662/(T)+6.72817-0.971307*math.log(T)
        C=40.5016/(T)-0.721404+0.103915*math.log(T)
        D=5321/(T)+233.76-0.9297*(T)+0.001417*(T)**2-0.0000008292*(T)**3
        S=1.17202*(sum(mizi_2)/sum(mi))*0.9982**0.5*(23375.556/(D*(T)))**1.5
        fi=1-S/(3.375*SI)*((1+1.5*SI**0.5)-2*math.log(1+1.5*SI**0.5)-1/(1+1.5*SI**0.5))+B*sum(mi)/2+C*(sum(mi)/2)**2
        sum_conc=sum(conc)
        p_osmo_psig=1.205*fi*(T)*sum(mi)
        p_osmo=p_osmo_psig/14.3 #bar
    return p_osmo

# Calculate number of cell-pairs
Acp_tot = np.ones(N) * Acp
Ncp = np.sum(Acp_tot) / Acp

# Initializations
Sc = np.zeros(N)
Sd = np.zeros(N)
Ns_c = np.zeros(N)
Ns_d = np.zeros(N)
Nw_c = np.zeros(N)
Nw_d = np.zeros(N)
Js = np.zeros(N)
Jw = np.zeros(N)
Qed_in=nf.Qconc/density_calc.density_calc(T-273, sum(nf.Cconc_out))


# Set initial values salt stream 
Sc[0] = Sc_i
Ns_c[0]=nf.Qconc*nf.Cconc_out[0]/MWs
Ms_in=nf.Qconc*sum(nf.Cconc_out)/1000 #kg salr/hr
Mw_in=nf.Qconc-Ms_in #kg water/hr
Nw_c[0]=Mw_in*1000/MWw

# Set initial values diluate stream 
Sd[0] = 35
Ns_d[0]=nf.Qconc*Sd[0]/MWs
Ms_in_d=nf.Qconc*Sd[0]/1000 #kg salr/hr
Mw_in_d=nf.Qconc-Ms_in_d #kg water/hr
Nw_d[0]=Mw_in_d*1000/MWw
Sd_o=20
# Iterate over cells
for j in range(1, N):
    Sc[j] = Sc[j-1] + (Sc_o - Sc_i) / (N - 1)
    Sd[j] = Sd[j-1] + (Sd_o - Sd[0]) / (N - 1)
    
    
    # Ns_c[j] = Acp_tot[j] * Sc[j] / (MWs)
    # Nw_c[j] = Acp_tot[j] * (1 - Sc[j] / 1000)
    
    Js[j] = Ts_cp(Sd[j]) * Ij / F - Ls_cp(Sc[j],Sd[j]) * (Sc[j] - Sd[j])
    Jw[j] = Tw_cp(Sc[j],Sd[j]) * Ij / F + Lw_cp(Sd[j]) * (p_osmo(Sc[j])-p_osmo(Sd[j]))
    
    Ns_c[j] = Ns_c[j-1] + Acp_tot[j] * Js[j]
    Ns_d[j] = Ns_d[j-1] - Acp_tot[j] * Js[j]
    
    Nw_c[j] = Nw_c[j-1] + Acp_tot[j] * Jw[j]
    Nw_d[j] = Nw_d[j-1] - Acp_tot[j] * Jw[j]
    

print("effluent concentration diluate is "+str(Sd))
print("effluent concentration brine is "+str(Sc))

Mc=Ns_c[N-1]*MWs/1000+Nw_c[N-1]*MWw/1000 #concnetrate stream mass flow rate (kg/hr)
print("Mass flowrate concentrate stream is "+str(Mc))
dc_out=density_calc.density_calc(T-273, Sc[N-1])/1000 #(kg/l)
Qc=Mc/dc_out #concnetrate stream volume flow rate (l/hr)

print("volume flowrate concentrate stream is "+str(Qc))

Md=Ns_d[N-1]*MWs/1000+Nw_d[N-1]*MWw/1000 #diluate stream mass flow rate (kg/hr)
Sd_f=Sd[N-1]
print("Mass flowrate diluate stream is "+str(Md))
dd_out=density_calc.density_calc(T-273, Sd[N-1])/1000
Qd=Md/dd_out #diluate stream volume flow rate (l/hr)
print("volume flowrate diluate stream is "+str(Qc))

print("in flow rate is "+str(nf.Qconc))
#power required 
bal=nf.Qconc*2-Md-Mc
print("balance is "+str(bal))

Ws = 0
for j in range(N):
    Ws += Ij * Acp * (Ncp * Vcp + Vel)

Ws=Ws/1000

Ppump=(2*Qed_in*1+(Qc*2+Qd*1))/3600*1e5/npump  #units: W ->  bar to J 1e5N/m2*1J/m ; hr to 3660s  
Eel_t=Ws+Ppump/1000
print("total energy consumption is "+str(Eel_t))


#specific cost
Amem_tot=2*np.sum(Acp_tot)/Mem_eff
Eq=600*Amem_tot
# # Calculate salt and water fluxes in each cell
# for j in range(N):
#     if j == 0:
#         ij = 0 # No current flows into first cell
#     else:
#         ij = (Vcp[j-1] - Vcp[j])/Rcp[j]
        
