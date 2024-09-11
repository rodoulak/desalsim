#costs calculations functions 
from desalsim import constants 
#%%
#symbols:
#hr-> hours
#lf -> lifetime 
#el -> electricity
#s-> steam
#pr -> price 
#c-> cost
#eq-> equipment 
#conc -> consumption 
#wat -> water 
#cw -> cooling water 
#E -> Energy 
#prd -> product mass rate 

#%%constants
class econom:
    def __init__(self, eq_c, el_conc, s_conc, chem1_conc, chem1_pr,chem2_conc, chem2_pr, cw_conc, wat_conc):
        """
        Initialize an instance of the economic costs.

        Args:
            equipment_cost (float): Cost of equipment in euros.
            electricity_consumption (float): Electricity consumption in kWh/year.
            steam_consumption (float): Steam consumption in kWh/year.
            chemical1_concentration (float): Concentration of chemical 1 in solution (e.g., 1M).
            chemical1_price (float): Price of chemical 1 per unit (e.g., euro/L).
            chemical2_concentration (float): Concentration of chemical 2 in solution (e.g., 1M).
            chemical2_price (float): Price of chemical 2 per unit (e.g., euro/L).
            cooling_water_consumption (float): Cooling water consumption in kg/year.
            water_consumption (float): Water consumption in kg/year.
        """
        self.eq_c= eq_c
        self.el_conc=el_conc
        self.s_conc=s_conc
        self.chem1_conc=chem1_conc
        self.chem1_pr=chem1_pr
        self.chem2_conc=chem2_conc
        self.chem2_pr=chem2_pr
        self.cw_conc=cw_conc
        self.wat_conc=wat_conc
        
    def capex_calc(self, capex_assumptions):
        """
        Calculate the capital expenditure (CAPEX) of the unit.
        
        Args:
            inst_c (float): Installation cost in euros 
            buil_c (float): Building, process and auxillary cost in euros
            land_c (float): Land cost in euros
            hard_c (float): Hardware costs are the sum of costs on purchased equipment and installation in euro 
            dir_c (float):  Direct costs in euro. Direct costs consist of purchased equipment, Purchased-equipment installation, Instrumentation and controls, 
                            Piping, Electrical systems, Buildings (including services), Yard improvements,Service facilities, Land 
            ind_c (float):  Indirect costs in euro. Indirect costs consist of engineering and supervision costs, legal expenses, Construction expenses, 
                            Constructor's fee, Contingency
            fix_c (float):  fixed-capital investment in euro. fixed-capital investment is the capital necessary for the in- stalled process equipment 
                            with all components that are needed for complete process operation. Expenses
            work_c (float): working capital in euro 
            t_capital_inv (float): total capital investment 
            
        """
        #Calculate Installation cost in euros
        self.inst_c=capex_assumptions[0]*self.eq_c  
        
        #Calculate Building, process and auxillary cost in euros
        self.buil_c=capex_assumptions[1]*self.eq_c   
        
        #Calculate Land cost in euros
        self.land_c=capex_assumptions[2]*self.eq_c   
        
        #Calculate Hardware costs in euro 
        self.hard_c=self.eq_c+self.inst_c 
        
        #Calculate Direct costs in euro
        self.dir_c= self.hard_c+self.buil_c+ self.land_c 
        
        #Calculate Indirect costs in euro
        self.ind_c=capex_assumptions[3]*self.dir_c 
        
        #Calculate fixed-capital investment in euro
        self.fix_c=self.dir_c+self.ind_c 
        
        #Calculate working capital in euro
        self.work_c=capex_assumptions[4]*self.fix_c
        
        #Calculate total capital investment in euro
        self.t_capital_inv=self.fix_c+self.work_c 
        return self.t_capital_inv
      
        
    
    def opex_calc(self, hr, el_pr, s_pr, cw_pr, w_pr, economic_assumptions):
        """
        Calculate the capital expenditure (OPEX) of the unit.
        
        Args:
            E_el (float): electricity cost in euro 
            E_th (float): thermal energy cost in euro 
            t_E_c total(float): energy costs in euro 
            chem_c (float): cost for chemicals in euro 
            cw_c (float): cost for cooling water in euro 
            wat_c (float): cost for water consumption in euro 
            main_c (float): maintenance cost  in euro 
            oper_sup_c (float): operating suppliers cost  in euro 
            oper_lab_c (float): operating labor in euro 
            super_c (float): direct supervisory and clerical labor in euro 
            lab_c (float): Laboratory charges in euro 
            pat_c (float): Patents and royalties in euro 
            fix_char (float): Fixed charges in euro 
            over_c (float): plant overhead costs in euro 
            OPEX (float): Operating costs in euro. It consists of utilities, maintenance, operating supplies, operating labor, direct supervisory and clerical labor, 
                         laboratory charges, patents and royalties, fixed charges, and plant overhead cost             
         """
         
        """utilities cost calculation """        


        #Calculate Electricity cost in euro 
        self.E_el=self.el_conc*hr*el_pr 
        
        #Calculate Thermal energy cost in euro 
        self.E_th=self.s_conc*hr*s_pr  
        
        #Calculate Total energy cost in euro 
        self.t_E_c=self.E_el+self.E_th     
        
        #Calculate cost for chemicals in euro 
        self.chem_c=self.chem1_conc*hr*self.chem1_pr+self.chem2_conc*hr*self.chem2_pr  
        
        #Calculate cost for cooling water in euro 
        self.cw_c=self.cw_conc*hr*cw_pr 
        
        #Calculate cost for water in euro 
        self.wat_c=self.wat_conc*hr*w_pr  
        
        """Operating costs calculation """
        main_c_percent, oper_sup_c_percent, oper_lab_c_percent, super_c_percent, lab_c_percent, pat_c_percent, fix_char_percent, over_c_percent, norm_factor = economic_assumptions
        
        #Calculate maintenance cost in euro
        self.main_c=main_c_percent*self.fix_c        
        
        #Calculate operating suppliers cost  in euro 
        self.oper_sup_c=oper_sup_c_percent*self.main_c  
        
        #Calculate operating labor in euro
        self.oper_lab_c=(self.t_E_c+self.chem_c+self.wat_c+self.cw_c+self.main_c+self.oper_sup_c)/norm_factor*oper_lab_c_percent     
        
        
        #Calculate direct supervisory and clerical labor in euro 
        self.super_c=super_c_percent*self.oper_lab_c 
        
        #Calculate Laboratory charges in euro 
        self.lab_c=lab_c_percent*self.oper_lab_c
        
        #Calculate Patents and royalties in euro
        self.pat_c=(self.t_E_c+self.chem_c+self.wat_c+self.cw_c+self.main_c+self.oper_sup_c)/norm_factor*pat_c_percent 
        
        #Calculate Fixed charges in euro
        self.fix_char=(self.t_E_c+self.chem_c+self.wat_c+self.cw_c+self.main_c+self.oper_sup_c)/norm_factor*fix_char_percent  
        
        #Calculate plant overhead costs in euro
        self.over_c=(self.t_E_c+self.chem_c+self.wat_c+self.cw_c+self.main_c+self.oper_sup_c)/norm_factor*over_c_percent
        
        #Calculate OPEX in euro 
        self.opex=self.t_E_c+self.chem_c+self.wat_c+self.cw_c+self.main_c+self.oper_sup_c+self.oper_lab_c+self.super_c+self.lab_c+self.pat_c+self.fix_char+self.over_c
        
        return self.opex
#%%
class revenue:
    """
    Calculate the revenues from selling products  of the unit.
    
    Args:
        prd_name (str): name of product 
        rev_prd (float): reveneues from product i in euro  
    """
    def __init__(self, prd, prd_name):
        self.prd=prd
        self.prd_name=prd_name
        
        
    def rev(self, hr, w_pr, nacl_pr, mgoh2_pr,na2so4_pr, naoh_pr, hcl_pr ):
        #density
        d_naoh=1.04 #kg/l for 1M solution 
        d_hcl=1.0585#kg/l for 1M solution 
        
        if self.prd_name=="Water":
            self.rev_prd=self.prd*w_pr*hr #euro/year
        elif self.prd_name=="NaCl":
            self.rev_prd=self.prd*nacl_pr*hr #euro/year
        elif self.prd_name=="Mg(OH)2":
            self.rev_prd=self.prd*mgoh2_pr*hr #euro/year
        elif self.prd_name=="Na2SO4":
            self.rev_prd=self.prd*na2so4_pr*hr #euro/year
        elif self.prd_name=="NaOH":
            self.prd=self.prd*d_naoh
            self.rev_prd=self.prd*naoh_pr*hr #euro/year
        elif self.prd_name=="HCl":
            self.prd=self.prd*d_hcl
            self.rev_prd=self.prd*hcl_pr*hr #euro/year

