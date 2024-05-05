import unittest
from thermal_cryst_f import thermal_calc, conc_cal, calculate_energy
import density_calc
import constants
import scaleup

class TestThermalCrystalUnit(unittest.TestCase):
    def setUp(self):
        # Define input parameters for testing
        self.T_op = 60  # Operating temperature (°C)
        self.Qf = 1000  # Inlet flow rate (kg/hr)
        self.Cf_s = 200  # Total ion concentration of solution (g/L)
        Cf_caso4 = 10  # Concentration of CaSO4 in the solution (g/L)
        self.T_in = 40  # Inlet feed temperature (°C)
        Cf_in = [80.42, 116.69, 2.29, 0.01, 0.04, 0.54]  # List of concentration of ions in the solution (g/L)
        self.salt_mois = 20  # Weight percent of NaCl at saturation
        self.LHV_v = 2199.7  # Boiling point elevation (°C)
        self.LHV_s = 2357.69  # Steam temperature (°C)
        self.T_cw_o = 40  # Cooling water inlet temperature (°C)
        self.T_cw_f = 25  # Cooling water outlet temperature (°C)
        self.d_sol=density_calc.density_calc(self.T_in, sum(Cf_in)) 
        
        # Create an instance of the thermal_calc class with the defined parameters
        self.thermal_calc_instance = thermal_calc(self.T_op, self.Qf, self.Cf_s, Cf_caso4,self.T_in, Cf_in, self.salt_mois, self.LHV_v, self.LHV_s,self.T_cw_o, self.T_cw_f)
        self.result_mass=self.thermal_calc_instance.mass_bal_cryst()
        self.result_heat=self.thermal_calc_instance.heat_bal_cryst()
        
        
    def test_mass_bal_cryst(self):
        # Call the mass_bal_cryst method
        self.thermal_calc_instance.mass_bal_cryst()
        self.result_mass=self.thermal_calc_instance.mass_bal_cryst()
        
        self.ev_mass= self.thermal_calc_instance.ev_mass
        self.solid_mass=self.thermal_calc_instance.solid_mass
        
        expected_ev_mass = self.Qf-self.solid_mass
        expected_solid_mass=(self.Qf*self.Cf_s/self.d_sol/1000)/(1-self.salt_mois/100)
        
        # Check if the calculated values are reasonable
        self.assertAlmostEqual(self.solid_mass, expected_solid_mass, delta=expected_solid_mass*0.01)
        self.assertAlmostEqual(self.ev_mass, expected_ev_mass, delta=expected_ev_mass*0.01)

    def test_heat_bal_cryst(self):
        result_mass=self.thermal_calc_instance.mass_bal_cryst()
        
        # Call the heat_bal_cryst method
        Cp_f=3.14 # Feed specific heat capacity (units: KJ* Kg*oC)
        CP_cw=4.187 # Water specific heat capacity (units: KJ* Kg*oC)
        UA=45990
        self.thermal_calc_instance.heat_bal_cryst()
        self.ev_mass= self.thermal_calc_instance.ev_mass
        self.result_heat=self.thermal_calc_instance.heat_bal_cryst()
        
        heat_req=self.thermal_calc_instance.heat_req
        cw_mass=self.thermal_calc_instance.cw_mass
        expected_heat_req=self.LHV_v*self.ev_mass+self.Qf*Cp_f*(self.T_op-self.T_in)
        expected_cw_mass=(self.ev_mass*self.LHV_v)/(CP_cw*(self.T_cw_o-self.T_cw_f)) 
        self.assertAlmostEqual(heat_req, expected_heat_req, delta=expected_heat_req*0.01)
        self.assertAlmostEqual(cw_mass, expected_cw_mass, delta=expected_cw_mass*0.01)

    def test_calculate_energy(self):
        # Define input parameters for calculate_energy function
        Qf = 1000  # Flow rate of the feed solution (m3/h)
        Q_evap_mass = 500  # Distillate flow rate (m3/h)
        Qcw = 200  # Flow rate of cooling water (m3/h)
        M_Nacl = 100  # Mass of NaCl (kg)
        heat_req = 5000  # Heat required (kJ/h)
        d_sol = 1.2  # Density of the feed solution (kg/m3)
        dp_f = 0.1  # Pressure drop for feed (bar)
        dp_w = 1  # Pressure drop for water streams (bar)
        dp_slurry = 3.5  # Pressure drop for slurry streams (bar)
        dp_cw = 2  # Pressure drop for cooling water (bar)
        npump = 0.8  # Pump efficiency

        # Call the calculate_energy function
        E_el_th_Cr, E_th_th_Cr, SEC_el_f, SEC_el_NaCl, SEC_el_w = calculate_energy(Qf, Q_evap_mass, Qcw, M_Nacl,
                                                                                          heat_req, d_sol, dp_f, dp_w,
                                                                                          dp_slurry, dp_cw, npump)
        expected_E_el_th_Cr=0.061
        # Check if the calculated energy values are reasonable
        self.assertAlmostEqual(E_el_th_Cr, expected_E_el_th_Cr, delta=expected_E_el_th_Cr*0.1)


if __name__ == '__main__':
    unittest.main()
