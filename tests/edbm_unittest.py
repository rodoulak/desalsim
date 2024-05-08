import unittest
from Desalsim.edbm_unit_f import EDBMCalc  

class TestEDBMCalc(unittest.TestCase):
    def setUp(self):
        # Define sample input data for initialization
        Qin = 1000
        A = 0.4
        I_d = 400
        F=96485.3 #Coulombs/mol
        N = 50 
        C_s_in = [13.44, 20.725, 1.146, 0, 0, 0.18, 0, 10**(-4.71), 3.01551E-11]
        C_b_in = [0, 0, 0, 0, 0, 0, 0, 10**(-7), 10**(-7)]
        C_a_in = [0, 0, 0, 0, 0, 0, 0, 10**(-7), 10**(-7)]
        T = 20 + 273.15
        self.ph_s = 4.71# Sample pH value for the salt channel
        self.I_ext=A*I_d
        self.JA=3.6*self.I_ext/F
        
        # Create an instance of the EDBMCalc class
        self.edbm_calc = EDBMCalc(Qin, A, I_d, N, C_s_in, C_b_in, C_a_in, T)
        self.results_flowrate=self.edbm_calc.flowrate()
        self.results_mass_flowrate=self.edbm_calc.in_mass_flow_rates(self.ph_s)
        self.results_acid_channel=self.edbm_calc.acid_channel()
        self.results_base_channel=self.edbm_calc.base_channel()
        
    def test_flowrate(self):
        # Call the flowrate method
        self.edbm_calc.flowrate()
        self.results_flowrate=self.edbm_calc.flowrate()
        
        # Check if the flow rates are calculated correctly
        expected_Q1_s_in = self.edbm_calc.Qin / self.edbm_calc.N_trip
        self.assertAlmostEqual(self.edbm_calc.Q1_s_in, expected_Q1_s_in, delta=0.001)
        
        expected_Q1_b_in = self.edbm_calc.Qin / self.edbm_calc.N_trip
        self.assertAlmostEqual(self.edbm_calc.Q1_b_in, expected_Q1_b_in, delta=0.001)
        
        expected_Q1_a_in = self.edbm_calc.Qin / self.edbm_calc.N_trip
        self.assertAlmostEqual(self.edbm_calc.Q1_a_in, expected_Q1_a_in, delta=0.001)


    def test_in_mass_flow_rates(self):
        # Call the in_mass_flow_rates method with a sample value for ph_s        
        results_flowrate=self.results_flowrate
        self.edbm_calc.in_mass_flow_rates(self.ph_s)
        self.results_mass_flowrate=self.edbm_calc.in_mass_flow_rates(self.ph_s)
        
        # Assert the correctness of the output
        # Check if the inlet mass flow rates of each ion are calculated correctly
        for i in range(9):
            # Assuming Ci_s_in, Ci_b_in, and Ci_a_in are properly initialized in setUp
            expected_M_s_in_i = self.edbm_calc.Q1_s_in * self.edbm_calc.Ci_s_in[i] * self.edbm_calc.PM_i[i] / 1000
            self.assertAlmostEqual(self.edbm_calc.M_s_in[i], expected_M_s_in_i, delta=0.001)
            
            expected_M_b_in_i = self.edbm_calc.Q1_b_in * self.edbm_calc.Ci_b_in[i] * self.edbm_calc.PM_i[i] / 1000
            self.assertAlmostEqual(self.edbm_calc.M_b_in[i], expected_M_b_in_i, delta=0.001)
            
            expected_M_a_in_i = self.edbm_calc.Q1_a_in * self.edbm_calc.Ci_a_in[i] * self.edbm_calc.PM_i[i] / 1000
            self.assertAlmostEqual(self.edbm_calc.M_a_in[i], expected_M_a_in_i, delta=0.001)

    def test_acid_channel(self):
        results_flowrate=self.results_flowrate
        results_mass_flowrate=self.edbm_calc.in_mass_flow_rates(self.ph_s)
        # Call the acid_channel method
        self.edbm_calc.acid_channel()
        self.results_acid_channel=self.edbm_calc.acid_channel()
        # Assert the correctness of the output
        # Check if the outlet mass flow rates for each ion are calculated correctly
        expected_M_a_out_i = [0.0, 0.21, 0.0, 0.0, 0.0, 0.0, 0.0, 0.006, 1.9e-09]
        for i in range(9):
            # Assuming M_a_out and Ci_a_out are properly calculated in the method
            #self.edbm_calc.M_a_in[i] + self.edbm_calc.JA * self.edbm_calc.PM_i[i]
            self.assertAlmostEqual(self.edbm_calc.M_a_out[i], expected_M_a_out_i[i], delta=expected_M_a_out_i[i]*0.1)
        
        # Check if the total outlet mass flow rate is calculated correctly
        expected_M_a_out_t = sum(self.edbm_calc.M_a_out)+self.edbm_calc.M_h2o_a_out
        self.assertAlmostEqual(self.edbm_calc.M_a_out_t, expected_M_a_out_t, delta=expected_M_a_out_t*0.01)
        
        # Check if the volumetric outlet flow rate is calculated correctly
        expected_Q1_a_out = self.edbm_calc.M_a_out_t / 1  # Assuming the density is 1 kg/l
        self.assertAlmostEqual(self.edbm_calc.Q1_a_out, expected_Q1_a_out, delta=expected_Q1_a_out*0.01)
        
        # Check if the outlet concentration of single ions in the channel is calculated correctly
        for i in range(9):
            expected_Ci_a_out_i = self.edbm_calc.M_a_out[i] / (self.edbm_calc.Q1_a_out * self.edbm_calc.PM_i[i] / 1000)
            self.assertAlmostEqual(self.edbm_calc.Ci_a_out[i], expected_Ci_a_out_i, delta=0.01)


    def test_base_channel(self):
        results_flowrate=self.results_flowrate
        results_mass_flowrate=self.edbm_calc.in_mass_flow_rates(self.ph_s)
        results_acid_channel=self.edbm_calc.acid_channel()
        # Call the acid_channel method
        self.edbm_calc.base_channel()
        self.results_base_channel=self.edbm_calc.base_channel()
        # Assert the correctness of the output
        # Check if the outlet mass flow rates for each ion are calculated correctly
        expected_M_b_out_i = [0.14, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.9e-09, 0.10]
        for i in range(9):
            # Assuming M_b_out and Ci_b_out are properly calculated in the method            
            self.assertAlmostEqual(self.edbm_calc.M_b_out[i], expected_M_b_out_i[i], delta=0.01)
        
        # Check if the total outlet mass flow rate is calculated correctly
        expected_M_b_out_t = sum(self.edbm_calc.M_b_out)+self.edbm_calc.M_h2o_b_out
        self.assertAlmostEqual(self.edbm_calc.M_b_out_t, expected_M_b_out_t, delta=expected_M_b_out_t*0.01)
        
        # Check if the volumetric outlet flow rate is calculated correctly
        expected_Q1_b_out = self.edbm_calc.M_b_out_t / 1  # Assuming the density is 1 kg/l
        self.assertAlmostEqual(self.edbm_calc.Q1_b_out, expected_Q1_b_out, delta=expected_Q1_b_out*0.01)
        
        # Check if the outlet concentration of single ions in the channel is calculated correctly
        for i in range(9):
            expected_Ci_b_out_i = self.edbm_calc.M_b_out[i] / (self.edbm_calc.Q1_b_out * self.edbm_calc.PM_i[i] / 1000)
            self.assertAlmostEqual(self.edbm_calc.Ci_b_out[i], expected_Ci_b_out_i, delta=0.01)

    def test_salt_channel(self):
        results_flowrate=self.results_flowrate
        results_mass_flowrate=self.edbm_calc.in_mass_flow_rates(self.ph_s)
        results_acid_channel=self.edbm_calc.acid_channel()
        results_base_channel=self.edbm_calc.base_channel()
        #Membrane characteristics
        Cm_bp_H= 0.0000001 #mol/l 
        Cm_bp_OH= 0.0000001 #mol/l 
        # Call the salt_channel method
        self.edbm_calc.salt_channel(Cm_bp_H, Cm_bp_OH)
        
        # Assert the correctness of the output
        # Check if the outlet mass flow rates for each ion are calculated correctly
        expected_M_s_out_i = [0.13, 0.20,0.02, 0.0, 0.0, 0.004,0.004, 3.89e-07, 6.03e-13]
        for i in range(9):
            # Assuming M_s_out and Ci_s_out are properly calculated in the method
            self.assertAlmostEqual(self.edbm_calc.M_s_out[i], expected_M_s_out_i[i], delta=0.01)
        
        # Check if the total outlet mass flow rate is calculated correctly
        expected_M_s_out_t = sum(self.edbm_calc.M_s_out)+self.edbm_calc.M_h2o_s_out
        self.assertAlmostEqual(self.edbm_calc.M_s_out_t, expected_M_s_out_t, delta=expected_M_s_out_t*0.01)
        
        # Check if the volumetric outlet flow rate is calculated correctly
        expected_Q1_s_out = self.edbm_calc.M_s_out_t / self.edbm_calc.d_s  # Assuming the density is d_s kg/l
        self.assertAlmostEqual(self.edbm_calc.Q1_s_out, expected_Q1_s_out, delta=expected_Q1_s_out*0.01)
        
        # Check if the outlet concentration of single ions in the channel is calculated correctly
        for i in range(9):
            expected_Ci_s_out_i = self.edbm_calc.M_s_out[i] / (self.edbm_calc.Q1_s_out * self.edbm_calc.PM_i[i] / 1000)
            self.assertAlmostEqual(self.edbm_calc.Ci_s_out[i], expected_Ci_s_out_i, delta=0.01)


if __name__ == '__main__':
    unittest.main()
