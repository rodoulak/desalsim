import unittest
from  Desalsim.economic_f import econom, revenue
from  Desalsim.density_calc import density_calc
from Desalsim import constants
from Desalsim import scaleup

class TestEconom(unittest.TestCase):
    def setUp(self):
        # Define sample input data for initialization
        eq_c = 1500000
        el_conc = 10
        s_conc = 0
        chem1_conc=10
        chem1_pr=0.002
        chem2_conc=0
        chem2_pr=0.002
        cw_conc = 0
        wat_conc = 0
        capex_assumptions=[0.25, 0.2, 0.06, 0.15, 0.2]
        
        prd=10  
        prd_name= "Water" 
        
        # Create an instance of the econom class
        self.econom_calc = econom(eq_c, el_conc, s_conc, chem1_conc, chem1_pr, chem2_conc, chem2_pr, cw_conc, wat_conc)
        self.result_capex = self.econom_calc.capex_calc(capex_assumptions)
        self.revenue_calc = revenue(prd, prd_name)
        
    def test_capex_calc(self):
        capex_result = self.result_capex
        
        # Set CAPEX assumptions 
        capex_assumptions=[0.25, 0.2, 0.06, 0.15, 0.2]

        # Call the capex_calc method
        self.result_capex = self.econom_calc.capex_calc(capex_assumptions)
        
        # Assert the correctness of the output
        expected_capex = 3125700  # Assuming this value based on your calculations
        self.assertAlmostEqual(self.result_capex, expected_capex, delta=expected_capex*0.01)


    def test_opex_calc(self):
        # Access the stored result from test_capex_calc
        capex_result = self.result_capex
        
        # Call the opex_calc method
        result_opex = self.econom_calc.opex_calc(7200, 0.253, 0, 0.000118, 0.001, [0.03, 0.05, 0.15, 0.15, 0.15, 0.03, 0.05, 0.05, 0.67])
        
        # Assert the correctness of the output
        expected_opex = 148541.66  # Assuming this value based on your calculations
        self.assertAlmostEqual(result_opex, expected_opex, delta=expected_opex*0.01)




if __name__ == '__main__':
    unittest.main()
