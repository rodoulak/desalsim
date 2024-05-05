import unittest
from  economic_f import econom, revenue
import density_calc
import constants
import scaleup

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
        
        # Create an instance of the econom class
        self.econom_calc = econom(eq_c, el_conc, s_conc, chem1_conc, chem1_pr, chem2_conc, chem2_pr, cw_conc, wat_conc)
        self.result_capex = self.econom_calc.capex_calc()
        
    def test_capex_calc(self):
        capex_result = self.result_capex
        
        # Call the capex_calc method
        self.result_capex = self.econom_calc.capex_calc()
        
        # Assert the correctness of the output
        expected_capex = 3125700  # Assuming this value based on your calculations
        self.assertAlmostEqual(self.result_capex[0], expected_capex, delta=1000)


    def test_opex_calc(self):
        # Access the stored result from test_capex_calc
        capex_result = self.result_capex
        
        # Call the opex_calc method
        result_opex = self.econom_calc.opex_calc(7200, 0.253, 0, 0.000118, 0.001, [0.03, 0.05, 0.15, 0.15, 0.15, 0.03, 0.05, 0.05, 0.67])
        
        # Assert the correctness of the output
        expected_opex = 148541.66  # Assuming this value based on your calculations
        self.assertAlmostEqual(result_opex, expected_opex, delta=1000)


class TestRevenue(unittest.TestCase):
    def setUp(self):
        # Define sample input data for initialization
        prd=10  
        prd_name= "Water"  
        
        # Create an instance of the revenue class
        self.revenue_calc = revenue(prd, prd_name)

    def test_rev_water(self):
        # Call the rev method with Water as product
        result_rev = self.revenue_calc.rev(7200, 0.001, 0.066, 1.0, 0.12, 7.2, 5.78)
        
        # Assert the correctness of the output
        expected_rev = 72.0  # Assuming this value based on your calculations
        self.assertAlmostEqual(result_rev, expected_rev, delta=1)


if __name__ == '__main__':
    unittest.main()
