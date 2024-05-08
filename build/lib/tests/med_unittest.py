import unittest
from Desalsim.med_unit_f import MEDCalculator, density_calc

class TestMEDCalculator(unittest.TestCase):
    def setUp(self):
        # Define the input parameters for testing
        Cin_med = [10.36, 15.39, 0.36, 0.028, 0.02, 0.07]
        Qf_med = 1000
        Mf_med = 174.475 # Assuming a density of 1 g/ml
        
        # Add any other necessary parameters for initialization
        T_in = 40
        
        # Create an instance of the MEDCalculator class with the defined parameters
        self.med_dat = MEDCalculator(Qf_med, Mf_med, Cin_med[0], Cin_med[1], Cin_med[2], Cin_med[3], Cin_med[4], Cin_med[5], T_in)

    def test_salinity_calc(self):
        # Call the salinity_calc method
        self.med_dat.salinity_calc()
        
        # Define the expected values
        expected_salinity_in = sum([10.36, 15.39, 0.36, 0.028, 0.02, 0.07])
        expected_xf = expected_salinity_in / density_calc(40, expected_salinity_in) * 1000
        expected_Mf = 174.475 / 3600 # Converting from kg/hr to kg/s
        
        # Assert the equality of the calculated values and expected values
        self.assertAlmostEqual(self.med_dat.salinity_in, expected_salinity_in, delta=expected_salinity_in*0.01)
        self.assertAlmostEqual(self.med_dat.xf, expected_xf, delta=expected_xf*0.01)
        self.assertAlmostEqual(self.med_dat.Mf, expected_Mf, delta=expected_Mf*0.01)

    # Add more test methods to cover other methods in the MEDCalculator class

if __name__ == '__main__':
    unittest.main()
