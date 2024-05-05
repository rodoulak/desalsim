import unittest
from mfpfr_unit_f import MFPFRCALC, HClAddition

class TestMFPFRCALC(unittest.TestCase):
    def setUp(self):
        # Define sample input data for initialization
        Qin = 1000
        Cin_mfpfr = [17.3, 38.6, 0.6, 4.3, 1.2, 9.9]
        C_NaOH_1 = 1
        C_NaOH_2 = 1
        conv_1 = 95
        conv_2 = 93
        
        # Create an instance of the MFPFRCALC class
        self.mfpfr_calc = MFPFRCALC(Qin, Cin_mfpfr, C_NaOH_1, C_NaOH_2, conv_1, conv_2)
        self.result = self.mfpfr_calc.calc_step1(5.61 * 0.000000000001, 2.34)
        
        
    def test_calc_step1(self):
        # Test the calc_step1 method
        
        # Call the calc_step1 method
        self.result = self.mfpfr_calc.calc_step1(5.61 * 0.000000000001, 2.34)
        
        # Assert the correctness of the output
        expected_value_QMgOH=9.80
        expected_value_QNaOH=336.14
        self.assertAlmostEqual(round(self.result[1],2), expected_value_QNaOH, delta=expected_value_QNaOH*0.01)
        self.assertAlmostEqual(round(self.result[2],2), expected_value_QMgOH, delta=expected_value_QMgOH*0.01)
        self.assertAlmostEqual(round(self.result[3],2), 1331.96, places=2)


    def test_calc_step2(self):
        # Test the calc_step2 method
        
        # Call the calc_step2 method
        self.result_2 = self.mfpfr_calc.calc_step2(2.34, 2.211)
        
        # Assert the correctness of the output

        expected_value_QNaOH_2=73.38
        self.assertAlmostEqual(round(self.result_2[1],2), expected_value_QNaOH_2, delta=expected_value_QNaOH_2*0.01)
        self.assertAlmostEqual(round(self.result_2[4],2), 1526.61, delta=1526.61*0.01)


class TestHClAddition(unittest.TestCase):
    def setUp(self):
        # Define sample input data for initialization
        Qout_2 = 1526.26
        Cout_all_m = [0.8, 0.7, 0.01, 0.0, 0.0, 0.07]
        MW_Cl = 35.45
        ph_2 = 13
        HCl_conc = 1
        
        # Create an instance of the HClAddition class
        self.hcl_addition = HClAddition(Qout_2, Cout_all_m, MW_Cl, ph_2, HCl_conc)

    def test_calculate_HCl_addition(self):
        # Test the calculate_HCl_addition method
        
        # Define additional input data
        Cout_mfpfr_g = [1, 2, 3, 4, 5, 6]
        
        # Call the calculate_HCl_addition method
        result = self.hcl_addition.calculate_HCl_addition(Cout_mfpfr_g)
        
        # Assert the correctness of the output
        expected_value_HCl=152.66
        self.assertAlmostEqual(round(result[0],2), expected_value_HCl, delta=expected_value_HCl*0.01)
        # Add more assertions as needed

if __name__ == '__main__':
    unittest.main()
