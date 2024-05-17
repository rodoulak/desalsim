import unittest
from desalsim import constants
from desalsim.ed_unit_f import ElectrodialysisCalc

# Define the input parameters for testing
C_feed = [10.36, 15.39, 0.36, 0.028, 0.02, 0.07]  # Feed concentrations in g/L
Q_feed = 1000  # Feed flow rate in L/hr
T_feed = 40  # Feed temperature in °C
T_dial = 25  # Dialysate temperature in °C
A_membrane = 10  # Membrane area in m^2
V_cell = 100  # Voltage across the cell in V
I_current = 50  # Current passing through the cell in A

class TestElectrodialysisCalc(unittest.TestCase):
    def setUp(self):        
        # Create an instance of the ElectrodialysisCalc class with the defined parameters
        self.ed_calc = ElectrodialysisCalc()
        
    def test_Ts_cp(self):
        # Call the Ts_cp method
        result = ElectrodialysisCalc.Ts_cp(10)
        
        # Define the expected value
        expected_result = 0.9676  # Expected result for S=10
        
        # Assert the equality of the calculated value and expected value
        self.assertAlmostEqual(result, expected_result, delta=expected_result*0.01)

    def test_Tw_cp(self):
        # Call the Tw_cp method
        Sc=43
        Sd=43
        result = ElectrodialysisCalc.Tw_cp(Sc, Sd)
        
        # Define the expected value
        expected_result = 10.3  
        
        # Assert the equality of the calculated value and expected value
        self.assertAlmostEqual(result, expected_result, delta=expected_result*0.01)
        
    def test_Lw_cp(self):
        # Call the Lw_cp method
        S = 43
        result = ElectrodialysisCalc.Lw_cp(S)
        
        # Define the expected value
        expected_result = 5 * S ** (-0.416)
        
        # Assert the equality of the calculated value and expected value
        self.assertAlmostEqual(result, expected_result, delta=expected_result*0.01)

    def test_p_osmo(self):
        # Call the p_osmo method
        S = 43
        T = 25
        MWs = constants.MW_NaCl
        
        result = ElectrodialysisCalc.p_osmo(S, T, MWs)
        
        # Define the expected value
        C1 = S / MWs * constants.MW_Na / constants.MW_Na
        C2 = S / MWs * constants.MW_cl / constants.MW_cl
        sum_Ci = C1 + C2
        expected_result = 0.0831446261815324 * sum_Ci * T
        
        # Assert the equality of the calculated value and expected value
        self.assertAlmostEqual(result, expected_result, delta=expected_result*0.01)

    def test_dC(self):
        # Define parameters
        Ts_cp = 0.95
        tcu = 0.5
        D = 1.61e-9 #Diffusion coefficient (m^2/s)
        Ij = 50
        h = 0.5
        Sh = 18

        
        # Call the dC method
        result = ElectrodialysisCalc.dC(Ts_cp, tcu, D, Ij, h, Sh)
        
        # Define the expected value
        F = 96485.3329  # Faraday constant (C/mol)
        Tcu = (Ts_cp + 1) / 2
        expected_result = -(Tcu - tcu) / D * (Ij / F) * (2 * h / 1000 / Sh)
        
        # Assert the equality of the calculated value and expected value
        self.assertAlmostEqual(result, expected_result, delta=expected_result*0.01)


if __name__ == '__main__':
    unittest.main()
