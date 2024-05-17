import unittest
from desalsim.nanofiltration_unit_f import NFMass, OsmoticPressure, NfEnergy

class TestNanofiltrationUnit(unittest.TestCase):

    def test_nfmass(self):
        # Test NFMass calculation
        components = ['Na', 'Cl', 'K', 'Mg', 'Ca', 'SO4']
        self.C_in = [12.33, 21.67, 0.45, 1.39, 0.45, 3.28]
        rjr_values = [0.16, 0.29, 0.21, 0.98, 0.95, 0.98] #Ions rejection rates based on membrane characteristics (units: -)
        Wrec = 0.7 # Water recovery based on membrane characteristics (units: -)
        Qf=1000
        
        nf_mass = NFMass(components[0], self.C_in[0], rjr_values[0], Wrec, Qf)  # Example values
        self.assertAlmostEqual(nf_mass.Qperm, 700, delta=10)  # Expected Qperm: 875 kg/hr
        self.assertAlmostEqual(nf_mass.Qconc, 300, delta=10)  # Expected Qconc: 2625 kg/hr
        # Add more tests for other attributes if needed

    def test_osmotic_pressure(self):
        # Test OsmoticPressure calculation
        self.C_in = [12.33, 21.67, 0.45, 1.39, 0.45, 3.28]
        self.z_values = [1, -1, 1, 2, 2, -2]
        osmotic_pressure = OsmoticPressure(self.C_in, self.z_values, 293)  # Example values
        self.assertAlmostEqual(osmotic_pressure.osmotic_pressure_calculation(), 32.86, delta=1)  # Expected osmotic pressure: 22.56 bar

    def test_nfenergy(self):
        # Test NfEnergy calculation
        nf_energy = NfEnergy(57.8, 32.86, 22.9, 2, 1018.1, 700.0, 1000,1028.28, 0.8)  # Example values
        result = nf_energy.calculate_energy_consumption()
        # Add assertions to test the result dictionary
        self.assertAlmostEqual(result['Applied pressure (Bar)'], 24.4, delta=2)  # Expected applied pressure: 12.56 bar


if __name__ == '__main__':
    unittest.main()
