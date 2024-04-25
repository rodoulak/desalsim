project_directory/
│
├── constants.py
├── density_calc.py
│
├── med/
│   ├── med_calculator.py
│   ├── __init__.py (can be empty)
│
└── nfmass/
    ├── molarity.py
    ├── nfmass.py
    ├── nfenergy.py
    ├── osmotic_pressure.py
    ├── __init__.py (can be empty)



from med.med_calculator import MEDCalculator
from nfmass.nfmass import NFMass
from nfmass.molarity import molarity
from nfmass.nfenergy import NfEnergy
from nfmass.osmotic_pressure import OsmoticPressure

