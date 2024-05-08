#setup.py

from setuptools import setup, find_packages

with open("README.md", "r") as f: 
    description =f.read()

setup(
      name='DesalSim',
      version='1.0.1',
      packages=find_packages(),
      url="https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-",
      author="rodoulak",
      author_email="r.ktori@tudelft.nl",
      license="MIT",
      install_requires=[
          "numpy>=1.0",
          "pandas>=1.0",  # Include pandas as a dependency
          "scipy>=1.0", # For figures 
          "matplotlib>=1.0", # For figures 
          "scikit-learn>=1.0",
          "plotly>=1.0", 
          "openpyxl>=1.0",  #read data from excel file
      ],
      long_description= description,
      long_description_content_type="text/markdown",
      )
