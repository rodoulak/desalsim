#setup.py

from setuptools import setup, find_packages

with open("README.md", "r") as f: 
    description =f.read()

setup(
      name='desalsim',
      version='1.0.2',
      packages=find_packages(),
      url="https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-",
      author="rodoulak",
      author_email="r.ktori@tudelft.nl",
      license="MIT",
      install_requires=[
          "numpy==1.21.5",
          "pandas==1.4.4",  # Include pandas as a dependency
          "scipy==1.9.1", # For figures 
          "matplotlib==3.5.2", # For figures 
          "scikit-learn==1.0.2",
          "plotly>=5.9.0", 
          "openpyxl==3.0.10",  #read data from excel file
      ],
      long_description= description,
      long_description_content_type="text/markdown",
      )
