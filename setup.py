#setup.py

from setuptools import setup, find_packages

with open("README.md", "r") as f: 
    description =f.read()

setup(
      name='desalsim',
      version='1.0.5',
      packages=find_packages(),
      url="https://github.com/rodoulak/desalsim.git",
      author="rodoulak",
      author_email="r.ktori@tudelft.nl",
      license="MIT",
      install_requires=[
          "numpy==1.26.4",
          "pandas==2.2.2",  # Include pandas as a dependency
          "scipy==1.13.0", # For figures 
          "matplotlib==3.8.4", # For figures 
          "scikit-learn==1.4.2",
          "plotly==5.22.0", 
          "openpyxl==3.1.2",  #read data from excel file
      ],
      long_description= description,
      long_description_content_type="text/markdown",
      )
