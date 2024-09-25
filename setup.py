#setup.py

from setuptools import setup, find_packages

with open("README.md", "r") as f: 
    description =f.read()

setup(
      name='desalsim',
      version='1.0.9.3',
      packages=find_packages(),
      url="https://github.com/rodoulak/desalsim.git",
      author="rodoulak",
      author_email="r.ktori@tudelft.nl",
      license="MIT",
      install_requires=[
          "numpy", #==1.26.4
          "pandas",  # Include pandas as a dependency ==2.2.2
          "scipy", # For figures ==1.13.0
          "matplotlib", # For figures ==3.8.4
          "scikit-learn", #==1.4.2
          "plotly", #==5.22.0
          "openpyxl",  #read data from excel file ==3.1.2
      ],
      long_description= description,
      long_description_content_type="text/markdown",
      )
