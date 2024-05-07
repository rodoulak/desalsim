#setup.py

from setuptools import setup, find_packages

with open("README.md", "r") as f: 
    description =f.read()

setup(
      name='desalination_and_brine_treatment',
      version='0.2',
      packages=find_packages(),
      url="https://github.com/rodoulak/Desalination-and-Brine-Treatment-Simulation-",
      author="rodoulak",
      author_email="r.ktori@tudelft.nl",
      license="MIT",
      install_requires=[
          "numpy>=1.0",
          "pandas>=1.0",  # Include pandas as a dependency
      ],
      long_description= description,
      long_description_content_type="text/markdown",
      )