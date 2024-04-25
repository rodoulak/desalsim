#setup.py

from setuptools import setup, find_packages

setup(
      name='desalination_and_brine_treatment',
      version='0.1',
      packages=find_packages(),
      install_requires=[
         'numpy>=1.0',  # Example dependency
        'other_package>=2.0',  # Another example dependency 
        ],
      )