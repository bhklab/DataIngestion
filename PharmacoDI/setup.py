from setuptools import setup, find_packages

setup(name='PharmacoDI',
      version='0.0.1',
      description="Tools for processing R PharmacoSet objects into .csv files of PharmacoDB database tables.",
      url='https://github.com/bhklab/DataIngestion/tree/master/PharmacoDI',
      author='Christopher Eeles, Benjamin Haibe-Kains',
      author_email='christopher.eeles@uhnresearch.ca, benjamin.haibe.kains@utoronto.ca',
      license='MIT',
      packages=find_packages(),
      zip_safe=False
      )