from setuptools import find_packages, setup

setup(
    name="nps-classifier",
    packages=find_packages(include=["nps"]),
    version="1.0",
    description="Module for substances classification as NPS",
    author='PD',
    install_requires=[
      "rdkit==2022.9.5",
      "matplotlib==3.7.2",
      "numpy==1.26.1",
    ]
)
