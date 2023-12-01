from setuptools import find_packages, setup

setup(
    name="nps-classifier",
    packages=find_packages(include=["nps"]),
    version="1.0",
    description="Module for substances classification as NPS",
    author='PD',
    install_requires=[
      "rdkit==2023.3.1",
      "matplotlib==3.8.0",
      "numpy==1.25.0rc1",
    ]
)
