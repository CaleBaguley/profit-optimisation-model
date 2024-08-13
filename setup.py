#!/usr/bin/env python

from setuptools import setup

packages = ['profit_optimisation_model',
            'profit_optimisation_model.src',
            'profit_optimisation_model.src.HydraulicConductanceModels',
            'profit_optimisation_model.src.HydraulicConductanceModels.DynamicModels',
            'profit_optimisation_model.src.PhotosynthesisModels',
            'profit_optimisation_model.src.ProfitModels',
            'profit_optimisation_model.src.ProfitModels.CO2GainModels',
            'profit_optimisation_model.src.ProfitModels.HydraulicCostModels',
            'profit_optimisation_model.src.TemperatureDependenceModels'
            ]

setup(name='Stomatal Optimisation Models',
      version='0.1',
      description='An implementation of different stomatal optimisation models.',
      author='Cale Baguley',
      url='https://github.com/CaleBaguley/profit-optimisation-model',
      packages=packages
      )

