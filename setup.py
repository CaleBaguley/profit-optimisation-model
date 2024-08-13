#!/usr/bin/env python

from setuptools import setup

packages = ['profit-optimisation-model',
            'profit-optimisation-model.src',
            'profit-optimisation-model.src.HydraulicConductanceModels',
            'profit-optimisation-model.src.HydraulicConductanceModels.DynamicModels',
            'profit-optimisation-model.src.PhotosynthesisModels',
            'profit-optimisation-model.src.ProfitModels',
            'profit-optimisation-model.src.ProfitModels.CO2GainModels',
            'profit-optimisation-model.src.ProfitModels.HydraulicCostModels',
            'profit-optimisation-model.src.ProfitModels.TemperatureModels'
            ]

setup(name='Stomatal Optimisation Models',
      version='0.1',
      description='An implementation of different stomatal optimisation models.',
      author='Cale Baguley',
      url='https://github.com/CaleBaguley/profit-optimisation-model',
      packages=packages
      )

