"""
-------------------------------------------------------------------------
Implementation of SOX model (C. B. Eller et al 2018)
-------------------------------------------------------------------------
"""

from src.ProfitModels.profit_optimisation_model import ProfitOptimisationModel


class SOXModel(ProfitOptimisationModel):

    def profit(self, CO2_gain, hydraulic_cost):
        return CO2_gain * (1 - hydraulic_cost)
