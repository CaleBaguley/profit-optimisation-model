"""
-------------------------------------------------------------------------
Implementation of Manon's Profit Max model (M. Sabot et al 2019)
-------------------------------------------------------------------------
"""


from src.ProfitModels.profit_optimisation_model import ProfitOptimisationModel

from numpy import zeros, linspace


class ProfitMaxModel(ProfitOptimisationModel):

    def profit(self, CO2_gain, hydraulic_cost):
        return CO2_gain - hydraulic_cost
