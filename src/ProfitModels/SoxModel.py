"""
-------------------------------------------------------------------------
Implementation of SOX model (C. B. Eller et al 2018)
-------------------------------------------------------------------------
"""

from src.ProfitModels.HydraulicCostModels.hydraulic_cost_model import HydraulicCostModel
from src.leaf_air_coupling_model import LeafAirCouplingModel
from src.ProfitModels.CO2GainModels.CO2_gain_model import CO2GainModelDummy
from src.ProfitModels.profit_optimisation_model import ProfitOptimisationModel

from numpy import zeros, linspace


class SOXModel(ProfitOptimisationModel):

    def profit(self, CO2_gain, hydraulic_cost):
        return CO2_gain * hydraulic_cost
