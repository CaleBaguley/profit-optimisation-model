"""
-------------------------------------------------------------------------

-------------------------------------------------------------------------
"""

from src.ProfitModels.CO2GainModels.CO2_gain_model import CO2GainModelDummy


class SOXCO2GainModel(CO2GainModelDummy):

    def gain_equation(self, net_CO2_uptake):
        return net_CO2_uptake
