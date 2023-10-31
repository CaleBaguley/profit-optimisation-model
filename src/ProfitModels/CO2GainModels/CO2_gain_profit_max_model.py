"""
-------------------------------------------------------------------------

-------------------------------------------------------------------------
"""

from numpy import zeros, nanmax
from src.ProfitModels.CO2GainModels.CO2_gain_model import CO2GainModelDummy

class ProfitMaxCO2GainModel(CO2GainModelDummy):

    def gain_equation(self, net_CO2_uptake):
        maximum_CO2_uptake = nanmax(net_CO2_uptake)

        if(maximum_CO2_uptake > 0.):
            return net_CO2_uptake/maximum_CO2_uptake

        return zeros(len(net_CO2_uptake))
