"""
-------------------------------------------------------------------------

-------------------------------------------------------------------------
"""

from numpy import zeros, max
from src.leaf_air_coupling_model import LeafAirCouplingModel
from src.photosynthesis_model import PhotosynthesisModel

class CO2GainModel:

    _leaf_air_coupling_model: LeafAirCouplingModel
    _photosynthesis_model: PhotosynthesisModel

    def __init__(self,
                 leaf_air_coupling_model,
                 photosynthesis_model):

        self._leaf_air_coupling_model = leaf_air_coupling_model
        self._photosynthesis_model = photosynthesis_model

    def CO2_gain(self,
                 transpiration_rates,
                 air_temperature,
                 air_vapour_pressure_deficit,
                 air_pressure,
                 atmospheric_CO2_concentration,
                 intercellular_O):

        """

        @param transpiration_rates: mmol m-2 s-1
        @param air_temperature: K
        @param air_vapour_pressure_deficit: kPa
        @param air_pressure: kPa
        @param atmospheric_CO2_concentration: umol mol-1
        @param intercellular_O: umol mol-1
        @return:
        """

        net_CO2_uptake = zeros(len(transpiration_rates))

        for i in range(len(transpiration_rates)):
            stomatal_conductance_to_CO2 = \
                self._leaf_air_coupling_model.stomatal_conductance_to_carbon(transpiration_rates[i],
                                                                             air_temperature,
                                                                             air_vapour_pressure_deficit,
                                                                             air_pressure)

            # Need to convert from mmol m-2 s-1 to mol m-2 s-1
            stomatal_conductance_to_CO2 /= 1000

            net_CO2_uptake[i] = \
                self._photosynthesis_model.net_rate_of_CO2_assimilation(stomatal_conductance_to_CO2,
                                                                        atmospheric_CO2_concentration,
                                                                        air_temperature,
                                                                        intercellular_O)

        maximum_CO2_uptake = max(net_CO2_uptake)

        return net_CO2_uptake / maximum_CO2_uptake
