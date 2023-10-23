"""
-------------------------------------------------------------------------

-------------------------------------------------------------------------
"""

from numpy import zeros, max, nan_to_num, where, nanmax
from src.leaf_air_coupling_model import LeafAirCouplingModel
from src.photosynthesis_model import PhotosynthesisModelDummy
from src.conversions import magnitude_conversion


class CO2GainModel:
    _leaf_air_coupling_model: LeafAirCouplingModel
    _photosynthesis_model: PhotosynthesisModelDummy

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
                 intercellular_O = None,
                 photosyntheticaly_active_radiation = None):
        """

        @param transpiration_rates: mmol m-2 s-1
        @param air_temperature: K
        @param air_vapour_pressure_deficit: kPa
        @param air_pressure: kPa
        @param atmospheric_CO2_concentration: umol mol-1
        @param intercellular_O: umol mol-1
        @param photosyntheticaly_active_radiation: umol m-2 s-1

        @return: normalised CO2 gain
        @return: maximum CO2 gain: umol m-2 s-1
        """

        net_CO2_uptake = zeros(len(transpiration_rates))
        intercellular_CO2_as_a_function_of_leaf_water_potential = zeros(len(transpiration_rates))
        stomatal_conductance_to_CO2_as_a_function_of_leaf_water_potential = zeros(len(transpiration_rates))

        for i in range(len(transpiration_rates)):
            stomatal_conductance_to_CO2 = \
                self._leaf_air_coupling_model.stomatal_conductance_to_carbon(transpiration_rates[i],
                                                                             air_temperature,
                                                                             air_vapour_pressure_deficit,
                                                                             air_pressure)

            # Need to convert from mmol m-2 s-1 to mol m-2 s-1
            stomatal_conductance_to_CO2 = magnitude_conversion(stomatal_conductance_to_CO2, 'm', '')

            (net_CO2_uptake[i],
            intercellular_CO2_as_a_function_of_leaf_water_potential[i]) = \
                self._photosynthesis_model.net_rate_of_CO2_assimilation(stomatal_conductance_to_CO2,
                                                                        atmospheric_CO2_concentration,
                                                                        air_temperature,
                                                                        intercellular_O,
                                                                        photosyntheticaly_active_radiation)

            stomatal_conductance_to_CO2_as_a_function_of_leaf_water_potential[i] = stomatal_conductance_to_CO2

        maximum_CO2_uptake = nanmax(net_CO2_uptake)

        #print("max CO2 uptake: ",maximum_CO2_uptake)

        if(maximum_CO2_uptake > 0.):
            CO2_gain = net_CO2_uptake / maximum_CO2_uptake

        else:
            CO2_gain = zeros(len(net_CO2_uptake))

        return (CO2_gain,
                maximum_CO2_uptake,
                intercellular_CO2_as_a_function_of_leaf_water_potential,
                stomatal_conductance_to_CO2_as_a_function_of_leaf_water_potential)
