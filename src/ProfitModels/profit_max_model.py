"""
-------------------------------------------------------------------------

-------------------------------------------------------------------------
"""

from src.ProfitModels.hydraulic_cost_model import HydraulicCostModel
from src.leaf_air_coupling_model import LeafAirCouplingModel
from src.ProfitModels.CO2_gain_model import CO2GainModel
from src.ProfitModels.profit_optimisation_model import ProfitOptimisationModel

from numpy import zeros, linspace
from numpy import argwhere, nanargmax


class ProfitMaxModel(ProfitOptimisationModel):
    _hydraulic_cost_model: HydraulicCostModel
    _leaf_air_coupling_model: LeafAirCouplingModel
    _CO2_gain_model: CO2GainModel

    def __init__(self,
                 hydraulic_cost_model,
                 leaf_air_coupling_model,
                 CO2_gain_model):

        super().__init__()
        self._hydraulic_cost_model = hydraulic_cost_model
        self._leaf_air_coupling_model = leaf_air_coupling_model
        self._CO2_gain_model = CO2_gain_model

    def profit_as_a_function_of_leaf_water_potential(self,
                                                     soil_water_potential,
                                                     air_temperature,
                                                     air_vapour_pressure_deficit,
                                                     air_pressure,
                                                     atmospheric_CO2_concentration,
                                                     intercellular_oxygen,
                                                     photosynthetically_active_radiation,
                                                     number_of_sample_points = 1000):
        """

        @param soil_water_potential: MPa
        @param air_temperature: K
        @param air_vapour_pressure_deficit: kPa
        @param air_pressure: kPa
        @param atmospheric_CO2_concentration: umol mol-1
        @param intercellular_oxygen: umol mol-1
        @param photosynthetically_active_radiation: umol m-2 s-1
        @param number_of_sample_points: 1000

        @return: profit
        @return: normalised CO2 gain
        @return: hydraulic cost
        @return: maximum CO2 uptake: umol m-2 s-1
        @return: transpiration: mmol m-2 s-1
        @return: leaf_water_potentials: kPa
        """

        critical_leaf_water_potential = self._hydraulic_cost_model.critical_leaf_water_potential

        leaf_water_potentials = linspace(soil_water_potential,
                                         critical_leaf_water_potential,
                                         num=number_of_sample_points)

        hydraulic_costs = \
            self._hydraulic_cost_model.hydraulic_cost_as_a_function_of_leaf_water_potential(leaf_water_potentials,
                                                                                            soil_water_potential)

        transpiration_as_a_function_of_leaf_water_potential = zeros(len(leaf_water_potentials))

        for i in range(len(leaf_water_potentials)):
            transpiration_as_a_function_of_leaf_water_potential[i] = \
                self._hydraulic_cost_model.transpiration(leaf_water_potentials[i],
                                                         soil_water_potential)

        (CO2_gain,
         maximum_CO2_uptake,
         intercellular_CO2_as_a_function_of_leaf_water_potential,
         stomatal_conductance_to_CO2_as_a_function_of_leaf_water_potential) = \
            self._CO2_gain_model.CO2_gain(transpiration_as_a_function_of_leaf_water_potential,
                                          air_temperature,
                                          air_vapour_pressure_deficit,
                                          air_pressure,
                                          atmospheric_CO2_concentration,
                                          intercellular_oxygen,
                                          photosynthetically_active_radiation)

        return (CO2_gain - hydraulic_costs,
                CO2_gain,
                hydraulic_costs,
                maximum_CO2_uptake,
                transpiration_as_a_function_of_leaf_water_potential,
                intercellular_CO2_as_a_function_of_leaf_water_potential,
                stomatal_conductance_to_CO2_as_a_function_of_leaf_water_potential,
                leaf_water_potentials)
