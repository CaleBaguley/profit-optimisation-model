"""
-------------------------------------------------------------------------

-------------------------------------------------------------------------
"""

from src.hydraulic_cost_model import HydraulicCostModel
from src.leaf_air_coupling_model import LeafAirCouplingModel
from src.CO2_gain_model import CO2GainModel

from numpy import zeros, linspace
from numpy import argmax


class ProfitOptimisationModel:
    _hydraulic_cost_model: HydraulicCostModel
    _leaf_air_coupling_model: LeafAirCouplingModel
    _CO2_gain_model: CO2GainModel

    def __init__(self,
                 hydraulic_cost_model,
                 leaf_air_coupling_model,
                 CO2_gain_model):
        self._hydraulic_cost_model = hydraulic_cost_model
        self._leaf_air_coupling_model = leaf_air_coupling_model
        self._CO2_gain_model = CO2_gain_model

    def profit_as_a_function_of_leaf_water_potential(self,
                                                     leaf_water_potentials,
                                                     soil_water_potential,
                                                     air_temperature,
                                                     air_vapour_pressure_deficit,
                                                     air_pressure,
                                                     atmospheric_CO2_concentration,
                                                     intercellular_oxygen):
        """

        @param leaf_water_potentials: MPa
        @param soil_water_potential: MPa
        @param air_temperature: K
        @param air_vapour_pressure_deficit: kPa
        @param air_pressure: kPa
        @param atmospheric_CO2_concentration: umol mol-1
        @param intercellular_oxygen: umol mol-1

        @return: profit
        @return: normalised CO2 gain
        @return: hydraulic cost
        @return: maximum CO2 uptake: umol m-2 s-1
        @return: transpiration: mmol m-2 s-1
        """

        hydraulic_costs = \
            self._hydraulic_cost_model.hydraulic_cost_as_a_function_of_leaf_water_potential(leaf_water_potentials,
                                                                                            soil_water_potential)

        transpiration_as_a_function_of_leaf_water_potential = zeros(len(leaf_water_potentials))

        for i in range(len(leaf_water_potentials)):
            transpiration_as_a_function_of_leaf_water_potential[i] = \
                self._hydraulic_cost_model.transpiration(leaf_water_potentials[i],
                                                         soil_water_potential)

        CO2_gain, maximum_CO2_uptake = \
            self._CO2_gain_model.CO2_gain(transpiration_as_a_function_of_leaf_water_potential,
                                          air_temperature,
                                          air_vapour_pressure_deficit,
                                          air_pressure,
                                          atmospheric_CO2_concentration,
                                          intercellular_oxygen)

        return (CO2_gain - hydraulic_costs,
                CO2_gain,
                hydraulic_costs,
                maximum_CO2_uptake,
                transpiration_as_a_function_of_leaf_water_potential)

    def optimal_state(self,
                      soil_water_potential,
                      air_temperature,
                      air_vapour_pressure_deficit,
                      air_pressure,
                      atmospheric_CO2_concentration,
                      intercellular_oxygen,
                      number_of_sample_points = 1000):
        """
        Uses profit optimisation to calculate the optimal leaf water potential.
        @param soil_water_potential: MPa
        @param air_temperature: K
        @param air_vapour_pressure_deficit: kPa
        @param air_pressure: kPa
        @param atmospheric_CO2_concentration: umol mol-1
        @param intercellular_oxygen: umol mol-1
        @param number_of_sample_points: Number of leaf water potentials to test

        @return: optimal leaf water potential(MPa)
        @return: net CO2 uptake: umol m-2 s-1
        @return: transpiration: mmol m-2 s-1
        """

        critical_leaf_water_potential = self._hydraulic_cost_model.critical_hydraulic_conductance

        leaf_water_potentials = linspace(soil_water_potential,
                                         critical_leaf_water_potential,
                                         num = number_of_sample_points)

        (profit, CO2_gain, hydraulic_costs, maximum_net_CO2_uptake,
         transpiration_as_a_function_of_leaf_water_potential) = \
            self.profit_as_a_function_of_leaf_water_potential(leaf_water_potentials,
                                                              soil_water_potential,
                                                              air_temperature,
                                                              air_vapour_pressure_deficit,
                                                              air_pressure,
                                                              atmospheric_CO2_concentration,
                                                              intercellular_oxygen)

        maximum_profit_id = argmax(profit)

        optimal_leaf_water_potential = leaf_water_potentials[maximum_profit_id]
        net_CO2_uptake = maximum_net_CO2_uptake * CO2_gain[maximum_profit_id]
        transpiration_rate = transpiration_as_a_function_of_leaf_water_potential[maximum_profit_id]

        return optimal_leaf_water_potential, net_CO2_uptake, transpiration_rate
