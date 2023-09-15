"""
-------------------------------------------------------------------------

-------------------------------------------------------------------------
"""

from src.hydraulic_cost_model import HydraulicCostModel
from src.leaf_air_coupling_model import LeafAirCouplingModel
from src.CO2_gain_model import CO2GainModel

from numpy import zeros


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
        hydraulic_costs = \
            self._hydraulic_cost_model.hydraulic_cost_as_a_function_of_leaf_water_potential(leaf_water_potentials,
                                                                                            soil_water_potential)

        transpiration_as_a_function_of_leaf_water_potential = zeros(len(leaf_water_potentials))

        for i in range(len(leaf_water_potentials)):
            transpiration_as_a_function_of_leaf_water_potential[i] = \
                self._hydraulic_cost_model.transpiration(leaf_water_potentials[i],
                                                         soil_water_potential)

        CO2_gain = self._CO2_gain_model.CO2_gain(transpiration_as_a_function_of_leaf_water_potential,
                                                 air_temperature,
                                                 air_vapour_pressure_deficit,
                                                 air_pressure,
                                                 atmospheric_CO2_concentration,
                                                 intercellular_oxygen)

        return CO2_gain - hydraulic_costs, CO2_gain, hydraulic_costs
