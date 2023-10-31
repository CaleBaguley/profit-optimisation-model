"""
-------------------------------------------------------------------------

-------------------------------------------------------------------------
"""

from numpy import zeros
from numpy import argwhere, nanargmax


class ProfitOptimisationModel:

    def __init__(self):
        pass

    def profit_as_a_function_of_leaf_water_potential(self,
                                                     soil_water_potential,
                                                     air_temperature,
                                                     air_vapour_pressure_deficit,
                                                     air_pressure,
                                                     atmospheric_CO2_concentration,
                                                     intercellular_oxygen,
                                                     photosynthetically_active_radiation,
                                                     number_of_sample_points=1000):
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

        raise Exception("profit_as_a_function_of_leaf_water_potential method not implemented in"
                        " ProfitOptimisationModel class")

    def optimal_state(self,
                      soil_water_potential,
                      air_temperature,
                      air_vapour_pressure_deficit,
                      air_pressure,
                      atmospheric_CO2_concentration,
                      intercellular_oxygen,
                      photosynthetically_active_radiation,
                      number_of_sample_points=1000):
        """
        Uses profit optimisation to calculate the optimal leaf water potential.
        @param soil_water_potential: MPa
        @param air_temperature: K
        @param air_vapour_pressure_deficit: kPa
        @param air_pressure: kPa
        @param atmospheric_CO2_concentration: umol mol-1
        @param intercellular_oxygen: umol mol-1
        @param photosynthetically_active_radiation: umol m-2 s-1
        @param number_of_sample_points: Number of leaf water potentials to test

        @return: optimal leaf water potential(MPa)
        @return: net CO2 uptake: umol m-2 s-1
        @return: transpiration: mmol m-2 s-1
        """

        (profit, CO2_gain,
         hydraulic_costs,
         maximum_net_CO2_uptake,
         transpiration_as_a_function_of_leaf_water_potential,
         intercellular_CO2_as_a_function_of_leaf_water_potential,
         stomatal_conductance_to_CO2_as_a_function_of_leaf_water_potential,
         leaf_water_potentials) = \
            self.profit_as_a_function_of_leaf_water_potential(soil_water_potential,
                                                              air_temperature,
                                                              air_vapour_pressure_deficit,
                                                              air_pressure,
                                                              atmospheric_CO2_concentration,
                                                              intercellular_oxygen,
                                                              photosynthetically_active_radiation,
                                                              number_of_sample_points)

        # Limit Ci/Ca to < 0.95
        realistic_intercellular_CO2_concentration_arg = argwhere(intercellular_CO2_as_a_function_of_leaf_water_potential
                                                                 < 0.95*atmospheric_CO2_concentration)
        maximum_profit_id = 0
        if(len(realistic_intercellular_CO2_concentration_arg) > 0):
            maximum_profit_id = nanargmax(profit[realistic_intercellular_CO2_concentration_arg])


        maximum_profit_id_old = nanargmax(profit)

        optimal_leaf_water_potential = leaf_water_potentials[maximum_profit_id]
        net_CO2_uptake = maximum_net_CO2_uptake * CO2_gain[maximum_profit_id]
        transpiration_rate = transpiration_as_a_function_of_leaf_water_potential[maximum_profit_id]
        intercellular_CO2 = intercellular_CO2_as_a_function_of_leaf_water_potential[maximum_profit_id]
        stomatal_conductance_to_CO2 = \
            stomatal_conductance_to_CO2_as_a_function_of_leaf_water_potential[maximum_profit_id]

        return (optimal_leaf_water_potential,
                net_CO2_uptake,
                transpiration_rate,
                intercellular_CO2,
                stomatal_conductance_to_CO2)


def run_optimisation_model_on_data(profit_optimisation_model: ProfitOptimisationModel,
                                   time_steps,
                                   soil_water_potential_values,
                                   air_temperature_values,
                                   air_vapour_pressure_deficit_values,
                                   air_pressure_values,
                                   atmospheric_CO2_concentration_values,
                                   intercellular_oxygen_values,
                                   photosynthetically_active_radiation_values,
                                   number_of_leaf_water_potential_sample_points=1000):
    """

    @param profit_optimisation_model:
    @param time_steps:
    @param soil_water_potential_values: MPa
    @param air_temperature_values: K
    @param air_vapour_pressure_deficit_values: kPa
    @param air_pressure_values: kPa
    @param atmospheric_CO2_concentration_values: umol mol-1
    @param intercellular_oxygen_values: umol mol-1
    @param photosynthetically_active_radiation_values: umol m-2 s-1
    @param number_of_leaf_water_potential_sample_points:

    @return: optimal leaf water potentials: MPa
    @return: net CO2 uptake values: umol m-2 s-1
    @return: transpiration rates: mmol m-2 s-1
    """

    # Setup output arrays
    optimal_leaf_water_potentials = zeros(len(time_steps))
    net_CO2_uptake_values = zeros(len(time_steps))
    transpiration_rate_values = zeros(len(time_steps))
    intercellular_CO2_values = zeros(len(time_steps))
    stomatal_conductance_to_CO2_values = zeros(len(time_steps))

    for i in range(len(time_steps)):
        (optimal_leaf_water_potentials[i],
         net_CO2_uptake_values[i],
         transpiration_rate_values[i],
         intercellular_CO2_values[i],
         stomatal_conductance_to_CO2_values[i]) = \
            profit_optimisation_model.optimal_state(soil_water_potential_values[i],
                                                    air_temperature_values[i],
                                                    air_vapour_pressure_deficit_values[i],
                                                    air_pressure_values[i],
                                                    atmospheric_CO2_concentration_values[i],
                                                    intercellular_oxygen_values[i],
                                                    photosynthetically_active_radiation_values[i],
                                                    number_of_leaf_water_potential_sample_points)

    return (optimal_leaf_water_potentials,
            net_CO2_uptake_values,
            transpiration_rate_values,
            intercellular_CO2_values,
            stomatal_conductance_to_CO2_values)
