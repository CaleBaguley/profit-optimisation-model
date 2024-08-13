"""
--------------------------------------------------------------------------------------------
Model for calculating the rate of carboxylation and oxygenation by Rubisco.
--------------------------------------------------------------------------------------------
"""

from profit_optimisation_model.src.TemperatureDependenceModels.arrhenius_and_peaked_arrhenius_function import(
    ArrheniusModel, PeakedArrheniusModel)
from profit_optimisation_model.src.michaelis_menten_response_function import (
    michaelis_menten_constant, michaelis_menten_response_function)


class RubiscoRates:
    def __init__(self,
                 maximum_carboxylation_rate_model = PeakedArrheniusModel(rate_at_25_centigrade=30.0,
                                                                         activation_energy=60000.,
                                                                         deactivation_energy=200000.,
                                                                         entropy_term=650.),
                 maximum_oxygenation_rate_model = PeakedArrheniusModel(rate_at_25_centigrade=60.0,
                                                                       activation_energy=3000.,
                                                                       deactivation_energy=200000.,
                                                                       entropy_term=650.),
                 michaelis_menten_constant_CO2_model = ArrheniusModel(rate_at_25_centigrade=404.9,
                                                                      activation_energy=79430.0),
                 michaelis_menten_constant_O_model = ArrheniusModel(rate_at_25_centigrade=278.4,
                                                                    activation_energy=36380.0)):
        """
        Model for calculating the RubiscoRates
        @param maximum_carboxylation_rate_model: Temperature dependence model
        @param maximum_oxygenation_rate_model: Temperature dependence model
        @param michaelis_menten_constant_CO2_model: Temperature dependence model
        @param michaelis_menten_constant_O_model: Temperature dependence model
        """

        self._maximum_carboxylation_rate_model = maximum_carboxylation_rate_model
        #self._maximum_oxygenation_rate_model = maximum_oxygenation_rate_model
        self._michaelis_menten_constant_CO2_model = michaelis_menten_constant_CO2_model

        if michaelis_menten_constant_O_model is None:
            michaelis_menten_constant_O_model = ArrheniusModel(rate_at_25_centigrade=278.4, activation_energy=36380.0)
        self._michaelis_menten_constant_O_model = michaelis_menten_constant_O_model

    def michaelis_menten_constant_carboxylation(self, leaf_temperature, inter_cellular_oxygen):
        """
        Calculates the Michaelis-Menten constant for carboxylation by Rubisco.
        Models the oxygen and carbon dioxide Michaelis-Menten constants temperature dependence by an Arrhenius function.
        @param leaf_temperature: (K)
        @param inter_cellular_oxygen: (umol mol-1)
        @return: Michaelis-Menten constant for carboxylation (umol mol-1)
        """
        michaelis_menten_constant_O = self.michaelis_menten_constant_O(leaf_temperature)
        michaelis_menten_constant_CO2 = self.michaelis_menten_constant_CO2(leaf_temperature)

        michaelis_menten_constant_carboxylation = michaelis_menten_constant(inter_cellular_oxygen,
                                                                            michaelis_menten_constant_O,
                                                                            michaelis_menten_constant_CO2)

        return michaelis_menten_constant_carboxylation

    def michaelis_menten_constant_oxygenation(self, leaf_temperature, inter_cellular_carbon):
        """
        Calculates the Michaelis-Menten constant for oxygenation by Rubisco.
        Models the oxygen and carbon dioxide Michaelis-Menten constants temperature dependence by an Arrhenius function.
        @param leaf_temperature: (K)
        @param inter_cellular_oxygen: (umol mol-1)
        @return: Michaelis-Menten constant for carboxylation (umol mol-1)
        """
        michaelis_menten_constant_O = self.michaelis_menten_constant_O(leaf_temperature)
        michaelis_menten_constant_CO2 = self.michaelis_menten_constant_CO2(leaf_temperature)

        michaelis_menten_constant_carboxylation = michaelis_menten_constant(inter_cellular_carbon,
                                                                            michaelis_menten_constant_CO2,
                                                                            michaelis_menten_constant_O)

        return michaelis_menten_constant_carboxylation

    def michaelis_menten_constant_CO2(self, leaf_temperature):
        """
        Calculates the Michaelis-Menten constant for carbon dioxide at a given temperature using an arrhenius function.
        @param leaf_temperature: (K)
        @return: Michaelis-Menten constant for carbon dioxide (umol mol-1)
        """
        return self._michaelis_menten_constant_CO2_model.get_value_at_temperature(leaf_temperature)

    def michaelis_menten_constant_O(self, leaf_temperature):
        """
        Calculates the Michaelis-Menten constant for oxygen at a given temperature  at a given temperature using an
        arrhenius function.
        @param leaf_temperature: (K)
        @return: Michaelis-Menten constant for oxygen (umol mol-1)
        """
        return self._michaelis_menten_constant_O_model.get_value_at_temperature(leaf_temperature)

    def maximum_carboxylation_rate(self, temperature):
        return self._maximum_carboxylation_rate_model.get_value_at_temperature(temperature)

    #def maximum_oxygenation_rate(self, temperature):
    #    return self._maximum_oxygenation_rate_model.get_value_at_temperature(temperature)

    def carboxylation_rate(self, intercellular_CO2_concentration, intercellular_O_concentration,
                           temperature):
        """
        Calculates the rate of carbon dioxide uptake during carboxylation by Rubisco.
        NOTE: Parameters denoted by (*) must share the same units.

        Variables
        @param intercellular_CO2_concentration: Intercellular concentration of carbon dioxide. (*)
        @param intercellular_O_concentration: Intercellular concentration of oxygen. (*)
        @param temperature: (K)

        @return: Rate of carbon dioxide uptake during carboxylation (same units as maximum_rate)
        """

        # Calculate shape parameters for given temperature
        maximum_rate = self._maximum_carboxylation_rate_model.get_value_at_temperature(temperature)
        michaelis_menten_constant_O = self._michaelis_menten_constant_O_model.get_value_at_temperature(temperature)
        michaelis_menten_constant_CO2 = self._michaelis_menten_constant_CO2_model.get_value_at_temperature(temperature)

        reaction_michaelis_menten_constant = michaelis_menten_constant(intercellular_O_concentration,
                                                                       michaelis_menten_constant_O,
                                                                       michaelis_menten_constant_CO2)

        return michaelis_menten_response_function(intercellular_CO2_concentration,
                                                  maximum_rate, reaction_michaelis_menten_constant)

    def oxygenation_rate(self, intercellular_O_concentration, intercellular_CO2_concentration,
                         temperature):
        """
        Calculates the rate of oxygen uptake during oxygenation by Rubisco.
        NOTE: Parameters denoted by (*) must share the same units.

        Variables
        @param intercellular_O_concentration: Intercellular concentration of oxygen. (*)
        @param intercellular_CO2_concentration: Intercellular concentration of carbon dioxide. (*)
        @param temperature: (K)

        @return: Rate of oxygen uptake during oxygenation (same units as maximum_rate)
        """

        # Calculate shape parameters for given temperature
        maximum_rate = self._maximum_carboxylation_rate_model.get_value_at_temperature(temperature)
        michaelis_menten_constant_O = self._michaelis_menten_constant_O_model.get_value_at_temperature(temperature)
        michaelis_menten_constant_CO2 = self._michaelis_menten_constant_CO2_model.get_value_at_temperature(temperature)

        reaction_michaelis_menten_constant = michaelis_menten_constant(intercellular_CO2_concentration,
                                                                       michaelis_menten_constant_CO2,
                                                                       michaelis_menten_constant_O)

        return michaelis_menten_response_function(intercellular_O_concentration,
                                                  maximum_rate, reaction_michaelis_menten_constant)
