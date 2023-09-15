"""
-------------------------------------------------------------------------

-------------------------------------------------------------------------
"""

from src.rubisco_CO2_and_O_model import RubiscoRates
from src.TemperatureDependenceModels.arrhenius_and_peaked_arrhenius_function import ArrheniusModel
from src.TemperatureDependenceModels.Q10_temperature_dependence_model import Q10TemperatureDependenceModel

from numpy import roots, max, min


class PhotosynthesisModel:

    def __init__(self,
                 rubisco_rates_model=RubiscoRates(),
                 CO2_compensation_point_model=ArrheniusModel(42.75, 37830.0),
                 mitochondrial_respiration_rate_model=Q10TemperatureDependenceModel(0.2, 2.)):
        self._rubisco_rates_model = rubisco_rates_model
        self._CO2_compensation_point_model = CO2_compensation_point_model
        self._mitochondrial_respiration_rate_model = mitochondrial_respiration_rate_model

    def intercellular_CO2_concentration(self,
                                        stomatal_conductance_to_CO2,
                                        atmospheric_CO2_concentration,
                                        leaf_temperature,
                                        intercellular_O):
        """

        @param stomatal_conductance_to_CO2: mol m-2 s-1
        @param atmospheric_CO2_concentration: umol mol-1
        @param leaf_temperature: K
        @param intercellular_O: umol mol-1
        @return: intercellular CO2 concentration (umol mol-1)
        """

        mitochondrial_respiration_rate = (
            self._mitochondrial_respiration_rate_model.get_value_at_temperature(leaf_temperature))

        michaelis_menten_constant_carboxylation = (
            self._rubisco_rates_model.michaelis_menten_constant_carboxylation(leaf_temperature, intercellular_O))

        maximum_carboxylation_rate = self._rubisco_rates_model.maximum_carboxylation_rate(leaf_temperature)

        CO2_compensation_point = self._CO2_compensation_point_model.get_value_at_temperature(leaf_temperature)

        # Quadratic equation components (Ax^2 + Bx + C = 0)
        A = -stomatal_conductance_to_CO2

        B = (atmospheric_CO2_concentration * stomatal_conductance_to_CO2
             + mitochondrial_respiration_rate
             - stomatal_conductance_to_CO2 * michaelis_menten_constant_carboxylation
             - maximum_carboxylation_rate)

        C = (atmospheric_CO2_concentration * stomatal_conductance_to_CO2 * michaelis_menten_constant_carboxylation
             + mitochondrial_respiration_rate * michaelis_menten_constant_carboxylation
             + CO2_compensation_point * maximum_carboxylation_rate)

        # Find the roots of the quadratic equation
        intercellular_CO2_concentratin = roots([A, B, C])

        return max(intercellular_CO2_concentratin)

    def net_rate_of_CO2_assimilation(self,
                                     stomatal_conductance_to_CO2,
                                     atmospheric_CO2_concentration,
                                     leaf_temperature,
                                     intercellular_O):
        """

        @param stomatal_conductance_to_CO2: mol m-2 s-1
        @param atmospheric_CO2_concentration: umol mol-1
        @param leaf_temperature: K
        @param intercellular_O: umol mol-1
        @return: net rate of CO2 assimilation umol m-2 s-1
        """

        intercellular_CO2_concentration = self.intercellular_CO2_concentration(stomatal_conductance_to_CO2,
                                                                               atmospheric_CO2_concentration,
                                                                               leaf_temperature,
                                                                               intercellular_O)

        return (atmospheric_CO2_concentration - intercellular_CO2_concentration) * stomatal_conductance_to_CO2

