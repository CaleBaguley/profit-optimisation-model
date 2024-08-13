"""
-------------------------------------------------------------------------
Rubisco and electron transport limited photosynthesis models calculated
using equations from "Climate Change and terestrial Ecosystem Modeling"
by Gordon Bonan
-------------------------------------------------------------------------
"""

import numpy as np
from numpy import roots, nanmax

from profit_optimisation_model.src.PhotosynthesisModels.photosynthesis_model \
    import PhotosynthesisModelDummy
from profit_optimisation_model.src.TemperatureDependenceModels.Q10_temperature_dependence_model \
    import Q10TemperatureDependenceModel
from profit_optimisation_model.src.TemperatureDependenceModels.arrhenius_and_peaked_arrhenius_function \
    import ArrheniusModel
from profit_optimisation_model.src.TemperatureDependenceModels.temperature_dependence_model \
    import TemperatureDependenceModel
from profit_optimisation_model.src.electron_transport_rate_model \
    import ElectronTransportRateModel
from profit_optimisation_model.src.rubisco_CO2_and_O_model import RubiscoRates


class PhotosynthesisModelRubiscoLimitedBonan(PhotosynthesisModelDummy):

    _rubisco_rates_model: RubiscoRates
    _CO2_compensation_point_model: TemperatureDependenceModel
    _mitochondrial_respiration_rate_model: TemperatureDependenceModel

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
                                        intercellular_O = None,
                                        utilized_photosynthetically_active_radiation = None):
        """


        @param stomatal_conductance_to_CO2: mol m-2 s-1
        @param atmospheric_CO2_concentration: umol mol-1
        @param leaf_temperature: K
        @param intercellular_O: umol mol-1
        @param utilized_photosynthetically_active_radiation: Not needed
        @return: intercellular CO2 concentration (umol mol-1)
        """

        mitochondrial_respiration_rate = (
            self._mitochondrial_respiration_rate_model.get_value_at_temperature(leaf_temperature))

        michaelis_menten_constant_carboxylation = (
            self._rubisco_rates_model.michaelis_menten_constant_carboxylation(leaf_temperature, intercellular_O))

        maximum_carboxylation_rate = self._rubisco_rates_model.maximum_carboxylation_rate(leaf_temperature)

        mitochondrial_respiration_rate = 0.015*maximum_carboxylation_rate

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
        intercellular_CO2_concentration = roots([A, B, C])

        intercellular_CO2_concentration = nanmax(intercellular_CO2_concentration)

        if(intercellular_CO2_concentration is None):
            return np.nan
        elif (intercellular_CO2_concentration < 0. or intercellular_CO2_concentration > atmospheric_CO2_concentration):
            return np.nan

        return intercellular_CO2_concentration


class PhotosynthesisModelElectronTransportLimitedBonan(PhotosynthesisModelDummy):

    _electron_transport_rate_model: ElectronTransportRateModel
    _CO2_compensation_point_model: TemperatureDependenceModel
    _mitochondrial_respiration_rate_model: TemperatureDependenceModel
    _rubisco_rates_model: RubiscoRates

    def __init__(self,
                 electron_transport_rate_model = ElectronTransportRateModel(),
                 CO2_compensation_point_model = ArrheniusModel(42.75, 37830.0),
                 mitochondrial_respiration_rate_model = Q10TemperatureDependenceModel(0.2, 2.),
                 rubisco_rates_model=RubiscoRates()
                 ):
        """

        @param electron_transport_rate_model:
        @param CO2_compensation_point_model:
        @param mitochondrial_respiration_rate_model:
        @param rubisco_rates_model:
        """

        self._electron_transport_rate_model = electron_transport_rate_model
        self._CO2_compensation_point_model = CO2_compensation_point_model
        self._mitochondrial_respiration_rate_model = mitochondrial_respiration_rate_model
        self._rubisco_rates_model = rubisco_rates_model

    def intercellular_CO2_concentration(self,
                                        stomatal_conductance_to_CO2,
                                        atmospheric_CO2_concentration,
                                        leaf_temperature,
                                        intercellular_O = None,
                                        utilized_photosynthetically_active_radiation = None):
        """

        @param stomatal_conductance_to_CO2: mol m-2 s-1
        @param atmospheric_CO2_concentration: umol mol-1
        @param leaf_temperature: K
        @param intercellular_O: umol mol-1
        @param utilized_photosynthetically_active_radiation: (umol m-2 unit time-1)
        @return: intercellular CO2 concentration (umol mol-1)
        """

        if(utilized_photosynthetically_active_radiation == 0.):
            return atmospheric_CO2_concentration

        maximum_carboxylation_rate = self._rubisco_rates_model.maximum_carboxylation_rate(leaf_temperature)

        mitochondrial_respiration_rate = 0.015 * maximum_carboxylation_rate

        #mitochondrial_respiration_rate = (
        #    self._mitochondrial_respiration_rate_model.get_value_at_temperature(leaf_temperature))

        electron_transport_rate = (
            self._electron_transport_rate_model.electron_transport_rate(leaf_temperature,
                                                                        utilized_photosynthetically_active_radiation))

        CO2_compensation_point = self._CO2_compensation_point_model.get_value_at_temperature(leaf_temperature)

        # Quadratic equation components (Ax^2 + Bx + C = 0)
        A = -stomatal_conductance_to_CO2

        B = (atmospheric_CO2_concentration * stomatal_conductance_to_CO2
             + mitochondrial_respiration_rate
             - 2 * stomatal_conductance_to_CO2 * CO2_compensation_point
             - electron_transport_rate / 4)

        C = (2 * atmospheric_CO2_concentration * stomatal_conductance_to_CO2 * CO2_compensation_point
             + 2 * mitochondrial_respiration_rate * CO2_compensation_point
             + CO2_compensation_point * CO2_compensation_point / 4)

        # Find the roots of the quadratic equation
        intercellular_CO2_concentration = roots([A, B, C])

        intercellular_CO2_concentration = nanmax(intercellular_CO2_concentration)

        if (intercellular_CO2_concentration is None):
            return np.nan
        elif (intercellular_CO2_concentration < 0. or intercellular_CO2_concentration > atmospheric_CO2_concentration):
            return np.nan

        return intercellular_CO2_concentration
        #return max(intercellular_CO2_concentration)
