"""
-------------------------------------------------------------------------

-------------------------------------------------------------------------
"""

from src.rubisco_CO2_and_O_model import RubiscoRates
from src.electron_transport_rate_model import ElectronTransportRateModel
from src.TemperatureDependenceModels.temperature_dependence_model import TemperatureDependenceModel
from src.TemperatureDependenceModels.arrhenius_and_peaked_arrhenius_function import ArrheniusModel
from src.TemperatureDependenceModels.Q10_temperature_dependence_model import Q10TemperatureDependenceModel
import math
import numpy as np

from numpy import roots, nanmax


class PhotosynthesisModelDummy:

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

        raise Exception("intercellular_CO2_concentration method not implemented in PhotosynthesisModelDummy class")

    def net_rate_of_CO2_assimilation(self,
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
        @return: net rate of CO2 assimilation umol m-2 s-1
        """

        intercellular_CO2_concentration = self.intercellular_CO2_concentration(
            stomatal_conductance_to_CO2,
            atmospheric_CO2_concentration,
            leaf_temperature,
            intercellular_O,
            utilized_photosynthetically_active_radiation)

        print("net CO2: ", (atmospheric_CO2_concentration - intercellular_CO2_concentration) * stomatal_conductance_to_CO2)

        return ((atmospheric_CO2_concentration - intercellular_CO2_concentration) * stomatal_conductance_to_CO2,
                intercellular_CO2_concentration)


class PhotosynthesisModelRubiscoLimited(PhotosynthesisModelDummy):

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

        print(stomatal_conductance_to_CO2)
        print(atmospheric_CO2_concentration)
        print(mitochondrial_respiration_rate)
        print(michaelis_menten_constant_carboxylation)
        print(maximum_carboxylation_rate)

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
        intercellular_CO2_concentration = max(intercellular_CO2_concentration)


        #intercellular_CO2_concentration = quadratic(a=A, b=B, c=C, large=True)

        print("roots: ", intercellular_CO2_concentration)
        if(intercellular_CO2_concentration is None):
            return np.nan
        elif (intercellular_CO2_concentration < 0. or intercellular_CO2_concentration > atmospheric_CO2_concentration):
            return np.nan

        return intercellular_CO2_concentration
        #return max(intercellular_CO2_concentration)


class PhotosynthesisModelElectronTransportLimited(PhotosynthesisModelDummy):

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

        if(len(intercellular_CO2_concentration) == 0):
            return np.nan

        intercellular_CO2_concentration = max(intercellular_CO2_concentration)

        #intercellular_CO2_concentration = quadratic(A, B, C, large=True)

        if (intercellular_CO2_concentration is None):
            return np.nan
        elif (intercellular_CO2_concentration < 0. or intercellular_CO2_concentration > atmospheric_CO2_concentration):
            return np.nan

        return intercellular_CO2_concentration
        #return max(intercellular_CO2_concentration)


class PhotosynthesisModel(PhotosynthesisModelDummy):
    _photosynthesis_rubisco_limited_model: PhotosynthesisModelRubiscoLimited
    _photosynthesis_electron_transport_limited_model: PhotosynthesisModelElectronTransportLimited

    def __init__(self,
                 photosynthesis_rubisco_limited_model=PhotosynthesisModelRubiscoLimited(),
                 photosynthesis_electron_transport_limited_model=PhotosynthesisModelElectronTransportLimited()):
        self._photosynthesis_rubisco_limited_model = photosynthesis_rubisco_limited_model
        self._photosynthesis_electron_transport_limited_model = photosynthesis_electron_transport_limited_model

    def intercellular_CO2_concentration(self,
                                        stomatal_conductance_to_CO2,
                                        atmospheric_CO2_concentration,
                                        leaf_temperature,
                                        intercellular_O=None,
                                        utilized_photosynthetically_active_radiation=None):

        intercellular_CO2_rubisco_limited = \
            (self._photosynthesis_rubisco_limited_model
             .intercellular_CO2_concentration(stomatal_conductance_to_CO2,
                                              atmospheric_CO2_concentration,
                                              leaf_temperature,
                                              intercellular_O,
                                              utilized_photosynthetically_active_radiation))

        intercellular_CO2_electron_transport_limited = \
            (self._photosynthesis_electron_transport_limited_model
             .intercellular_CO2_concentration(stomatal_conductance_to_CO2,
                                              atmospheric_CO2_concentration,
                                              leaf_temperature,
                                              intercellular_O,
                                              utilized_photosynthetically_active_radiation))

        print("Rubisco limited: ", intercellular_CO2_rubisco_limited)
        print("electron limited: ", intercellular_CO2_electron_transport_limited)

        if((intercellular_CO2_electron_transport_limited is None) & (intercellular_CO2_rubisco_limited is None)):
            return 0.
        elif(intercellular_CO2_rubisco_limited is None):
            return intercellular_CO2_electron_transport_limited
        elif(intercellular_CO2_electron_transport_limited is None):
            return intercellular_CO2_rubisco_limited

        return nanmax([intercellular_CO2_rubisco_limited, intercellular_CO2_electron_transport_limited])


def quadratic(a=None, b=None, c=None, large=False):
    """ minimilist quadratic solution as root for J solution should always
    be positive, so I have excluded other quadratic solution steps. I am
    only returning the smallest of the two roots

    Parameters:
    ----------
    a : float
        co-efficient
    b : float
        co-efficient
    c : float
        co-efficient

    Returns:
    -------
    val : float
        positive root
    """
    d = b ** 2.0 - 4.0 * a * c  # discriminant
    if d < 0.0:
        raise ValueError('imaginary root found')
    # root1 = np.where(d>0.0, (-b - np.sqrt(d)) / (2.0 * a), d)
    # root2 = np.where(d>0.0, (-b + np.sqrt(d)) / (2.0 * a), d)

    if large:
        if math.isclose(a, 0.0) and b > 0.0:
            root = -c / b
        elif math.isclose(a, 0.0) and math.isclose(b, 0.0):
            root = 0.0
            if c != 0.0:
                raise ValueError('Cant solve quadratic')
        else:
            root = (-b + np.sqrt(d)) / (2.0 * a)
    else:
        if math.isclose(a, 0.0) and b > 0.0:
            root = -c / b
        elif math.isclose(a, 0.0) and math.isclose(b, 0.0):
            root = 0.0
            if c != 0.0:
                raise ValueError('Cant solve quadratic')
        else:
            root = (-b - np.sqrt(d)) / (2.0 * a)

    return root


"""def get_root(a,b,c):
    ref_root = quad(a, b, c) * conv.ref_kPa * conv.FROM_MILI
    
    if (ref_root > Cs) or (ref_root < cst.zero):
    return quad(a, b, c, large_root=False) * conv.ref_kPa * conv.FROM_MILI
    else:
    return ref_root"""