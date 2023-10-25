"""
-------------------------------------------------------------------------

-------------------------------------------------------------------------
"""
from src.PhotosynthesisModels.Bonan_Model import PhotosynthesisModelRubiscoLimitedBonan, \
    PhotosynthesisModelElectronTransportLimitedBonan
import math
import numpy as np

from numpy import nanmax


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

        #print("net CO2: ", (atmospheric_CO2_concentration - intercellular_CO2_concentration) * stomatal_conductance_to_CO2)

        return ((atmospheric_CO2_concentration - intercellular_CO2_concentration) * stomatal_conductance_to_CO2,
                intercellular_CO2_concentration)


class PhotosynthesisModel(PhotosynthesisModelDummy):
    _photosynthesis_rubisco_limited_model: PhotosynthesisModelDummy
    _photosynthesis_electron_transport_limited_model: PhotosynthesisModelDummy

    def __init__(self,
                 photosynthesis_rubisco_limited_model=PhotosynthesisModelRubiscoLimitedBonan(),
                 photosynthesis_electron_transport_limited_model=PhotosynthesisModelElectronTransportLimitedBonan()):
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

        #print("Rubisco limited: ", intercellular_CO2_rubisco_limited)
        #print("electron limited: ", intercellular_CO2_electron_transport_limited)

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
