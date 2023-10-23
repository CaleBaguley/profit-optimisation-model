"""
---------------------------------------------------------------------------------
Class for modeling the electron transport rate during photosynthesis
---------------------------------------------------------------------------------
"""

from src.TemperatureDependenceModels.temperature_dependence_model import (TemperatureDependenceModel,
                                                                          LowTemperatureAdjustedModel)
from src.TemperatureDependenceModels.arrhenius_and_peaked_arrhenius_function import PeakedArrheniusModel
from numpy import roots, argwhere, amin
import math
import numpy as np


class ElectronTransportRateModel:

    _maximum_electron_transport_rate_model: TemperatureDependenceModel
    _curvature_parameter: float

    def __init__(self,
                 curvature_parameter = 0.85,
                 maximum_electron_transport_rate_model = LowTemperatureAdjustedModel(PeakedArrheniusModel(60., 30000., 200000., 650.))):

        self._curvature_parameter = curvature_parameter
        self._maximum_electron_transport_rate_model = maximum_electron_transport_rate_model

    def get_maximum_electron_transport_rate(self, temperature):
        return self._maximum_electron_transport_rate_model.get_value_at_temperature(temperature)

    def electron_transport_rate(self, temperature, utilized_photosynthetically_active_radiation = None):
        """
        Calculate the electron transport rate for a given temperature and amount of utilised photosynthetically
        active radiation
        @param temperature: (K)
        @param utilized_photosynthetically_active_radiation: (umol m-2 unit time-1). If set to None method returns
        the maximum electron transport rate for the given temperature.
        @return: electron transport rate (umol m-2 unit time-1)
        """
        maximum_transport_rate = self.get_maximum_electron_transport_rate(temperature)

        # If we are not supplied the amount of utilized light simply return the maximum electron transport rate
        if(utilized_photosynthetically_active_radiation is None):
            return maximum_transport_rate

        # to calculate electron transportation rate J we need to solve the following quadratic equation:
        # 0 = self._curvature_parameter * J^2
        #     - (utilized_photosynthetically_active_radiation + maximum_transport_rate) * J
        #     + utilized_photosynthetically_active_radiation * maximum_transport_rate

        a = self._curvature_parameter
        b = - (utilized_photosynthetically_active_radiation + maximum_transport_rate)
        c = utilized_photosynthetically_active_radiation * maximum_transport_rate

        # Uses numpy to calculate possible roots
        #solutions = roots([a, b, c])

        return quadratic(a, b, c, large=False)

        """
        # electron transport rate can only ever be positive
        solutions = solutions[argwhere(solutions >= 0)]

        # if there are no positive solutions return zero
        if(len(solutions) == 0):
            return 0.

        # Return the smallest possible positive solution
        return amin(solutions)"""

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
    d = b**2.0 - 4.0 * a * c # discriminant
    if d < 0.0:
        raise ValueError('imaginary root found')
    #root1 = np.where(d>0.0, (-b - np.sqrt(d)) / (2.0 * a), d)
    #root2 = np.where(d>0.0, (-b + np.sqrt(d)) / (2.0 * a), d)

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
