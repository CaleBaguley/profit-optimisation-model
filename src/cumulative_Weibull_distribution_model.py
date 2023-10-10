"""
-----------------------------------------------------------------------------------------

-----------------------------------------------------------------------------------------
"""

from numpy import exp, power, log, abs
from numpy import linspace, trapz


class CumulativeWeibullDistribution:

    _maximum_conductance: float
    _sensitivity_parameter: float
    _shape_parameter: float

    def __init__(self, maximum_conductance: float, sensitivity_parameter: float, shape_parameter: float):

        """

        @param maximum_conductance: (mmol m-2 s-1 MPa-1)
        @param sensitivity_parameter: (MPa)
        @param shape_parameter: (unitless)
        """

        self._maximum_conductance = maximum_conductance
        self._sensitivity_parameter = sensitivity_parameter
        self._shape_parameter = shape_parameter

    def conductance(self, water_potential):

        """

        @param water_potential: (MPa)
        @return: conductance (mmol m-2 s-1 MPa-1)
        """

        return cumulative_Weibull_distribution(water_potential,
                                               self.maximum_conductance,
                                               self.sensitivity_parameter,
                                               self.shape_parameter)

    def water_potential_from_conductivity_loss_fraction(self, conductivity_loss_fraction):
        """

        @param conductivity_loss_fraction: (unitless)
        @return: (MPa)
        """

        conductivity_fraction = 1 - conductivity_loss_fraction

        return self.sensitivity_parameter * power(- log(conductivity_fraction), 1/self.shape_parameter)

    def water_potential_from_conductance(self, conductance):
        """

        @param conductance: (mmol m-2 s-1 MPa-1)
        @return: (MPa)
        """

        conductivity_loss_fraction = 1 - conductance/self.maximum_conductance

        return self.water_potential_from_conductivity_loss_fraction(conductivity_loss_fraction)

    def transpiration(self, min_water_potential, max_water_potential, steps = 100):
        """
        Calculates the transpiration rate across the water potentials using the trapezium integral approximation
        @param min_water_potential: (MPa)
        @param max_water_potential: (MPa)
        @param steps: number of steps for integral approximation (unitless)
        @return: (mmol m-2 s-1)
        """

        water_potential_values = linspace(min_water_potential, max_water_potential, steps)

        conductance_values = self.conductance(water_potential_values)

        return trapz(conductance_values, water_potential_values)

    @property
    def maximum_conductance(self):
        """
        @return: (mmol m-2 s-1 MPa-1)
        """
        return self._maximum_conductance

    @property
    def sensitivity_parameter(self):
        """
        @return: (MPa)
        """
        return self._sensitivity_parameter

    @property
    def shape_parameter(self):
        """
        @return: (Unitless)
        """
        return self._shape_parameter


def cumulative_Weibull_distribution(water_potential, maximum_conductance, sensitivity_parameter, shape_parameter):

    """

    @param water_potentials: (MPa)
    @param maximum_conductance: (mmol m-2 s-1 MPa-1)
    @param sensitivity_parameter: (MPa)
    @param shape_parameter: (unitless)
    @return: (MPa)
    """

    exponent = - power(water_potential / sensitivity_parameter, shape_parameter)

    return maximum_conductance * exp(exponent)

def cumulative_weibull_distribution_from_conductance_loss_at_given_water_potentials(maximum_conductance,
                                                                                    water_potential_1,
                                                                                    water_potential_2,
                                                                                    conductance_loss_fraction_1,
                                                                                    conductance_loss_fraction_2):

    """
    Creates a cumulative Weibull distribution from the fractional conductive loss at two given water potentials

    @param maximum_conductance: (mmol m-2 s-1 MPa-1)
    @param water_potential_1: (MPa)
    @param water_potential_2: (MPa)
    @param conductance_loss_fraction_1: (unitless)
    @param conductance_loss_fraction_2: (unitless)

    @return: CumulativeWeibullDistribution
    """

    ln_of_conductance_fraction_1 = log(1-conductance_loss_fraction_1)
    ln_of_conductance_fraction_2 = log(1-conductance_loss_fraction_2)

    ln_of_water_potential_1 = log(abs(water_potential_1))
    ln_of_water_potential_2 = log(abs(water_potential_2))

    shape_parameter = ((log(ln_of_conductance_fraction_1 / ln_of_conductance_fraction_2))
                       / (ln_of_water_potential_1 - ln_of_water_potential_2))

    sensitivity_parameter = water_potential_1 / power(-ln_of_conductance_fraction_1, 1/shape_parameter)

    return CumulativeWeibullDistribution(maximum_conductance, sensitivity_parameter, shape_parameter)
