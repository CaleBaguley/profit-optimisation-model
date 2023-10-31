"""
-----------------------------------------------------------------------------------------

-----------------------------------------------------------------------------------------
"""

from numpy import exp, power, log, abs
from numpy import linspace, trapz


class HydraulicConductanceModel:

    _maximum_conductance: float

    def __init__(self, maximum_conductance: float):

        """

        @param maximum_conductance: (mmol m-2 s-1 MPa-1)
        """

        self._maximum_conductance = maximum_conductance

    def conductance(self, water_potential):

        """

        @param water_potential: (MPa)
        @return: conductance (mmol m-2 s-1 MPa-1)
        """

        raise Exception("conductance method not implemented in base class")

    def water_potential_from_conductivity_loss_fraction(self, conductivity_loss_fraction):
        """

        @param conductivity_loss_fraction: (unitless)
        @return: (MPa)
        """

        raise Exception(" water_potential_from_conductivity_loss_fraction method not implemented in base class")

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