"""
-----------------------------------------------------------------------------------------

-----------------------------------------------------------------------------------------
"""

from numpy import power
from src.HydraulicConductanceModels.hydraulic_conductance_model import HydraulicConductanceModel


class SOXHydraulicConductanceModel(HydraulicConductanceModel):
    _water_potential_at_half_conductance: float
    _shape_parameter: float

    def __init__(self, maximum_conductance: float, water_potential_at_half_conductance: float, shape_parameter: float):

        """

        @param maximum_conductance: (mmol m-2 s-1 MPa-1)
        @param water_potential_at_half_conductance: (MPa)
        @param shape_parameter: (unitless)
        """

        super().__init__(maximum_conductance)
        self._water_potential_at_half_conductance = water_potential_at_half_conductance
        self._shape_parameter = shape_parameter

    def conductance(self, water_potential):
        """

        @param water_potential: (MPa)
        @return: conductance (mmol m-2 s-1 MPa-1)
        """

        normalised_conductance = 1/(1+power(water_potential/self._water_potential_at_half_conductance,
                                            self._shape_parameter))

        return self._maximum_conductance * normalised_conductance

    def water_potential_from_conductivity_loss_fraction(self, conductivity_loss_fraction):
        """

        @param conductivity_loss_fraction: (unitless)
        @return: (MPa)
        """

        return self._water_potential_at_half_conductance * power(conductivity_loss_fraction - 1,
                                                                 1/self._shape_parameter)
