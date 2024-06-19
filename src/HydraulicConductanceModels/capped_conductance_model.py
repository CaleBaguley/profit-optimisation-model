"""
-----------------------------------------------------------------------------------------

-----------------------------------------------------------------------------------------
"""

from src.HydraulicConductanceModels.hydraulic_conductance_model import HydraulicConductanceModel
from numpy import clip

class CappedHydraulicConductanceModel(HydraulicConductanceModel):

    _conductance_cap: float
    _base_conductance_model: HydraulicConductanceModel

    def __init__(self, base_conductance_model: HydraulicConductanceModel, conductance_cap: float):

        """
        @param base_conductance_model: (HydraulicConductanceModel)
        @param conductance_cap: (mmol m-2 s-1 MPa-1)
        """

        self._base_conductance_model = base_conductance_model
        self._conductance_cap = conductance_cap
        super().__init__(base_conductance_model.maximum_conductance)

    def conductance(self, water_potential):

        """
        @param water_potential: (MPa)
        @return: conductance (mmol m-2 s-1 MPa-1)
        """

        return clip(self._base_conductance_model.conductance(water_potential), 0, self._conductance_cap)

    def water_potential_from_conductivity_loss_fraction(self, conductivity_loss_fraction):

        """
        @param conductivity_loss_fraction: (unitless)
        @return: (MPa)
        """

        raise Exception("water_potential_from_conductivity_loss_fraction not calculable in capped conductance model")

    def water_potential_from_conductance(self, conductance):

        """
        @param conductance: (mmol m-2 s-1 MPa-1)
        @return: (MPa)
        """

        raise Exception("water_potential_from_conductance not calculable in capped conductance model")

    def water_potential_at_cap_switch(self):

        """
        Function to get the water potential at which the conductance switches to the capped value
        @return: (MPa)
        """

        return self._base_conductance_model.water_potential_from_conductance(self._conductance_cap)

# -- Getters and Setters --
    def get_cap_conductance(self):

        """
        @return: (mmol m-2 s-1 MPa-1)
        """

        return self._conductance_cap

    def update_cap_cpnductance(self, conductance_cap):

        """
        @param conductance_cap: (mmol m-2 s-1 MPa-1)
        """

        self._conductance_cap = conductance_cap

    def get_base_conductance_model(self):

        """
        @return: (HydraulicConductanceModel)
        """

        return self._base_conductance_model
