"""
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
"""

from src.ProfitModels.HydraulicCostModels.hydraulic_cost_model import HydraulicCostModel
from src.HydraulicConductanceModels.hydraulic_conductance_model import HydraulicConductanceModel


class SOXHydraulicCostModel(HydraulicCostModel):

    def hydraulic_cost(self, leaf_water_potential, soil_water_potential):
        """

        @param leaf_water_potential: (MPa)
        @param soil_water_potential: (MPa)
        @return: hydraulic cost (unitless)
        """

        hydraulic_conductance = self.hydraulic_conductance(leaf_water_potential)

        return 1 - hydraulic_conductance / self.hydraulic_conductance_model.maximum_conductance
