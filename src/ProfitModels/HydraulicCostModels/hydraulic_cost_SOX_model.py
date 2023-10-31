"""
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
"""

from src.ProfitModels.HydraulicCostModels.hydraulic_cost_model import HydraulicCostModel
from src.HydraulicConductanceModels.hydraulic_conductance_model import HydraulicConductanceModel


class ProfitMaxHydraulicCostModel(HydraulicCostModel):
    _curvature_parameter: float

    def __init__(self,
                 hydraulic_conductance_model: HydraulicConductanceModel,
                 critical_leaf_water_potential: float,
                 curvature_parameter: float):

        super().__init__(hydraulic_conductance_model, critical_leaf_water_potential)
        self._curvature_parameter = curvature_parameter

    def hydraulic_cost(self, leaf_water_potential, soil_water_potential):
        """

        @param leaf_water_potential: (MPa)
        @param soil_water_potential: (MPa)
        @return: hydraulic cost (unitless)
        """

        hydraulic_conductance = self.hydraulic_conductance(leaf_water_potential)

        return hydraulic_conductance / self.hydraulic_conductance_model.maximum_conductance
