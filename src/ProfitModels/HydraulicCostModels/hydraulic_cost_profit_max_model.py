"""
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
"""

from src.ProfitModels.HydraulicCostModels.hydraulic_cost_model import HydraulicCostModel


class ProfitMaxHydraulicCostModel(HydraulicCostModel):

    def hydraulic_cost(self, leaf_water_potential, soil_water_potential):
        """

        @param leaf_water_potential: (MPa)
        @param soil_water_potential: (MPa)
        @return: hydraulic cost (unitless)
        """

        hydraulic_conductance = self.hydraulic_conductance(leaf_water_potential)

        instantaneous_maximum_hydraulic_conductance = \
            self.instantaneous_maximum_hydraulic_conductance(soil_water_potential)

        return ((instantaneous_maximum_hydraulic_conductance - hydraulic_conductance)
                / (instantaneous_maximum_hydraulic_conductance - self._critical_hydraulic_conductance))