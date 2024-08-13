"""
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
"""

from src.ProfitModels.HydraulicCostModels.hydraulic_cost_model import HydraulicCostModel


class SOXHydraulicCostModel(HydraulicCostModel):

    def hydraulic_cost(self, leaf_water_potential, soil_water_potential):
        """

        @param leaf_water_potential: (MPa)
        @param soil_water_potential: (MPa)
        @return: hydraulic cost (unitless)
        """

        # SOX uses the conductance calculated for the mean of the leaf and root zone
        # water potentials.
        mean_conductance = (leaf_water_potential + soil_water_potential)/2

        hydraulic_conductance = self.hydraulic_conductance(mean_conductance)

        return (1 - (hydraulic_conductance - self.critical_hydraulic_conductance)
                 / (self.hydraulic_conductance_model.maximum_conductance - self.critical_hydraulic_conductance))
