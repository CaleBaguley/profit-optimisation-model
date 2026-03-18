"""
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
"""

from profit_optimisation_model.src.ProfitModels.HydraulicCostModels.hydraulic_cost_model \
    import HydraulicCostModel


class SOXHydraulicCostModel(HydraulicCostModel):

    def hydraulic_cost(self, leaf_water_potential, soil_water_potential):
        """

        @param leaf_water_potential: (MPa)
        @param soil_water_potential: (MPa)
        @return: hydraulic cost (unitless)
        """

        # Note: original SOX used the mean of the leaf and root water potentials, here we use the leaf only.
        #       This matches the implementation in Sabot et al. 2022 and now JULES.

        hydraulic_conductance = self.hydraulic_conductance(leaf_water_potential,
                                                           leaf_water_potential,
                                                           soil_water_potential)

        return (1 - (hydraulic_conductance - self.critical_hydraulic_conductance)
                 / (self.hydraulic_conductance_model.maximum_conductance - self.critical_hydraulic_conductance))
