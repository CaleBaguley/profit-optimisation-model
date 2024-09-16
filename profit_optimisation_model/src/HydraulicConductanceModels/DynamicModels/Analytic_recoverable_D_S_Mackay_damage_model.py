"""
-----------------------------------------------------------------------------------------
Analytic implementation of the model of xylem damage by embolism as described by D.S.Mackay et al.
 (2015).
-----------------------------------------------------------------------------------------
"""

from profit_optimisation_model.src.HydraulicConductanceModels.DynamicModels.Analytic_D_S_Mackay_damage_model \
    import DSMackayXylemDamageModelAnalytic
from profit_optimisation_model.src.HydraulicConductanceModels.cumulative_Weibull_distribution_model \
    import cumulative_Weibull_distribution_parameters_from_conductance_loss


class DSMackayXylemDamageModelAnalyticRecovery(DSMackayXylemDamageModelAnalytic):

    _recovery_rate: float
    _damage_rate: float

    def __init__(self,
                 maximum_conductance,
                 sensitivity_parameter,
                 shape_parameter,
                 critical_conductance_loss_fraction=0.9,
                 recovery_rate = 0.01,
                 damage_rate = 0.01
                 ):

        self._recovery_rate = recovery_rate
        self._damage_rate = damage_rate

        super().__init__(maximum_conductance,
                         sensitivity_parameter,
                         shape_parameter,
                         critical_conductance_loss_fraction,
                         None,
                         None)

    def update_xylem_damage(self, water_potential, timestep, transpiration_rate, root_water_potential):
        """
        @param water_potential: (MPa)
        @param timestep: (s)
        @param transpiration_rate: (mmol m-2 s-1)
        @return: bool indicting if the model has changed
        """

        # Calculate the current conductance
        # NOTE: This is always less than or equal to the maximum conductance
        damage_target_maximum_conductance = self.conductance(water_potential)

        # Calculate the change due to damage
        damage_change = max(self._damage_rate * (self._k_max - damage_target_maximum_conductance),0)

        # Calculate the change due to recovery
        recovery_change = max(self._recovery_rate * (self._base_maximum_conductance - self._k_max),0)

        # Calculate the new maximum conductance given the current conductance and the recovery rate
        new_maximum_conductance = self._k_max + recovery_change - damage_change

        # Update the model given the new maximum conductance
        self._update_given_new_maximum_conductance(new_maximum_conductance)

        return False


def analytic_recoverable_D_S_Mackay_damage_model_from_conductance_loss(maximum_conductance,
                                                                       water_potential_1,
                                                                       water_potential_2,
                                                                       conductance_loss_fraction_1,
                                                                       conductance_loss_fraction_2,
                                                                       critical_conductance_loss_fraction = 0.9,
                                                                       recovery_rate = 0.01,
                                                                       damage_rate = 0.01):

    """
    @param maximum_conductance:
    @param water_potential_1: MPa
    @param water_potential_2: MPa
    @param conductance_loss_fraction_1: unitless
    @param conductance_loss_fraction_2: unitless
    @param critical_conductance_loss_fraction: unitless
    @param recovery_rate: unitless
    @param damage_rate: unitless
    @return: DSMackayXylemDamageModel
    """

    maximum_conductance, sensitivity_parameter, shape_parameter = \
        cumulative_Weibull_distribution_parameters_from_conductance_loss(maximum_conductance,
                                                                         water_potential_1,
                                                                         water_potential_2,
                                                                         conductance_loss_fraction_1,
                                                                         conductance_loss_fraction_2)

    return DSMackayXylemDamageModelAnalyticRecovery(maximum_conductance,
                                                    sensitivity_parameter,
                                                    shape_parameter,
                                                    critical_conductance_loss_fraction,
                                                    recovery_rate,
                                                    damage_rate)