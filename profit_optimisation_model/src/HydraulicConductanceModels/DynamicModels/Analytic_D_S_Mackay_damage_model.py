"""
-----------------------------------------------------------------------------------------
Analytic implementation of the model of xylem damage by embolism as described by D.S.Mackay et al.
 (2015).
-----------------------------------------------------------------------------------------
"""

from numpy import exp, power, log
from profit_optimisation_model.src.HydraulicConductanceModels.DynamicModels.D_S_Mackay_damage_model \
    import DSMackayXylemDamageModel
from profit_optimisation_model.src.HydraulicConductanceModels.cumulative_Weibull_distribution_model \
    import cumulative_Weibull_distribution_parameters_from_conductance_loss


class DSMackayXylemDamageModelAnalytic(DSMackayXylemDamageModel):

    def __init__(self,
                 maximum_conductance,
                 sensitivity_parameter,
                 shape_parameter,
                 critical_conductance_loss_fraction=0.9,
                 xylem_recovery_water_potnetial: float = 0.,
                 PLC_damage_threshold=0.05):

        super().__init__(maximum_conductance,
                         sensitivity_parameter,
                         shape_parameter,
                         None,
                         critical_conductance_loss_fraction,
                         xylem_recovery_water_potnetial,
                         PLC_damage_threshold)

    def _damage_xylem(self, water_potential, timestep, transpiration_rate, root_water_potential):

        """
        @param root_water_potential:
        @param water_potential: (MPa)
        @param timestep: (s)
        @param transpiration_rate: (mmol m-2 s-1)
        @return: bool indicting if the model has changed
        """

        new_k_max = self.conductance(water_potential)

        self._update_given_new_maximum_conductance(new_k_max)

        return True

    def _update_given_new_maximum_conductance(self, new_k_max):

        """
        Update the shape of the vulnerability curve given a new maximum conductance.
        @param new_k_max:
        @return: None
        """

        new_k_max = min(max(new_k_max,
                            (1. - self._critical_conductance_loss_fraction) * self._base_maximum_conductance),
                        self._base_maximum_conductance)

        # find the new sensitivity parameter. This is the water potential at which the
        # conductance of the healthy model is equal to the new k_max value times e to
        # the minus one.
        b_new = (self._base_sensitivity_parameter
                 * power(1 - log(new_k_max / self._base_maximum_conductance),
                         1 / self._base_shape_parameter))

        # Calculate the shape parameter c. This is done by forcing the new model to have
        # the same gradient as the previous one when the water potential is equal to the
        # new_b value.

        new_c = (self._base_shape_parameter
                 * power(b_new / self._base_sensitivity_parameter, self._base_shape_parameter))

        self._k_max = new_k_max
        self._sensitivity_parameter = b_new
        self._shape_parameter = new_c


def analytic_D_S_Mackay_damage_model_from_conductance_loss(maximum_conductance,
                                                           water_potential_1,
                                                           water_potential_2,
                                                           conductance_loss_fraction_1,
                                                           conductance_loss_fraction_2,
                                                           critical_conductance_loss_fraction = 0.9,
                                                           xylem_recovery_water_potnetial = 0.,
                                                           PLC_damage_threshold = 0.05):

    """
    @param maximum_conductance:
    @param water_potential_1: MPa
    @param water_potential_2: MPa
    @param conductance_loss_fraction_1: unitless
    @param conductance_loss_fraction_2: unitless
    @param N_sample_points_xylem_damage: int
    @param critical_conductance_loss_fraction: unitless
    @param xylem_recovery_water_potnetial: MPa
    @param PLC_damage_threshold: unitless
    @return: DSMackayXylemDamageModel
    """

    maximum_conductance, sensitivity_parameter, shape_parameter = \
        cumulative_Weibull_distribution_parameters_from_conductance_loss(maximum_conductance,
                                                                         water_potential_1,
                                                                         water_potential_2,
                                                                         conductance_loss_fraction_1,
                                                                         conductance_loss_fraction_2)

    return DSMackayXylemDamageModelAnalytic(maximum_conductance,
                                            sensitivity_parameter,
                                            shape_parameter,
                                            critical_conductance_loss_fraction,
                                            xylem_recovery_water_potnetial,
                                            PLC_damage_threshold)