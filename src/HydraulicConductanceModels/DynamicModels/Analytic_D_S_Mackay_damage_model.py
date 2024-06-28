"""
-----------------------------------------------------------------------------------------
Analytic implementation of the model of xylem damage by embolism as described by D.S.Mackay et al.
 (2015).
-----------------------------------------------------------------------------------------
"""

from numpy import exp, power
from src.HydraulicConductanceModels.DynamicModels.D_S_Mackay_damage_model \
    import DSMackayXylemDamageModel
from src.HydraulicConductanceModels.cumulative_Weibull_distribution_model \
    import cumulative_Weibull_distribution_parameters_from_conductance_loss


class DSMackayXylemDamageModelAnalytic(DSMackayXylemDamageModel):

    def __init__(self,
                 maximum_conductance,
                 sensitivity_parameter,
                 shape_parameter,
                 N_sample_points_xylem_damage=1000,
                 critical_conductance_loss_fraction=0.9,
                 xylem_recovery_water_potnetial: float = 0.,
                 PLC_damage_threshold=0.05):

        super().__init__(maximum_conductance,
                         sensitivity_parameter,
                         shape_parameter,
                         N_sample_points_xylem_damage,
                         critical_conductance_loss_fraction,
                         xylem_recovery_water_potnetial,
                         PLC_damage_threshold)

    def _damage_xylem(self, water_potential, timestep, transpiration_rate):

        """
        @param water_potential: (MPa)
        @param timestep: (s)
        @param transpiration_rate: (mmol m-2 s-1)
        @return: bool indicting if the model has changed
        """

        new_k_max = self.conductance(water_potential)

        # find the new sensitivity parameter. This is the water potential at which the
        # conductance of the current model is equal to the new k_max value times e to
        # the minus one.
        b_new = self.water_potential_from_conductance(new_k_max * exp(-1))

        # Calculate the shape parameter c. This is done by forcing the new model to have
        # the same grgadient as the previous one when the water potential is equal to the
        # new_b value.
        exponent = 1 - power(b_new/self._sensitivity_parameter, self._shape_parameter)

        new_c = (self._shape_parameter
                 * power(b_new/self._sensitivity_parameter, self._shape_parameter)
                 * (self._k_max / new_k_max)
                 * exp(exponent))

        self._k_max = new_k_max
        self._sensitivity_parameter = b_new
        self._shape_parameter = new_c

        return True


def analytic_D_S_Mackay_damage_model_from_conductance_loss(maximum_conductance,
                                                           water_potential_1,
                                                           water_potential_2,
                                                           conductance_loss_fraction_1,
                                                           conductance_loss_fraction_2,
                                                           N_sample_points_xylem_damage = 1000,
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
                                            N_sample_points_xylem_damage,
                                            critical_conductance_loss_fraction,
                                            xylem_recovery_water_potnetial,
                                            PLC_damage_threshold)