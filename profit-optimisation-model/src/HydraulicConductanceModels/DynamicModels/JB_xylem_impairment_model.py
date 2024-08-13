"""
-----------------------------------------------------------------------------------------

-----------------------------------------------------------------------------------------
"""

from src.HydraulicConductanceModels.DynamicModels.Analytic_D_S_Mackay_damage_model \
    import DSMackayXylemDamageModelAnalytic
from src.HydraulicConductanceModels.cumulative_Weibull_distribution_model \
    import (cumulative_Weibull_distribution_parameters_from_conductance_loss)


class JBDynamicXylemConductanceModel(DSMackayXylemDamageModelAnalytic):

    _sapwood_area: float
    _base_sapwood_area: float
    _recovery_rate: float
    _recovery_shape: float
    _impairment_rate: float
    _impairment_shape: float
    _growth_rate: float
    _death_rate: float
    _death_shape: float

    def __init__(self,
                 maximum_conductance,
                 sapwood_area,
                 sensitivity_parameter,
                 shape_parameter,
                 critical_conductance_loss_fraction=0.9,
                 recovery_rate = 0.01,
                 recovery_shape = 1.0,
                 impairment_rate = 0.01,
                 impairment_shape = 1.0,
                 growth_rate = 0.01,
                 death_rate = 0.01,
                 death_shape = 1.0
                 ):

        self._sapwood_area = sapwood_area
        self._base_sapwood_area = sapwood_area
        self._recovery_rate = recovery_rate
        self._recovery_shape = recovery_shape
        self._impairment_rate = impairment_rate
        self._impairment_shape = impairment_shape
        self._growth_rate = growth_rate
        self._death_rate = death_rate
        self._death_shape = death_shape

        super().__init__(maximum_conductance,
                         sensitivity_parameter,
                         shape_parameter,
                         critical_conductance_loss_fraction,
                         0.,
                         0.)

    def update_xylem_damage(self, water_potential, timestep, transpiration_rate):

        # Calculate the new sapwood area
        self._sapwood_area += (self._growth_rate - self._death_rate) * timestep

        # calculate the current leaf conductance
        k_leaf = self.conductance(water_potential)

        # Calculate the impairment and recovery rates
        recovery_rate = self._calc_recovery_rate(k_leaf)
        impairment_rate = self._calc_impairment_rate(k_leaf)

        # Calculate the change in impaired conductance
        recovery = (self._base_maximum_conductance - self._k_max) * recovery_rate * timestep
        impairment = self._k_max * impairment_rate * timestep
        growth = (self._base_maximum_conductance - self._k_max) * self._growth_rate * timestep
        death = ((self._base_maximum_conductance * self.healthy_fraction**self._death_shape
                  - self._k_max)
                 * self._death_rate * timestep)

        new_k_max = self._k_max + recovery - impairment + growth - death

        # Update b and c parameters
        #self._k_max = max(new_k_max, self.critical_conductance)
        self._update_given_new_maximum_conductance(new_k_max)

        return False

    def _calc_recovery_rate(self, k_leaf):
        return self._recovery_rate * (k_leaf / self.maximum_conductance)**self._recovery_shape

    def _calc_impairment_rate(self, k_leaf):
        return self._impairment_rate * (1 - k_leaf / self.maximum_conductance)**self._impairment_shape

    def reset_xylem_damage(self):
        """
        @return: None
        """
        self._k_max = self._base_maximum_conductance
        self._sensitivity_parameter = self._base_sensitivity_parameter
        self._shape_parameter = self._base_shape_parameter
        self._critical_conductance_loss_fraction = self._base_critical_conductance_loss_fraction
        self._sapwood_area = self._base_sapwood_area
        return None

    #def conductance(self, water_potential):
    #    return self._sapwood_area * super().conductance(water_potential)

    # -- Properties --
    @property
    def sapwood_area(self):
        return self._sapwood_area

    @property
    def recovery_rate(self):
        return self._recovery_rate

    @property
    def impairment_rate(self):
        return self._impairment_rate

    @property
    def growth_rate(self):
        return self._growth_rate

    @property
    def death_rate(self):
        return self._death_rate

    @property
    def death_shape(self):
        return self._death_shape

    @property
    def healthy_fraction(self):
        return self._k_max / self._base_maximum_conductance

    @property
    def maximum_conductance(self):
        return self._sapwood_area * self._k_max


def JB_xylem_damage_model_from_conductance_loss(maximum_conductance,
                                                sapwood_area,
                                                water_potential_1,
                                                water_potential_2,
                                                conductance_loss_fraction_1,
                                                conductance_loss_fraction_2,
                                                critical_conductance_loss_fraction = 0.9,
                                                recovery_rate = 0.01,
                                                recovery_shape = 1.0,
                                                impairment_rate = 0.01,
                                                impairment_shape = 1.0,
                                                growth_rate = 0.01,
                                                death_rate = 0.01,
                                                death_shape = 1.0
                                               ):

    """
    @param maximum_conductance:
    @param sapwood_area: m2
    @param water_potential_1: MPa
    @param water_potential_2: MPa
    @param conductance_loss_fraction_1: unitless
    @param conductance_loss_fraction_2: unitless
    @param critical_conductance_loss_fraction: unitless
    @param recovery_rate: unitless
    @param recovery_shape: unitless
    @param impairment_rate: unitless
    @param impairment_shape: unitless
    @param growth_rate: unitless
    @param death_rate: unitless
    @param death_shape: unitless
    @return: DSMackayXylemDamageModel
    """

    maximum_conductance, sensitivity_parameter, shape_parameter = \
        cumulative_Weibull_distribution_parameters_from_conductance_loss(maximum_conductance,
                                                                         water_potential_1,
                                                                         water_potential_2,
                                                                         conductance_loss_fraction_1,
                                                                         conductance_loss_fraction_2)

    return JBDynamicXylemConductanceModel(maximum_conductance,
                                          sapwood_area,
                                          sensitivity_parameter,
                                          shape_parameter,
                                          critical_conductance_loss_fraction,
                                          recovery_rate,
                                          recovery_shape,
                                          impairment_rate,
                                          impairment_shape,
                                          growth_rate,
                                          death_rate,
                                          death_shape)
