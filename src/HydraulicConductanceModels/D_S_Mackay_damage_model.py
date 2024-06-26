"""
-----------------------------------------------------------------------------------------
Implimentation of the model of xylem damage by embolism as described by D.S.Mackay et al.
 (2015).
-----------------------------------------------------------------------------------------
"""

from numpy import exp, linspace, asarray, clip, sum
from scipy.optimize import leastsq
from src.HydraulicConductanceModels.cumulative_Weibull_distribution_model \
    import (cumulative_Weibull_distribution,
            CumulativeWeibullDistribution,
            cumulative_Weibull_distribution_parameters_from_conductance_loss)


class DSMackayXylemDamageModel(CumulativeWeibullDistribution):

    _base_maximum_conductance: float
    _base_sensitivity_parameter: float
    _base_shape_parameter: float
    _base_critical_conductance_loss_fraction: float
    _N_sample_points_xylem_damage: int

    def __init__(self,
                 maximum_conductance,
                 sensitivity_parameter,
                 shape_parameter,
                 N_sample_points_xylem_damage = 1000,
                 critical_conductance_loss_fraction = 0.9,
                 xylem_recovery_water_potnetial: float = 0.,
                 PLC_damage_threshold = 0.05):

        self._base_maximum_conductance = maximum_conductance
        self._base_sensitivity_parameter = sensitivity_parameter
        self._base_shape_parameter = shape_parameter
        self._base_critical_conductance_loss_fraction = critical_conductance_loss_fraction
        self._N_sample_points_xylem_damage = N_sample_points_xylem_damage

        super().__init__(maximum_conductance,
                         sensitivity_parameter,
                         shape_parameter,
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

        # Setup the interim capped conductance model
        psi_array = linspace(0., self.critical_water_potential, self._N_sample_points_xylem_damage)
        capped_conductance_array = asarray([self.conductance(psi) for psi in psi_array])
        capped_conductance_array = clip(capped_conductance_array, 0, new_k_max)

        # Calculate the conductivity loss P of the capped conductance model at each
        # water potential.
        P_array = 1 - capped_conductance_array / new_k_max

        # Fit score to minimise as a function to pass to the curve_fit function
        def fit_score(c_current):
            k_prime = cumulative_Weibull_distribution(psi_array, new_k_max, b_new, c_current)

            numerator = sum(P_array * k_prime) ** 2
            denominator = (sum(P_array) ** 2) * (sum(k_prime) ** 2)

            return numerator / denominator

        # Fit the shape parameter c to the capped conductance model
        c_new = leastsq(func=fit_score,
                        x0=self.shape_parameter)[0][0]

        # Update model parameters
        self._k_max = new_k_max
        self._sensitivity_parameter = b_new
        self._shape_parameter = c_new
        self._critical_conductance_loss_fraction = self.critical_conductance / new_k_max

        return True

    def _recover_xylem(self, water_potential, timestep):
        """
        @param water_potential: (MPa)
        @param timestep: (s)
        @return: bool indicting if the model has changed
        """
        self.reset_xylem_damage()
        return True

    def reset_xylem_damage(self):
        """
        @return: None
        """
        self._k_max = self._base_maximum_conductance
        self._sensitivity_parameter = self._base_sensitivity_parameter
        self._shape_parameter = self._base_shape_parameter
        self._critical_conductance_loss_fraction = self._base_critical_conductance_loss_fraction
        return None

def D_S_Mackay_damage_model_from_conductance_loss(maximum_conductance,
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

    return DSMackayXylemDamageModel(maximum_conductance,
                                    sensitivity_parameter,
                                    shape_parameter,
                                    N_sample_points_xylem_damage,
                                    critical_conductance_loss_fraction,
                                    xylem_recovery_water_potnetial,
                                    PLC_damage_threshold)