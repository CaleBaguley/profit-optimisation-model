"""
-----------------------------------------------------------------------------------------
Implimentation of the model of xylem damage by embolism as described by D.S.Mackay et al.
 (2015).
-----------------------------------------------------------------------------------------
"""

from numpy import exp, linspace, asarray, clip, sum
from scipy.optimize import leastsq
from src.HydraulicConductanceModels.cumulative_Weibull_distribution_model import (cumulative_Weibull_distribution,
                                                                                  CumulativeWeibullDistribution)

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
                 PLC_damage_threshold = 0.5,):

        self._base_maximum_conductance = maximum_conductance
        self._base_sensitivity_parameter = sensitivity_parameter
        self._base_shape_parameter = shape_parameter
        self._base_critical_conductance_loss_fraction = critical_conductance_loss_fraction
        self._N_sample_points_xylem_damage = N_sample_points_xylem_damage


        super().__init__(maximum_conductance,
                         sensitivity_parameter,
                         shape_parameter,
                         critical_conductance_loss_fraction,
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
        psi_array = linspace(0, self.critical_water_potential, self._N_sample_points_xylem_damage)
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
        self._maximum_conductance = new_k_max
        self._sensitivity_parameter = b_new
        self._shape_parameter = c_new
        self._critical_conductance_loss_fraction = self.critical_conductance / new_k_max

        return True



def D_S_Mackay_damage_model(new_k_max, current_CW_model, sample_points = 1000):

    """
    Create a new cumulative Weibull distribution model from the existing
    model and the new k_max value.

    @param new_k_max: float
    @param current_CW_model: CumulativeWeibullDistribution
    @param sample_points: number of water potentials to use when fitting
                          the shape parameter c; int
    @return: new CumulativeWeibullDistribution
    """

    # find the new sensitivity parameter. This is the water potential at which the
    # conductance of the current model is equal to the new k_max value times e to
    # the minus one.
    b_new = current_CW_model.water_potential_from_conductance(new_k_max * exp(-1))

    # Setup the interim capped conductance model
    psi_array = linspace(0, current_CW_model.critical_water_potential, sample_points)
    capped_conductance_array = asarray([current_CW_model.conductance(psi) for psi in psi_array])
    capped_conductance_array = clip(capped_conductance_array, 0, new_k_max)

    # Calculate the conductivity loss P of the capped conductance model at each
    # water potential.
    P_array = 1 - capped_conductance_array / new_k_max

    # Fit score to minimise as a function to pass to the curve_fit function
    def fit_score(c_current):
        k_prime = cumulative_Weibull_distribution(psi_array, new_k_max, b_new, c_current)

        numerator = sum(P_array*k_prime)**2
        denominator = (sum(P_array)**2)*(sum(k_prime)**2)

        return numerator/denominator

    # Fit the shape parameter c to the capped conductance model
    c_new = leastsq(func = fit_score,
                    x0 = current_CW_model.shape_parameter)[0][0]

    # Update critical conductance loss fraction. This maintains the critical conductance value
    critical_conductance_loss_fraction = current_CW_model.critical_conductance / new_k_max

    return CumulativeWeibullDistribution(new_k_max,
                                         b_new,
                                         c_new,
                                         critical_conductance_loss_fraction)
