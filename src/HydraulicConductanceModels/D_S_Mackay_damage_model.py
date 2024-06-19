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
