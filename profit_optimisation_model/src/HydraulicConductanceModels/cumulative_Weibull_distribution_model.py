"""
-----------------------------------------------------------------------------------------

-----------------------------------------------------------------------------------------
"""

from src.HydraulicConductanceModels.hydraulic_conductance_model import HydraulicConductanceModel
from numpy import exp, power, log, abs


class CumulativeWeibullDistribution(HydraulicConductanceModel):

    _sensitivity_parameter: float
    _shape_parameter: float

    def __init__(self,
                 maximum_conductance: float,
                 sensitivity_parameter: float,
                 shape_parameter: float,
                 critical_conductance_loss_fraction: float = 0.9,
                 xylem_recovery_water_potnetial: float = 0.,
                 PLC_damage_threshold: float = 0.1):

        """

        @param maximum_conductance: (mmol m-2 s-1 MPa-1)
        @param sensitivity_parameter: (MPa)
        @param shape_parameter: (unitless)
        @param critical_conductance_loss_fraction: (unitless)
        @param xylem_recovery_water_potnetial: (MPa)
        """

        self._sensitivity_parameter = sensitivity_parameter
        self._shape_parameter = shape_parameter
        super().__init__(maximum_conductance,
                         critical_conductance_loss_fraction,
                         xylem_recovery_water_potnetial,
                         PLC_damage_threshold)

    def conductance(self, water_potential):

        """

        @param water_potential: (MPa)
        @return: conductance (mmol m-2 s-1 MPa-1)
        """

        return cumulative_Weibull_distribution(water_potential,
                                               self.maximum_conductance,
                                               self.sensitivity_parameter,
                                               self.shape_parameter)

    def water_potential_from_conductivity_loss_fraction(self, conductivity_loss_fraction):
        """

        @param conductivity_loss_fraction: (unitless)
        @return: (MPa)
        """

        conductivity_fraction = 1. - conductivity_loss_fraction

        return self.sensitivity_parameter * power(- log(conductivity_fraction), 1/self.shape_parameter)

    def water_potential_from_conductance(self, conductance):
        """

        @param conductance: (mmol m-2 s-1 MPa-1)
        @return: (MPa)
        """

        conductivity_loss_fraction = 1 - conductance/self.maximum_conductance

        return self.water_potential_from_conductivity_loss_fraction(conductivity_loss_fraction)

    @property
    def sensitivity_parameter(self):
        """
        @return: (MPa)
        """
        return self._sensitivity_parameter

    @property
    def shape_parameter(self):
        """
        @return: (Unitless)
        """
        return self._shape_parameter


def cumulative_Weibull_distribution(water_potentials, maximum_conductance, sensitivity_parameter, shape_parameter):

    """

    @param water_potentials: (MPa)
    @param maximum_conductance: (mmol m-2 s-1 MPa-1)
    @param sensitivity_parameter: (MPa)
    @param shape_parameter: (unitless)
    @return: (MPa)
    """

    exponent = - power(water_potentials / sensitivity_parameter, shape_parameter)

    return maximum_conductance * exp(exponent)


def cumulative_Weibull_distribution_from_conductance_loss_at_given_water_potentials(
        maximum_conductance,
        water_potential_1,
        water_potential_2,
        conductance_loss_fraction_1,
        conductance_loss_fraction_2,
        critical_conductance_loss_fraction = 0.9,
        xylem_recovery_water_potnetial = 0.,
        PLC_damage_threshold = 0.1
        ):

    """
    Creates a cumulative Weibull distribution from the fractional conductive loss at two given water potentials

    @param maximum_conductance: (mmol m-2 s-1 MPa-1)
    @param water_potential_1: (MPa)
    @param water_potential_2: (MPa)
    @param conductance_loss_fraction_1: (unitless)
    @param conductance_loss_fraction_2: (unitless)
    @param critical_conductance_loss_fraction: (unitless)
    @param xylem_recovery_water_potnetial: (MPa)
    @param PLC_damage_threshold: (unitless)

    @return: CumulativeWeibullDistribution
    """

    maximum_conductance, sensitivity_parameter, shape_parameter = \
        cumulative_Weibull_distribution_parameters_from_conductance_loss(maximum_conductance,
                                                                         water_potential_1,
                                                                         water_potential_2,
                                                                         conductance_loss_fraction_1,
                                                                         conductance_loss_fraction_2
                                                                         )

    return CumulativeWeibullDistribution(maximum_conductance,
                                         sensitivity_parameter,
                                         shape_parameter,
                                         critical_conductance_loss_fraction,
                                         xylem_recovery_water_potnetial,
                                         PLC_damage_threshold
                                         )

def cumulative_Weibull_distribution_parameters_from_conductance_loss(maximum_conductance,
                                                                     water_potential_1,
                                                                     water_potential_2,
                                                                     conductance_loss_fraction_1,
                                                                     conductance_loss_fraction_2):

    """

    @param maximum_conductance:
    @param water_potential_1:
    @param water_potential_2:
    @param conductance_loss_fraction_1:
    @param conductance_loss_fraction_2:
    @return: maximum conductance, sensitivity parameter, shape parameter
    """

    ln_of_conductance_fraction_1 = log(1 - conductance_loss_fraction_1)
    ln_of_conductance_fraction_2 = log(1 - conductance_loss_fraction_2)

    ln_of_water_potential_1 = log(abs(water_potential_1))
    ln_of_water_potential_2 = log(abs(water_potential_2))

    shape_parameter = ((log(ln_of_conductance_fraction_1 / ln_of_conductance_fraction_2))
                       / (ln_of_water_potential_1 - ln_of_water_potential_2))

    sensitivity_parameter = water_potential_1 / power(-ln_of_conductance_fraction_1, 1 / shape_parameter)

    return maximum_conductance, sensitivity_parameter, shape_parameter