"""
-----------------------------------------------------------------------------------------

-----------------------------------------------------------------------------------------
"""

from numpy import power, log
from profit_optimisation_model.src.HydraulicConductanceModels.hydraulic_conductance_model \
    import HydraulicConductanceModel


class SOXHydraulicConductanceModel(HydraulicConductanceModel):
    _water_potential_at_half_conductance: float
    _shape_parameter: float

    def __init__(self,
                 maximum_conductance: float,
                 water_potential_at_half_conductance: float,
                 shape_parameter: float,
                 critical_conductance_loss_fraction: float = 0.9,
                 xylen_recovery_water_potential: float = 0.):

        """

        @param maximum_conductance: (mmol m-2 s-1 MPa-1)
        @param water_potential_at_half_conductance: (MPa)
        @param shape_parameter: (unitless)
        @param critical_conductance_loss_fraction: (unitless)
        """

        super().__init__(maximum_conductance, critical_conductance_loss_fraction, xylen_recovery_water_potential)
        self._water_potential_at_half_conductance = water_potential_at_half_conductance
        self._shape_parameter = shape_parameter

    def conductance(self, water_potential, leaf_water_potential=None, soil_water_potential=None):
        """

        @param water_potential: (MPa)
        @return: conductance (mmol m-2 s-1 MPa-1)
        """

        normalised_conductance = 1/(1+power(water_potential/self._water_potential_at_half_conductance,
                                            self._shape_parameter))

        return self._k_max * normalised_conductance

    def water_potential_from_conductivity_loss_fraction(self, conductivity_loss_fraction):
        """

        @param conductivity_loss_fraction: (unitless)
        @return: (MPa)
        """

        conductance_fraction = 1 - conductivity_loss_fraction

        return self._water_potential_at_half_conductance * power(1/conductance_fraction - 1,
                                                                 1/self._shape_parameter)

    @property
    def sensitivity_parameter(self):
        """
        @return: (MPa)
        """
        return self._water_potential_at_half_conductance

    @property
    def shape_parameter(self):
        """
        @return: (Unitless)
        """
        return self._shape_parameter

# -- Model creation functions ------------------------------------------------


def SOX_conductance_model_from_conductance_loss_at_goven_water_potentials(
    maximum_conductance,
    water_potential_1,
    water_potential_2,
    conductance_loss_fraction_1,
    conductance_loss_fraction_2,
    critical_conductance_loss_fraction = 0.9,
    xylem_recovery_water_potnetial = 0.,
    PLC_damage_threshold = 0.1
    ):

    maximum_conductance, water_potential_at_half_conductance, shape_parameter = (
        SOX_conductance_parameters_from_conductance_loss_at_goven_water_potentials(
            maximum_conductance,
            water_potential_1,
            water_potential_2,
            conductance_loss_fraction_1,
            conductance_loss_fraction_2
            )
        )

    return SOXHydraulicConductanceModel(maximum_conductance,
                                        water_potential_at_half_conductance,
                                        shape_parameter,
                                        critical_conductance_loss_fraction,
                                        xylem_recovery_water_potnetial
                                        )


def SOX_conductance_parameters_from_conductance_loss_at_goven_water_potentials(
    maximum_conductance,
    water_potential_1,
    water_potential_2,
    conductance_loss_fraction_1,
    conductance_loss_fraction_2,
    ):

    conductance_fraction_1 = 1 - conductance_loss_fraction_1
    conductance_fraction_2 = 1 - conductance_loss_fraction_2

    shape_parameter = (log(1/conductance_fraction_1 - 1) - log(1/conductance_fraction_2 - 1)
                       / (log(water_potential_1 / water_potential_2)))

    water_potential_at_half_conductance = (water_potential_1
                                           * power(1/conductance_fraction_1 - 1, -1/shape_parameter))

    return maximum_conductance, water_potential_at_half_conductance, shape_parameter
