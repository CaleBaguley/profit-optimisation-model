"""
Code for the Arrhenius and peaked Arrhenius functions.
"""

from numpy import exp
from profit_optimisation_model.src.constants import MOLAR_GAS_CONSTANT
from profit_optimisation_model.src.conversions import TWENTY_FIVE_DEGREES_CENTIGRADE_IN_KELVIN
from profit_optimisation_model.src.TemperatureDependenceModels.temperature_dependence_model \
    import TemperatureDependenceModel


class ArrheniusModel(TemperatureDependenceModel):
    _activation_energy: float
    _rate_at_25_centigrade: float

    def __init__(self, rate_at_25_centigrade: float, activation_energy: float) -> object:
        super().__init__(rate_at_25_centigrade)
        self._rate_at_25_centigrade = rate_at_25_centigrade
        self._activation_energy = activation_energy

    def get_value_at_temperature(self, temperature):
        return arrhenius_function(temperature, self._rate_at_25_centigrade, self._activation_energy)



class PeakedArrheniusModel(ArrheniusModel):
    _deactivation_energy: float
    _entropy_term: float

    def __init__(self, rate_at_25_centigrade: float, activation_energy: float, deactivation_energy: float,
                 entropy_term: float) -> object:
        super().__init__(rate_at_25_centigrade, activation_energy)
        self._deactivation_energy = deactivation_energy
        self._entropy_term = entropy_term

    def get_value_at_temperature(self, temperature):
        return peaked_arrhenius_function(temperature, self._rate_at_25_centigrade, self._activation_energy,
                                         self._deactivation_energy, self._entropy_term)


def arrhenius_function(temperature_kelvin, rate_at_25_centigrade, activation_energy):
    """
    Uses the arrhenius function to calculate the rate of a biological reaction as a function of temperature.

    Variables:
    @param temperature_kelvin: Temperature (K)

    Shape parameters:
    @param rate_at_25_centigrade: Rate measured at 25 degrees centigrade.
    @param activation_energy: Activation energy (J mol^-1)

    @return: Reaction rate (same units as rate_at_25_centigrade)
    """

    exponent = 1 - TWENTY_FIVE_DEGREES_CENTIGRADE_IN_KELVIN / temperature_kelvin
    exponent *= activation_energy / (TWENTY_FIVE_DEGREES_CENTIGRADE_IN_KELVIN * MOLAR_GAS_CONSTANT)

    return rate_at_25_centigrade * exp(exponent)


def thermal_breakdown_normalised_to_twenty_five_degrees(temperature_kelvin, deactivation_energy, entropy_term):
    """
    Calculate the thermal breakdown of a biological reaction at a given temperature

    Variable
    @param temperature_kelvin: Temperature (K)

    Shape parameters
    @param deactivation_energy: Deactivation energy (J mol-1)
    @param entropy_term: Entropy term (J K-1 mol-1)

    @return: reaction rate normalised to 25 degrees centigrade (unitless)
    """
    numerator_exponent = ((entropy_term * TWENTY_FIVE_DEGREES_CENTIGRADE_IN_KELVIN - deactivation_energy)
                          / (MOLAR_GAS_CONSTANT * TWENTY_FIVE_DEGREES_CENTIGRADE_IN_KELVIN))

    numerator = 1 + exp(numerator_exponent)

    denominator_exponent = ((entropy_term * temperature_kelvin - deactivation_energy)
                            / (MOLAR_GAS_CONSTANT * temperature_kelvin))

    denominator = 1 + exp(denominator_exponent)

    return numerator / denominator


def peaked_arrhenius_function(temperature_kelvin, rate_at_25_centigrade, activation_energy,
                              deactivation_energy, entropy_term):
    """
    Calculates the rate as a function of temperature using the peaked arrhenius function

    variable
    @param temperature_kelvin: Temperature (K)

    shape parameters
    @param rate_at_25_centigrade: Rate measured at 25 degrees centigrade.
    @param activation_energy: Activation energy (J mol^-1)
    @param deactivation_energy:  Deactivation energy (J mol-1)
    @param entropy_term:  Entropy term (J K-1 mol-1)

    @return: Reaction rate (same units as rate_at_25_centigrade)
    """

    rate = arrhenius_function(temperature_kelvin, rate_at_25_centigrade, activation_energy)
    rate *= thermal_breakdown_normalised_to_twenty_five_degrees(temperature_kelvin, deactivation_energy, entropy_term)

    return rate
