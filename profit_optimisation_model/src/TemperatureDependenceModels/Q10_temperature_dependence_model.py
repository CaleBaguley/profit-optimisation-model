"""
-------------------------------------------------------------------------
Implementation of the Q10 temperature dependence model
-------------------------------------------------------------------------
"""

from profit_optimisation_model.src.TemperatureDependenceModels.temperature_dependence_model \
    import TemperatureDependenceModel
from profit_optimisation_model.src.conversions import TWENTY_FIVE_DEGREES_CENTIGRADE_IN_KELVIN
from numpy import power


class Q10TemperatureDependenceModel(TemperatureDependenceModel):
    _Q10_ratio: float

    def __init__(self, value_at_25C: float, Q10_parameter: float):
        """

        @param value_at_25C: measured value at 25C
        @param Q10_parameter: shape parameter
        """
        super().__init__(value_at_25C)
        self._Q10_ratio = Q10_parameter

    def get_value_at_temperature(self, temperature):
        power_component = (temperature - TWENTY_FIVE_DEGREES_CENTIGRADE_IN_KELVIN) / 10
        return self._value_at_25C * power(self._Q10_ratio, power_component)
