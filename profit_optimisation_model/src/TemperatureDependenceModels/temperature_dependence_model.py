"""
------------------------------------------------------------------------
Base class for temperature dependent models
------------------------------------------------------------------------
"""


from profit_optimisation_model.src.conversions import degrees_centigrade_to_kelvin

from numpy import full, ndarray, zeros, float64


class TemperatureDependenceModel:
    _value_at_25C: float

    def __init__(self, value_at_25C: float):
        self._value_at_25C = value_at_25C

    def get_value_at_temperature(self, temperature):

        if(type(temperature) is ndarray):
            return full(temperature.shape(), self._value_at_25C)

        return self._value_at_25C

    @property
    def get_value_at_25C(self):
        return self._value_at_25C

class LowTemperatureAdjustedModel(TemperatureDependenceModel):
    _base_temperature_dependent_model: TemperatureDependenceModel
    _lower_bound: float
    _upper_bound: float

    def __init__(self, base_temperature_dependent_model: TemperatureDependenceModel,
                 lower_bound_C: float = 0., upper_bound_C: float = 10.):
        """
        Temperature dependent model designed to force value to zero as temperature reduces from
         upper_bound_C to lower_bound_C. Value is zero for all temperature bellow lower_bound_C
        @param base_temperature_dependent_model: temperature dependent model above upper_bound_C
        @param lower_bound_C: temperature below which value is zero
        @param upper_bound_C: temperature above which value matches base_temperature_dependent_model
        """
        super().__init__(base_temperature_dependent_model.get_value_at_25C)
        self._base_temperature_dependent_model = base_temperature_dependent_model
        self._lower_bound = degrees_centigrade_to_kelvin(lower_bound_C)
        self._upper_bound = degrees_centigrade_to_kelvin(upper_bound_C)

    def get_value_at_temperature(self, temperature):

        if (type(temperature) is float or type(temperature) is float64):
            return self._get_value_at_single_temperature(temperature)

        elif (type(temperature) is ndarray):
            values = zeros(len(temperature))

            for i in range(len(temperature)):
                values[i] = self._get_value_at_single_temperature(temperature[i])

            return values

        return None

    def _get_value_at_single_temperature(self, temperature):
        # value is zero for temperature bellow the lower bound
        if temperature <= self._lower_bound:
            return 0.

        # get the unmodified value from the base temperature model
        base_value = self._base_temperature_dependent_model.get_value_at_temperature(temperature)

        # if the temperature is above the upper bound the value from the base model is unchanged
        if temperature >= self._upper_bound:
            return base_value

        # scale the value from the base temperature model if temperature is within the upper and lower bounds.
        return base_value * (temperature - self._lower_bound) / (self._upper_bound - self._lower_bound)
