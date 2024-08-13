"""
-------------------------------------------------------------------------

-------------------------------------------------------------------------
"""

from src.conversions import magnitude_conversion, degrees_kelvin_to_centigrade
from numpy import exp
from numpy import where


def saturation_vapour_pressure(air_temperature_K):

    """

    @param air_temperature_K: K
    @return: kPa
    """

    air_temperature_C = degrees_kelvin_to_centigrade(air_temperature_K)

    return 0.61078 * exp((17.27 * air_temperature_C) / (237.3 + air_temperature_C))


def vapour_pressure(specific_humidity, air_pressure):

    """

    @param specific_humidity: kg kg-1
    @param air_pressure: Pa
    @return: Pa
    """

    return (specific_humidity * air_pressure) / (0.622 + (1.0 - 0.622) * specific_humidity)


def vapour_pressure_deficit(air_temperature, specific_humidity, air_pressure, minimum = 0.05):
    """

    @param air_temperature: K
    @param specific_humidity: kg kg-1
    @param air_pressure: Pa
    @param minimum: minimum vapour pressure deficit kPa
    @return: kPa
    """

    # calculate saturated vapour pressure and current vapour pressure
    saturation_vapour_pressure_values = saturation_vapour_pressure(air_temperature)  # kPa
    vapour_pressure_values = vapour_pressure(specific_humidity, air_pressure)  # Pa

    # Calculate vapour pressure deficit values and convert to kPa
    vapour_pressure_values = magnitude_conversion(vapour_pressure_values, '', 'k')  # kPa
    vapour_pressure_deficit_values = saturation_vapour_pressure_values - vapour_pressure_values  # kPa

    # Impose minimum vapour pressure deficit
    vapour_pressure_deficit_values = where(vapour_pressure_deficit_values < minimum,
                                           minimum,
                                           vapour_pressure_deficit_values)

    return vapour_pressure_deficit_values
