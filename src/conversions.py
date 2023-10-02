"""
A set of functions for converting between various units
"""

from src.constants import MOLAR_MASS_OF_WATER
from numpy import power

# -- Conversion constants: --
# grams per unit moles
GRAMS_PER_MOLE_OF_WATER: float = 18.02  # g mol-1
GRAMS_PER_MOLE_OF_CARBON: float = 12.  # g mol-1

# radiation conversions
PHOTOSYNTHETICALLY_ACTIVE_RADIATION_PER_UNIT_SHORT_WAVE_RADIATION: float = 0.5
JOULES_PER_MICRO_MOLE_OF_LIGHT: float = 1/4.57  # J umol-1

# Seconds per unit time
SECONDS_PER_HALF_HOUR: float = 1800.  # s
SECONDS_PER_HOUR: float = 3600.  # s
SECONDS_PER_DAY: float = 86400.  # s

# Temperature conversions
ZERO_DEGREES_CENTIGRADE_IN_KELVIN: float = 273.15  # K
TWENTY_FIVE_DEGREES_CENTIGRADE_IN_KELVIN: float = 298.15  # K

# Stomatal conductance ratios
RATIO_OF_STOMATAL_CONDUCTANCE_OF_WATER_TO_CARBON: float = 1.57

# Leaf boundary conductance ratios
RATIO_OF_LEAF_BOUNDARY_CONDUCTANCE_OF_CARBON_TO_HEAT: float = 1.32
RATIO_OF_LEAF_BOUNDARY_CONDUCTANCE_OF_HEAT_TO_WATER: float = 1.075
RATIO_OF_LEAF_BOUNDARY_CONDUCTANCE_OF_CARBON_TO_WATER: float = (RATIO_OF_LEAF_BOUNDARY_CONDUCTANCE_OF_CARBON_TO_HEAT
                                                                * RATIO_OF_LEAF_BOUNDARY_CONDUCTANCE_OF_HEAT_TO_WATER)

# Latent heat of water
LATENT_HEAT_OF_VAPORISATION_OF_WATER: float = 44200  # J mol-1
LATENT_HEAT_OF_WATER: float = 2.501e6  # J kg-1

# Dict matching si unit prefixes to their scale. Used in magnitude_conversion().
symbol_magnitude_dict = {'Y': power(10., 24),
                         'Z': power(10., 21),
                         'E': power(10., 18),
                         'P': power(10., 15),
                         'T': power(10., 12),
                         'G': power(10., 9),
                         'M': power(10., 6),
                         'k': 1000.,
                         'h': 100.,
                         'da': 10.,
                         '': 1.,
                         'd': 0.1,
                         'c': 0.01,
                         'm': 0.001,
                         'u': power(10., -6),
                         'n': power(10., -9),
                         'p': power(10., -12),
                         'f': power(10., -15),
                         'a': power(10., -18),
                         'z': power(10., -21),
                         'y': power(10., -24)
                         }


def magnitude_conversion(value, symbol_current_magnitude='', symbol_target_magnitude=''):
    """
    Function to convert between si magnitudes
    @param value: value to be converted
    @param symbol_current_magnitude: current si magnitude symbol
    @param symbol_target_magnitude: target si magnitude symbol
    @return: converted value
    """

    scale: float = symbol_magnitude_dict[symbol_current_magnitude] / symbol_magnitude_dict[symbol_target_magnitude]

    return value * scale


# ---- Converting between moles and grams ----

# -- water --
def mole_water_to_gram(value_moles):
    """
    Converts from moles of water to grams
    @param value_moles: amount of water in moles
    @return: amount of water in grams
    """
    return value_moles * GRAMS_PER_MOLE_OF_WATER


def gram_water_to_mole(value_grams):
    """
    Converts from grams of water to moles
    @param value_grams: amount of water in grams
    @return: amount of water in moles
    """
    return value_grams / GRAMS_PER_MOLE_OF_WATER


# -- carbon --
def mole_carbon_to_grams(value_moles):
    """
    converts moles of carbon to grams.
    @param value_moles: moles of carbon
    @return: grams of carbon
    """
    return value_moles * GRAMS_PER_MOLE_OF_CARBON


def grams_carbon_to_moles(value_grams):
    """
    converts grams of carbon to moles.
    @param value_grams: grams of carbon
    @return: moles of carbon
    """
    return value_grams / GRAMS_PER_MOLE_OF_CARBON


# ---- Radiation conversions ----
# -- Photosynthetically active radiation and short wave conversions --
def short_wave_to_photosynthetically_active_radiation(short_wave_radiation):
    """
    convert short wave radiation to photosynthetically active radiation
    @param short_wave_radiation:
    @return: photosynthetically active radiation
    """
    return short_wave_radiation * PHOTOSYNTHETICALLY_ACTIVE_RADIATION_PER_UNIT_SHORT_WAVE_RADIATION


def photosynthetically_active_radiation_to_short_wave(photosynthetically_active_radiation):
    """
    convert short wave radiation to photosynthetically active radiation
    @param photosynthetically_active_radiation:
    @return: short wave radiation
    """
    return photosynthetically_active_radiation / PHOTOSYNTHETICALLY_ACTIVE_RADIATION_PER_UNIT_SHORT_WAVE_RADIATION


# -- micromoles of light and joules --
def micro_moles_of_light_to_joules(micro_moles_of_light):
    """
    Converts from photons of light in micromoles to energy in joules
    @param micro_moles_of_light: umol
    @return: energy in joules: J
    """
    return micro_moles_of_light * JOULES_PER_MICRO_MOLE_OF_LIGHT


def light_energy_in_joules_to_micro_moles_of_light(light_energy_joules):
    """
        Converts from energy in joules to photons of light in micromoles
        @param light_energy_joules: J
        @return: photons in micromoles: umol
        """
    return light_energy_joules / JOULES_PER_MICRO_MOLE_OF_LIGHT


# ---- Time conversions ----
# -- seconds and half an hour --
def half_hours_to_seconds(time_in_half_hours):
    """
    Convert from half hour units to seconds
    @param time_in_half_hours:
    @return: time in seconds
    """
    return time_in_half_hours * SECONDS_PER_HALF_HOUR


def seconds_to_half_hours(time_in_seconds):
    """
    Convert from seconds to half hour units
    @param time_in_seconds:
    @return:
    """
    return time_in_seconds / SECONDS_PER_HALF_HOUR


def per_second_to_per_half_hour(per_second):
    """

    @param per_second: s-1
    @return: (0.5h)-1
    """
    return per_second * SECONDS_PER_HALF_HOUR


def per_half_hour_to_per_second(per_half_hour):
    """

    @param per_half_hour: (0.5h)-1
    @return: s-1
    """
    return per_half_hour / SECONDS_PER_HALF_HOUR


# -- seconds and hours --
def hours_to_seconds(time_in_hours):
    """
    Convert from hours to seconds
    @param time_in_hours:
    @return: time in seconds
    """
    return time_in_hours * SECONDS_PER_HOUR


def seconds_to_hours(time_in_seconds):
    """
    Convert from seconds to hours
    @param time_in_seconds:
    @return: time in hours
    """
    return time_in_seconds / SECONDS_PER_HOUR


# -- seconds and days --
def days_to_seconds(time_days):
    """
    Convert time from days to seconds
    @param time_days:
    @return: time in seconds
    """
    return time_days * SECONDS_PER_DAY


def seconds_to_days(time_seconds):
    """
    Convert time from seconds to days
    @param time_seconds:
    @return: time in days
    """
    return time_seconds / SECONDS_PER_DAY


def per_day_to_per_second(per_day):
    """

    @param per_day: d-1
    @return: s-1
    """

    return per_day / SECONDS_PER_DAY


def per_second_to_per_day(per_second):
    """

    @param per_second: s-1
    @return: d-1
    """

    return per_second * SECONDS_PER_DAY


# ---- Temperature conversions ----
# kelvin and centigrade
def degrees_centigrade_to_kelvin(temperature_centigrade):
    """
    Convert from centigrade to kelvin
    @param temperature_centigrade:
    @return: temperature in kelvin
    """
    return temperature_centigrade + ZERO_DEGREES_CENTIGRADE_IN_KELVIN


def degrees_kelvin_to_centigrade(temperature_kelvin):
    """
    Convert from kelvin to centigrade
    @param temperature_kelvin:
    @return: temperature in centigrade
    """
    return temperature_kelvin - ZERO_DEGREES_CENTIGRADE_IN_KELVIN


# ---- Stomatal Conductance ----
# -- water and carbon --
def convert_stomatal_conductance_of_carbon_to_water(stomatal_conductance_of_carbon):
    """
    Converts from the stomatal conductance of carbon to that of water
    @param stomatal_conductance_of_carbon:
    @return: stomatal conductance of water
    """
    return stomatal_conductance_of_carbon * RATIO_OF_STOMATAL_CONDUCTANCE_OF_WATER_TO_CARBON


def convert_stomatal_conductance_of_water_to_carbon(stomatal_conductance_of_water):
    """
    Converts from the stomatal conductance of water to that of carbon
    @param stomatal_conductance_of_water:
    @return: stomatal conductance of carbon
    """
    return stomatal_conductance_of_water / RATIO_OF_STOMATAL_CONDUCTANCE_OF_WATER_TO_CARBON


# ---- Leaf boundary conductance ----
# -- Carbon and heat --
def convert_leaf_boundary_conductance_of_carbon_to_that_of_heat(conductance_of_carbon):
    """
    Convert the leaf boundary conductance of carbon to that of heat
    @param conductance_of_carbon: leaf boundary conductance of carbon
    @return: leaf boundary conductance of heat
    """
    return conductance_of_carbon * RATIO_OF_LEAF_BOUNDARY_CONDUCTANCE_OF_CARBON_TO_HEAT


def convert_leaf_boundary_conductance_of_heat_to_that_of_carbon(conductance_of_heat):
    """
    Convert the leaf boundary conductance of heat to that of carbon
    @param conductance_of_heat: leaf boundary conductance of heat
    @return: leaf boundary conductance of carbon
    """
    return conductance_of_heat / RATIO_OF_LEAF_BOUNDARY_CONDUCTANCE_OF_CARBON_TO_HEAT


# -- Carbon and water --
def convert_leaf_boundary_conductance_of_carbon_to_that_of_water(conductance_of_carbon):
    """
    Convert the leaf boundary conductance of carbon to that of water
    @param conductance_of_carbon: leaf boundary conductance of carbon
    @return: leaf boundary conductance of water
    """
    return conductance_of_carbon * RATIO_OF_LEAF_BOUNDARY_CONDUCTANCE_OF_CARBON_TO_WATER


def convert_leaf_boundary_conductance_of_water_to_that_of_carbon(conductance_of_water):
    """
    Convert the leaf boundary conductance of water to that of carbon
    @param conductance_of_water: leaf boundary conductance of water
    @return: leaf boundary conductance of carbon
    """
    return conductance_of_water / RATIO_OF_LEAF_BOUNDARY_CONDUCTANCE_OF_CARBON_TO_WATER


# -- Heat and water --
def convert_leaf_boundary_conductance_of_heat_to_that_of_water(conductance_of_heat):
    """
    Convert the leaf boundary conductance of heat to that of water
    @param conductance_of_heat: leaf boundary conductance of heat
    @return: leaf boundary conductance of water
    """
    return conductance_of_heat * RATIO_OF_LEAF_BOUNDARY_CONDUCTANCE_OF_HEAT_TO_WATER


def convert_leaf_boundary_conductance_of_water_to_that_of_heat(conductance_of_water):
    """
    Convert the leaf boundary conductance of water to that of heat
    @param conductance_of_water: leaf boundary conductance of water
    @return: leaf boundary conductance of heat
    """
    return conductance_of_water / RATIO_OF_LEAF_BOUNDARY_CONDUCTANCE_OF_HEAT_TO_WATER


# ---- Stomatal conductance ----
# -- Carbon and water --
def convert_stomatal_conductance_of_carbon_to_that_of_water(conductance_of_carbon):
    """
    Convert the leaf boundary conductance of carbon to that of water
    @param conductance_of_carbon: leaf boundary conductance of carbon
    @return: leaf boundary conductance of water
    """
    return conductance_of_carbon * RATIO_OF_STOMATAL_CONDUCTANCE_OF_WATER_TO_CARBON


def convert_stomatal_conductance_of_water_to_that_of_carbon(conductance_of_water):
    """
    Convert the leaf boundary conductance of water to that of carbon
    @param conductance_of_water: leaf boundary conductance of water
    @return: leaf boundary conductance of carbon
    """
    return conductance_of_water / RATIO_OF_STOMATAL_CONDUCTANCE_OF_WATER_TO_CARBON


# ---- Transpiration and latent heat of water --
def convert_transpiration_rate_to_latent_energy(transpiration_rate, air_temperature):
    """

    @param transpiration_rate: mol m-2 s-1
    @param air_temperature: K
    @return: J m-2 s-1 or W m-2
    """

    air_temperature_C = degrees_kelvin_to_centigrade(air_temperature)

    latent_heat_of_vaporisation_of_water = (LATENT_HEAT_OF_WATER - 2.365E3 * air_temperature_C) * MOLAR_MASS_OF_WATER

    return transpiration_rate * latent_heat_of_vaporisation_of_water


def convert_latent_energy_to_transpiration(latent_energy, air_temperature):
    """

    @param latent_energy: J m-2 s-1 or W m-2
    @param air_temperature: K
    @return: mol m-2 s-1
    """

    #TODO: temperature should be in C

    air_temperature_C = degrees_kelvin_to_centigrade(air_temperature)

    latent_heat_of_vaporisation_of_water = (LATENT_HEAT_OF_WATER - 2.365E3 * air_temperature_C) * MOLAR_MASS_OF_WATER

    return latent_energy / LATENT_HEAT_OF_VAPORISATION_OF_WATER
