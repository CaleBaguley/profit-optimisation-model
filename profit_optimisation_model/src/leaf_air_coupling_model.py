"""
-------------------------------------------------------------------------

-------------------------------------------------------------------------
"""

from profit_optimisation_model.src.conversions import convert_stomatal_conductance_of_water_to_carbon


class LeafAirCouplingModel:

    def __init__(self):
        """
        Leaf air coupling model that assumes perfect leaf atmosphere coupling and that the temperature of the leaf is
        equal to that of the air.
        """
        return

    def stomatal_conductance_to_water(self,
                                      transpiration,
                                      air_temperature,
                                      vapour_pressure_deficit_of_the_air,
                                      air_pressure):

        """
        Calculates stomatal conductance to water assuming perfect leaf atmosphere coupling and that the temperature of
        the is equal to that of the air.
        Note: Potentially not static in child classes.
        :param transpiration: (mmol m-2 s-1)
        :param air_temperature: (k)
        :param vapour_pressure_deficit_of_the_air: (kPa)
        :param air_pressure: (kPa)
        :return: stomatal conductance to water (mmol m-2 s-1)
        """

        return transpiration * air_pressure / vapour_pressure_deficit_of_the_air

    def stomatal_conductance_to_carbon(self,
                                       transpiration,
                                       air_temperature,
                                       vapour_pressure_deficit_of_the_air,
                                       air_pressure):

        """
        Calculates stomatal conductance to carbon assuming perfect leaf atmosphere coupling and that the temperature of
        the is equal to that of the air.
        :param transpiration: (mmol m-2 s-1)
        :param air_temperature: (k)
        :param vapour_pressure_deficit_of_the_air: (kPa)
        :param air_pressure: (kPa)
        :return: stomatal conductance to water (mmol m-2 s-1)
        """

        stomatal_conductance_to_water = self.stomatal_conductance_to_water(transpiration,
                                                                           air_temperature,
                                                                           vapour_pressure_deficit_of_the_air,
                                                                           air_pressure)

        return convert_stomatal_conductance_of_water_to_carbon(stomatal_conductance_to_water)
