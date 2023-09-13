"""
-------------------------------------------------------------------------

-------------------------------------------------------------------------
"""

from src.saturated_vapour_pressure_model import SaturatedVapourPressureModel

class LeafAirCouplingModel:

    def __init__(self, saturated_vapour_pressure_model: SaturatedVapourPressureModel):
        """
        Leaf air coupling model that assumes perfect leaf atmosphere coupling and that the temperature of the leaf is
        equal to that of the air.
        :type saturated_vapour_pressure_model: SaturatedVapourPressureModel
        """
        self._saturated_vapour_pressure_model = saturated_vapour_pressure_model
        return

    def stomatal_conductance_to_water(self,
                                      transpiration,
                                      air_temperature,
                                      vapour_pressure_deficit_of_the_air,
                                      air_pressure):

        """
        Calculates stomatal conductance to water assuming perfect leaf atmosphere coupling and that the temperature of
        the is equal to that of the air.
        :param transpiration: (mmol m-2 s-1)
        :param air_temperature: (k)
        :param vapour_pressure_deficit_of_the_air: (kPa)
        :param air_pressure: (kPa)
        :return: stomatal conductance to water (mmol m-2 s-1)
        """

        return transpiration * air_pressure / vapour_pressure_deficit_of_the_air