"""
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
"""

from src.cumulative_Weibull_distribution_model import CumulativeWeibullDistribution


class HydraulicCostModel:

    _hydraulic_conductance_model: CumulativeWeibullDistribution
    _saturated_water_potential: float
    _critical_leaf_water_potential: float
    _critical_hydraulic_conductance: float

    def __init__(self,
                 hydraulic_conductance_model: CumulativeWeibullDistribution,
                 saturated_water_potential: float,
                 critical_leaf_water_potential: float):

        """

        @param hydraulic_conductance_model:
        @param saturated_water_potential:
        @param critical_leaf_water_potential:
        """
        self._hydraulic_conductance_model = hydraulic_conductance_model
        self._saturated_water_potential = saturated_water_potential
        self._critical_leaf_water_potential = critical_leaf_water_potential

        self._critical_hydraulic_conductance = (
            hydraulic_conductance_model.conductance(critical_leaf_water_potential))

    # ----- Hydraulic cost ------
    def hydraulic_cost_as_a_function_of_leaf_water_potential(self, leaf_water_potentials, soil_water_potential):

        """

        @param leaf_water_potentials: 1d numpy array (Pa)
        @param soil_water_potential: float (Pa)
        @return: array of hydraulic costs
        """

        return self.hydraulic_cost(leaf_water_potentials, soil_water_potential)

    def hydraulic_cost(self, leaf_water_potential, soil_water_potential):
        """

        @param leaf_water_potential: (Pa)
        @param soil_water_potential: (Pa)
        @return: hydraulic cost (unitless)
        """

        hydraulic_conductance = self.hydraulic_conductance(leaf_water_potential)

        instantaneous_maximum_hydraulic_conductance = \
            self.instantaneous_maximum_hydraulic_conductance(soil_water_potential)

        return ((instantaneous_maximum_hydraulic_conductance - hydraulic_conductance)
                / (instantaneous_maximum_hydraulic_conductance - self._critical_hydraulic_conductance))

    def hydraulic_conductance(self, water_potential):

        """

        @param water_potential: (Pa)
        @return: hydraulic conductance (mol m-2 s-1 Pa-1)
        """

        return self._hydraulic_conductance_model.conductance(water_potential)

    def instantaneous_maximum_hydraulic_conductance(self, soil_water_potential):

        """

        @param soil_water_potential: (Pa)
        @return: instantaneous maximum hydraulic conductance (mol m-2 s-1 Pa-1)
        """

        return self.hydraulic_conductance(soil_water_potential)

    def transpiration(self, min_water_potential, max_water_potential, steps = 100):
        """
        Calculates the transpiration rate across the water potentials using the trapezium integral approximation.
        Simply calls the same method from the hydraulic conductance model
        @param min_water_potential: (Pa)
        @param max_water_potential: (Pa)
        @param steps: number of steps for integral approximation (unitless)
        @return: (mol m-2 s-1)
        """
        return self._hydraulic_conductance_model.transpiration(min_water_potential, max_water_potential, steps)

    @property
    def hydraulic_conductance_model(self):
        return self._hydraulic_conductance_model

    @property
    def saturated_water_potential(self):
        return self._saturated_water_potential

    @property
    def critical_leaf_water_potential(self):
        return self._critical_leaf_water_potential

    @property
    def critical_hydraulic_conductance(self):
        return self._critical_hydraulic_conductance
