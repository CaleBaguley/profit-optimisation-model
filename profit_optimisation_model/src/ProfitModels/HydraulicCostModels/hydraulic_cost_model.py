"""
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
"""

from profit_optimisation_model.src.HydraulicConductanceModels.hydraulic_conductance_model \
    import HydraulicConductanceModel


class HydraulicCostModel:

    _hydraulic_conductance_model: HydraulicConductanceModel
    _critical_leaf_water_potential: float
    _critical_hydraulic_conductance: float

    def __init__(self,
                 hydraulic_conductance_model: HydraulicConductanceModel,
                 critical_leaf_water_potential: float):

        """

        @param hydraulic_conductance_model:
        @param critical_leaf_water_potential: (MPa)
        """
        self._hydraulic_conductance_model = hydraulic_conductance_model
        self._critical_leaf_water_potential = critical_leaf_water_potential

        self._critical_hydraulic_conductance = \
            hydraulic_conductance_model.conductance(critical_leaf_water_potential, 0.0, 0.0)

    # ----- Hydraulic cost ------
    def hydraulic_cost_as_a_function_of_leaf_water_potential(self, leaf_water_potentials, soil_water_potential):

        """

        @param leaf_water_potentials: 1d numpy array (MPa)
        @param soil_water_potential: float (MPa)
        @return: array of hydraulic costs
        """

        return self.hydraulic_cost(leaf_water_potentials, soil_water_potential)

    def hydraulic_cost(self, leaf_water_potential, soil_water_potential):
        """

        @param leaf_water_potential: (MPa)
        @param soil_water_potential: (MPa)
        @return: hydraulic cost (unitless)
        """

        raise Exception("hydraulic_cost method not implemented in base class.")

    def hydraulic_conductance(self, water_potential, leaf_water_potential, soil_water_potential):

        """

        @param water_potential: (MPa)
        @return: hydraulic conductance (mmol m-2 s-1 MPa-1)
        """

        return self._hydraulic_conductance_model.conductance(water_potential, leaf_water_potential,
                                                             soil_water_potential)

    def instantaneous_maximum_hydraulic_conductance(self, soil_water_potential, leaf_water_potential = None):

        """

        @param soil_water_potential: (MPa)
        @return: instantaneous maximum hydraulic conductance (mmol m-2 s-1 MPa-1)
        """

        return self.hydraulic_conductance(soil_water_potential, leaf_water_potential, soil_water_potential)

    def transpiration(self, min_water_potential, max_water_potential, steps = 100):
        """
        Calculates the transpiration rate across the water potentials using the trapezium integral approximation.
        Simply calls the same method from the hydraulic conductance model
        @param min_water_potential: (MPa)
        @param max_water_potential: (MPa)
        @param steps: number of steps for integral approximation (unitless)
        @return: (mmol m-2 s-1)
        """
        return self._hydraulic_conductance_model.transpiration(min_water_potential, max_water_potential, steps)

    def update_xylem_damage(self, water_potential, timestep, transpiration_rate, root_water_potential):
        """
        @param water_potential: (MPa)
        @param timestep: (s)
        @param transpiration_rate: (mmol m-2 s-1)
        @param root_water_potential: (MPa)
        @return: bool indicting if the conductance model has changed
        """
        return self._hydraulic_conductance_model.update_xylem_damage(water_potential,
                                                                     timestep,
                                                                     transpiration_rate,
                                                                     root_water_potential)

    @property
    def hydraulic_conductance_model(self):
        return self._hydraulic_conductance_model

    @property
    def critical_leaf_water_potential(self):
        return self._critical_leaf_water_potential

    @property
    def critical_hydraulic_conductance(self):
        return self._critical_hydraulic_conductance
