"""
-----------------------------------------------------------------------------------------

-----------------------------------------------------------------------------------------
"""

from numpy import exp, power, log, abs
from numpy import linspace, trapz


class HydraulicConductanceModel:

    _maximum_conductance: float
    _critical_conductance_loss_fraction: float
    _xylem_recovery_water_potnetial: float
    _PLC_damage_threshold: float

    def __init__(self,
                 maximum_conductance: float,
                 critical_conductance_loss_fraction: float = 0.9,
                 xylem_recovery_water_potnetial: float = 0.,
                 PLC_damage_threshold: float = 0.1):

        """

        @param maximum_conductance: (mmol m-2 s-1 MPa-1)
        @param critical_conductance_loss_fraction: (unitless)
        """

        self._maximum_conductance = maximum_conductance
        self._critical_conductance_loss_fraction = critical_conductance_loss_fraction
        self._xylem_recovery_water_potnetial = xylem_recovery_water_potnetial
        self._PLC_damage_threshold = PLC_damage_threshold

    def conductance(self, water_potential):

        """

        @param water_potential: (MPa)
        @return: conductance (mmol m-2 s-1 MPa-1)
        """

        raise Exception("conductance method not implemented in base class")

    def water_potential_from_conductivity_loss_fraction(self, conductivity_loss_fraction):
        """

        @param conductivity_loss_fraction: (unitless)
        @return: (MPa)
        """

        raise Exception(" water_potential_from_conductivity_loss_fraction method not implemented in base class")

    def water_potential_from_conductance(self, conductance):
        """

        @param conductance: (mmol m-2 s-1 MPa-1)
        @return: (MPa)
        """

        conductivity_loss_fraction = 1 - conductance/self.maximum_conductance

        return self.water_potential_from_conductivity_loss_fraction(conductivity_loss_fraction)

    def transpiration(self, min_water_potential, max_water_potential, steps = 100):
        """
        Calculates the transpiration rate across the water potentials using the trapezium integral approximation
        @param min_water_potential: (MPa)
        @param max_water_potential: (MPa)
        @param steps: number of steps for integral approximation (unitless)
        @return: (mmol m-2 s-1)
        """

        water_potential_values = linspace(min_water_potential, max_water_potential, steps)

        conductance_values = self.conductance(water_potential_values)

        return trapz(conductance_values, water_potential_values)

    def update_xylem_damage(self, water_potential, timestep, transpiration_rate):
        """
        @param water_potential: (MPa)
        @param timestep: (s)
        @param transpiration_rate: (mmol m-2 s-1)
        @return: bool indicting if the model has changed
        """

        if water_potential >= self._xylem_recovery_water_potnetial:
            return self._recover_xylem(water_potential, timestep)

        conductance = self.conductance(water_potential)

        if conductance <= self.maximum_conductance * (1 - self._PLC_damage_threshold):
            return self._damage_xylem(water_potential, timestep, transpiration_rate)

        return False

    def _damage_xylem(self, water_potential, timestep, transpiration_rate):
        """
        @param water_potential: (MPa)
        @param timestep: (s)
        @param transpiration_rate: (mmol m-2 s-1)
        @return: bool indicting if the model has changed
        """
        return False

    def _recover_xylem(self, water_potential, timestep):
        """
        @param water_potential: (MPa)
        @param timestep: (s)
        @return: bool indicting if the model has changed
        """
        return False

    def reset_xylem_damage(self):
        """
        @return: None
        """
        return None

    @property
    def maximum_conductance(self):
        """
        @return: (mmol m-2 s-1 MPa-1)
        """
        return self._maximum_conductance

    @property
    def critical_conductance_loss_fraction(self):
        """
        @return: (unitless)
        """
        return self._critical_conductance_loss_fraction

    @property
    def critical_conductance(self):
        """
        @return: (mmol m-2 s-1 MPa-1)
        """
        return self.maximum_conductance * self.critical_conductance_loss_fraction

    @property
    def critical_water_potential(self):
        """
        @return: (MPa)
        """
        return self.water_potential_from_conductivity_loss_fraction(self.critical_conductance_loss_fraction)
    
    @property
    def xylem_recovery_water_potnetial(self):
        """
        @return: (MPa)
        """
        return self._xylem_recovery_water_potnetial

    @property
    def PLC_damage_threshold(self):
        """
        @return: (unitless)
        """
        return self._PLC_damage_threshold