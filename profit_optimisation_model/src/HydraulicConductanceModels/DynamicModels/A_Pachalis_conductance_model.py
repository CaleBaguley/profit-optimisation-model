"""
-----------------------------------------------------------------------------------------
Implementation of xylem impairment, as described by A. Pachalis et al. (2023). The
specifics of this model are described in subsection 2.6 of the paper. Source:
https://pubmed.ncbi.nlm.nih.gov/37962234/
-----------------------------------------------------------------------------------------
"""

from src.HydraulicConductanceModels.hydraulic_conductance_model import HydraulicConductanceModel
from numpy import array as np_array
from numpy import ndarray, float64
from numpy import linspace, ones, zeros
from numpy import clip

class APachalisConductanceModel(HydraulicConductanceModel):

    _base_vulnerability_curve: HydraulicConductanceModel
    _num_ages: int
    _time_step_size: float
    _growth_rate: float
    _turnover_rate: float

    # Arrays to stor information on the xylem as a function of their age
    _xylem_conductance : np_array
    _xylem_population = np_array
    _xylem_age = np_array

    def __init__(self,
                 base_vulnerability_curve: HydraulicConductanceModel,
                 num_ages: int,
                 time_step_size: float,
                 growth_rate: float = 1.0,
                 turnover_rate: float = 0.01):

        """
        @param base_vulnerability_curve: (HydraulicConductanceModel)
        @param num_ages: (int) Number of age groups to consider
        @param time_step_size: (float) (s) Time step size for the model
        @param growth_rate: (float) (s-1) Growth rate of new xylem
        @param turnover_rate: (float) (s-1) Turnover rate of the xylem population
        """

        self._base_vulnerability_curve = base_vulnerability_curve
        self._num_ages = num_ages
        self._time_step_size = time_step_size
        self._growth_rate = growth_rate
        self._turnover_rate = turnover_rate

        # Setup arrays to contain age information
        self._xylem_age = linspace(0, time_step_size*num_ages, num_ages)
        self._xylem_conductance = ones(num_ages) * base_vulnerability_curve.maximum_conductance
        self._xylem_population = zeros(num_ages)

        self.initialise_xylem_population(growth_rate, turnover_rate)

        super().__init__(base_vulnerability_curve.maximum_conductance,
                         base_vulnerability_curve.critical_conductance_loss_fraction,
                         base_vulnerability_curve.xylem_recovery_water_potnetial,
                         base_vulnerability_curve.PLC_damage_threshold)

    def conductance(self, water_potential):

        """
        Calculate the conductance of the xylem given a water potential.
        @param water_potential: (MPa)
        @return: conductance (mmol m-2 s-1 MPa-1)
        """

        if(type(water_potential) is float
           or type(water_potential) is float64):

            # Get the healthy xylem conductance imposed by the base vulnerability curve
            current_maximum_conductance = self._base_vulnerability_curve.conductance(water_potential)

            # Produce a temporary xylem conductance array that is clipped to the current maximum conductance
            tmp_xylem_conductance = clip(self._xylem_conductance, 0, current_maximum_conductance)

            # Return the sum of the conductance of all xylem
            return sum(tmp_xylem_conductance * self._xylem_population)

        elif(isinstance(water_potential, ndarray)):
            conductances = zeros(water_potential.shape)

            for i in range(len(water_potential)):
                conductances[i] = self.conductance(water_potential[i])

            return conductances

        else:
            raise Exception("Invalid input type for water_potential, {}. Must be float or np.ndarray"
                            .format(type(water_potential)))

    def water_potential_from_conductivity_loss_fraction(self, conductivity_loss_fraction):

        """
        Calculate the water potential given a conductivity loss fraction.
        @param conductivity_loss_fraction: (unitless)
        @return: (MPa)
        """

        # Calculate the conductance of the xylem given the conductivity loss fraction
        target_conductance = conductivity_loss_fraction * self.maximum_conductance

        # Water potential if all xylem were healthy
        water_potential_healthy = (self._base_vulnerability_curve.water_potential_from_conductance(target_conductance)
                       * self.total_xylem_popuation)

        # Step through water potentials, starting from psi_healthy, until the conductance is larger than the
        # target conductance. The water potential at this point is the water potential at the target conductance.
        N = 1000
        for i in range(1, N):
            # Calculate the current water potential
            current_water_potential = water_potential_healthy * (1 - i/N)

            # calculate the current conductance
            current_conductance = self.conductance(current_water_potential)

            # If the current conductance is larger than the target conductance, return the current water potential
            if current_conductance > target_conductance + 1e-6:
                return current_water_potential

        return 0.0

    def _update_xylem_conductance(self, current_psi, timestep):

        """
        Update the conductance of the xylem as a function of age given a current water potential and timestep.
        @param current_psi:
        @param timestep:
        """

        if(timestep != self._time_step_size):
            raise Exception("Timestep size does not match the model timestep size")

        # Get the current conductance from the base vulnerability curve and current water potential
        current_conductance = self._base_vulnerability_curve.conductance(current_psi)

        # Clip those xylem ages that have a conductance higher than the current conductance
        self._xylem_conductance = clip(self._xylem_conductance, 0, current_conductance)

        # Age the xylem by shifting the conductance values one age group
        self._xylem_conductance[1:] = self._xylem_conductance[:-1]
        self._xylem_conductance[0] = current_conductance

    def _update_xylem_population(self, growth_rate, turnover_rate,  timestep):

        """
        Update the xylem population as a function of age given a growth rate and turnover rate.
        @param growth_rate:
        @param turnover_rate:
        @param timestep:
        """

        if(timestep != self._time_step_size):
            raise Exception("Timestep size does not match the model timestep size")

        # Age the xylem population by shifting the population values one age group
        self._xylem_population[1:] = self._xylem_population[:-1]

        # Apply turnover to the entire population
        self._xylem_population -= self._xylem_population * turnover_rate

        # Set the new population of the youngest xylem age group
        self._xylem_population[0] = growth_rate * timestep

    def update_xylem_damage(self, water_potential, timestep, transpiration_rate):

        """
        Update the xylem conductance and population given a current water potential, timestep and transpiration rate.
        @param water_potential: (MPa)
        @param timestep: (s)
        @param transpiration_rate: (mmol m-2 s-1)
        """

        # Update the xylem conductance and population
        self._update_xylem_conductance(water_potential, timestep)
        self._update_xylem_population(self._growth_rate, self._turnover_rate, timestep)

        return True

    def initialise_xylem_population(self, initial_growth_rate, initial_turnover_rate):

        """
        Reset the xylem population as a function of age given a fixed growth rate and turnover rate.
        @param initial_growth_rate: (float) Initial growth rate of the xylem population
        @param initial_turnover_rate: (float) Initial turnover rate of the xylem population
        """

        # Reset the xylem population and conductance arrays
        self._xylem_population = zeros(self._num_ages)

        # Update the population _num_ages times to reach a steady state
        for i in range(self._num_ages):
            self._update_xylem_population(initial_growth_rate, initial_turnover_rate, self._time_step_size)

        return self._xylem_population

    def reset_xylem_damage(self):

        """
        Reset the xylem population and conductance arrays to their starting form.
        """

        # reset the xylem population
        self.initialise_xylem_population(self._growth_rate, self._turnover_rate)

        # reset the xylem conductance
        self._xylem_conductance = ones(self._num_ages) * self._base_vulnerability_curve.maximum_conductance

        return None

    @property
    def maximum_conductance(self):
        """
        @return: (mmol m-2 s-1 MPa-1)
        """

        # Calculates the maximum conductance of all xylem. Essentially the sum of all xylem conductance's
        # if each xylem had maximum conductance.
        return self._base_vulnerability_curve.maximum_conductance * self.total_xylem_popuation

    @property
    def xylem_conductivity(self):
        """
        @return: (mmol m-2 s-1 MPa-1)
        """
        return self._xylem_conductance

    @property
    def xylem_population(self):
        """
        @return: (unitless)
        """
        return self._xylem_population

    @property
    def xylem_age(self):
        """
        @return: (s)
        """
        return self._xylem_age

    @property
    def PLC(self):
        """
        @return: (unitless)
        """

        return self.conductance(0.) / self.maximum_conductance

    @property
    def total_xylem_popuation(self):
        """
        @return: (unitless)
        """

        return sum(self._xylem_population)
