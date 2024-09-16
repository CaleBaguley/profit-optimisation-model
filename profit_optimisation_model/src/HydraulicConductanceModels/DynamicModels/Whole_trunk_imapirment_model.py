"""
-----------------------------------------------------------------------------------------

-----------------------------------------------------------------------------------------
"""
from jupyter_server.transutils import trans

from profit_optimisation_model.src.HydraulicConductanceModels.hydraulic_conductance_model import (
    HydraulicConductanceModel)

class WholeTrunkImapirmentModel(HydraulicConductanceModel):

    _base_conductance_model : HydraulicConductanceModel
    _psi_leaf_extreme : float
    _psi_root_extreme : float

    def __init__(self,
                 base_conductance_model : HydraulicConductanceModel,
                 xylem_recovery_water_potnetial: float = 0.0,
                 PLC_damage_threshold: float = 0.1):

        self._base_conductance_model = base_conductance_model
        self._psi_leaf_extreme = 0.
        self._psi_root_extreme = 0.

        super().__init__(base_conductance_model.maximum_conductance,
                         base_conductance_model.critical_conductance_loss_fraction,
                         xylem_recovery_water_potnetial,
                         PLC_damage_threshold)

    def conductance(self, water_potential, leaf_water_potential=None, soil_water_potential=None):

        if(leaf_water_potential is None):
            raise ValueError("Leaf water potential must be provided for the whole trunk impairment model")
        elif(soil_water_potential is None):
            raise ValueError("Soil water potential must be provided for the whole trunk impairment model")

        tmp_psi_leaf_extreme = min(self._psi_leaf_extreme, leaf_water_potential)
        tmp_psi_root_extreme = min(self._psi_root_extreme, soil_water_potential)

        print("tmp_psi_leaf_extreme: ", tmp_psi_leaf_extreme)
        print("tmp_psi_root_extreme: ", tmp_psi_root_extreme)

        #            (Psi - Psi_soil)
        # Psi' = ----------------------- * (Psi_leaf_extreme - Psi_root_extreme) + Psi_root_extreme
        #         (Psi_leaf - Psi_soil)

        modified_psi = (water_potential - soil_water_potential)
        modified_psi /= (leaf_water_potential - soil_water_potential)
        modified_psi *= (tmp_psi_leaf_extreme - tmp_psi_root_extreme)
        modified_psi += tmp_psi_root_extreme

        return self._base_conductance_model.conductance(modified_psi)

    def PLC(self, water_potential):

        transpiration_over_limits = self._base_conductance_model.transpiration(self._psi_leaf_extreme,
                                                                               self._psi_root_extreme)

        maximum_transpiration = self._k_max * (self._psi_leaf_extreme - self._psi_root_extreme)

        return 100*(1 - transpiration_over_limits/maximum_transpiration)

    def water_potential_from_conductivity_loss_fraction(self, conductivity_loss_fraction):
        return self._base_conductance_model.water_potential_from_conductivity_loss_fraction(conductivity_loss_fraction)

    def transpiration(self, min_water_potential, max_water_potential, steps = 100):

        # Calculate the extreme water potentials for the current calculation
        tmp_psi_leaf_extreme = min(self._psi_leaf_extreme, min_water_potential)
        tmp_psi_root_extreme = min(self._psi_root_extreme, max_water_potential)

        # Calculate the transpiration over the extreme water potentials
        transpiration = self._base_conductance_model.transpiration(tmp_psi_leaf_extreme,
                                                                   tmp_psi_root_extreme)

        # Scale the transpiration to the current water potential limits
        transpiration *= (max_water_potential - min_water_potential) / (tmp_psi_root_extreme - tmp_psi_leaf_extreme)

    def _damage_xylem(self, water_potential, timestep, transpiration_rate, root_water_potential):
        """
        @param water_potential: (MPa)
        @param timestep: (s)
        @param transpiration_rate: (mmol m-2 s-1)
        @param root_water_potential:
        @return: bool indicting if the model has changed
        """

        self._psi_leaf_extreme = min(self._psi_leaf_extreme, water_potential)
        self._psi_root_extreme = min(self._psi_root_extreme, root_water_potential)

        return False