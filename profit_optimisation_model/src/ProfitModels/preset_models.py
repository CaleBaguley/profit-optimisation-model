"""
-------------------------------------------------------------------------
The contents of this file are designed to generate preset optimisation
models. This is useful for quick testing as it avoids the need to write
out the full model creation code.
-------------------------------------------------------------------------
"""

# -- Import the required libraries --
from profit_optimisation_model.src.HydraulicConductanceModels.cumulative_Weibull_distribution_model import (
    cumulative_Weibull_distribution_from_conductance_loss_at_given_water_potentials)
from profit_optimisation_model.src.HydraulicConductanceModels.SOX_hydraulic_conductance_model import (
    SOX_conductance_model_from_conductance_loss_at_goven_water_potentials)

from profit_optimisation_model.src.ProfitModels.HydraulicCostModels.hydraulic_cost_profit_max_model import (
    ProfitMaxHydraulicCostModel)
from profit_optimisation_model.src.ProfitModels.HydraulicCostModels.hydraulic_cost_SOX_model import (
    SOXHydraulicCostModel)

from profit_optimisation_model.src.leaf_air_coupling_model import LeafAirCouplingModel

from profit_optimisation_model.src.PhotosynthesisModels.photosynthesis_model import PhotosynthesisModel

from profit_optimisation_model.src.ProfitModels.CO2GainModels.CO2_gain_profit_max_model import ProfitMaxCO2GainModel
from profit_optimisation_model.src.ProfitModels.CO2GainModels.CO2_gain_SOX_model import SOXCO2GainModel

from profit_optimisation_model.src.ProfitModels.profit_max_model import ProfitMaxModel
from profit_optimisation_model.src.ProfitModels.SoxModel import SOXModel

# -- Define the preset models -------------------------------------------


def build_profit_max_model():

    # conductance model
    P50 = -4.739999999999999  # MPa
    P88 = -5.476604474129552  # MPa
    k_max = 0.2  # mmol m-2 s-1 MPa-1

    conductance_model = cumulative_Weibull_distribution_from_conductance_loss_at_given_water_potentials( P50,
                                                                                                         P88,
                                                                                                         k_max)

    # hydraulic cost model
    critical_water_potential = conductance_model.water_potential_from_conductivity_loss_fraction(0.95)  # MPa
    hydraulic_cost_model = ProfitMaxHydraulicCostModel(conductance_model,
                                                       critical_water_potential)

    # leaf to air conductance model
    leaf_air_coupling_model = LeafAirCouplingModel()

    # photosynthesis model
    photosynthesis_model = PhotosynthesisModel()

    # CO2 gain model
    CO2_gain_model = ProfitMaxCO2GainModel(leaf_air_coupling_model,
                                           photosynthesis_model)

    profit_optimisation_model = ProfitMaxModel(hydraulic_cost_model, leaf_air_coupling_model, CO2_gain_model)

    return profit_optimisation_model


def build_SOX_model():

    # conductance model
    P50 = -4.739999999999999  # MPa
    P88 = -5.476604474129552  # MPa
    k_max = 0.2  # mmol m-2 s-1 MPa-1
    conductance_model = SOX_conductance_model_from_conductance_loss_at_goven_water_potentials(P50,
                                                                                              P88,
                                                                                              k_max)

    # hydraulic cost model
    critical_water_potential = conductance_model.water_potential_from_conductivity_loss_fraction(0.95)  # MPa
    hydraulic_cost_model = SOXHydraulicCostModel(conductance_model,
                                                 critical_water_potential)

    # leaf to air conductance model
    leaf_air_coupling_model = LeafAirCouplingModel()

    # photosynthesis model
    photosynthesis_model = PhotosynthesisModel()

    # CO2 gain model
    CO2_gain_model = SOXCO2GainModel(leaf_air_coupling_model,
                                     photosynthesis_model)

    profit_optimisation_model = SOXModel(hydraulic_cost_model, leaf_air_coupling_model, CO2_gain_model)

    return profit_optimisation_model
