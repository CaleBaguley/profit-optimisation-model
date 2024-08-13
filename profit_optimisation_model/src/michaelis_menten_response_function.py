"""
--------------------------------------------------------------------------------------
    This file contains the Michaelis-Menten response function and specific uses
including CO2 uptake during carboxylation and O uptake during oxygenation by Rubisco.
--------------------------------------------------------------------------------------
"""


def michaelis_menten_response_function(substrate_concentration, maximum_rate, michaelis_menten_constant):
    """
    Returns the enzyme reaction rate as a function of substrate concentration.
    @param substrate_concentration:
    @param maximum_rate:
    @param michaelis_menten_constant:
    @return: Enzyme reaction rate (same units as maximum rate)
    """
    return maximum_rate * substrate_concentration / (substrate_concentration + michaelis_menten_constant)


def michaelis_menten_constant(intercellular_concentration_a, michaelis_menten_constant_a, michaelis_menten_constant_b):
    """
    Calculates the Michaelis-Menten constant for the reaction transforming a to b
    @param intercellular_concentration_a: Concentration of initial material
    @param michaelis_menten_constant_a: Michaelis_Menten constant of initial material
    @param michaelis_menten_constant_b: Michaelis_Menten constant of produced material
    @return: Mechaelis-Menten constant for the reaction
    """
    return michaelis_menten_constant_b * (1 + intercellular_concentration_a / michaelis_menten_constant_a)

