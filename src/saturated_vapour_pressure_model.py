"""
-------------------------------------------------------------------------

-------------------------------------------------------------------------
"""


class SaturatedVapourPressureModel:

    _polynomial_constants: list

    def __init__(self, polynomial_constants: list = [6.11213476,
                                                     4.44007856E-1,
                                                     1.43064234E-2,
                                                     2.64461437E-4,
                                                     3.05903558E-6,
                                                     1.96237241E-8,
                                                     8.92344772E-11,
                                                    -3.73208410E-13,
                                                     2.09339997E-16]):
        """
        Models vapour pressure of air as a function of temperature using a polynomial.

        e_{sat}(T) = 100[a_0 + T(a_1 + T(...(a_{N-2} + T(a_{N-1} + Ta_N))...))]

        :param polynomial_constants: i.e. a_0, ..., a_N
        """

        self._polynomial_constants = polynomial_constants

    def vapour_pressure(self, temperature):

        """
        Calculate the vapour pressure of air at a given temperature using the polynomial:

        e_{sat}(T) = 100[a_0 + T(a_1 + T(...(a_{N-1} + Ta_N)...))]

        :param temperature: (K)
        :return:
        """

        result = 0

        for current in reversed(self._polynomial_constants):
            result *= temperature
            result += current

        return 100*result

