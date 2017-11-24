"""
Functions/classes for handling atomic masses.
"""

import re
import numpy as np
import pandas as pd
from mofi import configure
import pandas.core.series


def formstring_to_composition(formstring):
    """
    Converts an elemental composition string to a Series.

    :param str formstring: collection of elements followed by their counts
                           (example: "C50 H100 N-3 Cl")
    :return: a Series labeled by the element
             (example: C: 50, H: 100, N: -3, Cl: 1)
    :rtype: pd.Series
    :raises ValueError: if the composition string is invalid
    """

    pattern = re.compile(r"([A-Z][a-z]?)(-?\d*)$")
    composition = {}
    for formula_part in formstring.split():
        match = pattern.match(formula_part)
        if match:
            atom, count = match.groups()
            if count:
                composition[atom] = int(count)
            else:
                composition[atom] = 1
        else:
            raise ValueError("Invalid formula part: {}".format(formula_part))
    return pd.Series(composition)


class Formula:
    """
    A molecular formula comprising the elements C, H, N, O, P and S.

    :ivar pd.Series _composition: atomic composition
    :ivar float _mass: atomic mass

    .. automethod:: __init__
    .. automethod:: __add__
    .. automethod:: __sub__
    .. automethod:: __mul__
    .. automethod:: __rmul__
    .. automethod:: __bool__
    .. automethod:: __repr__
    """

    def __init__(self, composition=None):
        """
        Create a new Formula.

        :param composition: one of the following:

               * a pd.Series like ``pd.Series(dict(C=6, H=12, O=6))``
               * a dict like ``{"C": 6; "H": 12; "O": 6}``
               * a string like ``"C6 H12 O6".``
        :raises ValueError: if composition is of a different type
        """
        if isinstance(composition, dict):
            self._composition = pd.Series(composition)
        elif isinstance(composition, pd.core.series.Series):
            self._composition = composition
        elif isinstance(composition, str):
            self._composition = formstring_to_composition(composition)
        elif composition is None:
            self._composition = pd.Series()
        else:
            raise ValueError("Could not create formula from "
                             + str(composition))
        self._composition = self._composition.replace(np.nan, 0).astype(int)
        self._update_masses()

    def change_atom_count(self, atom, count):
        """
        Change the number of a given atom type

        :param str atom: The element whose number should be changed
        :param int count: Increase the element's count by this value
        :return: nothing
        """

        self._composition[atom] = self._composition[atom] + count

    def _update_masses(self):
        """
        Update the mass from the composition.
        Changes :attr:`~Formula._mass`.

        :return: nothing
        :raises ValueError: if an unknown atom symbol is found.
        """

        self._mass = 0
        try:
            for a in self._composition.index:
                self._mass += (self._composition[a]
                               * configure.current_mass_set[a])
        except KeyError as e:
            raise ValueError("Atom symbol '" + e.args[0] + "' unknown.")

    @property
    def mass(self):
        """
        Property, calculated when the Formula is generated or the atom count
        is changed using the current settings from the configure module.

        :return: the mass of the Formula
        :rtype: float
        """
        self._update_masses()
        return self._mass

    @property
    def composition(self):
        """
        Property, calculated when the Formula is generated or the atom count
        is changed.

        :return: the composition of the Formula
        :rtype: pd.Series
        """
        return self._composition

    def __add__(self, other):
        """
        Sum two formulas.

        :param Formula other: formula to be added
        :return: the sum of self and other
        :rtype: Formula
        """
        return Formula(self._composition
                       .add(other.composition, fill_value=0)
                       .astype(int))

    def __sub__(self, other):
        """
        Subtract two formulas.

        :param Formula other: formula to be subtracted
        :return: self minus other
        :rtype: Formula
        """
        return Formula(self._composition
                       .sub(other.composition, fill_value=0)
                       .astype(int))

    def __mul__(self, factor):
        """
        Multiply each atom by a factor.

        :param int factor: multiplicative factor
        :return: the multiplied formula
        :rtype: Formula
        """
        if type(factor) == int:
            return Formula(self._composition * factor)


    def __rmul__(self, factor):
        """
        Allow the formula to be the right operand in a multiplication.

        :param int factor: multiplicative factor
        :return: the multiplied formula
        :rtype: Formula
        """
        return self.__mul__(factor)

    def __bool__(self):
        """
        Implement truth value testing.

        :return: True if the formula contains at least one atom;
                 False otherwise
        :rtype: bool
        """

        return bool(self._composition.sum())

    def __repr__(self):
        """
        Convert the formula to a human-readable string.

        :return: a string representing the formula
        :rtype: str
        """

        result = []
        for element, count in self._composition.items():
            if count == 1:
                result.append(element)
            elif count != 0:
                result.append("{:s}{:d}".format(element, count))
        return " ".join(result)
