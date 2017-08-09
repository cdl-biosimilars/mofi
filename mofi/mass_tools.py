"""
mass_tools.py

Helper functions/classes for atomic masses.

Authors: Stefan Senn, Wolfgang Skala

(c) 2017 Christian Doppler Laboratory for Biosimilar Characterization
"""

import re
import os
import numpy as np
import pandas as pd
from mofi import configure
import pandas.core.series


def read_massfile(massfile):
    """
    Reads a list of masses (peaks) from an external Excel or CSV file.

    The file must contain at least one column labeled "Average Mass"
    or "Average mass (mean)". If there is a column labeled
    "Relative Abundance", its values will be used for the peak heights.
    Otherwise, all peaks will have the same height.

    :param massfile: name of the file containing the mass list
    :return: a dataframe if the import was successful,
             otherwise None
    """

    filename, ext = os.path.splitext(massfile)
    try:
        if ext in [".xls", ".xlsx"]:
            inputframe = pd.read_excel(massfile)
        elif ext in [".txt", ".csv"]:
            inputframe = pd.read_csv(massfile, sep=None, engine="python")
        else:
            return None
    except (TypeError, OSError):
        return None

    # find column with average masses
    try:
        massframe = pd.DataFrame(inputframe["Average Mass"])
    except KeyError:
        try:
            massframe = pd.DataFrame(inputframe["Average Mass (mean)"])
        except KeyError:
            return None

    # find column with relative intensities; if absent provide default values
    try:
        massframe["Relative Abundance"] = inputframe["Relative Abundance"]
    except KeyError:
        massframe["Relative Abundance"] = 100

    # convert to float
    try:
        massframe = (massframe
                     .astype(float)
                     .sort_values("Average Mass", ascending=True))
    except ValueError:
        return None

    massframe.reset_index(drop=True, inplace=True)
    return massframe


def formstring_to_composition(formstring):
    """
    Converts an elemental composition string to a pandas series.
    Raises a ValueError if the composition string is invalid.

    :param formstring: collection of elements followed by their counts
                       (example: "C50 H100 N-3 Cl")
    :return: pd.Series labelled by the element
             (example: C: 50, H: 100, N: -3, Cl: 1)
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

    Properties:
        mass: calculated when the Formula is generated or the atom count is
              changed using the current settings from the configure module.

    Protected members:
        self._composition: A pandas.Series
        self._mass
    """

    def __init__(self, composition=None):
        """
        Create a new Formula.

        :param composition: one of the following:
                   - a pandas Series, like pd.Series(dict(C=6, H=12, O=6))
                   - a dict, like {"C": 6; "H": 12; "O": 6}
                   - a string, like "C6 H12 O6".
                   Raises a ValueError if composition is of a different type
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

        :param atom: The element whose number should be changed
        :param count: Increase the element's count by this value
        :return: nothing
        """

        self._composition[atom] = self._composition[atom] + count

    def _update_masses(self):
        """
        Update the mass from the composition.
        Raises a ValueError if an unknown atom symbol is found.

        Changes:
            self._mass

        :return: nothing
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
        self._update_masses()
        return self._mass

    @property
    def composition(self):
        return self._composition

    def __add__(self, other):
        return Formula(self._composition + other.composition)

    def __mul__(self, factor):
        if type(factor) == int:
            return Formula(self._composition * factor)

    def __repr__(self):
        result = []
        for element, count in self._composition.items():
            if count == 1:
                result.append(element)
            elif count > 1:
                result.append("{:s}{:d}".format(element, count))
        return " ".join(result)
