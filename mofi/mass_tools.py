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
from . import configure


def read_massfile(massfile, sort_by=None):
    """
    Reads a list of masses (peaks) from an external file as generated,
    e.g., by Thermo BioPharma Finder.
    Supported file types: Excel, CSV.
    These files should contain at least two columns with labels "Average Mass"
    and "Relative Abundance", containing the mass and relative abundance
    of a peak, respectively.

    :param massfile: Name of the file containing the mass list.
    :param sort_by: Sort the generated dataframe by this column
    :return: a pandas dataframe if the import was successful, otherwise None.
    """

    filename, filext = os.path.splitext(massfile)
    if filext in [".xls", ".xlsx"]:
        massframe = pd.read_excel(massfile)
    elif filext in [".txt", ".csv"]:
        massframe = pd.read_table(massfile,
                                  sep=None,
                                  header=0,
                                  engine="python")
    else:
        return None

    # Only select the two required column and abort if they are missing
    massframe = massframe[["Average Mass", "Relative Abundance"]]
    if len(massframe.columns) != 2:
        return None

    if sort_by:
        massframe.sort_values(sort_by, ascending=True, inplace=True)
        massframe.index = range(len(massframe))

    return massframe


def find_delta_masses(df_masses, mass_difference, tolerance=2.5):
    """
    Finds all pairs of peaks in massframe whose masses
    differ by mass_difference +/- tolerance.

    :param df_masses: A dataframe containing MS peaks and associated
                      information, as returned by read_massfile.
    :param mass_difference: The mass difference to search for.
    :param tolerance: Accept mass differences even if they deviate
                      from mass_difference by this value.
    :return: A list of tuples, each of which describes
             one of the found mass differences:
             (1) Mass of peak 1
             (2) Mass of peak 2
             (3) The larger of the relative abundances of peak 1 or peak 2
    """

    masses = np.array(df_masses["Average Mass"])
    hits = []
    for i in range(1, len(masses)):
        delta = np.abs(masses[i:] - masses[:-i])
        indices = np.where((delta >= mass_difference - tolerance)
                           & (delta <= mass_difference + tolerance))[0]
        for j in indices:
            height = min(df_masses.ix[j]["Relative Abundance"],
                         df_masses.ix[j+i]["Relative Abundance"])
            hits.append((masses[j], masses[j+i], height))
    return hits


def combine_formulas(formulas):
    """
    Calculate a compound Formula from a list of Formulas.

    :param formulas: An iterable containing Formula objects.
    :return: the combined Formula
    """

    result = Formula()
    for f in formulas:
        result += f
    return result


def formstring_to_composition(formstring):
    """
    Converts an elemental composition string to a pandas series.

    :param formstring: collection of elements followed by their counts
                       (example: "C50 H100 N20")
    :return: pd.Series labelled by the element (example: C: 50, H: 100, N: 20)
    """

    pattern = re.compile(r"([A-Z][a-z]?)(\d*)")
    composition = {}
    for atom, count in pattern.findall(formstring):
        if count:
            composition[atom] = int(count)
        else:
            composition[atom] = 1
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
        # noinspection PyUnresolvedReferences
        if type(composition) == dict:
            self._composition = pd.Series(composition)
        elif type(composition) == pd.core.series.Series:
            self._composition = composition
        elif type(composition) == str:
            self._composition = formstring_to_composition(composition)
        elif composition is None:
            self._composition = pd.Series()
        else:
            raise ValueError("Could not create formula from "
                             + str(composition))
        self._composition = self._composition.replace(np.nan, 0).astype(int)
        self._composition = self._composition

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
