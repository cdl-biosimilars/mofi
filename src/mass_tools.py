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
import configure


def read_massfile(massfile, sort_by=None):
    """
    Reads a list of masses (peaks) from an external file as generated, e.g., by Thermo BioPharma Finder.
    Supported file types: Excel, CSV
    These files should contain at least two columns with labels "Average Mass" and "Relative Abundance",
     containing the mass and relative abundance of a peak, respectively.

    :param massfile: Name of the file containing the mass list.
    :param sort_by: Sort the generated dataframe by this column
    :return: a pandas dataframe if the import was successful, otherwise None.
    """

    filename, filext = os.path.splitext(massfile)
    if filext in [".xls", ".xlsx"]:
        massframe = pd.read_excel(massfile)
    elif filext in [".txt", ".csv"]:
        massframe = pd.read_table(massfile, sep=None, header=0, engine="python")
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


def find_delta_masses(massframe, mass_difference, tolerance=2.5):
    """
    Finds all pairs of peaks in massframe whose masses differ by mass_difference +/- tolerance.

    :param massframe: A pandas dataframe containing MS peaks and associated information, as returned by read_massfile.
    :param mass_difference: The mass difference to search for.
    :param tolerance: Accept mass differences even if they deviate from mass_difference by this value.
    :return: A list of tuples, each of which describes one of the found mass differences:
             (1) index of peak 1
             (2) index of peak 2
             (3) Mass difference between the two peaks
             (4) Mass of peak 1
             (5) Mass of peak 2
             (6) The larger of the relative abundances of peak 1 or peak 2
    """

    # indices will be an array of [peak index 1, peak index 2] for all peaks
    # whose mass difference is mass_difference +- tolerance
    masses = massframe["Average Mass"]
    hits = []
    for i in range(len(masses)):
        for j in range(i + 1, len(masses)):
            delta = abs(masses[i] - masses[j])
            if mass_difference - tolerance <= delta <= mass_difference + tolerance:
                height = min(massframe.ix[i]["Relative Abundance"], massframe.ix[j]["Relative Abundance"])
                hits.append((i, j, delta, masses[i], masses[j], height))
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

    :param formstring: collection of elements followed by their counts (example: "C50 H100 N20")
    :return: pd.Series labelled by the element (example: C: 50, H: 100, N: 20)
    """

    pattern = re.compile(r'[A-Z][a-z]?\d*')
    composition = {}
    for i in pattern.findall(formstring):
        if len(i) == 1:
            composition[i] = 1
        else:
            composition[i[0]] = int(i[1:])
    return pd.Series(composition)


# def compare_masses(masslist_a, masslist_b, tolerance=0.01):
#     hits = {}
#     for i in range(len(masslist_a)):
#         for j in range(len(masslist_b)):
#             delta = abs(masslist_a[i]-masslist_b[j])
#             if delta <= tolerance:
#                 hits[(i,j)] = (masslist_a[i], masslist_b[j], delta)
#     return hits


class Formula:
    """
    A molecular formula comprising the elements C, H, N, O, P and S.

    Members:
        average_mass
        monoisotopic_mass: Those are calculated when the Formula is generated or the atom count is changed
                           using the current settings from the configure module.

    Protected members:
        self._composition: A pandas.Series
    """

    def __init__(self, input=None):
        if type(input) == dict:
            self._composition = pd.Series(input)
        elif type(input) == pd.core.series.Series:
            self._composition = input
        elif type(input) == str:
            self._composition = formstring_to_composition(input)
        else:
            self._composition = pd.Series()
        self._composition = self._composition.replace(np.nan, 0)
        for element in "CHNOPS":
            if element not in self._composition.index:
                self._composition.set_value(element, 0)

        self._composition = self._composition.astype(int)
    
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
        Update the average and monoisotopic mass from the composition

        Changes:
            self._average_mass
            self._monoisotopic_mass

        :return: nothing
        """

        self._average_mass = 0
        self._monoisotopic_mass = 0
        for a in self._composition.index:
            self._average_mass += self._composition[a] * configure.average_masses[a]
            self._monoisotopic_mass += self._composition[a] * configure.monoisotopic_masses[a]

    @property
    def average_mass(self):
        self._update_masses()
        return self._average_mass

    @property
    def monoisotopic_mass(self):
        self._update_masses()
        return self._monoisotopic_mass

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
        for element in "CHNOPS":
            if element in self._composition.index:
                if self._composition[element] == 1:
                    result.append(element)
                elif self._composition[element] > 1:
                    result.append("%s%d" % (element, self._composition[element]))
        return " ".join(result)


# if __name__ == "__main__":
#     import sequence_tools
#     print(sum(
#         [Formula(sequence_tools.amino_acid_compositions["G"]),
#          Formula(sequence_tools.amino_acid_compositions["A"])]))
