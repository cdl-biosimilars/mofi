"""
    glyco_tools.py
    
    Helper functions/classes for glycan compositions and masses.
    
    Authors: Stefan Senn, Wolfgang Skala
    
    (c) 2017 Christian Doppler Laboratory for Biosimilar Characterization
    
"""
import os
import re
import copy
import pickle
import pandas as pd

import mass_tools

# list of monosaccharides
monosaccharides = ["Hex", "HexNAc", "Neu5Ac", "Neu5Gc", "Fuc", "Pent", "N-core", "O-core"]

# mapping from monosaccharide abbreviations to single-letter codes
mono_single_letter = {"Hex": "H",
                      "Fuc": "F",
                      "HexNAc": "N",
                      "Pent": "P",
                      "Neu5Ac": "S",
                      "Neu5Gc": "G"}

# mapping from monosaccharide single-letter codes to abbreviations
single_letter_mono = {v: k for k, v in mono_single_letter.items()}

# # composition of typical mAB glycans  TODO remove
# mab_glycan_atoms = {"Man5": {"C": 46, "H": 76, "N": 2, "O": 35},
#                     "G0": {"C": 50, "H": 82, "N": 4, "O": 35},
#                     "G0F": {"C": 56, "H": 92, "N": 4, "O": 39},
#                     "G1": {"C": 56, "H": 92, "N": 4, "O": 40},
#                     "G1F": {"C": 62, "H": 102, "N": 4, "O": 44},
#                     "G2": {"C": 62, "H": 202, "N": 4, "O": 45},
#                     "G2F": {"C": 68, "H": 112, "N": 4, "O": 49},
#                     "G2FSA": {"C": 79, "H": 129, "N": 5, "O": 57},
#                     "G2FSA2": {"C": 90, "H": 146, "N": 6, "O": 65}}

# composition of monosaccharides
mono_glycan_atoms = {"Hex": {"C": 6, "H": 10, "O": 5},
                     "Fuc": {"C": 6, "H": 10, "O": 4},
                     "HexNAc": {"C": 8, "H": 13, "O": 5, "N": 1},
                     "Pent": {"C": 5, "H": 8, "O": 4},
                     "Neu5Ac": {"C": 11, "H": 17, "O": 8, "N": 1},
                     "Neu5Gc": {"C": 11, "H": 17, "O": 9, "N": 1},
                     "O-core": {"C": 14, "H": 23, "O": 10, "N": 1},
                     "N-core": {"C": 34, "H": 56, "O": 25, "N": 2}}

# Formula objects for monosaccharides in mono_glycan_atoms
# Example: "Hex": Formula("C6 H10 O5")
glycan_formula = {mg: mass_tools.Formula(mono_glycan_atoms[mg]) for mg in mono_glycan_atoms}

fc_glycans = {"Man5": [("Hex", 2), ("N-core", 1)],
              "G0": [("HexNAc", 2), ("N-core", 1)],
              "G1": [("Hex", 1), ("HexNAc", 2), ("N-core", 1)],
              "G2": [("Hex", 2), ("HexNAc", 2), ("N-core", 1)],
              "G0F": [("HexNAc", 2), ("Fuc", 1), ("N-core", 1)],
              "G1F": [("Hex", 1), ("HexNAc", 2), ("Fuc", 1), ("N-core", 1)],
              "G2F": [("Hex", 2), ("HexNAc", 2), ("Fuc", 1), ("N-core", 1)],
              "G2FSA1": [("Hex", 2), ("HexNAc", 2), ("Fuc", 1), ("Neu5Ac", 1), ("N-core", 1)],
              "G2FSA2": [("Hex", 2), ("HexNAc", 2), ("Fuc", 1), ("Neu5Ac", 2), ("N-core", 1)]}


def glycanlist_to_modlist(glycans, max_counts=2):
    """
    Convert a list of glycans to a list of modifications.

    :param glycans: Dict {name: composition} of glycans; format like fc_glycans
    :param max_counts: maximum count for the modification search
    :return: A list of triples: (1) Name of the glycan, (2) its mass, (3) max count
    """
    result = []
    for name, composition in glycans.items():
        formula = mass_tools.combine_formulas([glycan_formula[name] * count for name, count in composition])
        result.append((name, formula.mass, max_counts))
    return result


def glycanlist_generator(glycans):
    """
    Generator that convert a list of glycans to a list of modifications

    :param glycans: Dict {name: composition} of glycans; format like fc_glycans
    :return: name, composition
    """
    for name, monomers in glycans.items():
        composition = []
        for monomer, count in monomers:
            if count == 1:
                composition.append(monomer)
            else:
                composition.append(str(count) + " " + monomer)
        yield name, ", ".join(composition)


# class Nglycan(object):
#
#     def __init__(self, ginput, name=None):
#         try:
#             if type(ginput) == dict:
#                 self._gdict = ginput
#             elif type(ginput) == str:
#                 pattern = re.compile(r"([NHFPSG]\d+)")
#                 self._gdict = {single_letter_mono[g[0]]: int(g[1:]) for g in pattern.findall(ginput)}
#         except:
#             print("Unable to recognize input")
#         if name:
#             self._tname = name
#         else:
#             self._tname = None
#         self.calculate()
#
#     def remove_tname(self):
#         self._tname = None
#
#     def make_name(self):
#         if self._gdict:
#             self._name = "".join(["{}{}".format(*(mono_single_letter[i[0]], i[1])) for i in self._gtuple])
#         else:
#             self._name = "Unglycosylated"
#
#     def calculate(self):
#         self._gdict = {k:v for k,v in self._gdict.items() if v}
#         self._form_dict = {}
#         for i in self._gdict:
#             self._form_dict[(i, self._gdict[i])]= glycan_formula[i]*self._gdict[i]
#         self._formula = mass_tools.Formula("")
#         for i in self._form_dict.values():
#             self._formula += i
#         self.average_mass = self._formula.average_mass
#         self.monoisotopic_mass = self._formula.monoisotopic_mass
#         self._gtuple = tuple(sorted(self._gdict.items()))
#         self.make_name()
#
#     def sialylated(self):
#         self.sialidase()
#         return self
#
#     def sialidase(self):
#         try:
#             del self._gdict["Neu5Ac"]
#             self.calculate()
#         except:
#             pass
#
#     def remove_core(self):
#         try:
#             self._gdict["Hex"] = self._gdict["Hex"]-3
#             self._gdict["HexNAc"] = self._gdict["HexNAc"]-2
#             self.calculate()
#         except:
#             print("No core")
#
#     def endoS(self):
#         try:
#             if self._gdict["Fuc"] != 0:
#                 self.fucosidase()
#             self._gdict["HexNAc"] = self._gdict["HexNAc"]-1
#             self.calculate()
#         except:
#             pass
#
#     def fucosidase(self):
#         try:
#             del self._gdict["Fuc"]
#             self.calculate()
#         except:
#             pass
#
#     def strip(self):
#         self.fucosidase()
#         self.sialidase()
#         self.remove_core()
#
#     def __hash__(self):
#         return(hash(self._gtuple))
#
#     def __eq__(self, other):
#         return(self._gtuple == other._gtuple)
#
#     def __ne__(self, other):
#         return(not(self == other))
#
#     def __lt__(self, other):
#         if self.average_mass < other.average_mass:
#             return -1
#         elif self.average_mass > other.average_mass:
#             return 1
#         else:
#             return 0
#
#     def __repr__(self):
#         # if hasattr(self, "_name"):
#         restr = "{:s}".format(self._name)
#         # else:
#             # restr = ""
#         # restr += "".join(["{}({})".format(*i)for i in self._gtuple])
#         restr += "|m{:.2f}".format(self.monoisotopic_mass)
#         if self._tname:
#             restr += "|{:s}".format(self._tname)
#         return(restr)
#
#
# def dataframe_to_modlist(df):
#     """
#     Function turns N-glycan DataFrame into a list of modifications.
#
#     :param df: Dataframe with the following format:
#                   Name	  N149  N171  Nr.
#                0  F1H3N4     0     2    2
#                1  F1H4N4     0     2    2
#
#     :returns: A list of tuples. Each tuple describes a modification:
#               (1) name, (2) mass, (3) maximum occurrences
#     """
#
#     result = []
#     for i, r in df.iterrows():
#         mass = Nglycan(r["Name"]).average_mass
#         result.append((r["Name"], mass, r["Nr."]))
#     return result
