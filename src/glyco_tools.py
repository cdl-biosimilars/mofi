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
# monosaccharides = ["Fuc", "Hex", "HexNAc", "Pent", "Neu5Ac", "Neu5Gc"]

# mapping from monosaccharide abbreviations to single-letter codes
mono_single_letter = {"Hex": "H",
                      "Fuc": "F",
                      "HexNAc": "N",
                      "Pent": "P",
                      "Neu5Ac": "S",
                      "Neu5Gc": "G"}

# mapping from monosaccharide single-letter codes to abbreviations
single_letter_mono = {v: k for k, v in mono_single_letter.items()}

# composition of typical mAB glycans
mab_glycan_atoms = {"Man5": {"C": 46, "H": 76, "N": 2, "O": 35},
                    "G0": {"C": 50, "H": 82, "N": 4, "O": 35},
                    "G0F": {"C": 56, "H": 92, "N": 4, "O": 39},
                    "G1": {"C": 56, "H": 92, "N": 4, "O": 40},
                    "G1F": {"C": 62, "H": 102, "N": 4, "O": 44},
                    "G2": {"C": 62, "H": 202, "N": 4, "O": 45},
                    "G2F": {"C": 68, "H": 112, "N": 4, "O": 49},
                    "G2FSA": {"C": 79, "H": 129, "N": 5, "O": 57},
                    "G2FSA2": {"C": 90, "H": 146, "N": 6, "O": 65}}

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

fc_glycans = {"Man5": [("Hex", 5), ("HexNAc", 2)],
              "G0": [("Hex", 3), ("HexNAc", 4)],
              "G1": [("Hex", 4), ("HexNAc", 4)],
              "G2": [("Hex", 5), ("HexNAc", 4)],
              "G0F": [("Hex", 3), ("HexNAc", 4), ("Fuc", 1)],
              "G1F": [("Hex", 4), ("HexNAc", 4), ("Fuc", 1)],
              "G2F": [("Hex", 5), ("HexNAc", 4), ("Fuc", 1)],
              "G2FSA1": [("Hex", 5), ("HexNAc", 4), ("Fuc", 1), ("Neu5Ac", 1)],
              "G2FSA2": [("Hex", 5), ("HexNAc", 4), ("Fuc", 1), ("Neu5Ac", 2)]}

# typical combinations of fc_glycans found on the two N-glycosylation sites of mABs
mab_glycans = {"(Man5)/2": [("Hex", 10), ("HexNAc", 4)],
               "(G0F)/2": [("Hex", 6), ("HexNAc", 8), ("Fuc", 2)],
               "G0F/G1F": [("Hex", 7), ("HexNAc", 8), ("Fuc", 2)],
               "G0F/G2F or (G1F)2": [("Hex", 8), ("HexNAc", 8), ("Fuc", 2)],
               "G1F/G2F": [("Hex", 9), ("HexNAc", 8), ("Fuc", 2)],
               "(G2F)2": [("Hex", 10), ("HexNAc", 8), ("Fuc", 2)],
               "G1F/G2F SA": [("Hex", 9), ("HexNAc", 8), ("Fuc", 2), ("Neu5Ac", 1)],
               "G1F/G2F (SA)2": [("Hex", 9), ("HexNAc", 8), ("Fuc", 2), ("Neu5Ac", 2)],
               "G2F/G2F SA": [("Hex", 10), ("HexNAc", 8), ("Fuc", 2), ("Neu5Ac", 1)],
               "G2F/G2F (SA)2": [("Hex", 10), ("HexNAc", 8), ("Fuc", 2), ("Neu5Ac", 2)]}


def glycanlist_to_modlist(glycans, max_counts=2, use_average_masses=True):
    """
    Convert a list of glycans to a list of modifications.

    :param glycans: Dict {name: composition} of glycans; format like fc_glycans
    :param max_counts: maximum count for the modification search
    :param use_average_masses: if true, use average masses for calculation, otherwise monoisotopic masses
    :return: A list of triples: (1) Name of the glycan, (2) its mass, (3) max count
    """
    result = []
    for name, composition in glycans.items():
        # mg[0] is the monosaccharide name, mg[1] the number of monosaccharides
        formula = mass_tools.combine_formulas([glycan_formula[mg[0]] * mg[1] for mg in composition])
        if use_average_masses:
            result.append((name, formula.average_mass, max_counts))
        else:
            result.append((name, formula.monoisotopic_mass, max_counts))
    return result


# def glycopattern_formula(pattern, name=None):
#     mg,c = pattern[0]
#     gp_atoms = glycan_formula[mg]*c
#     for mg,c in pattern[1:]:
#         gp_atoms += glycan_formula[mg]*c
#     if name:
#         gp_atoms.name = name
#     return mass_tools.Formula(gp_atoms)
#
#
# def glycan_composition_formula(composition):
#     gformula = []
#     for g,c in composition:
#         gformula.append(glycan_formula[g]*c)
#     gformula = mass_tools.combine_formulas(gformula)
#     return gformula
#
#
# def glycan_string_formula(gstring):
#     #Byonic Naming conventions catch.
#     gstring = gstring.replace("NeuAc", "Neu5Ac")
#     gstring = gstring.replace("NeuGc", "Neu5Gc")
#     pattern = re.compile(r"(\w+)")
#     hits = pattern.findall(gstring)
#     composition = zip(hits[::2],[int(i) for i in hits[1::2]])
#     return glycan_composition_formula(composition)
#
#
# def read_byonic_library(filename):
#     fh = open(filename)
#     glycans = {}
#     for g in [i.split(" % ")[0] for i in fh.readlines()[1:]]:
#         glycans[g] = glycan_string_formula(g)
#     return glycans
#
#
# def complete_mods(glycan_library, positions=2, min_mass=0):
#     mods = []
#     for g in glycan_library:
#         m = glycan_library[g]._average_mass
#         if m >= min_mass:
#             mods.append((g,m,positions))
#     return mods
#
#
class Nglycan(object):

    def __init__(self, ginput, name=None):
        try:
            if type(ginput) == dict:
                self._gdict = ginput
            elif type(ginput) == str:
                pattern = re.compile(r"([NHFPSG]\d+)")
                self._gdict = {single_letter_mono[g[0]]:int(g[1:]) for g in pattern.findall(ginput)}
        except:
            print("Unable to recognize input")
        # self._gtuple = tuple(sorted(self._gdict.items()))
        # self.calculate()
        if name:
            self._tname = name
        else:
            self._tname = None
        self.calculate()

    def remove_tname(self):
        self._tname = None

    def make_name(self):
        # self._name = "".join(["{}({})".format(*i) for i in self._gtuple])
        if self._gdict:
            self._name = "".join(["{}{}".format(*(mono_single_letter[i[0]], i[1])) for i in self._gtuple])
        else:
            self._name = "Unglycosylated"

    def calculate(self):
        self._gdict = {k:v for k,v in self._gdict.items() if v}
        self._form_dict = {}
        for i in self._gdict:
            self._form_dict[(i, self._gdict[i])]= glycan_formula[i]*self._gdict[i]
        self._formula = mass_tools.Formula("")
        for i in self._form_dict.values():
            self._formula += i
        self._average_mass = self._formula.average_mass
        self._monoisotopic_mass = self._formula.monoisotopic_mass
        self._gtuple = tuple(sorted(self._gdict.items()))
        self.make_name()

    def sialylated(self):
        self.sialidase()
        return self

    def sialidase(self):
        try:
            del self._gdict["Neu5Ac"]
            self.calculate()
        except:
            pass

    def remove_core(self):
        try:
            self._gdict["Hex"] = self._gdict["Hex"]-3
            self._gdict["HexNAc"] = self._gdict["HexNAc"]-2
            self.calculate()
        except:
            print("No core")

    def endoS(self):
        try:
            if self._gdict["Fuc"] != 0:
                self.fucosidase()
            self._gdict["HexNAc"] = self._gdict["HexNAc"]-1
            self.calculate()
        except:
            pass

    def fucosidase(self):
        try:
            del self._gdict["Fuc"]
            self.calculate()
        except:
            pass

    def strip(self):
        self.fucosidase()
        self.sialidase()
        self.remove_core()

    def __hash__(self):
        return(hash(self._gtuple))

    def __eq__(self, other):
        return(self._gtuple == other._gtuple)

    def __ne__(self, other):
        return(not(self == other))

    def __lt__(self, other):
        if self._average_mass < other._average_mass:
            return -1
        elif self._average_mass > other._average_mass:
            return 1
        else:
            return 0

    def __repr__(self):
        # if hasattr(self, "_name"):
        restr = "{:s}".format(self._name)
        # else:
            # restr = ""
        # restr += "".join(["{}({})".format(*i)for i in self._gtuple])
        restr += "|m{:.2f}".format(self._monoisotopic_mass)
        if self._tname:
            restr += "|{:s}".format(self._tname)
        return(restr)
#
#
# class NglycanLibrary(object):
#
#     def __init__(self, site_specific_dict, site_specifi_abundance=False, dimeric=False):
#         self._ssd = {j:[Nglycan(i[1], name=i[0]) for i in site_specific_dict[j]]
#                                       for j in site_specific_dict}
#         self._all_nglycans = list(set([i for j in self._ssd
#                                          for i in self._ssd[j]]))
#         self._all_nglycans.sort(key=lambda t:t._average_mass, reverse=True)
#         self._library = pd.DataFrame({"Nglycan":[i for i in self._all_nglycans]})
#         self._sites = []
#         for site in self._ssd:
#             if dimeric:
#                 self._df = 2
#             else:
#                 self._df = 1
#             self._library[site] = [1*self._df if i in self._ssd[site] else 0
#                                    for i in self._all_nglycans]
#             self._sites.append(site)
#         self._sites.sort()
#         self._library["Nr."] = self._library[self._sites].T.sum()
#         self._library["Name"] = self._library.Nglycan.apply(lambda x:x._tname)
#         self._max_ngs = len(self._ssd)*self._df
#         self._cols = ["Name", "Nglycan"]
#         self._cols.extend(self._sites)
#         self._cols.append("Nr.")
#         self._library = self._library[self._cols]
#
#     def sialidase(self):
#         self._library.Nglycan.apply(lambda x: x.sialidase())
#         self._library.Name = self._library.Nglycan.apply(lambda x:x._name)
#         lgrp = self._library.groupby("Name")
#         n = []
#         for k,v in lgrp:
#             if len(v) > 1:
#                 v[self._sites]
#                 sitesum = v[self._sites].sum()
#                 if any(sitesum[sitesum>self._df]):
#                     for site, count in sitesum.iteritems():
#                         if count > self._df:
#                             sitesum[site] = self._df
#                 ns = v.iloc[0][["Name", "Nglycan"]]
#                 ns = ns.append(sitesum)
#                 ns["Nr."] = ns[self._sites].sum()
#                 n.append(pd.DataFrame(ns).T)
#             else:
#                 n.append(pd.DataFrame(v))
#         n = pd.concat(n)
#         n.index = range(len(n))
#         n["Nglycan"].apply(lambda x:x.remove_tname())
#         self._library = n
#
#     def extract_sites(self, nsites):
#         dropsites = set(self._sites)-set(nsites)
#         nlib = copy.deepcopy(self)
#         for site in dropsites:
#             nlib._library.drop(site, axis=1, inplace=True)
#         nlib._library["Nr."] = nlib._library[nsites].T.sum()
#         nlib._library = nlib._library[nlib._library["Nr."]!=0]
#         nlib._library["Nr."] = nlib._library["Nr."].astype(int)
#         nlib._sites = nsites
#         nlib._cols = ["Name", "Nglycan"]
#         nlib._cols.extend(nlib._sites)
#         nlib._cols.append("Nr.")
#         nlib._max_ngs = len(nsites)*nlib._df
#         return nlib
#
#     def __add__(self, other):
#         nlib = copy.deepcopy(self)
#         new = pd.concat([self._library, other._library])
#         ngrp = new.groupby("Name")
#         n = []
#         for k,v in ngrp:
#             if len(v) > 1:
#                 v[self._sites]
#                 sitesum = v[self._sites].sum()
#                 if any(sitesum[sitesum>self._df]):
#                     for site, count in sitesum.iteritems():
#                         if count > self._df:
#                             sitesum[site] = self._df
#                 ns = v.iloc[0][["Name", "Nglycan"]]
#                 ns = ns.append(sitesum)
#                 ns["Nr."] = ns[self._sites].sum()
#                 n.append(pd.DataFrame(ns).T)
#             else:
#                 n.append(pd.DataFrame(v))
#         n = pd.concat(n)
#         n.index = range(len(n))
#         nlib._library = n
#         return nlib
#
#     def to_excel(self):
#         pass
#
#     def to_modlist(self):
#         modlist = []
#         for i,row in self._library.iterrows():
#             modlist.append((row["Name"], row["Nglycan"]._average_mass, row["Nr."]))
#         return modlist
#
#     def __repr__(self):
#         return(self._library.to_string())
#
#
# def glycan_composition():
#     pass
#
#
# def make_byonic_glycan(byonic_string):
#     composition_list = [i.split("(") for i in byonic_string.split(")")[:-1]]
#     gdict = {}
#     for msac, count in composition_list:
#         if msac =="NeuAc":
#             gdict["Neu5Ac"] = int(count)
#         elif msac == "NeuGc":
#             gdict["Neu5Gc"] = int(count)
#         else:
#             gdict[msac] = int(count)
#     return Nglycan(gdict)
#
#
# def make_nglycan(bp_series):
#     #building glycan from information in BioPharma library, starting with core
#     glycdict  = {"Hex":3, "HexNAc":2}
#     #add fucose
#     if bp_series["Fuc"] != 0:
#         glycdict["Fuc"] = 1
#     #going through antennas
#     antennae = {"Hex":0, "HexNAc":0, "Neu5Ac":0, "Neu5Gc":0}
#     for i in range(1,5):
#         a = bp_series["Antenna %d" %i]
#         if not a == "None":
#             asp = a.split("-")[:-1]
#             for ms in asp:
#                 if ms == "S":
#                     antennae["Neu5Ac"]+=1
#                 if ms == "Sg":
#                     antennae["Neu5Gc"]+=1
#                 if ms == "Gn":
#                     antennae["HexNAc"]+=1
#                 if ms == "G":
#                     antennae["Hex"]+=1
#                 if ms == "Ga":
#                     antennae["Hex"]+=1
#     #adding mannoses
#     glycdict["Hex"] += int(bp_series["Man"]-3)
#     glycdict = { k: glycdict.get(k, 0) + antennae.get(k, 0) for k in set(glycdict) | set(antennae) }
#     glycdict = {k:v for k,v in glycdict.items() if v}
#     return Nglycan(glycdict, name=bp_series["Glycan Name"])
#
#
# def read_nglycan_library(filename, sheetname="CHO glycans"):
#     nglycans_cho = pd.read_excel(filename, sheetname=sheetname)
#     nglycans_cho = nglycans_cho[:-7] #cutting of Gn, GnF since the don"t parse
#     nglycans = []
#     for i in nglycans_cho.index:
#         nglycans.append(make_nglycan(nglycans_cho.ix[i]))
#     return nglycans
#
#
# def pickle_bpf_glycans(filename, sheetname="CHO glycans"):
#     nglycans_raw = pd.read_excel(filename, sheetname=sheetname)
#     nglycans_raw = nglycans_raw[:-7] #cutting of Gn, GnF since the don"t parse
#     nglycans = {}
#     for i, line in nglycans_raw.iterrows():
#         nglycan = make_nglycan(line)
#         nglycans[nglycan._tname] = nglycan._name
#     dumpname = "".join([os.path.splitext(filename)[0], ".dump"])
#     fh = open(dumpname, "wb")
#     pickle.dump(nglycans, fh)
#     fh.close()
#     print("%s pickled to %s"% (filename, dumpname))
#
#
# def bpf_lookup(bpf_glycan, filename):
#     try:
#         dumpfh = open(filename, "rb")
#         ngdict = pickle.load(dumpfh)
#         dumpfh.close()
#         return Nglycan(ngdict[bpf_glycan], bpf_glycan)
#     except:
#         return Nglycan("H0",)


def dataframe_to_modlist(df):
    """
    Function turns site-specific N-glycan DataFrame into a modlist.
    DataFrame ought to look like so and Name has to be a composition string.


    :param df: Dataframe with the following format:
                  Name	  N149  N171  Nr.
               0  F1H3N4     0     2    2
               1  F1H4N4     0     2    2

    :returns: A list of tuples. Each tuple describes a modification:
              (1) name, (2) mass, (3) maximum occurrences
    """

    modlist = []
    for i, r in df.iterrows():
        mass = Nglycan(r["Name"])._average_mass
        modlist.append((r["Name"], mass, r["Nr."]))
    return modlist
