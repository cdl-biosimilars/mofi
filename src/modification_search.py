"""
    modification_search.py

    Helper functions for searching modifications.

    Authors: Stefan Senn, Wolfgang Skala, Gabriel Huber

    (c) 2017 Christian Doppler Laboratory for Biosimilar Characterization
"""

import numpy as np
import pandas as pd
from findmods_source import findmods
import re
from itertools import product


def find_monomers(mods, unexplained_masses, mass_tolerance=5.0, explained_mass=0, progress_bar=None):
    """
    Wrapper function that runs the C function findmods.examine_modifications on a list of target masses.

    :param mods: list of tuples (name, mass, maxcount)
    :param unexplained_masses: iterable of unexplained masses (floats)
    :param mass_tolerance: tolerance for the unexplained mass in Da; either a single value, which applies to all
                           unexplained_masses, or a list of tolerances (one value per unexplained mass)
    :param explained_mass: mass explained by the protein sequence and known modifications
    :param progress_bar: a QProgressBar that gets updated during the search
    :return: a dataframe with index (1) Massindex (corresponds to index of experimental mass in input list of peaks),
                                    (2) Isobar (0-based consecutive numbering of found masses) and
                                    (3) Hit (0-based consecutive numbering of hits per Massindex
             columns: one column for each modification
                      Exp. Mass
                      Theo. Mass
                      Da.
                      ppm
                      Modstring
    """

    # max-first sorting decreases the time of the algorithm by a factor of approx 1.5
    sorted_mods = sorted(mods, key=lambda t: t[1], reverse=True)
    mod_names = [m[0] for m in sorted_mods]
    mod_masses = np.array([m[1] for m in sorted_mods])
    sorted_mods = [(m[1], m[2]) for m in sorted_mods]
    combinations = {}

    # run the search on each peak
    for mass_index, unexplained_mass in enumerate(unexplained_masses):
        try:
            current_tolerance = mass_tolerance[mass_index]
        except TypeError:
            current_tolerance = mass_tolerance

        if progress_bar is not None:
            progress_bar.setValue(int(mass_index / (len(unexplained_masses) - 1) * 100))
        result = findmods.examine_modifications(sorted_mods, unexplained_mass, current_tolerance)

        # (a) transform result to combs_per_mass.
        # combs_per_mass will be a list of dicts with the following keys:
        # [mod_name]: number of occurrences ... for each type of modification used in the search
        # Exp. Mass: mass measured in the experiment
        # Theo. Mass: mass of protein and found modifications
        # Da.: remaining unexplained mass difference
        # abs(Delta): absolute remaining mass difference
        # ppm: mass difference in ppm
        combs_per_mass = []
        if result:
            for r in result:
                theoretical_mass = explained_mass + sum(np.array(r) * mod_masses)  # mass of protein + found mods
                experimental_mass = explained_mass + unexplained_mass  # mass measured by the mass spectrometer
                r = dict(zip(mod_names, r))
                r["Theo. Mass"] = theoretical_mass
                r["Exp. Mass"] = experimental_mass
                r["Da."] = experimental_mass - theoretical_mass
                r["abs(Delta)"] = abs(r["Da."])
                r["ppm"] = (experimental_mass - theoretical_mass) / theoretical_mass * 1_000_000
                combs_per_mass.append(r)
        else:  # what to do if examine_modifications didn"t find any modifications
            r = np.zeros(len(sorted_mods), dtype="int")
            theoretical_mass = 0
            experimental_mass = explained_mass + unexplained_mass
            r = dict(zip(mod_names, r))
            r["Theo. Mass"] = theoretical_mass
            r["Exp. Mass"] = experimental_mass
            r["Da."] = 0.0
            r["abs(Delta)"] = 0.0
            r["ppm"] = 0.0
            combs_per_mass.append(r)

        # (b) create a dataframe from combs_per_mass
        combs_frame = pd.DataFrame(combs_per_mass)

        # only accept solutions with at most 2 Neu5Ac per O-core  TODO remove?
        if "O-core" in mod_names and "Neu5Ac" in mod_names:
            combs_frame = combs_frame[2 * combs_frame["O-core"] >= combs_frame["Neu5Ac"]]

        # sort by abs(Delta) and reindex the frame starting from 1 instead of 0
        # thereby, alternative combninations will be numbered 1, 2, ...
        combs_frame = combs_frame.sort_values("abs(Delta)")
        combs_frame.index = range(1, len(combs_frame) + 1)

        # (c) finally, append the solutions for the current peak to the overall list of combinations
        combinations[mass_index] = combs_frame.set_index("Theo. Mass", append=True, drop=False).swaplevel(1, 0)

    # after processing all peaks, combinations is a dict of dataframes
    # each dataframe contains information on all found combinations for a single peak
    # now, combine all those dataframes into a single dataframe
    # there will be three indices: (1) consecutive numbers
    #                              (2) theoretical mass
    #                              (3) 1-based counter for hits per experimental mass
    combinations: pd.DataFrame = pd.concat(combinations)

    # sort columns so that they have the original order from the monomer table
    sortlist = [m[0] for m in mods] + ["Exp. Mass", "Theo. Mass", "Da.", "ppm"]
    combinations = combinations.reindex_axis(sortlist, axis="columns")

    # amend the index
    # (a) index "Isobar", which is currently a list of (float) theoretical masses,
    #                     but should be a consecutive (integer) numbering
    combinations.index.names = ["Massindex", "Isobar", "Hit"]
    isodict = {v: i for i, v in enumerate(combinations["Theo. Mass"].unique())}  # a {mass: running counter} dict
    combinations.reset_index("Isobar", drop=True, inplace=True)  # delete index "Isobar"
    combinations["Isobar"] = combinations["Theo. Mass"].map(isodict)
    combinations.set_index("Isobar", append=True, inplace=True)

    # (b) index "Hit", which is wrong at this stage, since theo. masses from different searches (peaks)
    #                  may have been grouped under a single theo. mass index
    combinations.reset_index("Hit", drop=True, inplace=True)  # delete index "Hit"
    combinations["Hit"] = combinations.groupby(level="Isobar").cumcount()  # create column "Hit", a counter per isobar
    combinations.set_index("Hit", inplace=True, append=True)  # move column "Hit" to index
    combinations.reorder_levels(["Massindex", "Isobar", "Hit"])  # organize indices

    # final structure of combinations:
    # - multiindex (1) Massindex: consecutive numbers for peaks (=experimental masses)
    #              (2) Isobar: consecutive numbers for combinations (=theoretical masses)
    #              (3) Hit: consecutive numbers for all combinations within an isobar
    # - one column for each modification (Hex, HexNAc, N-core, ...): number of each modification
    # - Exp. Mass: peak mass
    # - Theo. Mass: mass of found combination
    # - Da.: Theo. Mass - Exp. Mass
    # - ppm: relative difference wrt Theo. Mass
    # - Modstring: string representation of modifications
    return combinations


_re_monomer_list = re.compile(r"(\d*)\s*([\w-]+)(?:,|$)")  # regex for decomposing glycan compositions

def _calc_monomer_counts(value, monomers=None):
    """
    Calculate monomer counts for a given complex glycan

    :param value: string indicating the composition of the glycan (like "1 N-core, 2 Hex")
    :param monomers: list of monomers in the library
    :return: numpy array containing the monomer counts (including absent monomers)
    """
    composition = {}  # generate a {name: count} dict
    for (count, monomer) in _re_monomer_list.findall(value):
        try:
            composition[monomer] = int(count)
        except ValueError:
            composition[monomer] = 1

    result = np.zeros(len(monomers), dtype=np.uint16)
    for i, m in enumerate(monomers):
        result[i] = composition.setdefault(m, 0)  # convert to a numpy array
    return result


def _sum_glycan_monomer_counts(row, glycan_composition=None, monomers=None):
    """
    Sum the monomer counts of several glycans.

    :param row: row from a dataframe containing the name of a glycan per cell
    :param glycan_composition: dataframe indication the glycan composition
    :param monomers: list of monomers in the library
    :return: a tuple containing the numbers of each atom
    """
    result = np.zeros(len(monomers), dtype=np.uint16)
    for glycan in row:
        result += glycan_composition["Monomers"][glycan]
    return tuple(result)


def get_monomers_from_library(glycan_library):
    """
    Find all monomers that appear in the glycan library

    :param glycan_library: dataframe with glycan library (see find_polymers)
    :return: alphabetically ordered list of monomers
    """
    rows = []
    for name, data in glycan_library.iterrows():
        for (count, monomer) in _re_monomer_list.findall(data["Composition"]):
            try:
                count = int(count)
            except ValueError:
                count = 1
            rows.append([name, monomer, count])
    df_monomers = pd.DataFrame(rows, columns=["Glycan", "Monomer", "Count"])
    return list(df_monomers["Monomer"].unique())


def find_polymers(combinations, glycan_library, monomers, progress_bar=None):
    """
    Search a polymer library with the results from a combinatorial monomer search.

    :param combinations: dataframe with results from combinatorial search
    :param glycan_library: dataframe containing a glycan library
                           index: glycan names
                           columns: Composition, Sites, Abundance
    :param monomers: list of monomers in the library as returned by get_monomers_from_library
    :param progress_bar: a QProgressBar that gets updated during the search
    :return: nothing
    """

    if progress_bar is not None:
        progress_bar.setValue(0)

    # mods_per_site: a Series with the glycosylation sites as index and lists of possible glycans as values
    rows = []
    for name, data in glycan_library.iterrows():
        for site in data["Sites"].split(","):
            rows.append([site.strip(), name])
    mods_per_site = pd.DataFrame(rows, columns=["Sites", "Name"]) \
        .groupby("Sites")["Name"] \
        .apply(list)
    if progress_bar is not None:
        progress_bar.setValue(20)

    # df_glycan_composition: dataframe with glycans as index and columns "Composition" and "Monomers" (list of counts)
    df_glycan_composition = pd.DataFrame(index=glycan_library.index)
    df_glycan_composition["Monomers"] = glycan_library["Composition"] \
        .apply(_calc_monomer_counts, monomers=monomers)
    if progress_bar is not None:
        progress_bar.setValue(40)

    # df_glycan_combinations: a dataframe with monomers (Hex, HexNAc, ...) as multiindex and one column per site,
    # indicating which glycan composition belongs to a given elemental composition
    df_glycan_combinations = pd.DataFrame(list(product(*list(mods_per_site))), columns=list(mods_per_site.index))
    df_glycan_combinations["Monomers"] = df_glycan_combinations \
        .apply(_sum_glycan_monomer_counts, axis=1, glycan_composition=df_glycan_composition, monomers=monomers)
    df_glycan_combinations = df_glycan_combinations \
        .set_index(pd.MultiIndex.from_tuples(df_glycan_combinations["Monomers"], names=monomers)) \
        .sort_index() \
        .drop("Monomers", 1)
    df_glycan_combinations["Abundance"] = df_glycan_combinations \
        .apply(lambda row:
               np.prod(glycan_library.loc[row]["Abundance"]) / (100 ** (len(df_glycan_combinations.columns) - 1)),
               axis=1)
    if progress_bar is not None:
        progress_bar.setValue(60)

    # df_found_polymers: combinations with monomer composition as multiindex,
    # with additional columns for the polymer sites
    old_index = combinations.index.names
    df_found_polymers = combinations \
        .reset_index(old_index) \
        .set_index(monomers) \
        .sort_index() \
        .join(df_glycan_combinations, how="inner")
    if progress_bar is not None:
        progress_bar.setValue(80)

    df_found_polymers = df_found_polymers \
        .reset_index(df_found_polymers.index.names) \
        .set_index(old_index) \
        .sort_index()
    if progress_bar is not None:
        progress_bar.setValue(100)

    return df_found_polymers








#
# The following functions are merely for performance testing
#
# def fast_find_modifications_python(mods, unexplained_masses, mass_tolerance=5.0):
#     """
#     test function to implement the search in pure python
#     same as C++ function
#     """
#
#     def next_mod(used_mass=0, index=0):
#         next_index = index + 1
#         for x in range(mods[index][2] + 1):
#             counts[index] = x
#             remaining = unexplained_mass - used_mass
#             if abs(remaining) <= mass_tolerance:
#                 result.append(counts[:])
#             elif next_index < len(counts) and remaining >= mass_tolerance:
#                 next_mod(used_mass, next_index)
#             used_mass += mods[index][1]
#         counts[index] = 0
#
#     mods.sort(key=lambda t: t[1], reverse=True)
#
#     # run the search on each peak
#     for mass_index, unexplained_mass in enumerate(unexplained_masses):
#         result = []
#         counts = [0 for _ in range(len(mods))]
#         next_mod()
#
#
# def fast_find_modifications_cpp(mods, unexplained_masses, mass_tolerance=5.0):
#     """
#     Wrapper function that runs the C function findmods.examine_modifications on a list of target masses.
#
#     """
#
#     # max-first sorting decreases the time of the algorithm by a factor of approx 1.5
#     mods.sort(key=lambda t: t[1], reverse=True)
#     mods = [(m[1], m[2]) for m in mods]
#
#     # run the search on each peak
#     for mass_index, unexplained_mass in enumerate(unexplained_masses):
#         result = findmods.examine_modifications(mods, unexplained_mass, mass_tolerance)
#
#
# def fast_find_modifications_db(mods, unexplained_masses, mass_tolerance=5.0):
#     """
#     test function: find all possible combinations and store in database
#     for each peak, find the possible database entries
#     """
#
#     def explore(used_mass=0, index=0):
#         next_index = index + 1
#         for x in range(mods[index][2] + 1):
#             counts[index] = x
#             if next_index < len(counts) and largest_mass - used_mass >= mass_tolerance:
#                 explore(used_mass, next_index)
#             else:
#                 database.append((used_mass, counts[:]))
#             used_mass += mods[index][1]
#         counts[index] = 0
#
#     mods.sort(key=lambda t: t[1], reverse=True)
#     counts = [0 for i in range(len(mods))]
#
#     # find all combinations that may produce any mass leq than the largest peak
#     database = []
#     largest_mass = unexplained_masses[-1]
#     explore()
#     database.sort()
#     database = list(zip(*database))
#     database = pd.Series(database[1], index=database[0])
#
#     # extract solutions for each peak from the database
#     for mass_index, unexplained_mass in enumerate(unexplained_masses):
#         # current_time = time.time()
#         # print(mass_index, "{:.4f}".format(current_time - last_time), sep="\t", end="\t")
#         # last_time = current_time
#         result = database[unexplained_mass - mass_tolerance:unexplained_mass + mass_tolerance]
#         # print(mass_index, len(result), list(result), sep="\t")
#
#
#
#
#
#
#
# def speed_test():
#     import time
#
#     # test data from Kadcyla.mofi
#     modlist_small = [
#         ("DM1", 957.00, 10),
#         ("G0F", 1445.34, 2),
#         ("G1F", 1607.48, 2),
#         ("G2FSA2", 2352.13, 2),
#         ("G0", 1299.19, 2),
#         ("G2F", 1769.62, 2),
#         ("Man5", 1217.09, 2),
#         ("G2", 1623.48, 2),
#         ("G2FSA1", 2060.87, 2),
#         ("G1", 1461.33, 2)]
#
#     modlist_large = [
#         ("Hex", 162.14, 10),
#         ("Fuc", 146.14, 10),
#         ("Neu5Ac", 291.26, 10),
#         ("Neu5Gc", 307.25, 10),
#         ("Pent", 132.11, 10),
#         ("N-core", 892.81, 10),
#         ("O-core", 365.33, 10),
#         ("HexNAc", 203.19, 10)]
#
#     unexpl_masses = [2891.59, 2956.68, 3054.58, 3203.98, 3215.83, 3232.54, 3303.90, 3343.96, 3376.25, 3403.45, 3521.36,
#                      3650.07, 3714.38, 3719.19, 3798.27, 3850.46, 3878.96, 4013.87, 4067.96, 4089.39, 4173.50, 4233.29,
#                      4287.07, 4330.80, 4470.43, 4551.04, 4588.41, 4615.67, 4660.94, 4716.18, 4808.26, 4877.73, 4880.69,
#                      4904.34, 4969.80, 5027.53, 5052.26, 5127.60, 5131.48, 5189.57, 5213.72, 5251.64, 5294.28, 5354.30,
#                      5420.16, 5507.31, 5519.00, 5561.94, 5577.87, 5618.95, 5663.94, 5766.53, 5800.32, 5836.69, 5868.01,
#                      5926.92, 5987.60, 5999.96, 6031.49, 6090.08, 6149.26, 6250.56, 6252.27, 6308.23, 6421.27, 6473.93,
#                      6531.91, 6574.23, 6665.33, 6724.26, 6747.70, 6818.92, 6884.06, 6941.44, 6944.23, 7047.95, 7105.60,
#                      7209.33, 7266.45, 7339.15, 7371.41, 7435.66, 7491.26, 7535.30, 7591.28, 7610.76, 7681.71, 7841.83,
#                      7902.60, 7915.18, 7933.71, 7961.13, 8004.98, 8063.93, 8112.17, 8163.52, 8223.90, 8226.63, 8257.46,
#                      8294.31, 8449.78, 8496.19, 8637.39, 8649.98, 8668.11, 8711.36, 8800.37, 8858.87, 8862.39, 8901.50,
#                      8959.10, 8962.09, 9020.06, 9057.90, 9122.87, 9177.56, 9183.76, 9186.83, 9284.89, 9340.71, 9372.35,
#                      9593.22, 9596.52, 9612.99, 9670.15, 9755.44, 9812.47, 9817.96, 9825.00, 9917.68, 9973.48, 9977.10,
#                      10078.55, 10138.86, 10148.74, 10209.43, 10263.88, 10308.55, 10312.59, 10423.32, 10454.34, 10529.87,
#                      10666.29, 10712.98, 10879.41, 10910.56, 10973.47, 11036.43, 11079.21, 11093.55, 11224.69, 11314.54,
#                      11375.35, 11562.29, 11720.59, 11739.16, 11793.62, 11802.98]
#
#     test_functions = [fast_find_modifications_cpp, fast_find_modifications_python, fast_find_modifications_db]
#     test_modlists = [modlist_small]
#
#     for f in test_functions:
#         for m in test_modlists:
#             print("Testing function", f, "on modlist", m)
#             for i in range(5):
#                 start_time = time.time()
#                 f(m, unexpl_masses)
#                 print("\t", "{:.4f}".format(time.time() - start_time))
#
#
# if __name__ == "__main__":
#     speed_test()
