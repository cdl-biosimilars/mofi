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
    :return: None if the search completely failed, i.e., no single combination was found; otherwise:
             a dataframe with index (1) Mass_ID (corresponds to index of experimental mass in input list of peaks),
                                    (2) Isobar (0-based consecutive numbering of found masses) and
                                    (3) Stage1_hit (0-based consecutive numbering of hits per Mass_ID
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
    any_combination_found = False  # is set to true as soon as a combination is found

    # run the search on each peak
    for mass_index, unexplained_mass in enumerate(unexplained_masses):
        try:
            current_tolerance = mass_tolerance[mass_index]
        except TypeError:
            current_tolerance = mass_tolerance

        if progress_bar is not None:
            progress_bar.setValue(int((mass_index + 1) / len(unexplained_masses) * 100))
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
            any_combination_found = True
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
        else:  # no appropriate combination was found
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

        # only accept solutions with at most 2 Neu5Ac per O-core  TODO remove this filter? -> Therese
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
    if not any_combination_found:
        return

    # sort columns so that they have the original order from the monomer table
    sortlist = [m[0] for m in mods] + ["Exp. Mass", "Theo. Mass", "Da.", "ppm"]
    combinations = combinations.reindex_axis(sortlist, axis="columns")
    combinations.index.names = ["Mass_ID", "Isobar", "Stage1_hit"]

    # amend the index
    # (a) index "Isobar", which is currently a list of (float) theoretical masses,
    #                     but should be a consecutive (integer) numbering
    isodict = {v: i for i, v in enumerate(combinations["Theo. Mass"].unique())}  # a {mass: running counter} dict
    combinations.reset_index("Isobar", drop=True, inplace=True)  # delete index "Isobar"
    combinations["Isobar"] = combinations["Theo. Mass"].map(isodict)
    combinations.set_index("Isobar", append=True, inplace=True)

    # (b) index "Stage1_hit", which is wrong at this stage, since theo. masses from different searches (peaks)
    #                  may have been grouped under a single theo. mass index
    combinations.reset_index("Stage1_hit", drop=True, inplace=True)  # delete index "Stage1_hit"
    combinations["Stage1_hit"] = combinations \
        .groupby(level="Isobar")\
        .cumcount()  # create column "Stage1_hit", a counter per isobar
    combinations.set_index("Stage1_hit", inplace=True, append=True)  # move column "Stage1_hit" to index
    combinations.reorder_levels(["Mass_ID", "Isobar", "Stage1_hit"])  # organize indices

    # final structure of combinations:
    # - multiindex (1) Mass_ID: consecutive numbers for peaks (=experimental masses)
    #              (2) Isobar: consecutive numbers for combinations (=theoretical masses)
    #              (3) Stage1_hit: consecutive numbers for all combinations within an isobar
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


def _calc_glycan_composition_and_abundance(row, glycan_composition=None, monomers=None):
    """
    Sum the monomer counts of several glycans and calculate their abundance.

    :param row: row from a dataframe containing the name of a glycan per cell
    :param glycan_composition: dataframe indication the glycan composition
    :param monomers: list of monomers in the library
    :return: a tuple containing the numbers of each atom
    """
    composition = np.zeros(len(monomers), dtype=np.uint16)
    abundance = 1
    for site, glycan in row.iteritems():
        composition += glycan_composition["Monomers"][glycan, site]
        abundance *= glycan_composition["Abundance"][glycan, site]
    abundance /= 100 ** (len(row) - 1)
    return tuple(composition), abundance


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
    :return: a datframe that
    """

    if progress_bar is not None:
        progress_bar.setValue(0)

    # mods_per_site: series, like
    # Sites
    # N149_A                         [A2G2, A2G2F, M5]
    # N149_B                         [A2G2, A2G2F, M5]
    # N171_A    [A2G2F, A2G0F, A2G1F, A3G3F, A2Ga1G2F]
    # N171_B    [A2G2F, A2G0F, A2G1F, A3G3F, A2Ga1G2F]
    rows = []
    for name, data in glycan_library.iterrows():
        for site in data["Sites"].split(","):
            rows.append([site.strip(), name])
    mods_per_site = pd.DataFrame(rows, columns=["Sites", "Name"]) \
        .groupby("Sites")["Name"] \
        .apply(list)

    if progress_bar is not None:
        progress_bar.setValue(20)

    # df_glycan_composition: dataframe, like
    #                                       Composition  Abundance      Monomers
    # Name     Sites
    # A2G0F    N171_A         1 N-core, 2 HexNAc, 1 Fuc      10.00  [0, 2, 1, 1]
    #          N171_B         1 N-core, 2 HexNAc, 1 Fuc      10.00  [0, 2, 1, 1]
    # A2G1F    N171_A  1 N-core, 2 HexNAc, 1 Fuc, 1 Hex       6.20  [1, 2, 1, 1]
    #          N171_B  1 N-core, 2 HexNAc, 1 Fuc, 1 Hex       6.20  [1, 2, 1, 1]
    df_glycan_composition = glycan_library.set_index(["Composition", "Abundance"], append=True)
    df_glycan_composition["Sites"] = df_glycan_composition["Sites"] \
        .apply(lambda x: [x.strip() for x in x.split(",")])  # convert sites to list
    df_glycan_composition = pd.melt(
            df_glycan_composition["Sites"].apply(pd.Series)
                                          .reset_index(),
            id_vars=["Name", "Composition", "Abundance"],
            value_name="Site") \
        .dropna() \
        .set_index(["Name", "Site"]) \
        .drop("variable", axis=1) \
        .sort_index()  # explode sites and move to index
    df_glycan_composition["Monomers"] = df_glycan_composition["Composition"] \
        .apply(_calc_monomer_counts, monomers=monomers)

    df_glycan_combinations = pd.DataFrame(list(product(*list(mods_per_site))),
                                          columns=list(mods_per_site.index))
    df_glycan_combinations["Monomers"], df_glycan_combinations["Abundance"] = zip(
        *df_glycan_combinations.apply(_calc_glycan_composition_and_abundance,
                                      axis=1,
                                      glycan_composition=df_glycan_composition,
                                      monomers=monomers))
    df_glycan_combinations = df_glycan_combinations \
        .set_index(pd.MultiIndex.from_tuples(df_glycan_combinations["Monomers"], names=monomers)) \
        .sort_index() \
        .drop("Monomers", axis=1)

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

    # create an additional index "Stage2_hit", which counts different possible combinations per monomer hit
    df_found_polymers["Stage2_hit"] = df_found_polymers.groupby(["Isobar", "Stage1_hit"]).cumcount()
    df_found_polymers.set_index("Stage2_hit", inplace=True, append=True)
    df_found_polymers.reorder_levels(["Mass_ID", "Isobar", "Stage1_hit", "Stage2_hit"])

    if progress_bar is not None:
        progress_bar.setValue(100)

    # df_found_polymers.to_excel("df_found_polymers.xls")
    return df_found_polymers
