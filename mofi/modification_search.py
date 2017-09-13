"""
Helper functions for searching modifications.
"""

import numpy as np
import pandas as pd
from mofi import findmods
import re
from itertools import product


def find_monomers(mods, unexplained_masses, mass_tolerance=5.0,
                  explained_mass=0, progress_bar=None):
    """
    Wrapper function that runs the C function
    :func:`findmods.examine_modifications()` on a list of target masses.

    :param mods: list of tuples (name, mass, maxcount)
    :param unexplained_masses: iterable of unexplained masses (floats)
    :param mass_tolerance: tolerance for the unexplained mass in Da;
                           either a single value, which applies to all
                           unexplained_masses, or a list of tolerances
                           (one value per unexplained mass)
    :param explained_mass: mass explained by the protein sequence
                           and known modifications
    :param QProgressBar progress_bar: a progress bar that gets updated
                                      during the search
    :return: None if the search completely failed, i.e., no single combination
             was found; otherwise: a dataframe with multiindex

             (1) Mass_ID (corresponds to index of experimental mass
                 in input list of peaks),
             (2) Isobar (0-based consecutive numbering of found masses) and
             (3) Stage1_hit (0-based numbering of hits per Mass_ID)
             (4) Stage2_hit (0-based numbering of hits per Stage1_hit)

             columns:

             * one column for each modification
             * Exp. Mass
             * Theo. Mass
             * Da
             * ppm
    """

    sorted_mods = sorted(mods, key=lambda t: t[1], reverse=True)
    mod_names = [m[0] for m in sorted_mods]
    mod_masses = np.array([m[1] for m in sorted_mods])
    sorted_mods = [(m[1], m[2]) for m in sorted_mods]
    combinations = {}
    any_combination_found = False  # true as soon as a combination is found

    print("mass index\tsearch space size")

    # run the search on each peak
    for mass_index, unexplained_mass in enumerate(unexplained_masses):
        try:
            current_tolerance = mass_tolerance[mass_index]
        except TypeError:
            current_tolerance = mass_tolerance

        if progress_bar is not None:
            progress_bar.setValue(int((mass_index + 1)
                                      / len(unexplained_masses) * 100))
        print(mass_index, end="\t", flush=True)
        result = findmods.examine_modifications(sorted_mods,
                                                unexplained_mass,
                                                current_tolerance)
        print()

        # (a) transform result to combs_per_mass.
        # combs_per_mass will be a list of dicts with the following keys:
        # [mod_name]: number of occurrences ... for each type of modification
        # Exp. Mass: mass measured in the experiment
        # Theo. Mass: mass of protein and found modifications
        # Da: remaining unexplained mass difference
        # abs(Delta): absolute remaining mass difference
        # ppm: mass difference in ppm
        combs_per_mass = []
        if result:
            any_combination_found = True
            for r in result:
                theoretical_mass = (explained_mass
                                    + sum(np.array(r) * mod_masses))
                experimental_mass = explained_mass + unexplained_mass
                r = dict(zip(mod_names, r))
                r["Theo. Mass"] = theoretical_mass
                r["Exp. Mass"] = experimental_mass
                r["Da"] = experimental_mass - theoretical_mass
                r["abs(Delta)"] = abs(r["Da"])
                r["ppm"] = ((experimental_mass - theoretical_mass)
                            / theoretical_mass * 1000000)
                combs_per_mass.append(r)
        else:  # no appropriate combination was found
            combs_per_mass.append({"Theo. Mass": 0.0,
                                   "abs(Delta)": 0.0})

        # (b) create a dataframe from combs_per_mass
        # sort by abs(Delta) and reindex the frame starting from 1 instead of 0
        # thereby, alternative combninations will be numbered 1, 2, ...
        combs_frame = pd.DataFrame(combs_per_mass).sort_values("abs(Delta)")
        combs_frame.index = range(1, len(combs_frame) + 1)

        # (c) finally, append the solutions for the current peak
        # to the overall list of combinations
        combinations[mass_index] = (
            combs_frame
            .set_index("Theo. Mass", append=True, drop=False)
            .swaplevel(1, 0))

    # after processing all peaks, combinations is a dict of dataframes
    # each dataframe contains information on all found combinations
    # for a single peak
    # now, combine all those dataframes into a single dataframe
    # there will be three indices: (1) consecutive numbers
    #                              (2) theoretical mass
    #                              (3) 1-based counter for hits per exp. mass
    # drop all rows with NaN values (i.e., corresponding to failed searches)
    combinations = (
        pd.concat(combinations)
        .dropna()
        .astype({m: int for m in mod_names}))
    if not any_combination_found:
        return

    # sort columns so that they have the original order from the monomer table
    sortlist = [m[0] for m in mods] + ["Exp. Mass", "Theo. Mass", "Da", "ppm"]
    combinations = combinations.reindex_axis(sortlist, axis="columns")
    combinations.index.names = ["Mass_ID", "Isobar", "Stage1_hit"]

    # amend the index
    # (a) index "Isobar", which is currently a list of (float) theo. masses,
    #     but should be a zero-based consecutive (integer) numbering
    # isodict is a {mass: running counter} dict
    isodict = {v: i for i, v in enumerate(combinations["Theo. Mass"].unique())}
    combinations.reset_index("Isobar", drop=True, inplace=True)
    combinations["Isobar"] = combinations["Theo. Mass"].map(isodict)
    combinations.set_index("Isobar", append=True, inplace=True)

    # (b) index "Stage1_hit", which is wrong at this stage, since theo. masses
    #     from different searches (peaks) may have been grouped
    #     under a single theo. mass index
    combinations.reset_index("Stage1_hit", drop=True, inplace=True)
    combinations["Stage1_hit"] = (
        combinations
        .groupby(level="Isobar")
        .cumcount())  # create column "Stage1_hit", a counter per isobar
    combinations.set_index("Stage1_hit", inplace=True, append=True)

    combinations.reorder_levels(["Mass_ID", "Isobar", "Stage1_hit"])

    # final structure of combinations:
    # - multiindex with three levels:
    #   (1) Mass_ID: consecutive numbers for peaks (=exp. masses)
    #   (2) Isobar: consecutive numbers for combinations (=theo. masses)
    #   (3) Stage1_hit: consecutive numbers for all
    #                  combinations within an isobar
    # - one column for each modification (Hex, HexNAc, N-core, ...):
    #   number of each modification
    # - Exp. Mass: peak mass
    # - Theo. Mass: mass of found combination
    # - Da: Theo. Mass - Exp. Mass
    # - ppm: relative difference wrt Theo. Mass
    return combinations


# regex for decomposing glycan compositions
_re_monomer_list = re.compile(r"(\d*)\s*([\w-]+)(?:,|$)")


def _calc_monomer_counts(value, monomers=None):
    """
    Calculate monomer counts for a given complex glycan

    :param value: string indicating the composition of the glycan
                  (like "1 N-core, 2 Hex")
    :param monomers: list of monomers in the library
    :return: numpy array containing the monomer counts
             (including absent monomers)
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


def _calc_glycan_composition_abundance(row, glycan_composition=None,
                                       monomers=None):
    """
    Sum the monomer counts of several glycans and calculate their abundance.

    :param row: row from a dataframe containing the name of a glycan per cell
    :param glycan_composition: dataframe indicating the glycan composition
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

    :param glycan_library: dataframe with glycan library
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


def calc_polymer_combinations(glycan_library, monomers, progress_bar=None):
    """
    Calculate all possible polymer combinations based on a library.

    :param glycan_library: dataframe containing a glycan library;
                           index: glycan names;
                           columns: Composition, Sites, Abundance
    :param monomers: list of monomers in the library
                     as returned by :func:`get_monomers_from_library()`
    :param QProgressBar progress_bar: a progress bar that gets updated
                                      during the search
    :return: a dataframe like::

        #                         N149_A N149_B    N171_A    N171_B  Abundance
        #   Hex HexNAc Fuc N-core
        #   4   4      2   4          M5     M5     A2G0F     A2G0F        0.0
        #       6      2   4        A2G2     M5     A2G0F     A2G0F        0.0
        #                  4          M5   A2G2     A2G0F     A2G0F        0.0
        #              3   4       A2G2F     M5     A2G0F     A2G0F        0.0
        #                  4          M5  A2G2F     A2G0F     A2G0F        0.0

    :raises ValueError: if the library contains duplicate glycan/site pairs
    """

    if not monomers:  # we need a nonempty list of monomers to continue
        return

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
    mods_per_site = (
        pd.DataFrame(rows, columns=["Sites", "Name"])
        .groupby("Sites")["Name"]
        .apply(list))

    if progress_bar is not None:
        progress_bar.setValue(50)

    # df_glycan_composition: dataframe, like
    #                                       Composition  Abundance     Monomers
    # Name     Sites
    # A2G0F    N171_A         1 N-core, 2 HexNAc, 1 Fuc      10.00 [0, 2, 1, 1]
    #          N171_B         1 N-core, 2 HexNAc, 1 Fuc      10.00 [0, 2, 1, 1]
    # A2G1F    N171_A  1 N-core, 2 HexNAc, 1 Fuc, 1 Hex       6.20 [1, 2, 1, 1]
    #          N171_B  1 N-core, 2 HexNAc, 1 Fuc, 1 Hex       6.20 [1, 2, 1, 1]
    df_glycan_composition = (
        glycan_library
        .set_index(["Composition", "Abundance"], append=True))
    df_glycan_composition["Sites"] = (
        df_glycan_composition["Sites"]
        .apply(lambda x: [x.strip() for x in x.split(",")]))  # sites to list
    df_glycan_composition = (
        pd.melt(
            df_glycan_composition["Sites"].apply(pd.Series)
                                          .reset_index(),
            id_vars=["Name", "Composition", "Abundance"],
            value_name="Site")
        .dropna()
        .set_index(["Name", "Site"])
        .drop("variable", axis=1)
        .sort_index())  # explode sites and move to index
    df_glycan_composition["Monomers"] = (
        df_glycan_composition["Composition"]
        .apply(_calc_monomer_counts, monomers=monomers))

    # check for duplicate glycan/site combinations in the list of glycans
    df_duplicates = df_glycan_composition[
        (df_glycan_composition
         .groupby(level=["Name", "Site"])
         .size() > 1)
        .reindex(df_glycan_composition.index)]
    if not df_duplicates.empty:
        df_duplicates = (df_duplicates[~df_duplicates.index.duplicated()]
                         .reset_index())
        duplicates = ["{} at site {}".format(row["Name"], row["Site"])
                      for _, row in df_duplicates.iterrows()]
        raise ValueError("Duplicate glycan/site pairs detected: "
                         + ", ".join(duplicates))

    # df_glycan_combinations: see above
    df_glycan_combinations = pd.DataFrame(list(product(*list(mods_per_site))),
                                          columns=list(mods_per_site.index))
    df_glycan_combinations["Monomers"], df_glycan_combinations["Abundance"] = (
        zip(*df_glycan_combinations.apply(
            _calc_glycan_composition_abundance,
            axis=1,
            glycan_composition=df_glycan_composition,
            monomers=monomers)))
    df_glycan_combinations = (
        df_glycan_combinations
        .set_index(
            pd.MultiIndex.from_tuples(df_glycan_combinations["Monomers"],
                                      names=monomers))
        .sort_index()
        .drop("Monomers", axis=1))

    if progress_bar is not None:
        progress_bar.setValue(100)

    return df_glycan_combinations


def find_polymers(stage_1_results, polymer_combinations,
                  monomers, progress_bar=None):
    """
    Search a polymer library with the results from search stage 1.

    :param stage_1_results: dataframe with results from search stage 1
    :param polymer_combinations: dataframe with all possible
                                 polymer combinations as returned by
                                 :func:`calc_polymer_combinations`
    :param monomers: list of monomers in the library
                     as returned by get_monomers_from_library
    :param QProgressBar progress_bar: a progress bar that gets updated
                                      during the search
    :return pd.DataFrame: a dataframe like::

        #                                       Hex  HexNAc  Neu5Ac  Fuc  ...
        # Mass_ID Isobar Stage1_hit Stage2_hit
        # 0       6      2          0             0       4       0    2  ...
        # 2       137    1          0             1       4       0    2  ...
        #                           1             1       4       0    2  ...
        # 4       288    0          0             2       4       0    2  ...
        #                           1             2       4       0    2  ...
        #
        #
        #     Exp. Mass    Theo. Mass        Da        ppm    ch_A    ch_B
        #
        # 148057.122228  148056.20272  0.919508   6.210536     G0F     G0F
        # 148220.112210  148218.34356  1.768650  11.932736     G0F     G1F
        # 148220.112210  148218.34356  1.768650  11.932736     G1F     G0F
        # 148381.360467  148380.48440  0.876067   5.904193     G0F     G2F
        # 148381.360467  148380.48440  0.876067   5.904193     G1F     G1F
        #
        #
        # Abundance                  Hash  Permutations  Permutation abundance
        #       0.0    329177915115358981             1                    nan
        #       0.0  -4442633047439962303             2                    nan
        #       0.0  -4442633047439962303             2                    nan
        #       0.0    335986899926952527             1                    nan
        #       0.0  -4554265004199314174             1                    nan
    """

    if progress_bar is not None:
        progress_bar.setValue(0)

    old_index = stage_1_results.index.names
    try:
        df_found_polymers = (
            stage_1_results
            .reset_index(old_index)
            .set_index(monomers)
            .sort_index()
            .join(polymer_combinations, how="inner"))
    except TypeError:
        return None

    if progress_bar is not None:
        progress_bar.setValue(33)

    df_found_polymers = (
        df_found_polymers
        .reset_index(df_found_polymers.index.names)
        .set_index(old_index)
        .sort_index())

    # create an additional index "Stage2_hit", which counts
    # different possible combinations per monomer hit,
    # and rename column "Abundance" to "Score"
    df_found_polymers["Stage2_hit"] = (
        df_found_polymers
        .groupby(["Isobar", "Stage1_hit"])
        .cumcount())
    df_found_polymers = (
        df_found_polymers
        .set_index("Stage2_hit", append=True)
        .reorder_levels(["Mass_ID", "Isobar",
                         "Stage1_hit", "Stage2_hit"])
        .rename(columns={"Abundance": "Score"}))

    # recalculate the scores to represent the contribution
    # of each annotation to the peak height
    df_found_polymers["Score"] = (
        df_found_polymers["Score"]
        * 100
        / df_found_polymers.groupby("Mass_ID")["Score"].sum())

    if progress_bar is not None:
        progress_bar.setValue(67)

    # calculate glycan hash and counts/abundance per hash/peak/isobar
    df_found_polymers["Hash"] = (
        df_found_polymers
        .iloc[:, df_found_polymers.columns.get_loc("ppm") + 1: -1]
        .apply(lambda x: hash(frozenset(x)), axis=1))

    df_found_polymers = (
        df_found_polymers
        .reset_index()
        .set_index(["Mass_ID", "Isobar", "Hash"]))

    hash_group = df_found_polymers.groupby(df_found_polymers.index)
    df_found_polymers["Permutations"] = hash_group.size()
    df_found_polymers["Permutation score"] = hash_group["Score"].sum()

    df_found_polymers = (
        df_found_polymers
        .reset_index()
        .set_index(["Mass_ID", "Isobar", "Stage1_hit", "Stage2_hit"]))

    if progress_bar is not None:
        progress_bar.setValue(100)

    return df_found_polymers


def drop_glycan_permutations(df):
    """
    Drop duplicate permutations of glycans
    at different glycosylation sites.

    :param pd.DataFrame df: stage 2 search results
    :return: the deduplicated dataframe
    """

    return (df.reset_index()
              .drop_duplicates(["Mass_ID", "Isobar", "Hash"])
              .set_index(["Mass_ID", "Isobar", "Stage1_hit", "Stage2_hit"]))
