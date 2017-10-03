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
             was found; otherwise two dataframes:
             (1) a dataframe describing the search results;

             multiindex:

             a. Mass_ID (corresponds to index of experimental mass
                in input list of peaks),
             b. Isobar (0-based consecutive numbering of found masses) and
             c. Stage1_hit (0-based numbering of hits per Mass_ID)

             columns:

             * one column for each modification (e.g., Hex, HexNAc, ...)
             * Exp_Mass
             * Theo_Mass
             * Da
             * ppm

             If no hit was found, all columns contain the value 0

             (2) a dataframe describting the search statistics
                 with one row per unexplained mass
                 columns: `search_space_size`, `stage1_results`
    """

    sorted_mods = sorted(mods, key=lambda t: t[1], reverse=True)
    mod_names = [m[0] for m in sorted_mods]
    mod_masses = np.array([m[1] for m in sorted_mods])
    sorted_mods = [(m[1], m[2]) for m in sorted_mods]
    combinations = {}
    statistics = []
    any_combination_found = False  # true as soon as a combination is found

    # run the search on each peak
    for mass_index, unexplained_mass in enumerate(unexplained_masses):
        try:
            current_tolerance = mass_tolerance[mass_index]
        except TypeError:
            current_tolerance = mass_tolerance

        if progress_bar is not None:
            progress_bar.setValue(int((mass_index + 1)
                                      / len(unexplained_masses) * 100))
        search_space_size, result = findmods.examine_modifications(
            sorted_mods, unexplained_mass, current_tolerance)
        statistics.append(search_space_size)

        # (a) transform result to combs_per_mass, a list of dicts
        combs_per_mass = []
        if result:
            any_combination_found = True
            for r in result:
                theoretical_mass = (explained_mass
                                    + sum(np.array(r) * mod_masses))
                experimental_mass = explained_mass + unexplained_mass
                delta = experimental_mass - theoretical_mass
                r = dict(
                    zip(mod_names, r),  # number of occurrences for each mod
                    Theo_Mass=theoretical_mass,  # mass of protein + found mods
                    Exp_Mass=experimental_mass,  # experimentally measured mass
                    Da=delta,  # residual unexplained mass
                    abs_Da=abs(delta),  # absolute residual unexplained mass
                    ppm=(delta / theoretical_mass * 1000000))  # Da in ppm
                combs_per_mass.append(r)
        else:  # no appropriate combination was found
            combs_per_mass.append(dict(Theo_Mass=0.0, abs_Da=0.0))

        # (b) create a dataframe from combs_per_mass
        # sort by abs_Da and reindex the frame starting from 1 instead of 0
        # thereby, alternative combinations will be numbered 1, 2, ...
        combs_frame = pd.DataFrame(combs_per_mass).sort_values("abs_Da")
        combs_frame.index = range(1, len(combs_frame) + 1)

        # (c) finally, append the solutions for the current peak
        # to the overall list of combinations
        combinations[mass_index] = (
            combs_frame
            .set_index("Theo_Mass", append=True, drop=False)
            .swaplevel(1, 0))

    # after processing all peaks, combinations is a dict of dataframes
    # each dataframe contains information on all found combinations
    # for a single peak
    # now, combine all those dataframes into a single dataframe
    # there will be three indices: (1) consecutive numbers
    #                              (2) theoretical mass
    #                              (3) 1-based counter for hits per exp. mass
    # replace NaNs in mod columns by -1 and convert those columns to int
    if not any_combination_found:
        return None, None
    combinations = (
        pd.concat(combinations)
        .fillna(0)
        .astype({m: int for m in mod_names}))

    # sort columns so that they have the original order from the monomer table
    sortlist = [m[0] for m in mods] + ["Exp_Mass", "Theo_Mass", "Da", "ppm"]
    combinations = combinations.reindex_axis(sortlist, axis="columns")
    combinations.index.names = ["Mass_ID", "Isobar", "Stage1_hit"]

    # amend the index
    combinations.reset_index(["Isobar", "Stage1_hit"],
                             drop=True, inplace=True)
    # (a) "Isobar" should be a zero-based consecutive (integer) numbering
    #     except for unannotated peaks, where it should be -1
    isodict = {v: i for i, v in enumerate(combinations["Theo_Mass"].unique())}
    combinations["Isobar"] = combinations["Theo_Mass"].map(isodict)
    if isodict.get(0.0) is not None:
        combinations["Isobar"] -= 1

    # (b) "Stage1_hit" should be a counter per isobar
    #     except for unannotated peaks, where it should be -1
    combinations["Stage1_hit"] = (
        combinations
        .groupby("Isobar")
        .cumcount())
    combinations.loc[combinations["Isobar"] == -1, "Stage1_hit"] = -1

    # prepare search statistics
    df_statistics = pd.DataFrame(statistics, columns=["search_space_size"])
    df_statistics.index.name = "Mass_ID"
    df_statistics["stage1_results"] = (
        combinations[combinations["Stage1_hit"] != -1]
        .groupby(level="Mass_ID")
        .size())
    return (combinations.set_index(["Isobar", "Stage1_hit"], append=True),
            df_statistics.fillna(0).astype(int))


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


# noinspection PyUnresolvedReferences
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
    :return: two dataframes:
            (1) a dataframe describing the search results, like::

                                             Hex  HexNAc  Neu5Ac  Fuc  [...]
        # Mass_ID Isobar Stage2_hit Perm_ID
        # 0       5      0          0          6       8       0    2  [...]
        # 1       27     0          0          8       6       1    1  [...]
        #                           1          8       6       1    1  [...]
        #                1          0          5       4       1    1  [...]
        #                           1          5       4       1    1  [...]
        # 2       37     0          0          7       8       0    2  [...]
        #                           1          7       8       0    2  [...]
        #
        #                                     Exp_Mass    Theo_Mass        Da
        # Mass_ID Isobar Stage2_hit Perm_ID
        # 0       5      0          0        148057.12  148057.8620 -0.739772
        # 1       27     0          0        148122.21  148120.8700  1.341018
        #                           1        148122.21  148120.8700  1.341018
        #                1          0        148122.21  148120.8700  1.341018
        #                           1        148122.21  148120.8700  1.341018
        # 2       37     0          0        148220.11  148220.0030  0.109210
        #                           1        148220.11  148220.0030  0.109210
        #
        #                                          ppm     ch_A     ch_B
        # Mass_ID Isobar Stage2_hit Perm_ID
        # 0       5      0          0        -4.996506     2G0F    A2G0F
        # 1       27     0          0         9.053539  A2S1G0F       M4
        #                           1         9.053539       M4  A2S1G0F
        #                1          0         9.053539   Unglyc  A2S1G1F
        #                           1         9.053539  A2S1G1F   Unglyc
        # 2       37     0          0         0.736810    A2G0F    A2G1F
        #                           1         0.736810    A2G1F    A2G0F
        #
        #                                        Score  Permut  Permutation
        # Mass_ID Isobar Stage2_hit Perm_ID             ations  score
        # 0       5      0          0        100.000000      1  100.000000
        # 1       27     0          0          4.403297      2    8.806594
        #                           1          4.403297      2    8.806594
        #                1          0         45.596703      2   91.193406
        #                           1         45.596703      2   91.193406
        # 2       37     0          0         50.000000      2  100.000000
        #                           1         50.000000      2  100.000000

            (2) a dataframe describing the search statistics
                one row per mass
                columns: `stage2_results`, `stage2_uniques`
    """

    if progress_bar is not None:
        progress_bar.setValue(0)

    try:
        df_found_polymers = (
            stage_1_results
            .join(polymer_combinations, on=monomers, how="inner")
            .sort_index()
            .rename(columns={"Abundance": "Permutation score"}))
    except TypeError:
        return None

    if progress_bar is not None:
        progress_bar.setValue(25)

    # recalculate the scores to represent the contribution
    # of each annotation to the peak height
    df_found_polymers["Permutation score"] = (
        df_found_polymers["Permutation score"]
        * 100
        / df_found_polymers.groupby("Mass_ID")["Permutation score"].sum())

    if progress_bar is not None:
        progress_bar.setValue(50)

    # calculate glycan hash and counts/abundance per hash/peak/isobar
    df_found_polymers["Hash"] = (
        df_found_polymers
        .iloc[:, df_found_polymers.columns.get_loc("ppm") + 1: -1]
        .apply(lambda x: hash(frozenset(x)), axis=1))

    # column "Stage2_hit" numbers unique permutations per mass index
    df_found_polymers.reset_index(inplace=True)
    unique_hashes = df_found_polymers.groupby("Mass_ID")["Hash"].unique()
    df_found_polymers["Stage2_hit"] = (
        df_found_polymers
        .apply(lambda x: list(unique_hashes[x["Mass_ID"]]).index(x["Hash"]),
               axis=1))

    if progress_bar is not None:
        progress_bar.setValue(75)

    # create columns Perm_ID (index of permutation), Permutations (count)
    # and Permutation score
    df_found_polymers.set_index(["Isobar", "Stage2_hit"], inplace=True)
    hash_group = df_found_polymers.groupby(df_found_polymers.index.names)
    df_found_polymers["Perm_ID"] = hash_group.cumcount()
    df_found_polymers["Permutations"] = hash_group.size()
    df_found_polymers["Hit score"] = hash_group["Permutation score"].sum()

    df_found_polymers = (
        df_found_polymers
        .drop(["Stage1_hit", "Hash"], axis=1)
        .set_index(["Mass_ID", "Perm_ID"], append=True)
        .reorder_levels(["Mass_ID", "Isobar", "Stage2_hit", "Perm_ID"])
        .sort_index()
        .fillna({"Permutation score": 0, "Hit score": 0}))

    if progress_bar is not None:
        progress_bar.setValue(100)

    # add rows for unannotated peaks
    missing_rows = pd.DataFrame(
        0,
        columns=df_found_polymers.columns,
        index=stage_1_results.index.levels[0].difference(
            df_found_polymers.index.levels[0]))
    missing_rows.index.name = "Mass_ID"

    missing_rows = (
        missing_rows
        .assign(Isobar=-1, Stage2_hit=-1, Perm_ID=-1)
        .set_index(["Isobar", "Stage2_hit", "Perm_ID"], append=True))
    df_stage2_results = (pd.concat([df_found_polymers, missing_rows])
                         .sort_index())

    # calculate search statistics
    df_statistics = pd.DataFrame(
        df_found_polymers.groupby(level="Mass_ID").size(),
        columns=["stage2_results"])
    df_statistics["stage2_uniques"] = (
        df_found_polymers
        .reset_index("Stage2_hit")
        .groupby(level="Mass_ID")
        ["Stage2_hit"]
        .nunique())

    return df_stage2_results, df_statistics


def drop_glycan_permutations(df):
    """
    Drop duplicate permutations of glycans
    at different glycosylation sites.

    :param pd.DataFrame df: stage 2 search results
    :return: the deduplicated dataframe
    """

    old_columns = df.columns
    return (df.reset_index()
              .drop_duplicates(["Mass_ID", "Isobar", "Hash"])
              .set_index(old_columns))
