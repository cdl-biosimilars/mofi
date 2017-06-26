import numpy as np
import pandas as pd
import glyco_tools
import re
from itertools import product


def calc_monomer_counts(value):
    """
    Calculate monomer counts for a given complex glycan

    :param value: list of (count, monomer) tuples
    :return: numpy array containing the monomer counts (including absent monomers)
    """
    composition = dict([(v, int(k)) for (k, v) in value])
    result = np.zeros(len(glyco_tools.monosaccharides), dtype=np.uint16)
    for i, m in enumerate(glyco_tools.monosaccharides):
        result[i] = composition.setdefault(m, 0)
    return result


def sum_glycan_monomer_counts(row):
    """
    Sum the monomer counts of several glycans.

    :param row: row from a dataframe containing the name of a glycan per cell
    :return: a tuple containing the numbers of each atom
    """
    result = np.zeros(len(glyco_tools.monosaccharides), dtype=np.uint16)
    for glycan in row:
        result += df_glycans_mono["Monomers"][glycan]
    return tuple(result)


# dataframe that contains the original glycan library
df_library = pd.read_csv("../data/tnfr_fig_3a/library_large_ocore.csv", sep="\t")

# a Series with the glycosylation sites as index and lists of possible glycans as values
rows = []
for _, row in df_library.iterrows():
    for site in row["Site"].split(","):
        rows.append([site.strip(), row["Name"]])
mods_per_site = pd.DataFrame(rows, columns=["Site", "Name"])\
    .groupby("Site")["Name"]\
    .apply(list)

# dataframe with glycans as index and monomer composition as rows
prog = re.compile(r"(\d+)\s+([\w-]*)")
df_glycans_mono = df_library[["Name", "Composition"]]\
    .drop_duplicates(subset=["Name"])\
    .set_index("Name")\
    .fillna("")
df_glycans_mono["Monomers"] = df_glycans_mono["Composition"]\
        .apply(prog.findall)\
        .apply(calc_monomer_counts)

# a dataframe with monomers (Hex, HexNAc, ...) as multiindex and one column per site,
# indicating which glycan composition belongs to a given elemental composition
df_glycan_combinations = pd.DataFrame(list(product(*list(mods_per_site))), columns=list(mods_per_site.index))
df_glycan_combinations["Monomers"] = df_glycan_combinations\
    .apply(sum_glycan_monomer_counts, axis=1)
df_glycan_combinations = df_glycan_combinations\
    .set_index(
        pd.MultiIndex.from_tuples(df_glycan_combinations["Monomers"], names=glyco_tools.monosaccharides)
    )\
    .sort_index()\
    .drop("Monomers", 1)

# print(df_library)
# print(mods_per_site)
print(df_glycans_mono)
# df_glycan_combinations.to_csv("../data/tnfr (Fig 3a)/glycancombs.csv")
# print(df_glycan_combinations)
# print(df_glycan_combinations.loc[10, 10, 0, 0, 2, 0, 4, 0])


# df_combinations = pd.read_csv("../data/tnfr (Fig 3a)/annotation.csv")
# print(df_combinations)
