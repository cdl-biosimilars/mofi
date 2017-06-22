import numpy as np
import pandas as pd
import glyco_tools
import re
from itertools import product

prog = re.compile(r"(\d+)\s+([\w-]*)")

# dataframe with monosaccharides as columns, atoms as index, and atom counts as values
df_monosaccharides = pd.DataFrame(glyco_tools.mono_glycan_atoms)\
    .fillna(0)\
    .astype(np.uint16)

# dataframe that contains the original glycan library
df_library = pd.read_csv("../data/tnfr (Fig 3a)/library_large.csv", sep="\t")


def calc_atom_counts(value):
    """
    Calculare atom counts for glycans from their monosaccharide composition

    :param value: List of (count, monosaccharide) tuples
                  as returned from the regex search on column "Composition" in df_library
    :return: a list with atom counts
    """
    result = None
    for v in value:
        try:
            result = result.add(df_monosaccharides[v[1]] * int(v[0]))
        except AttributeError:
            result = df_monosaccharides[v[1]] * int(v[0])
    return list(result)


# dataframe that contains the glycans from the library as index and the elemental composition (amongst others) as column
df_glycans = df_library[["Name", "Composition"]]\
    .drop_duplicates(subset=["Name"])\
    .set_index("Name")
df_glycans["Elements"] = df_glycans["Composition"]\
    .apply(prog.findall)\
    .apply(calc_atom_counts)

# a Series with the glycosylation sites as index and lists of possible glycans as values
se_mods_per_site = df_library\
    .groupby("Site")["Name"]\
    .apply(list)


def calc_total(row):
    result = None
    for glycan in row:
        try:
            result += np.array(df_glycans["Elements"][glycan])
        except TypeError:
            result = np.array(df_glycans["Elements"][glycan])
    return tuple(result)


df_glycan_combinations = pd.DataFrame(list(product(*list(se_mods_per_site))), columns=list(se_mods_per_site.index))
df_glycan_combinations["Formula"] = df_glycan_combinations\
    .apply(calc_total, axis=1)
df_glycan_combinations = df_glycan_combinations\
    .set_index(
        pd.MultiIndex.from_tuples(df_glycan_combinations["Formula"], names=["C", "H", "N", "O"])
    )\
    .sort_index()
df_glycan_combinations.drop("Formula", 1, inplace=True)



# print(df_monosaccharides)
# print(df_glycans)
# print(se_mods_per_site)
print(df_glycan_combinations.loc[226, 370, 14, 164])

# TODO: creating a dataframe with an element (CHNO) multiindex in reasonably short time works;
# TODO: next: search with results from combinatoric search in this dataframe and output results! i.e., integrate with MoFi
