import pandas as pd
import glyco_tools
import mass_tools
import re

df = pd.read_csv("../data/tnfr (Fig 3a)/library.csv", sep="\t")
# print(df)
prog = re.compile(r"(\d+)\s+([\w-]*)")

for _, row in df.iterrows():
    composition = prog.findall(row["Composition"])
    formula = mass_tools.combine_formulas([glyco_tools.glycan_formula[glycan] * int(count) for count, glycan in composition])
    # print(row["Name"], dict(formula.composition))
se = df.groupby("Site")["Name"].apply(list)
for i, v in enumerate(se):
