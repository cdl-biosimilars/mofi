import re
import pandas as pd

def formstring_to_composition(formstring):
    """
    Converts an elemental composition string to a pandas series.

    :param formstring: collection of elements followed by their counts (example: "C50 H100 N20")
    :return: pd.Series labelled by the element (example: C: 50, H: 100, N: 20)
    """

    pattern = re.compile(r"([A-Z][a-z]?)(\d*)")
    composition = {}
    for atom, count in pattern.findall(formstring):
        if count:
            composition[atom] = int(count)
        else:
            composition[atom] = 1
    return pd.Series(composition)

print(formstring_to_composition("C6 H12 O6 N Na2"))