"""
    output_tools.py

    Helper functions/classes for output.

    Authors: Stefan Senn, Wolfgang Skala

    (c) 2017 Christian Doppler Laboratory for Biosimilar Characterization
"""


import time


def write_hits_to_csv(hitdf, outfilename, parameters=None):
    """
    Write the results from a combinatorial search to a CSV file.

    :param hitdf: Dataframe that should be saved
    :param outfilename: name of the CSV file
    :param parameters: list of extra (parameter name, value) tuples that should be stored in the CSV file
    :return: nothing
    """
    if parameters:
        parameter_string = ["{}: {}".format(k, v) for k, v in parameters]
        with open(outfilename, "w") as f:
            f.write("# Combinatorial search results by ModFinder\n")
            f.write("# Date: " + time.strftime("%c") + "\n")
            f.write("# Parameters: ")
            f.write(", ".join(parameter_string))
            f.write("\n")
            hitdf.to_csv(f)
    else:
        hitdf.to_csv(outfilename)
