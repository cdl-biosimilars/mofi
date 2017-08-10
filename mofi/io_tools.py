"""
io_tools.py

Helper functions for input and output.

Author: Wolfgang Skala

(c) 2017 Christian Doppler Laboratory for Biosimilar Characterization
"""

import re
import xml.etree.ElementTree as ETree
import pandas as pd

# regular expression for parsing glycans according to the nomenclature
# of Thermo Fisher BioPharma Finder
_re_bpf_glycan = re.compile(r"""
        ^(?:(A)(\d)+)?  # g[1]:  antennas
        (?:(Sg)(\d)+)?  # g[3]:  Neu5Gc
        (?:(S)(\d)+)?   # g[5]:  Neu5Ac
        (?:(Ga)(\d)+)?  # g[7]:  Ga
        (?:(G)(\d)+)?   # g[9]:  Gal
        (?:(M)(\d)+)?   # g[11]: Man
        (F)?            # g[12]: Fuc
        (B)?            # g[13]: GlcNAc
        """, re.VERBOSE)


def _parse_bpf_glycan(glycan):
    """
    Convert a glycan abbreviation (e.g., "A2G1F") to a composition string
    as used in the polymer table (e.g., "1 Hex, 1 HexNAc, 1 Fuc, 1 N-core").

    :param glycan: String representing a glycan as described
                   in the BioPharma Finder manual

    :return: a composition string
    """

    g = list(_re_bpf_glycan.findall(glycan)[0])
    if not "".join(g):
        # deal with two non-standard abbreviations
        if glycan == "Gn":
            counts = [0, 1, 0, 0, 0]
        elif glycan == "GnF":
            counts = [0, 1, 0, 0, 1]
        else:
            counts = [0, 0, 0, 0, 0]
    else:
        if g[11] == "":
            g[11] = 3
        for i in range(1, 13, 2):
            try:
                g[i] = int(g[i])
            except ValueError:
                g[i] = 0
        counts = [
            g[3] + g[5] + 2 * g[7] + g[9] + g[11],  # Hex
            g[1] + 2 + (1 if g[13] else 0),  # HexNAc
            g[5],  # Neu5Ac
            g[3],  # Neu5Gc
            1 if g[12] else 0]  # Fuc

    monosaccharides = ["Hex", "HexNAc", "Neu5Ac", "Neu5Gc", "Fuc"]
    return ", ".join("{} {}".format(c, m)
                     for m, c in zip(monosaccharides, counts)
                     if c > 0)


def prettify_xml(elem, level=0):
    """
    Prettify an XML tree inplace.

    :param elem: a xml.etree.ElementTree Element
    :param level: which level to prettify; only required for recirsive call
    :return: nothing (changes the ElementTree inplace)
    """
    i = "\n" + level * "  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            prettify_xml(elem, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


def dataframe_to_xml(df):
    """
    Convert a pandas DataFrame to an XML etree with the following format:

    <dataframe>
        <columns>
            <column dtype="[data dype]">[column name]</column>
            ...
        </columns>
        <rows>
            <row id=[id]>
                <cell>[value]</cell>
                ...
            </row>
            ...
        </rows>
     </dataframe>

    :param df: a dataframe
    :return: a xml.etree.ElementTree Element
    """

    root = ETree.Element("dataframe")
    if df is not None:
        columns = ETree.SubElement(root, "columns")
        for i, column in enumerate(df.columns.values):
            item = ETree.SubElement(columns, "column", dtype=str(df.dtypes[i]))
            item.text = column
        rows = ETree.SubElement(root, "rows")
        for row_id, row in df.iterrows():
            item = ETree.SubElement(rows, "row", id=str(row_id))
            for cell_id, cell in row.iteritems():
                ETree.SubElement(item, "cell").text = str(cell)
    return root


def dataframe_from_xml(root):
    """
    Convert an XML etree to a pandas dataframe.

    :param root: an xml.etree.ElementTree Element
    :return: a dataframe
    """

    root = root.find("dataframe")
    columns = root.find("columns")
    if columns is None:
        return pd.DataFrame()
    else:
        column_names = []
        dtypes = []
        for column in columns.findall("column"):
            column_names.append(column.text)
            dtypes.append(column.attrib["dtype"])

        rows = root.find("rows")
        rows_data = []
        for row in rows.findall("row"):
            row_data = []
            for cell in row.findall("cell"):
                row_data.append(cell.text)
            rows_data.append(row_data)

        return (pd.DataFrame(rows_data, columns=column_names)
                .replace({"True": True, "False": False})
                .astype(dict(zip(column_names, dtypes))))


def read_bpf_library(filename):
    """
    Read an N-glycan library as returned by Thermo BioPharma finder.
    A column labeled "Modification" is required; sites and composition
    will be deduced from this column.

    :param filename: name of the library file (Excel format)
    :return: a dataframe as required by ModFinder.table_from df()
    """

    df = pd.read_excel(filename)
    df["Sites"] = df["Modification"].apply(lambda x: x.split("+")[0])
    df["Name"] = df["Modification"].apply(lambda x: x.split("+")[1])
    df["Composition"] = df["Name"].apply(_parse_bpf_glycan)
    return df
