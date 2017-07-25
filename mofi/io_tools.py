"""
io_tools.py

Helper functions for input and output.

Authors: Stefan Senn, Wolfgang Skala

(c) 2017 Christian Doppler Laboratory for Biosimilar Characterization
"""

import xml.etree.ElementTree as ETree

import pandas as pd


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
    df = pd.read_excel(filename)

