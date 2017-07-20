import pandas as pd
import xml.etree.ElementTree as ETree
from io import StringIO


def indent(elem, level=0):
    i = "\n" + level * "  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


df = pd.DataFrame([(1, 2, "foo"), (3, 4, "bar")], columns=["eggs", "spam", "lobster"])

root = ETree.Element("dataframe")
columns = ETree.SubElement(root, "columns")
for name in df.columns.values:
    ETree.SubElement(columns, "name").text = name
rows = ETree.SubElement(root, "rows")
for row_id, row in df.iterrows():
    item = ETree.SubElement(rows, "row", attrib=dict(id=str(row_id)))
    for cell_id, cell in row.iteritems():
        ETree.SubElement(item, "cell").text = str(cell)

indent(root)
# etree.dump(root)

buf = StringIO(df.to_csv(index=False))
print(pd.read_csv(buf))
