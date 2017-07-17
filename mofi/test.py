import pandas as pd
import xml.etree.ElementTree as etree


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

root = etree.Element("dataframe")
etree.SubElement(root, "test").text = "123"
for row_id, row in df.iterrows():
    item = etree.SubElement(root, "row", attrib=dict(id=str(row_id)))
    for cell_id, cell in row.iteritems():
        subitem = etree.SubElement(item, "cell", attrib=dict(col=cell_id))
        subitem.text = str(cell)

indent(root)
etree.dump(root)