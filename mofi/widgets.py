"""
Custom widgets and GUI functions.
"""

import os

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import (QHeaderView, QLineEdit, QFileDialog,
                             QTableWidgetItem, QTreeWidgetItem)
from PyQt5.QtGui import QColor, QBrush

from matplotlib.widgets import RectangleSelector

# noinspection PyPep8Naming,PyUnresolvedReferences
class FilterHeader(QHeaderView):
    """
    A :class:`QHeaderView` with :class:`QLineEdit` widgets below the columns,
    whose values filter the table's entries.

    :ivar list(QLineEdit) _editors: (column index, line edit widget) tuples
    :ivar int _vertical_padding: vertical padding of the editors
    :ivar int _horizontal_padding: horizontal padding of the editors

    .. automethod:: __init__
    """
    filterChanged = pyqtSignal()


    def __init__(self, parent):
        """
        Create a new filter header.

        :param QWidget parent: parent widget
        """

        super().__init__(Qt.Horizontal, parent)
        self._editors = []
        self._vertical_padding = 4
        self._horizontal_padding = 4
        self.sectionResized.connect(self.adjustPositions)
        parent.horizontalScrollBar().valueChanged.connect(
            self.adjustPositions)


    def setFilterBoxes(self, sections):
        """
        Add :class:`QLineEdit` widgets below the column headers.

        :param list(int) sections: column indices for which to create a filter
        :return: nothing
        """
        while self._editors:
            editor = self._editors.pop()
            editor[1].deleteLater()
        for index in sections:
            editor = QLineEdit(self.parent())
            editor.setPlaceholderText("Filter")
            editor.returnPressed.connect(self.filterChanged.emit)
            self._editors.append((index, editor))
        self.adjustPositions()


    def sizeHint(self):
        """
        Calculate a size hinz which takes the filters into account.

        :return: a size hint
        :rtype: QSizeHint
        """
        size = super().sizeHint()
        if self._editors:
            height = self._editors[0][1].sizeHint().height()
            size.setHeight(size.height() + height + self._vertical_padding)
        return size


    def updateGeometries(self):
        """
        Takes editors into account for the geometry.

        :return: nothing
        """
        if self._editors:
            height = self._editors[0][1].sizeHint().height()
            self.setViewportMargins(0, 0, 0, height + self._vertical_padding)
        else:
            self.setViewportMargins(0, 0, 0, 0)
        super().updateGeometries()
        self.adjustPositions()


    def adjustPositions(self):
        """
        Adjust the position of the filter widgets.

        :return: nothing
        """

        for index, editor in self._editors:
            header_height = super().sizeHint().height()
            editor_height = editor.sizeHint().height()
            editor.move(
                self.sectionPosition(index)
                - self.offset()
                + self._horizontal_padding // 2,
                header_height + self._vertical_padding // 2)
            editor.resize(
                self.sectionSize(index)
                - self._horizontal_padding, editor_height)


    def filterText(self, col_index):
        """
        Return the text of a line edit.

        :param int col_index: index of the line edit to be queried
        :return: the line edit's text
        :rtype: str
        :raises KeyError: if no filter widget is at ``col_index``
        """

        for index, editor in self._editors:
            if index == col_index:
                return editor.text()
            raise KeyError("No filter at column index '{}'".format(col_index))


    def allFilters(self):
        """
        A generator that returns the contents
        of :attr:`~FilterHeader._editors`.

        :return: a (col_index, :class:`QLineEdit`) tuple generator
        :rtype: generator
        """
        yield from self._editors


    def setFilterText(self, col_index, text):
        """
        Set the text of a line edit.

        :param int col_index: the line edit's index
        :param str text: text to be set
        :return: nothing
        :raises KeyError: if no filter widget is at ``col_index``
        """
        for index, editor in self._editors:
            if index == col_index:
                editor.setText(text)
                return
            raise KeyError("No filter at column index '{}'".format(col_index))


    def clearFilters(self):
        """
        Clear contents of all line edits.

        :return: nothing
        """
        for _, editor in self._editors:
            editor.clear()
        self.filterChanged.emit()


class SortableTreeWidgetItem(QTreeWidgetItem):
    """
    A :class:`QTreeWidgetItem` which supports numerical sorting
    and implements custom methods.

    :ivar function default_key: default function used in sorting;
      applied if :attr:`~SortableTreeWidgetItem.col_key` does not
      provide a mapping for a column
    :ivar dict col_key: {column index: key function} mapping

    .. automethod:: __init__
    .. automethod:: __lt__
    """

    def __init__(self, parent=None, default_key=None, col_key=None):
        """
        Create a new sortable tree widget item.

        :param QWidget parent: parent widget
        :param func default_key: default function used in sorting
        :param dict col_key: {column index: key function} mapping
        """
        super().__init__(parent)
        if default_key is None:
            self.default_key = lambda x: x
        else:
            self.default_key = default_key
        self.col_key = col_key


    def __lt__(self, other):
        """
        Used for sorting. Apply the appropriate key function.

        :param SortableTreeWidgetItem other: object to which self is compared
        :return: True if self is less than other
        :rtype: bool
        """

        column = self.treeWidget().sortColumn()
        key_func = self.col_key.get(column, self.default_key)

        key1 = self.text(column)
        key2 = other.text(column)
        try:
            return key_func(key1) < key_func(key2)
        except ValueError:
            return key1 < key2


    # noinspection PyPep8Naming
    def setTotalBackground(self, color=None, column_count=None,
                           alternating=None, even_color=None, odd_color=None):
        """
        Set the background color for all columns, possibly depending on
        whether the widget has an even or odd line number.

        :param str color: the color to use
        :param int column_count: the intended number of columns
        :param int alternating: If this parameter is an even (odd) int,
          the color specified by ``even_color`` (``odd_color``) will be used.
        :param str even_color: Color for even widgets.
                               If None, fall back to the value of ``color``.
        :param str odd_color: Color for odd widgets.
                              If None, fall back to the value of ``color``.
        :return: nothing
        """

        if column_count is not None and self.columnCount() != column_count:
            self.setText(column_count - 1, "")
        if alternating is not None:
            if alternating % 2 == 0:
                if even_color is not None:
                    color = even_color
            else:
                if odd_color is not None:
                    color = odd_color

        for i in range(self.columnCount()):
            self.setBackground(i, QBrush(QColor(color)))


    # noinspection PyPep8Naming
    def getTopParent(self):
        """
        Find the top-level ancestor of self.

        :return: the top-level parent of self
        :rtype: SortableTreeWidgetItem
        """
        node = self
        while node.parent():
            node = node.parent()
        return node


class SortableTableWidgetItem(QTableWidgetItem):
    """
    A :class:`QTableWidgetItem` which supports numerical sorting.

    .. automethod:: __init__
    .. automethod:: __lt__
    """

    def __init__(self, parent=None):
        """
        Create a new sortable table widget item.

        :param QWidget parent: parent widget
        """
        super().__init__(parent)

    def __lt__(self, other):
        """
        Compare two items.

        :param SortableTableWidgetItem other: item to which self is compared
        :return: True if self is less than other
        :rtype: bool
        """
        key1 = self.text()
        key2 = other.text()
        try:
            return float(key1) < float(key2)
        except ValueError:
            return key1 < key2


class CollapsingRectangleSelector(RectangleSelector):
    """
    Select a rectangular region of an axes.

    :ivar float collapsex: The rectangle collapses to a line
                           if its x-dimension is less than this ...
    :ivar float collapsey: ... and its y-dimension is less than this variable.

    .. automethod:: __init__
    """

    def __init__(self, *args, collapsex=0, collapsey=0, **kwargs):
        """
        Create a new selector.

        :param args: positional arguments passed to the superclass
        :param float collapsex: minimum height of the rectangle (in data space)
        :param float collapsey: minimum width of the rectangle (in data space)
        :param kwargs: optional arguments passed to the superclass
        """
        super().__init__(*args, **kwargs)
        self.collapsex = collapsex
        self.collapsey = collapsey


    def draw_shape(self, extents):
        """
        Calculate the coordinates of the drawn rectangle
        (overrides method from parent).

        :param tuple(float) extents: Coordinates of the rectangle selector
        :return: nothing
        """
        xmin, xmax, ymin, ymax = extents

        # Collapse coordinates if their distance is too small
        if abs(xmax - xmin) < self.collapsex:
            xmax = xmin
        if abs(ymax - ymin) < self.collapsey:
            ymax = ymin

        self.to_draw.set_x(xmin)
        self.to_draw.set_y(ymin)
        self.to_draw.set_width(xmax - xmin)
        self.to_draw.set_height(ymax - ymin)


_default_file_types = {
    "xls": "Excel files",
    "xlsx": "Excel files",
    "csv": "CSV files",
    "fasta": "Sequence files",
    "xml": "ModFinder XML settings",
    "": ""
}

def get_filename(parent, kind="save", caption="", directory="",
                 extensions=None, file_types=None):
    """
    Get a filename by a :class:`QFileDialog`
    and automatically add extensions.

    :param QWidget parent: parent of the dialog
    :param str kind: ``"save"`` or ``"open"``, chooses the dialog type
    :param str caption: caption of the dialog
    :param str directory: initial directory
    :param list(str) extensions: file extensions
    :param dict file_types: {extension: description} dict
    :return: the file name with extension and the last path
    :rtype: tuple(str, str)
    :raise ValueError: if an invalid value was supplied to ``kind``
    """

    if extensions is None:
        extensions = [""]
    if file_types is None:
        file_types = _default_file_types
    ext_list = {k: "{0} [{1}] (*.{1})".format(v, k)
                for k, v in file_types.items()}
    reverse_extensions = {v: k for k, v in ext_list.items()}

    if kind == "save":
        dialog = QFileDialog.getSaveFileName
    elif kind == "open":
        dialog = QFileDialog.getOpenFileName
    else:
        raise ValueError("Unknown value for 'kind': " + kind)
    filename, filefilter = dialog(
        parent,
        caption,
        directory,
        ";;".join([ext_list[ext] for ext in extensions]))
    new_path = os.path.split(filename)[0]

    if not filename:
        return None, new_path

    if not filename.endswith(reverse_extensions[filefilter]):
        filename += "." + reverse_extensions[filefilter]
    return filename, new_path
