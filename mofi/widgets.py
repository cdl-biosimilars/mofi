#  from https://stackoverflow.com/questions/44343738/how-to-inject-
# widgets-between-qheaderview-and-qtableview

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import QHeaderView, QLineEdit, QTreeWidgetItem
from PyQt5.QtGui import QColor, QBrush

from matplotlib.widgets import RectangleSelector

# noinspection PyPep8Naming,PyUnresolvedReferences
class FilterHeader(QHeaderView):
    filterChanged = pyqtSignal()


    def __init__(self, parent):
        super().__init__(Qt.Horizontal, parent)
        self._editors = []  # list of (col index, QLineEdit) tuples
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
        size = super().sizeHint()
        if self._editors:
            height = self._editors[0][1].sizeHint().height()
            size.setHeight(size.height() + height + self._vertical_padding)
        return size


    def updateGeometries(self):
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
        :raises KeyError: if no filter widget is at :arg:`col_index`
        """

        for index, editor in self._editors:
            if index == col_index:
                return editor.text()
            raise KeyError("No filter at column index '{}'".format(col_index))


    def allFilters(self):
        """
        A generator that returns the contents of :var:`self._editors`.

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
        :raises KeyError: if no filter widget is at :arg:`col_index`
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
    A QTreeWidget which supports numerical sorting.
    :meth:`setTotalBackground` set the background of all columns.
    """

    def __init__(self, parent=None):
        super().__init__(parent)

    def __lt__(self, other):
        column = self.treeWidget().sortColumn()
        key1 = self.text(column)
        key2 = other.text(column)
        try:
            return float(key1) < float(key2)
        except ValueError:
            return key1 < key2

    # noinspection PyPep8Naming
    def setTotalBackground(self, color, column_count=None):
        """
        Set the background color for all columns.

        :param QColor color: the color to use
        :param int column_count: the intended number of columns
        :return: nothing
        """

        if column_count and self.columnCount() != column_count:
            self.setText(column_count - 1, "")
        for i in range(self.columnCount()):
            self.setBackground(i, QBrush(color))

    # noinspection PyPep8Naming
    def getTopParent(self):
        """
        Find the top-level ancestor of self.

        :return: the top-level parent of self
        """
        node = self
        while node.parent():
            node = node.parent()
        return node


class CollapsingRectangleSelector(RectangleSelector):
    """
    Select a rectangular region of an axes.
    The rectangle collapses to a line if a dimension
    is less than minspanx and minspany.

    .. automethod:: __init__
    """

    def __init__(self, *args, collapsex=0, collapsey=0, **kwargs):
        """
        Introduces the keys ``collapsex`` and ``collapsey``.

        :param args: Positional arguments passed to the superclass.
        :param collapsex: Minimum height of the rectangle (in data space).
        :param collapsey: Minimum width of the rectangle (in data space).
        :param kwargs: Optional arguments passed to the superclass.
        """
        super().__init__(*args, **kwargs)
        self.collapsex = collapsex
        self.collapsey = collapsey


    def draw_shape(self, extents):
        """
        Calculate the coordinates of the drawn rectangle
        (overrides method from parent).

        :param extents: Coordinates of the rectangle selector.
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
