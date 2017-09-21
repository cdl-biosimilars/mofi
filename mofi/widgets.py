from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import QHeaderView, QLineEdit


# noinspection PyPep8Naming,PyUnresolvedReferences
class FilterHeader(QHeaderView):
    filterActivated = pyqtSignal()


    def __init__(self, parent):
        super().__init__(Qt.Horizontal, parent)
        self._editors = []
        self._vertical_padding = 4
        self._horizontal_padding = 4
        self.sectionResized.connect(self.adjustPositions)
        parent.horizontalScrollBar().valueChanged.connect(
            self.adjustPositions)


    def setFilterBoxes(self, count):
        """
        Add :class:`QLineEdit` widgets below the column headers.

        :param int count: number of line edits
        :return: nothing
        """
        while self._editors:
            editor = self._editors.pop()
            editor.deleteLater()
        for index in range(count):
            editor = QLineEdit(self.parent())
            editor.setPlaceholderText("Filter")
            editor.returnPressed.connect(self.filterActivated.emit)
            self._editors.append(editor)
        self.adjustPositions()


    def sizeHint(self):
        size = super().sizeHint()
        if self._editors:
            height = self._editors[0].sizeHint().height()
            size.setHeight(size.height() + height + self._vertical_padding)
        return size


    def updateGeometries(self):
        if self._editors:
            height = self._editors[0].sizeHint().height()
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

        for index, editor in enumerate(self._editors):
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


    def filterText(self, index):
        """
        Return the text of a lne edit.

        :param index: index of the line edit
        :return: the line edit's text
        """
        if 0 <= index < len(self._editors):
            return self._editors[index].text()
        return ""


    def setFilterText(self, index, text):
        """
        Set the text of a line edit.

        :param index: the line edit's index
        :param text: text to be set
        :return: nothing
        """
        if 0 <= index < len(self._editors):
            self._editors[index].setText(text)


    def clearFilters(self):
        """
        Clear contents of all line edits.

        :return: nothing
        """
        for editor in self._editors:
            editor.clear()
