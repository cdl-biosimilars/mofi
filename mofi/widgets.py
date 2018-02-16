"""
Custom widgets and GUI functions.
"""

from collections import Counter
import os
import re
from urllib.request import pathname2url
import webbrowser

import pandas as pd

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import (QHeaderView, QLineEdit, QFileDialog,
                             QTableWidgetItem, QTreeWidgetItem, QDialog,
                             QDialogButtonBox, QGridLayout, QLabel,
                             QComboBox, QMessageBox)
from PyQt5.QtGui import QColor, QBrush

from matplotlib.widgets import RectangleSelector

from mofi.configure import dec_places
from mofi.createpointmutation_ui import Ui_CreatePointMutation
from mofi.createtruncation_ui import Ui_CreateTruncation
from mofi.importtabdata_ui import Ui_ImportTabData
from mofi.mass_tools import Formula
from mofi.paths import docs_dir
from mofi.sequence_tools import (amino_acid_compositions, amino_acid_names,
                                 read_fasta_string)

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


class ImportTabDataDialog(QDialog, Ui_ImportTabData):
    """
    A dialog for importing tabluar data represented by a CSV or XLS file.

    .. automethod:: __init__
    """

    def __init__(self, parent=None):
        """
        Initialize the dialog.

        :param QWidget parent: parent widget
        """

        # initialize the GUI
        super().__init__(parent)
        self.setupUi(self)

        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        self.buttonBox.button(QDialogButtonBox.Apply).clicked.connect(
            self.apply_options)
        self.buttonBox.button(QDialogButtonBox.Help).clicked.connect(
            self.show_help)

        self.cb_columns = []  # data on the comboboxes for selecting columns
        self.csv_mode = True  # whether the dialog is in csv or xls mode
        self.df_in = None  # dataframe that stores imported data
        self.filename = None  # name of the import file

    @staticmethod
    def show_help():
        """
        Open the manual in a browser and jump to the import dialog section.

        :return: nothing
        """

        site = "workflow.html"
        anchor = "import-dialog"
        path = os.path.abspath(os.path.join(docs_dir, "html", site))
        url = "#".join([pathname2url(path), anchor])
        webbrowser.open("file:{}".format(url))


    def set_options(self, filename, cols, mode):
        """
        Set up the dialog to show the contents of the opened file
        and the applicable options.

        :param str filename: input file
        :param cols: see parameter ``cols``
                     in :meth:`mofi.MainWindow.table_from_df()`
        :param str mode: selects "csv" or "xls" mode
        :return: nothing
        """

        self.filename = filename
        self.csv_mode = mode == "csv"

        # set the window title, load the file contents and enable line edits
        # depending whether we're opening a csv or xls file
        if self.csv_mode:
            self.setWindowTitle("Import CSV")
            self.gbCsvFile.setTitle(
                "Contents of " + os.path.basename(filename))
            with open(filename) as f:
                self.teFileContents.setText("".join(f.readlines()))
            for w in (self.lbSheetName, self.cbSheetName):
                w.setEnabled(False)
            self.cbSheetName.setEnabled(False)
        else:
            self.setWindowTitle("Import XLS")
            self.gbCsvFile.setTitle("")
            for w in (self.leSep, self.leComment, self.leDecimal,
                      self.leThousands, self.lbSep, self.lbComment,
                      self.lbDecimal, self.lbThousands, self.lbQuote,
                      self.leQuote):
                w.setEnabled(False)
            with pd.ExcelFile(filename) as f:
                self.cbSheetName.addItems(f.sheet_names)


        # create comboboxes for selecting required/optional columns
        self.cb_columns = []
        layout = QGridLayout(self.gbColumns)
        for i, col in enumerate(cols):
            row_id = i // 2
            col_id = i % 2 * 3
            cb = QComboBox()
            cb.setMinimumHeight(22)
            layout.addWidget(QLabel(col[2] + ":"), row_id, col_id, 1, 1)
            layout.addWidget(cb, row_id, col_id + 1, 1, 1)
            self.cb_columns.append((cb, col[0], col[1]))
        layout.setColumnMinimumWidth(2, 50)


    def apply_options(self):
        """
        Generate the preview table from the data to be imported
        using the current option settings.

        :return: nothing
        """

        # calculate arguments for pandas reader
        if self.cbHeader.isChecked():
            header = self.sbSkipRows.value()
        else:
            header = None
        decimal = self.leDecimal.text()
        if decimal == "":
            decimal = "."
        thousands = self.leThousands.text()
        if thousands == "":
            thousands = None
        quotechar = self.leQuote.text()
        if quotechar == "":
            quotechar = '"'

        # import data
        try:
            if self.csv_mode:
                self.df_in = pd.read_csv(
                    self.filename,
                    sep=self.leSep.text(),
                    comment=self.leComment.text(),
                    decimal=decimal,
                    thousands=thousands,
                    quotechar=quotechar,
                    skiprows=self.sbSkipRows.value(),
                    header=header,
                    keep_default_na=False)
            else:
                self.df_in = pd.read_excel(
                    self.filename,
                    sheetname=self.cbSheetName.currentText(),
                    skiprows=self.sbSkipRows.value(),
                    header=header,
                    keep_default_na=False)
        except (ValueError, OSError) as e:
            QMessageBox.critical(self, "Error", e.args[0])
            return

        # fill the preview table
        try:
            self.twPreview.clearContents()
            self.twPreview.setColumnCount(self.df_in.columns.size)
            self.twPreview.setRowCount(self.df_in.index.size)
            self.twPreview.setHorizontalHeaderLabels(
                [str(c) for c in self.df_in.columns])
            self.twPreview.horizontalHeader().setSectionResizeMode(
                QHeaderView.ResizeToContents)
            self.twPreview.verticalHeader().setSectionResizeMode(
                QHeaderView.ResizeToContents)
            for row_id, row in self.df_in.iterrows():
                for col_id, value in enumerate(row):
                    self.twPreview.setItem(row_id, col_id,
                                           QTableWidgetItem(str(value)))
        except TypeError:
            QMessageBox.critical(self, "Error", "Invalid quote character.")
            return

        # fill the column selectors
        for cb, label, default_value in self.cb_columns:
            cb.clear()
            if default_value is not None:
                cb.addItem("(Use default: {})".format(default_value))
            for column in self.df_in.columns:
                cb.addItem(str(column))
            cb.setCurrentText(label)


    def get_df(self):
        """
        Return the input data filtered according to the column selectors.

        :return: a dataframe if available; None otherwise
        :rtype: pd.DataFrame or NoneType
        """

        if self.df_in is None:
            return None
        df = pd.DataFrame(index=self.df_in.index)
        for cb, label, default_value in self.cb_columns:
            if default_value is not None and cb.currentIndex() == 0:
                df[label] = default_value
            else:
                df[label] = self.df_in[cb.currentText()]
        return df


    @staticmethod
    def get_data(parent=None, filename=None, cols=None,
                 mode="csv", force_accept=False):
        """
        Opens an Import CSV/XLS dialog and returns the imported data.

        :param QWidget parent: parent widget
        :param str filename: input file
        :param cols: see parameter ``cols``
                     in :meth:`mofi.MainWindow.table_from_df()`
        :param str mode: selects "csv" or "xls" mode
        :param boolean force_accept: True returns the imported data without
                                     allowing the user to modify import
                                     parameters
        :return: a dataframe with the imported data if Ok is clicked;
                 ``None`` otherwise
        :rtype: pd.DataFrame or NoneType
        """

        dialog = ImportTabDataDialog(parent)
        dialog.set_options(filename, cols, mode)
        dialog.apply_options()

        if force_accept:
            return dialog.get_df()
        else:
            result = dialog.exec_()
            if result == QDialog.Accepted:
                return dialog.get_df()
            else:
                return None


class FileTypes:
    """
    A class for managing file types in :class:`~PyQt5.QtWidgets.QFileDialog` .

    :cvar dict _default_file_types: default file types
    :ivar list types: actually used file types
    :ivar list filters: filter strings

    .. automethod:: __init__
    """

    _default_file_types = {
        "xls": ("xls", "Excel 97-2003"),
        "xls-bpf": ("xls", "BioPharma Finder PepMap results"),
        "xlsx": ("xlsx", "Excel"),
        "csv": ("csv", "Comma-separated value"),
        "fasta": ("fasta", "Protein sequence"),
        "xml": ("xml", "MoFi XML settings"),
        "": ("", "all files")
    }

    def __init__(self, type_list=None):
        """
        Create a new :class:`FileTypes` object.

        :param type_list: list of (extension, description) tuples
                          or {ext: description} dict
                          or list of extensions, which are then filtered
                          from ``_default_file_types``
        :type type_list: list(tuple) or dict or list(str)
        :return: nothing
        :raise TypeError: if an invalid type was specified for ``ext_list``
        """

        if type_list is None:
            self.types = self._default_file_types
            return

        self.types = []
        try:  # dict
            for ext, description in sorted(type_list.items()):
                self.types.append((ext, description))
        except AttributeError:
            try:  # list of tuples
                for ext, description in type_list:
                    self.types.append((ext, description))
            except (TypeError, ValueError):
                try:  # list of strings
                    for ext in type_list:
                        try:
                            self.types.append(self._default_file_types[ext])
                        except KeyError:
                            pass
                except (TypeError, ValueError):
                    raise TypeError("Argument for 'ext_list' has wrong type")
        self.filters = ["{1} [{0}] (*.{0})".format(k, v)
                        for k, v in self.types]


    def get_filter(self):
        """
        Returns  a file filter suitable for the ``filter`` parameter of
        :class:`~PyQt5.QtWidgets.QFileDialog`.

        :return: a file filter
        :rtype: str
        """

        return ";;".join(self.filters)


    def get_type(self, filefilter):
        """
        Returns the extension and description associated with a file filter.

        :param str filefilter: filefilter as returned by the dialog
        :return: the extension and description of the file type
        :rtype: tuple(str, str)
        """

        return self.types[self.filters.index(filefilter)]


def get_filename(parent, kind="save", caption="", directory="",
                 file_types=None):
    """
    Get a filename by a :class:`QFileDialog`
    and automatically add extensions.

    :param QWidget parent: parent of the dialog
    :param str kind: ``"save"`` or ``"open"``, chooses the dialog type
    :param str caption: caption of the dialog
    :param str directory: initial directory
    :param FileTypes file_types: file extensions
    :return: the file name with extension, the description
             of the used filetype and the last path
    :rtype: tuple(str, str, str)
    :raise ValueError: if an invalid value was supplied to ``kind``
    """

    if file_types is None:
        file_types = FileTypes()

    if kind == "save":
        dialog = QFileDialog.getSaveFileName
    elif kind == "open":
        dialog = QFileDialog.getOpenFileName
    else:
        raise ValueError("Unknown value for 'kind': " + kind)

    filename, used_filter = dialog(
        parent,
        caption,
        directory,
        file_types.get_filter())
    new_path = os.path.split(filename)[0]

    if not filename:
        return None, None, new_path

    ext, desc = file_types.get_type(used_filter)
    if not filename.endswith(ext):
        filename += "." + ext
    return filename, desc, new_path


class CreatePointMutationDialog(QDialog, Ui_CreatePointMutation):
    """
    A dialog for creating a modification describing a point mutation.

    :ivar str mutation_name: The name of the mutation
    :ivar Formula mutation_formula: The formula corresponding to the mutation

    .. automethod:: __init__
    """

    def __init__(self, parent=None):
        """
        Initialize the dialog.

        :param QWidget parent: parent widget
        """

        # initialize the GUI
        super().__init__(parent)
        self.setupUi(self)

        self.mutation_name = None
        self.mutation_formula = None

        self.leOldResidue.textChanged.connect(self._update_info)
        self.leNewResidue.textChanged.connect(self._update_info)
        self._update_info()

    def _update_info(self):
        """
        Calculate the mutation name/formula and update the info text
        whenever the contents of the line edits change.

        :return: nothing
        """

        # clear mutation data and info text
        self.lbMutation.setText("")
        self.lbFormula.setText("")
        self.lbMass.setText("")
        self.mutation_name = None
        self.mutation_formula = None

        # get formulas of residues
        try:
            old_res = self.leOldResidue.text().upper()
            new_res = self.leNewResidue.text().upper()
            old_formula = Formula(amino_acid_compositions[old_res])
            new_formula = Formula(amino_acid_compositions[new_res])
        except KeyError:
            return

        if old_res == new_res:
            return

        self.mutation_name = old_res + " → " + new_res
        self.mutation_formula = new_formula - old_formula

        # update info text
        self.lbMutation.setText(
            amino_acid_names[old_res][1]
            + " → "
            + amino_acid_names[new_res][1])
        self.lbFormula.setText(str(self.mutation_formula))
        self.lbMass.setText(dec_places().format(self.mutation_formula.mass))

    @staticmethod
    def get_mutation(parent=None):
        """
        Opens an Create point mutation dialog
        and returns the created mutation.

        :param QWidget parent: parent widget
        :return: the name and formula of the mutation
        :rtype: tuple(str, Formula)
        """

        dialog = CreatePointMutationDialog(parent)
        result = dialog.exec_()
        if result == QDialog.Accepted:
            return dialog.mutation_name, dialog.mutation_formula


def extract_positions(pos_list):
    """
    Extract positions from a string like "1, 3-5, 9".

    :param str pos_list: list of single positions and ranges
    :return: list of individual positions (e.g., (1, 3, 4, 5, 9))
    :rtype: list(int)
    :raises ValueError: if parsing fails
    """

    re_range = re.compile("(?P<start>\d+)-(?P<stop>\d+)")
    pos = []
    for s in pos_list.split(","):
        s = s.strip()
        if not s:
            continue
        match = re_range.match(s)
        if match:
            g = match.groupdict()
            pos += [i for i in range(int(g["start"]), int(g["stop"]) + 1)]
        else:
            pos.append(int(s))  # may raise ValueError
    return pos


class CreateTruncationDialog(QDialog, Ui_CreateTruncation):
    """
    A dialog for creating an N-or C-terminal truncation.

    :ivar list chains: chain names
    :ivar str sequences: chain sequences
    :ivar list modifications: data for the table of modifications
    :ivar list structures: data for the table of structures
    :ivar bool stage2: indicates whether data for stage 2 should be generated

    .. automethod:: __init__
    """

    def __init__(self, parent=None, main_sequence=None, stage2=True):
        """
        Initialize the dialog.

        :param QWidget parent: parent widget
        :param str main_sequence: sequence to be loaded
                                  if the user clicks 'From parameters'
        :param bool stage2: True if truncations should be searched in stage 2;
                            otherwise False
        """

        # initialize the GUI
        super().__init__(parent)
        self.setupUi(self)

        # signal-slot connections
        self.btFromParameters.clicked.connect(self.sequence_from_parameters)
        self.btPreview.clicked.connect(self.preview)
        self.rbResidueList.clicked.connect(self.switch_residue_input)
        self.rbResidueNth.clicked.connect(self.switch_residue_input)

        # initialize instance variables
        self.chains, self.sequences = read_fasta_string(main_sequence,
                                                        join_sequences=False)
        self.modifications = None
        self.structures = None
        self.stage2 = stage2

        # fill the chain combobox
        if self.sequences:
            self.cbChain.addItem("(all)")
            for chain in self.chains:
                self.cbChain.addItem(chain)
        else:
            self.btFromParameters.setEnabled(False)
            self.lbChain.setEnabled(False)
            self.cbChain.setEnabled(False)

    def sequence_from_parameters(self):
        """
        Fill the sequence text edit with (parts of) the sequence
        of the main window.

        :return: nothing
        """

        index = self.cbChain.currentIndex() - 1
        if index == -1:
            self.teSequence.setText("".join(self.sequences))
        else:
            self.teSequence.setText(self.sequences[index])

    def switch_residue_input(self):
        """
        Enable the appropriate entry fields for truncation positions.
        :return: nothing
        """

        if self.rbResidueNth.isChecked():
            self.sbResidueNth.setEnabled(True)
            self.leResidueList.setEnabled(False)
        else:
            self.sbResidueNth.setEnabled(False)
            self.leResidueList.setEnabled(True)

    def get_subsequences(self):
        """
        Calculate all subsequences from the current parameters

        :return: a list of subsequences
        :rtype: list(str)
        """
        sequence = self.teSequence.toPlainText()
        if not sequence:
            return []

        # get all truncation positions
        if self.rbResidueList.isChecked():
            try:
                raw_pos = extract_positions(self.leResidueList.text())
            except ValueError:
                return []
            if self.rbNterminus.isChecked():
                pos = [0] + raw_pos
            else:
                pos = [-i for i in raw_pos[::-1]] + [len(sequence)]
        else:
            step = self.sbResidueNth.value()
            if self.rbNterminus.isChecked():
                pos = range(0, len(sequence) + 1, step)
            else:
                pos = range(len(sequence), -1, -step)

        # create all subsequences
        if self.rbNterminus.isChecked():
            return [sequence[i:] for i in pos]
        else:
            return [sequence[:i] for i in pos]

    def get_modification_data(self, subsequences):
        """
        Calculate data for modifications and structures
        from a given list of subsequences.

        :param list(str) subsequences: list of subsequences
        :return: two lists of dicts containing data for creating rows in the
                 table of modifications and table of structures, respectively
        :rtype: tuple(list(dict), list(dict))
        """

        sequence = self.teSequence.toPlainText()
        if self.rbNterminus.isChecked():
            prefix = "N_"
        else:
            prefix = "C_"

        if self.stage2:
            # create data for the required stage 1 modifications
            modifications = []
            for aa, count in Counter(sequence).items():
                modifications.append(dict(
                    active=True,
                    name=prefix + amino_acid_names[aa][1],
                    composition=str(Formula(amino_acid_compositions[aa])),
                    min_count=0,
                    max_count=count))

            # create data for the corresponding stage 2 structures
            structures = []
            for s in subsequences:
                composition = []
                for aa, count in Counter(s).items():
                    composition.append("{} {}{}".format(
                        count,
                        prefix,
                        amino_acid_names[aa][1]))
                structures.append(dict(
                    active=True,
                    name=prefix + s,
                    composition=", ".join(composition),
                    sites=prefix + "trunc",
                    abundance=0.0))
        else:
            # create data for the required stage 1 modifications
            modifications = []
            for s in subsequences:
                if not s:
                    continue
                formula = Formula()
                for aa, count in Counter(s).items():
                    formula += Formula(amino_acid_compositions[aa]) * count
                modifications.append(dict(
                    active=True,
                    name=prefix + s,
                    composition=str(formula),
                    min_count=0,
                    max_count=1))

            # stage 2 structures are not required
            structures = []

        return modifications, structures

    def preview(self):
        """
        Fill the preview text edit with the generated truncations.

        :return: nothing
        """

        self.tePreview.clear()
        subsequences = self.get_subsequences()
        if subsequences:
            preview_text = []
            for s in subsequences:
                if not s:
                    continue
                formula = Formula()
                for aa, count in Counter(s).items():
                    try:
                        formula += Formula(amino_acid_compositions[aa]) * count
                    except KeyError:
                        formula = Formula()
                        break
                preview_text.append(
                    "<b>{}</b>: {} <font color='#808080'>({})</font>".format(
                        s,
                        str(formula),
                        dec_places().format(formula.mass)))
            self.tePreview.setText("<br/>".join(preview_text))

    def done(self, r):
        """
        Check whether there are any truncations when Ok is clicked.
        Ask if the user really wants to quit if there are none.

        :param int r: return code
        :return: nothing
        """

        if r == QDialog.Accepted:
            self.modifications, self.structures = self.get_modification_data(
                self.get_subsequences())
            if self.modifications:
                super().done(r)
            else:
                answer = QMessageBox.question(
                    self,
                    "MoFi",
                    "The current parameters do not create any truncations.\n"
                    "Quit nevertheless?",
                    QMessageBox.No | QMessageBox.Yes,
                    QMessageBox.No)
                if answer == QMessageBox.Yes:
                    super().done(QDialog.Rejected)
                else:
                    pass
        else:
            super().done(r)

    @staticmethod
    def get_truncation(parent=None, main_sequence=None, stage2=True):
        """
        Opens an Create truncation dialog
        and returns the created truncation.

        :param QWidget parent: parent widget
        :param str main_sequence: sequence to be loaded
                                  if the user clicks 'From parameters'
        :param bool stage2: True if truncations should be searched in stage 2;
                            otherwise False
        :return: two lists of dicts containing data for creating rows in the
                 table of modifications and table of structures, respectively
        :rtype: tuple(list(dict)), list(dict)) or None
        """

        dialog = CreateTruncationDialog(parent, main_sequence, stage2)
        result = dialog.exec_()
        if result == QDialog.Accepted:
            return dialog.modifications, dialog.structures
