#! python

# Important shell commands:
# Converting the GUI:    pyuic5 ModFinder_UI.ui > ModFinder_UI.py
# Converting resources:  pyrcc5 mofi.qrc > mofi_rc.py

# Freezing to .exe:
# change self._path in config.ini to something useful ('C:/'?).
# Finalize directory and copy to new freeze directory (versioned)
# cd Freeze
# Change Config.ini location in configure.py!!!
# Make sure to enable/disable logging in ModFinder.py
# pyinstaller -w ModFinder.py -i mofi.ico
# custom hook-numpy.py in Lib\site-packages\PyInstaller\hooks for importing mkl
# binaries. Anaconda installs numpy slightly differently ... *sigh*.
# Move config directory to the correct path
# Move docs directory to correct path
# If I was motivated, I'd write a freeze wrapper that does those things
# automatically, but I can't be bothered anymore. ;)

import math
import os
import re
import sys
import time
import webbrowser
import xml.etree.ElementTree as ETree

from PyQt5.QtWidgets import (QApplication, QMainWindow, QMenu, QActionGroup,
                             QTableWidgetItem, QCheckBox, QMessageBox,
                             QTreeWidgetItem, QHeaderView, QButtonGroup,
                             QSpinBox, QDoubleSpinBox, QWidget, QHBoxLayout,
                             QAction, QProgressBar, QLabel, QSizePolicy,
                             QFileDialog)
from PyQt5.QtGui import QColor, QBrush
from PyQt5.QtCore import Qt

import numpy as np
import pandas as pd

import matplotlib
import matplotlib.text
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.widgets import SpanSelector, RectangleSelector
from matplotlib.figure import Figure

from mofi import configure
from mofi import mass_tools
from mofi import modification_search
from mofi import sequence_tools
from mofi.paths import data_dir
from mofi.modfinder_ui import Ui_ModFinder

pd.set_option('display.max_rows', 5000)
matplotlib.use("Qt5Agg")

# default values for columns in the monomer/polymer table
_monomer_table_columns = [
    ("Checked", False),
    ("Name", None),
    ("Composition", None),
    ("Min", 0),
    ("Max", -1)
]

_polymer_table_columns = [
    ("Checked", True),
    ("Name", None),
    ("Composition", None),
    ("Sites", ""),
    ("Abundance", 0.0)
]


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


class CollapsingRectangleSelector(RectangleSelector):
    """
    Select a rectangular region of an axes.
    The rectangle collapses to a line if a dimension
    is less than minspanx and minspany.
    """

    def __init__(self, *args, collapsex=0, collapsey=0, **kwargs):
        """
        Introduces the keys "collapsex" and "collapsey".

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


class SortableTreeWidgetItem(QTreeWidgetItem):
    """
    A QTreeWidget which supports numerical sorting.
    Additionally, the attribute mass_index stores the mass index
    of a top level widget.
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        self.mass_index = None

    def __lt__(self, other):
        column = self.treeWidget().sortColumn()
        key1 = self.text(column)
        key2 = other.text(column)
        return self.natural_sort_key(key1) < self.natural_sort_key(key2)

    @staticmethod
    def natural_sort_key(key):
        regex = "(\d*\.\d+|\d+)"
        parts = re.split(regex, key)
        return ((e if i % 2 == 0 else float(e)) for i, e in enumerate(parts))


class MainWindow(QMainWindow, Ui_ModFinder):
    def __init__(self, parent=None):

        # initialize the GUI
        super(MainWindow, self).__init__(parent)
        self.setupUi(self)

        # connect signals to slots
        self.acAbout.triggered.connect(self.show_about)
        self.acHelp.triggered.connect(self.show_help)
        self.acLoadSettings.triggered.connect(self.load_settings)
        self.acOpenFasta.triggered.connect(self.load_fasta_file)
        self.acOpenPeaks.triggered.connect(self.load_mass_file)
        self.acQuit.triggered.connect(QApplication.instance().quit)
        self.acSaveSettings.triggered.connect(self.save_settings)

        self.btClearMonomers.clicked.connect(
            lambda: self.table_clear(self.tbMonomers))
        self.btClearPolymers.clicked.connect(
            lambda: self.table_clear(self.tbPolymers))
        self.btDeleteRowMonomers.clicked.connect(
            lambda: self.table_delete_row(self.tbMonomers))
        self.btDeleteRowPolymers.clicked.connect(
            lambda: self.table_delete_row(self.tbPolymers))
        self.btSaveResults.clicked.connect(self.save_search_results)
        self.btFindModifications.clicked.connect(self.sample_modifications)
        self.btInsertRowAboveMonomers.clicked.connect(
            lambda: self.table_insert_row(self.tbMonomers, above=True))
        self.btInsertRowAbovePolymers.clicked.connect(
            lambda: self.table_insert_row(self.tbPolymers, above=True))
        self.btInsertRowBelowMonomers.clicked.connect(
            lambda: self.table_insert_row(self.tbMonomers, above=False))
        self.btInsertRowBelowPolymers.clicked.connect(
            lambda: self.table_insert_row(self.tbPolymers, above=False))
        self.btLabelPeaks.clicked.connect(lambda: self.show_results())
        self.btLoadMonomers.clicked.connect(
            lambda: self.load_table(
                default=False,
                dialog_title="Import modifications",
                table_widget=self.tbMonomers,
                cols=_monomer_table_columns
            )
        )
        self.btLoadPolymers.clicked.connect(
            lambda: self.load_table(
                default=False,
                dialog_title="Import glycans",
                table_widget=self.tbPolymers,
                cols=_polymer_table_columns
            )
        )
        self.btResetZoom.clicked.connect(self.reset_zoom)
        self.btSaveMonomers.clicked.connect(
            lambda: self.save_table("Export modifications", "monomers"))
        self.btSavePolymers.clicked.connect(
            lambda: self.save_table("Export glycans", "polymers"))
        self.btSaveSpectrum.clicked.connect(self.save_spectrum)
        self.btUpdateMass.clicked.connect(self.calculate_protein_mass)

        self.cbTolerance.activated.connect(self.choose_tolerance_units)

        self.chDelta1.clicked.connect(self.toggle_delta_series)
        self.chDelta2.clicked.connect(self.toggle_delta_series)
        self.chFilterStructureHits.clicked.connect(self.filter_structure_hits)
        self.chPngase.clicked.connect(self.calculate_protein_mass)

        self.lwPeaks.itemSelectionChanged.connect(self.select_peaks_in_list)

        self.sbDeltaRepetition1.valueChanged.connect(self.show_results)
        self.sbDeltaRepetition2.valueChanged.connect(self.show_results)
        self.sbDeltaTolerance1.valueChanged.connect(self.show_results)
        self.sbDeltaTolerance2.valueChanged.connect(self.show_results)
        self.sbDeltaValue1.valueChanged.connect(self.show_results)
        self.sbDeltaValue2.valueChanged.connect(self.show_results)
        self.sbDisulfides.valueChanged.connect(self.calculate_protein_mass)

        self.tbMonomers.cellChanged.connect(self.calculate_mod_mass)

        self.teSequence.textChanged.connect(
            lambda: self.teSequence.setStyleSheet(
                "QTextEdit { background-color: rgb(255, 225, 225) }"))

        self.twResults.itemClicked.connect(self.select_peaks_in_tree)

        # generate mass set selectors from config file
        self.agSelectMassSet = QActionGroup(self.menuAtomicMasses)
        for set_id, mass_set in enumerate(configure.mass_sets):
            ac_select_set = QAction(self)
            ac_select_set.setCheckable(True)
            if set_id == 0:
                ac_select_set.setChecked(True)
            ac_select_set.setObjectName("acSelectMassSet{:d}".format(set_id))
            ac_select_set.setText(mass_set)
            ac_select_set.setToolTip(
                configure.mass_sets[mass_set].get("description", ""))
            self.menuAtomicMasses.addAction(ac_select_set)
            # self.tlMassSet.addAction(ac_select_set)
            self.agSelectMassSet.addAction(ac_select_set)
        # noinspection PyUnresolvedReferences
        self.agSelectMassSet.triggered.connect(self.choose_mass_set)

        # status bar
        sb_widget = QWidget(self)
        sb_layout = QHBoxLayout(sb_widget)
        sb_layout.setContentsMargins(0, 0, 0, 0)
        sb_widget.setLayout(sb_layout)
        sb_layout.addWidget(QLabel("Progress:"))
        self.pbSearchProgress = QProgressBar(self)
        self.pbSearchProgress.setValue(0)
        self.pbSearchProgress.setMinimumWidth(300)
        self.pbSearchProgress.setSizePolicy(
            QSizePolicy.Minimum, QSizePolicy.Fixed)
        sb_layout.addWidget(self.pbSearchProgress)
        self.statusbar.addPermanentWidget(sb_widget)

        # mass spectrum plot
        self.spectrum_fig = Figure(dpi=100, frameon=False,
                                   tight_layout={"pad": 0}, edgecolor="white")
        self.spectrum_canvas = FigureCanvas(self.spectrum_fig)
        self.spectrumViewLayout.addWidget(self.spectrum_canvas)
        # single Axes of the figure
        self.spectrum_axes = None
        # LineCollection representing the peaks in the spectrum
        self.spectrum_peak_lines = None
        # SpanSelector to select multiple peaks
        self.spectrum_span_selector = None
        # CollapsingRectangleSelector to zoom the spectrum
        self.spectrum_rectangle_selector = None
        # original limits of the x axis when the plot is drawn
        self.spectrum_x_limits = None
        # the single peak that was picked, or a representative singe peak
        self.spectrum_picked_peak = None

        # button group for specrum interaction mode
        self.bgSpectrum = QButtonGroup()
        self.bgSpectrum.addButton(self.btModeDelta)
        self.bgSpectrum.addButton(self.btModeSelection)
        # noinspection PyUnresolvedReferences
        self.bgSpectrum.buttonClicked.connect(self.toggle_spectrum_mode)

        # monomer table and associated buttons
        self.tbMonomers.horizontalHeader().setSectionResizeMode(
            1, QHeaderView.Stretch)
        for col, width in [(0, 40), (2, 130), (3, 45), (4, 45)]:
            self.tbMonomers.setColumnWidth(col, width)
        self.tbMonomers.verticalHeader().setSectionResizeMode(
            QHeaderView.Fixed)
        self.tbMonomers.verticalHeader().setDefaultSectionSize(22)
        self.tbMonomers.create_row = self._monomer_table_create_row

        menu = QMenu()
        for files in os.listdir(os.path.join(data_dir, "modifications")):
            library = os.path.splitext(files)[0]
            menu.addAction(
                library,
                lambda: self.load_table(
                    default=True,
                    subdir="modifications",
                    table_widget=self.tbMonomers,
                    cols=_monomer_table_columns
                )
            )
        self.btDefaultModsMonomers.setMenu(menu)

        # polymer table and associated buttons
        self.tbPolymers.horizontalHeader().setSectionResizeMode(
            1, QHeaderView.Stretch)
        for col, width in [(0, 40), (2, 130), (3, 80), (4, 50)]:
            self.tbPolymers.setColumnWidth(col, width)
        self.tbPolymers.verticalHeader().setSectionResizeMode(
            QHeaderView.Fixed)
        self.tbPolymers.verticalHeader().setDefaultSectionSize(22)
        self.tbPolymers.create_row = self._polymer_table_create_row

        menu = QMenu()
        for files in os.listdir(os.path.join(data_dir, "glycans")):
            library = os.path.splitext(files)[0]
            menu.addAction(
                library,
                lambda: self.load_table(
                    default=True,
                    subdir="glycans",
                    table_widget=self.tbPolymers,
                    cols=_polymer_table_columns
                )
            )
        self.btDefaultModsPolymers.setMenu(menu)

        # private members
        self._monomer_hits = None  # results from the monomer search
        self._exp_mass_data = None  # peak list (mass + relative abundance)
        self._known_mods_mass = 0  # mass of known modification
        self._path = configure.path  # last path selected in a file dialog
        self._polymer_hits = None  # results from the polymer search
        self._protein_mass = 0  # mass of the current Protein object


    def _monomer_table_create_row(self, row_id, active=False, name="",
                                  composition="", min_count=0, max_count=-1):
        """
        Create a new row in the table of modifications.
        This row describes a modification with the given parameters.

        Private function, to be called by general row insertion functions and
        loaders for modification libraries.

        :param row_id: Row index passed to QTableWidget.insertRow()
        :param active: true if the monomer should be used in search stage 1
        :param name: name of the modification (str)
        :param composition: composition or mass of the modification (str)
        :param min_count: minimum number of modifications (int)
        :param max_count: maximum number of modifications (int)
        :return: nothing
        """

        self.tbMonomers.blockSignals(True)
        self.tbMonomers.insertRow(row_id)

        active_checkbox = QCheckBox()
        active_checkbox.setChecked(active)
        # noinspection PyUnresolvedReferences
        active_checkbox.clicked.connect(self.calculate_mod_mass)
        active_ch_widget = QWidget()
        active_ch_layout = QHBoxLayout(active_ch_widget)
        active_ch_layout.addWidget(active_checkbox)
        active_ch_layout.setAlignment(Qt.AlignCenter)
        active_ch_layout.setContentsMargins(0, 0, 0, 0)
        self.tbMonomers.setCellWidget(row_id, 0, active_ch_widget)

        self.tbMonomers.setItem(row_id, 1, QTableWidgetItem(name))

        self.tbMonomers.setItem(row_id, 2, QTableWidgetItem(composition))

        min_spinbox = QSpinBox()
        min_spinbox.setMinimum(0)
        min_spinbox.setFrame(False)
        min_spinbox.setValue(min_count)
        min_spinbox.setStyleSheet(configure.spin_box_flat_style)
        # noinspection PyUnresolvedReferences
        min_spinbox.valueChanged.connect(self.calculate_mod_mass)
        self.tbMonomers.setCellWidget(row_id, 3, min_spinbox)

        max_spinbox = QSpinBox()
        max_spinbox.setMinimum(-1)
        max_spinbox.setSpecialValueText("inf")
        max_spinbox.setFrame(False)
        max_spinbox.setValue(max_count)
        max_spinbox.setStyleSheet(configure.spin_box_flat_style)
        self.tbMonomers.setCellWidget(row_id, 4, max_spinbox)

        self.tbMonomers.blockSignals(False)


    def _polymer_table_create_row(self, row_id, active=True, name="",
                                  composition="", sites="", abundance=0.0):
        """
        Create a new row in the table of modifications.
        This row describes a modification with the given parameters.

        Private function, to be called by general row insertion functions and
        loaders for modification libraries.

        :param row_id: Row index passed to QTableWidget.insertRow()
        :param active: true if the polymer should be used in ssearch stage 1
        :param name: name of the modification (str)
        :param composition: composition or mass of the modification (str)
        :param sites: identifier for the modification site;
                      empty means any site (str)
        :param abundance: relative abundance (float)
        :return: nothing
        """

        self.tbPolymers.blockSignals(True)
        self.tbPolymers.insertRow(row_id)

        active_checkbox = QCheckBox()
        active_checkbox.setChecked(active)
        active_ch_widget = QWidget()
        active_ch_layout = QHBoxLayout(active_ch_widget)
        active_ch_layout.addWidget(active_checkbox)
        active_ch_layout.setAlignment(Qt.AlignCenter)
        active_ch_layout.setContentsMargins(0, 0, 0, 0)
        self.tbPolymers.setCellWidget(row_id, 0, active_ch_widget)

        self.tbPolymers.setItem(row_id, 1, QTableWidgetItem(name))

        self.tbPolymers.setItem(row_id, 2, QTableWidgetItem(composition))

        self.tbPolymers.setItem(row_id, 3, QTableWidgetItem(sites))

        abundance_spinbox = QDoubleSpinBox()
        abundance_spinbox.setMinimum(0)
        abundance_spinbox.setMaximum(100)
        abundance_spinbox.setSingleStep(.1)
        abundance_spinbox.setFrame(False)
        abundance_spinbox.setValue(abundance)
        abundance_spinbox.setStyleSheet(configure.double_spin_box_flat_style)
        self.tbPolymers.setCellWidget(row_id, 4, abundance_spinbox)

        self.tbPolymers.blockSignals(False)


    @staticmethod
    def table_insert_row(table_widget, above=True):
        """
        Insert a row into the table of modifications.
        The row will be inserted relative to the current selection
        (if one exists) or to all rows otherwise.

        :param table_widget: the QTableWidget to modify
        :param above: True if the row should be inserted
                      above the current selection
        :return: nothing
        """

        selected_rows = table_widget.selectionModel().selectedRows()
        if selected_rows:
            if above:
                last_row = selected_rows[0].row()
            else:
                last_row = selected_rows[-1].row() + 1
        else:
            if above:
                last_row = 0
            else:
                last_row = table_widget.rowCount()

        table_widget.create_row(last_row)


    @staticmethod
    def table_clear(table_widget):
        """
        Delete all rows in a QTableWidget.

        :param table_widget: the QTableWidget to modify
        :return: nothing
        """

        table_widget.clearContents()
        table_widget.setRowCount(0)


    @staticmethod
    def table_delete_row(table_widget):
        """
        Delete selected rows in the table of modifications.

        :param table_widget: the QTableWidget to modify
        :return: nothing
        """
        if table_widget.selectionModel().selectedRows():
            for i in table_widget.selectionModel().selectedRows()[::-1]:
                table_widget.removeRow(i.row())


    def table_to_df(self, which="monomers"):
        """
        Create a dataframe from the contents of the monomer or polymer table.

        :param which: table to convert, "monomers" or "polymers"
        :return: dataframe
        """
        if which == "monomers":
            columns = ["Checked", "Name", "Composition", "Mass", "Min", "Max"]
            dtypes = ["bool", "str", "str", "float64", "int64", "int64"]
            df = pd.DataFrame(self.calculate_mod_mass(), columns=columns)
        elif which == "polymers":
            columns = ["Checked", "Name", "Composition", "Sites", "Abundance"]
            dtypes = ["bool", "str", "str", "str", "float64"]
            df = pd.DataFrame(self.get_polymers(), columns=columns)

        else:
            raise ValueError("Invalid value for 'which': {}".format(which))

        return df.astype(dict(zip(columns, dtypes)))


    def save_table(self, dialog_title=None, which="monomers"):
        """
        Export the contents of the monomer or polymer table.

        :param dialog_title: title of the QFileDialog
        :param which: table to export; "monomers" or "polymers"
        :return: nothing
        """
        filename = QFileDialog.getSaveFileName(
            self,
            dialog_title,
            self._path,
            "Comma-separated value (*.csv)")[0]
        self._path = os.path.split(filename)[0]
        if filename:
            if not filename.endswith(".csv"):
                filename += ".csv"
            try:
                df = self.table_to_df(which=which)  # type: pd.DataFrame
                df.to_csv(filename, index=False)
            except OSError:
                QMessageBox.critical(
                    self,
                    "Error",
                    "Error when writing to " + filename + OSError.args)


    def table_from_df(self, df, table_widget=None, cols=None):
        """
        Read monomers/polymers from a dataframe.

        :param df: a dataframe
        :param table_widget: QTableWidget to fill with values
        :param cols: list of (column header, default value) tuples,
                     sorted according to the order of the arguments to
                     _[monomer/polymer]_table_create_row.
                     If the default value is None, a column must
                     exist in the input file; otherwise, it will be filled
                     with the default value if missing in the input file.
        :return: nothing
        """

        for label, default in cols:
            if label not in df.columns:
                if default is None:
                    QMessageBox.warning(
                        self,
                        "Warning",
                        "Column '{}' missing in input. ".format(label)
                        + "No data imported.")
                    return
                else:
                    df[label] = default

        self.table_clear(table_widget)
        for row_id, data in df.iterrows():
            table_widget.create_row(row_id, *[data[c[0]] for c in cols])
        self.calculate_mod_mass()


    def load_table(self, default=False, subdir=None, dialog_title=None,
                   table_widget=None, cols=None):
        """
        Import the contents of the monomer/polymer table.

        :param default: true if a default monomer library should be loaded
                        false if the user should choose a monomer library file
        :param subdir: directory in config containing the default libraries
        :param dialog_title: title of the QFileDialog
        :param table_widget: QTableWidget to fill with values
        :param cols: see cols in self.table_from_df
        :return: nothing
        """

        if default:
            filename = os.path.join(data_dir, subdir,
                                    self.sender().text() + ".csv")
            file_format = "csv"
        else:
            filename = QFileDialog.getOpenFileName(
                self,
                dialog_title,
                self._path,
                "Excel files (*.xlsx *.xls);; CSV files (*.csv *.txt)")[0]
            self._path = os.path.split(filename)[0]
            ext = os.path.splitext(filename)[1]

            if ext == "":
                return
            elif ext in [".xls", ".xlsx"]:
                file_format = "xls"
            elif ext in [".txt", ".csv"]:
                file_format = "csv"
            else:
                QMessageBox.warning(
                    self,
                    "Warning",
                    "Unknown file format: {}.".format(ext))
                return

        if file_format == "csv":
            df = pd.read_csv(filename)
        else:
            df = pd.read_excel(filename)
        self.table_from_df(df, table_widget, cols)


    def show_about(self):
        """
        Show the about dialog.

        :return: nothing
        """

        QMessageBox.about(
            self,
            "About ModFinder",
            "ModFinder version: {}\n{}\nContact: {}".format(
                configure.version,
                configure.rights,
                configure.contact))


    @staticmethod
    def show_help():
        """
        Open the manual in the default webbrowser.

        :return: nothing
        """

        webbrowser.open("docs/MoFi_manual.pdf")


    def choose_tolerance_units(self):
        """
        Adjust the settings of the tolerance spin box
        when PPM or Da are selected.

        :return: nothing
        """
        if self.cbTolerance.currentIndex() == 0:  # that is, Da.
            self.sbTolerance.setDecimals(2)
            self.sbTolerance.setMinimum(0.0)
            self.sbTolerance.setMaximum(50.0)
            self.sbTolerance.setSingleStep(0.1)
            self.sbTolerance.setValue(configure.default_da)
        else:
            self.sbTolerance.setDecimals(0)
            self.sbTolerance.setMinimum(0)
            self.sbTolerance.setMaximum(150)
            self.sbTolerance.setSingleStep(1)
            self.sbTolerance.setValue(configure.default_ppm)


    def load_fasta_file(self):
        """
        Opens a FASTA file and displays its contents in self.teSequence

        Changes:
            self._path to directory of selected file
            text of self.teSequence to contents of input file

        :return: nothing
        """

        filename = QFileDialog.getOpenFileName(
            self,
            "Open FASTA file",
            self._path,
            "Sequence files (*.fasta *.txt)")[0]
        self._path = os.path.split(filename)[0]
        if filename:
            try:
                with open(filename) as f:
                    fasta_input = f.read()
            except OSError:
                QMessageBox.critical(self,
                                     "Error",
                                     "Could not open sequence file.")
                return
            self.teSequence.setText(fasta_input)
            self.calculate_protein_mass()


    def fill_peak_list(self, masses):
        """
        Fill the peak list (self.lwPeaks).

        :param masses: an iterable of numbers representing masses
        :return: nothing
        """
        self.lwPeaks.blockSignals(True)
        self.lwPeaks.clear()
        self.lwPeaks.addItems(["{:.2f}".format(i) for i in masses])
        self.lwPeaks.setCurrentRow(0)
        self.lwPeaks.blockSignals(False)


    def load_mass_file(self):
        """
        Open a mass list and display its contents.

        :return: nothing
        """

        filename = QFileDialog.getOpenFileName(
            self,
            "Open mass list",
            self._path,
            "Excel files (*.xlsx *.xls);; CSV files (*.csv *.txt)")[0]
        self._path = os.path.split(filename)[0]
        if filename:
            mass_data = mass_tools.read_massfile(filename)
            if mass_data is None:
                QMessageBox.critical(self,
                                     "Error",
                                     "Could not load mass file.")
                return

            self._exp_mass_data = mass_data
            self.twResults.clear()
            self._monomer_hits = None
            self._polymer_hits = None
            self.chFilterStructureHits.setEnabled(False)
            self.fill_peak_list(self._exp_mass_data["Average Mass"])
            self.draw_spectrum()


    def save_search_results(self):
        """
        Write the dearch results to a CSV file.

        :return: nothing
        """
        if self._monomer_hits is None:
            return

        filename = QFileDialog.getSaveFileName(
            self,
            "Save results",
            self._path,
            "Comma-separated value (*.csv)")[0]
        self._path = os.path.split(filename)[0]

        if not filename:
            return
        if not filename.endswith(".csv"):
            filename += ".csv"

        try:
            with open(filename, "w") as f:
                f.write("# Combinatorial search results by ModFinder\n")
                f.write("# Date: " + time.strftime("%c") + "\n")
                f.write("#   Tolerance: {:.2f} {}\n".format(
                    self.sbTolerance.value(),
                    self.cbTolerance.currentText()))

                f.write("#   Composition:\n")
                for (is_checked, name, _, mass,
                     min_count, max_count) in self.calculate_mod_mass():
                    if is_checked:
                        if max_count == -1:
                            max_count = "inf"
                        f.write("#     {} ({:.2f} Da), ".format(name, mass))
                        f.write("min {}, ".format(min_count))
                        f.write("max {}\n".format(max_count))

                f.write("#   Structures:\n")
                for (is_checked, name, composition,
                     sites, abundance) in self.get_polymers():
                    if is_checked:
                        f.write("#     {} ({}); ".format(name, composition))
                        f.write("sites: {}; ".format(sites))
                        f.write("abundance: {:.2f}\n".format(abundance))

                if self.chFilterStructureHits.isChecked():
                    f.write("# Results from structure search (stage 2).\n")
                    df_hits = self._polymer_hits
                else:
                    f.write("# Results from composition search (stage 1).\n")
                    df_hits = self._monomer_hits

                # add a column "Relative abundance" to the hits dataframe
                df_relative_abundance = (
                    self._exp_mass_data["Relative Abundance"].to_frame())
                df_relative_abundance.index.names = ["Mass_ID"]
                df_hits = df_relative_abundance.join(df_hits, how="inner")
                df_hits.to_csv(f)
        except OSError:
            QMessageBox.critical(
                self,
                "Error",
                "Error while writing to {}. No output saved.".format(filename))


    def choose_mass_set(self):
        """
        Chooses the set of atomic weights to use for mass calculations

        :return: nothing
        """
        configure.select_mass_set(self.agSelectMassSet.checkedAction().text())
        if self.teSequence.toPlainText():
            self.calculate_protein_mass()
            self.calculate_mod_mass()


    def calculate_protein_mass(self):
        """
        Calculates the mass of the protein from the current data.

        Changes self._protein_mass
        Updates the value of self.lbMassProtein, self.lbMassMods
                and self.lbMassTotal

        :return: nothing
        """

        protein_sequence = self.teSequence.toPlainText()
        chains, sequence = sequence_tools.read_fasta_string(protein_sequence)
        try:
            protein = sequence_tools.Protein(sequence,
                                             chains,
                                             self.sbDisulfides.value(),
                                             self.chPngase.isChecked())
        except KeyError as e:
            QMessageBox.critical(
                self,
                "Error",
                "Error when parsing sequence: "
                + "{} is not a valid symbol".format(e.args[0]))
            return

        self.sbDisulfides.setEnabled(True)
        self.chPngase.setEnabled(True)
        self.sbDisulfides.setMaximum(
            protein.amino_acid_composition["C"] / 2)
        self.teSequence.setStyleSheet(
            "QTextEdit { background-color: rgb(240, 251, 240) }")
        self._protein_mass = protein.mass
        self.lbMassProtein.setText("{:,.2f}".format(self._protein_mass))
        self.lbMassMods.setText("{:,.2f}".format(self._known_mods_mass))
        self.lbMassTotal.setText("{:,.2f}".format(self._protein_mass
                                                  + self._known_mods_mass))


    def calculate_mod_mass(self):
        """
        Calculate the mass of known modifications.

        Changes:
            self._known_mods_mass to the mass of known modifications
            updates the value of self.lbMassProtein, self.lbMassMods
            and self.lbMassTotal

        :return: list of (checked, name, composition,
                          mass,min count, max count) tuples
        """
        self._known_mods_mass = 0
        result = []

        # add min counts for monomers to the theoretical mass
        for row_id in range(self.tbMonomers.rowCount()):
            ch = self.tbMonomers.cellWidget(row_id, 0).findChild(QCheckBox)
            name = self.tbMonomers.item(row_id, 1).text()
            composition = self.tbMonomers.item(row_id, 2).text().strip()
            min_count = self.tbMonomers.cellWidget(row_id, 3).value()
            max_count = self.tbMonomers.cellWidget(row_id, 4).value()
            mass = 0
            error_in_formula = False
            try:  # composition could be a mass ...
                mass = float(composition)
            except ValueError:  # ... or a formula
                try:
                    formula = mass_tools.Formula(composition)
                except ValueError:
                    error_in_formula = True
                else:
                    mass = formula.mass
            if error_in_formula or mass == 0:
                self.tbMonomers.item(row_id, 2).setBackground(
                    QColor(255, 225, 225))
                self.tbMonomers.item(row_id, 2).setToolTip("")
            else:  # set the tooltip
                self.tbMonomers.item(row_id, 2).setBackground(
                    QColor(255, 255, 255))
                self.tbMonomers.item(row_id, 2).setToolTip(
                    "{:.2f} Da".format(mass))

            if ch.isChecked():
                self._known_mods_mass += mass * min_count
            result.append((ch.isChecked(), name, composition, mass,
                           min_count, max_count))

        self.lbMassProtein.setText("{:,.2f}".format(self._protein_mass))
        self.lbMassMods.setText("{:,.2f}".format(self._known_mods_mass))
        self.lbMassTotal.setText("{:,.2f}".format(self._protein_mass
                                                  + self._known_mods_mass))
        return result


    def get_polymers(self):
        """
        Read the current polymer library from the table and return as list.

        :return: list of (is_checked, name, compos., sites, abundance) tuples
        """

        result = []
        for row_id in range(self.tbPolymers.rowCount()):
            is_checked = (self.tbPolymers.cellWidget(row_id, 0)
                          .findChild(QCheckBox)
                          .isChecked())
            name = self.tbPolymers.item(row_id, 1).text()
            composition = self.tbPolymers.item(row_id, 2).text()
            sites = self.tbPolymers.item(row_id, 3).text()
            abundance = self.tbPolymers.cellWidget(row_id, 4).value()
            result.append((is_checked, name, composition, sites, abundance))
        return result


    def sample_modifications(self):
        """
        Prepare data for the two-stage search and process its results.

        :return: nothing
        """

        self.calculate_protein_mass()

        # calculate required input if a single mass was entered
        # (i.e., no peak list was loaded)
        if self.rbSingleMass.isChecked():
            self._exp_mass_data = pd.DataFrame(
                {"Average Mass": self.sbSingleMass.value(),
                 "Relative Abundance": 100.0},
                index=[0])
            self.fill_peak_list(self._exp_mass_data["Average Mass"])
            self.draw_spectrum()

        if self._exp_mass_data is None:
            QMessageBox.critical(
                self, "Error", "No mass list loaded. Aborting search.")
            return

        monomers = [(m[1], m[3], m[4], m[5])
                    for m in self.calculate_mod_mass()
                    if m[0]]
        modifications = []  # list of modifications for search stage 1
        explained_mass = self._protein_mass + self._known_mods_mass
        unknown_masses = (self._exp_mass_data["Average Mass"]
                          - explained_mass)  # type: pd.DataFrame

        if self.cbTolerance.currentIndex() == 0:  # that is, "Da."
            # calculate largest mass plus tolerance
            max_tol_mass = (max(self._exp_mass_data["Average Mass"])
                            + self.sbTolerance.value())
            mass_tolerance = self.sbTolerance.value()
        else:
            # calculate a mass tolerance for each peak
            # if we're working with ppm tolerance
            max_tol_mass = (max(self._exp_mass_data["Average Mass"])
                            * (1 + self.sbTolerance.value() / 1000000))
            mass_tolerance = []
            for _, m in self._exp_mass_data["Average Mass"].iteritems():
                mass_tolerance.append(m * self.sbTolerance.value() / 1000000)

        # prepare polymer combinations for search stage 2
        available_monomers = [m[0] for m in monomers]
        polymers = self.get_polymers()

        if polymers:  # list contains at least one entry
            df_polymers = pd.DataFrame.from_records(
                polymers,
                columns=["ch", "Name", "Composition", "Sites", "Abundance"])
            df_polymers = (df_polymers[df_polymers["ch"]]
                           .set_index("Name", drop=True))
            monomers_in_library = set(
                modification_search.get_monomers_from_library(df_polymers))
            monomers_for_polymer_search = [m for m in available_monomers
                                           if m in monomers_in_library]
            polymer_combs = modification_search.calc_polymer_combinations(
                df_polymers,
                monomers_for_polymer_search,
                self.pbSearchProgress)
        else:  # list remained empty, i.e., there are no polymers
            monomers_in_library = set()
            monomers_for_polymer_search = None
            polymer_combs = None

        # calculate the upper limit of glycans that may appear
        # add all checked single glycans to the list of modifications
        for name, mass, min_count, max_count in monomers:
            if max_count == -1:
                if (polymer_combs is not None
                        and name in polymer_combs.index.names):
                    max_count = max(polymer_combs.index.get_level_values(name))
                else:  # only search stage 1
                    max_count = min(
                        int((max_tol_mass - self._protein_mass) / mass),
                        configure.maxmods)
            modifications.append((name, mass, max_count - min_count))

        if not modifications:
            QMessageBox.critical(
                self, "Error", "List of monomers is empty. Aborting search.")
            return

        # check whether all monomers required in search stage 2 are available
        missing_monomers = monomers_in_library - set(available_monomers)
        if missing_monomers:
            error_message = [
                "The following monosaccharides appear in the library ",
                "but do not appear in the composition list: "]
            error_message += " ".join(sorted(missing_monomers))
            error_message.append(".\n")
            QMessageBox.critical(self,
                                 "Error",
                                 "".join(error_message))
            return

        self.statusbar.showMessage("Starting monomer search (stage 1) ...")
        print("Experimental Masses:",
              self._exp_mass_data["Average Mass"].head(),
              sep="\n")
        print("Explained mass (protein + known modifications):",
              explained_mass)
        print("Unknown masses searched:",
              unknown_masses.head(),
              sep="\n")
        print("Mass tolerance: {:f} {}".format(self.sbTolerance.value(),
                                               self.cbTolerance.currentText()))
        print("Modifications used in search state 1:\nName\tMass\tmax")
        for m in modifications:
            print("{}\t{:.2f}\t{:d}".format(*m))
        if monomers_for_polymer_search:
            print("Modifications used in search stage 2:")
            print(", ".join(monomers_for_polymer_search))

        # stage 1: monomer search
        self._monomer_hits = modification_search.find_monomers(
            modifications,
            list(unknown_masses),
            mass_tolerance=mass_tolerance,
            explained_mass=explained_mass,
            progress_bar=self.pbSearchProgress)

        self.statusbar.showMessage(
            "Monomer search done! Starting polymer search (stage 2) ...")

        if self._monomer_hits is None:
            # the modification search was not successful
            QMessageBox.critical(
                self, "Error", "Combinatorial search was unsuccessful.")
            return

        # add the minimum monomer counts to the result data frame
        for name, _, min_count, _ in monomers:
            self._monomer_hits[name] += min_count

        # stage 2: polymer search
        if polymers:
            self._polymer_hits = modification_search.find_polymers(
                self._monomer_hits,
                polymer_combinations=polymer_combs,
                monomers=monomers_for_polymer_search,
                progress_bar=self.pbSearchProgress)
            self.statusbar.showMessage("Polymer search DONE!", 5000)
        self.chFilterStructureHits.setEnabled(True)
        self.show_results()


    def draw_spectrum(self):
        """
        Show a mass spectrum after loading a peak list.

        :return: nothing
        """

        self.spectrum_fig.clear()
        self.spectrum_axes = self.spectrum_fig.add_subplot(111)
        self.spectrum_axes.set_xmargin(.02)
        self.spectrum_peak_lines = self.spectrum_axes.vlines(
            x=self._exp_mass_data["Average Mass"],
            ymin=0,
            ymax=self._exp_mass_data["Relative Abundance"],
            linewidth=1,
            picker=5)
        self.spectrum_axes.set_ylim(0, 110)
        self.spectrum_axes.set_xlabel("Mass (Da)")
        self.spectrum_axes.set_ylabel("Relative Abundance (%)")
        self.spectrum_axes.yaxis.set_ticks_position("left")
        self.spectrum_axes.xaxis.set_ticks_position("bottom")
        self.spectrum_axes.tick_params(direction="out")
        self.spectrum_x_limits = self.spectrum_axes.get_xlim()
        self.spectrum_canvas.draw()

        # set up the pick and span selectors.
        self.spectrum_fig.canvas.mpl_connect(
            "pick_event", self.select_peaks_by_pick)
        self.spectrum_span_selector = SpanSelector(
            self.spectrum_axes,
            self.select_peaks_by_span,
            "horizontal",
            minspan=10,  # in order not to interfere with the picker
            rectprops=dict(alpha=.1),
            button=1)  # only the left mouse button spans
        self.spectrum_rectangle_selector = CollapsingRectangleSelector(
            self.spectrum_axes,
            self.zoom_spectrum,
            collapsex=50,
            collapsey=5,
            button=3,
            rectprops=dict(
                facecolor=(1, 0, 0, .1),
                edgecolor="black",
                fill=True,
                linewidth=1.5))


    def zoom_spectrum(self, start, stop):
        """
        Set the new limits of the x and y axis
        after drawing the zoom rectangle.

        :param start: MouseEvent that describes the coordinates
                      where the mouse button was pressed.
        :param stop: MouseEvent that describes the coordinates
                     where the mouse button was released.
        :return: nothing
        """

        if math.isclose(start.ydata, stop.ydata):  # only zoom x axis
            self.spectrum_axes.set_xlim(start.xdata, stop.xdata)
        elif math.isclose(start.xdata, stop.xdata):  # only zoom y axis
            self.spectrum_axes.set_ylim(start.ydata, stop.ydata)
        else:  # zoom to rectangle
            self.spectrum_axes.set_xlim(start.xdata, stop.xdata)
            self.spectrum_axes.set_ylim(start.ydata, stop.ydata)


    def reset_zoom(self):
        """
        Restore the view of the complete spectrum.

        :return: nothing
        """
        try:
            self.spectrum_axes.set_xlim(*self.spectrum_x_limits)
            self.spectrum_axes.set_ylim(0, 110)
            self.spectrum_canvas.draw()
        except AttributeError:  # i.e., there is currently no spectrum
            pass


    def toggle_spectrum_mode(self):
        """
        Toggle between peak selection for monomer/polymer analysis
        and peak selection for mass difference analysis.

        :return: nothing
        """
        try:
            if self.bgSpectrum.checkedButton() == self.btModeSelection:
                self.gbDeltaSeries.setEnabled(False)
                self.spectrum_span_selector.active = True
            else:
                self.gbDeltaSeries.setEnabled(True)
                self.toggle_delta_series()
                self.spectrum_span_selector.active = False
        except AttributeError:  # i.e., there is currently no spectrum
            pass


    def toggle_delta_series(self):
        """
        Enable the delta series spin boxes
        whose corresponding checkbox is checked.

        :return: nothing
        """
        for ch, sb_delta, sb_tolerance, sb_repetitions in [
                (self.chDelta1, self.sbDeltaValue1,
                 self.sbDeltaTolerance1, self.sbDeltaRepetition1),
                (self.chDelta2, self.sbDeltaValue2,
                 self.sbDeltaTolerance2, self.sbDeltaRepetition2)]:
            if ch.isChecked():
                sb_delta.setEnabled(True)
                sb_tolerance.setEnabled(True)
                sb_repetitions.setEnabled(True)
            else:
                sb_delta.setEnabled(False)
                sb_tolerance.setEnabled(False)
                sb_repetitions.setEnabled(False)
        self.highlight_delta_series()

    def filter_structure_hits(self):
        """
        Called when the checkbox "Filter structure hits" is (un)checked.
        Ensure that only a single mass is displayed in the result tree
        if this box is unchecked.

        :return: nothing
        """
        if not self.chFilterStructureHits.isChecked():
            i = 0
            for i in range(self.lwPeaks.count()):
                if self.lwPeaks.item(i).isSelected():
                    break
            self.show_results([i])
        else:
            self.show_results()



    def select_peaks_in_list(self):
        """
        Select one or several peaks in the peak list.

        :return: nothing
        """
        self.spectrum_picked_peak = self.lwPeaks.currentRow()
        self.show_results()


    def select_peaks_in_tree(self):
        """
        Select the peak corresponding to an entry in the result tree.

        :return: nothing
        """

        toplevel_item = self.twResults.currentItem()
        while toplevel_item.parent():
            toplevel_item = toplevel_item.parent()
        self.spectrum_picked_peak = toplevel_item.mass_index
        self.show_results([toplevel_item.mass_index])


    def select_peaks_by_pick(self, event):
        """
        Select a peak picked by a mouseclick on the spectrum.

        :param event: PickEvent from the canvas
        :return: nothing
        """

        self.spectrum_picked_peak = event.ind[0]
        self.show_results(event.ind)


    def select_peaks_by_span(self, min_mass, max_mass):
        """
        Select all peaks that fall within a selected span.

        :param min_mass: lower end of the span
        :param max_mass: upper end of the span
        :return: nothing
        """

        peak_indices = []
        for i in range(self.lwPeaks.count()):
            if min_mass <= float(self.lwPeaks.item(i).text()) <= max_mass:
                peak_indices.append(i)
        self.chFilterStructureHits.blockSignals(True)
        self.chFilterStructureHits.setChecked(True)
        self.chFilterStructureHits.blockSignals(False)
        self.spectrum_picked_peak = peak_indices[0]
        self.show_results(peak_indices)


    @staticmethod
    def concat_interval_names(row):
        """
        Concatenate rows from a dataframe indicating interval names.
        This is an auxiliary function used by highlight_delta_series().
        Examples: ["", -1] -> "-1"
                  [2, -1] -> "2/-1"

        :param row: Row from dataframe.apply(axis=1)
        :return: a string
        """
        result = []
        for value in row:
            if value:
                result.append(str(value))
        return "/".join(result)


    def highlight_delta_series(self):
        """
        Highlights a series of peaks that differ by a given mass.

        :return: list of peak indices in the delta series
        """

        main_mass = float(self._exp_mass_data
                          .iloc[self.spectrum_picked_peak]["Average Mass"])
        min_mass = float(min(self._exp_mass_data["Average Mass"]))
        max_mass = float(max(self._exp_mass_data["Average Mass"]))

        # dataframe with the mass index as index,
        # the delta series index as columns,
        # and the intervall names als values
        df_delta_distances = pd.DataFrame(
            index=self._exp_mass_data.index,
            dtype=str)

        # dataframe with the following values:
        # 0 - peaks not in any delta mass series
        # 1 - peaks in series 1
        # 2 - peaks in series 2
        # 3 - peaks in both series
        # 4 - selected (central) peak
        df_delta_peaks = pd.DataFrame(
            index=self._exp_mass_data.index,
            dtype=int)

        # the algorithm is the same for both delta series; hence, use this loop
        for delta_id, ch_delta, sb_value, sb_tolerance, sb_repetitions in [
                (1, self.chDelta1, self.sbDeltaValue1,
                 self.sbDeltaTolerance1, self.sbDeltaRepetition1),
                (2, self.chDelta2, self.sbDeltaValue2,
                 self.sbDeltaTolerance2, self.sbDeltaRepetition2)]:
            if not ch_delta.isChecked():
                continue

            if sb_repetitions.value() == -1:
                max_iterations = int(sb_value.value()
                                     / sb_tolerance.value()
                                     / 2)
            else:
                max_iterations = sb_repetitions.value()
            intervals = {}  # a {number of differences: (start, end)} dict

            # calculate possible intervals for peaks separated
            # from the selected peaks by an integer number of mass differences;
            # for each added difference, increase the interval size
            # by the original tolerance times two
            current_mass = main_mass
            tolerance = sb_tolerance.value()
            i = 1
            while i <= max_iterations and current_mass > min_mass:
                current_mass -= sb_value.value()
                intervals[-i] = (current_mass - tolerance,
                                 current_mass + tolerance)
                tolerance += sb_tolerance.value()
                i += 1

            current_mass = main_mass
            tolerance = sb_tolerance.value()
            i = 1
            while i <= max_iterations and current_mass < max_mass:
                current_mass += sb_value.value()
                intervals[i] = (current_mass - tolerance,
                                current_mass + tolerance)
                tolerance += sb_tolerance.value()
                i += 1

            df_delta_distances[delta_id] = (
                self._exp_mass_data["Average Mass"]
                .apply(self.find_in_intervals, intervals=intervals))
            df_delta_distances[delta_id][self.spectrum_picked_peak] = "0"

            df_delta_peaks[delta_id] = (
                df_delta_distances[delta_id]
                .map(lambda x: delta_id if x else 0))

        # calculate "total" columns describing the results
        # from both delta series searches
        df_delta_peaks["total"] = df_delta_peaks.sum(axis=1).astype(int)
        df_delta_peaks["total"][self.spectrum_picked_peak] = 4
        if len(df_delta_distances.columns) == 0:
            df_delta_distances["total"] = ""
        else:
            df_delta_distances["total"] = (
                df_delta_distances
                .apply(self.concat_interval_names, axis=1))

        # color the peaks in the delta series and increase their line width
        color_set = np.array(["#b3b3b3",   # light gray
                              "#aa0088",   # violet
                              "#2ca05a",   # greenish
                              "#6b5071",   # mixture of the prevoius two
                              "#ff0000"])  # red

        linewidth_set = np.array([1, 2, 2, 2, 3])
        self.spectrum_peak_lines.set_color(color_set[df_delta_peaks["total"]])
        self.spectrum_peak_lines.set_linewidth(
            linewidth_set[df_delta_peaks["total"]])

        # annotate the peaks in the delta series
        for annotation in self.spectrum_axes.findobj(
                matplotlib.text.Annotation):
            annotation.remove()
        for peak_id in np.where(df_delta_peaks["total"] > 0)[0]:
            label = df_delta_distances["total"][peak_id]
            self.spectrum_axes.annotate(
                s=label,
                xy=(self._exp_mass_data.iloc[peak_id]["Average Mass"],
                    self._exp_mass_data.iloc[peak_id]["Relative Abundance"]),
                xytext=(0, 5),
                textcoords="offset pixels",
                horizontalalignment="center")
        self.spectrum_canvas.draw()

        return list(np.where(df_delta_peaks["total"] > 0)[0])


    @staticmethod
    def find_in_intervals(value, intervals=None):
        """
        Simple O(n) algorithm to determine whether a value
        falls into a set of intervals.
        Example: value=12, intervals={"a": (1, 6), "b": (9, 14)} -> "b"
                 value=8,  intervals={"a": (1, 6), "b": (9, 14)} -> ""

        :param value: Value to search
        :param intervals: {interval name: (lower interval boundary,
                                           upper interval boundary)} dict
        :return: Name of the interval containing the value;
                 empty string if no such interval exists
        """
        for name, (lower, upper) in intervals.items():
            if lower <= value <= upper:
                return name
        return ""


    def highlight_selected_peaks(self, peak_indices):
        """
        Highlight selected peaks in the spectrum.

        :param peak_indices: list of indices of peaks
                             that should be highlighted
        :return: nothing
        """

        peaks_with_result = np.zeros(self.lwPeaks.count(), dtype=int)
        try:
            peaks_with_result[self._polymer_hits.index.levels[0]] = 2
        except AttributeError:
            try:
                peaks_with_result[self._monomer_hits.index.levels[0]] = 2
            except AttributeError:
                pass


        selected_peaks = np.zeros(self.lwPeaks.count(), dtype=int)
        selected_peaks[peak_indices] = 1

        # peak colors will be an array with one entry per peak:
        # no polymers: 0 - not selected, 1 - selected
        # polymers:    2 - not selected, 3 - selected
        # alternative colors: orange [1, .49, .16, 1.0],
        #                     light red [1, .66, .66, 1.0]
        peak_colors = selected_peaks + peaks_with_result
        color_set = np.array(["#000000",   # black: unsel. peaks, no glycans
                              "#ffcc00",   # yellow: selected peaks, no glycans
                              "#00c000",   # green: unsel. peaks, glycans
                              "#ff0000"])  # red: selected peaks, glycans
        self.spectrum_peak_lines.set_color(color_set[peak_colors])
        self.spectrum_peak_lines.set_linewidth(1)

        # label the selection by masses
        for annotation in self.spectrum_axes.findobj(
                matplotlib.text.Annotation):
            annotation.remove()
        if self.btLabelPeaks.isChecked():
            for peak_id in peak_indices:
                self.spectrum_axes.annotate(
                    s="{:.2f}".format(self._exp_mass_data
                                      .iloc[peak_id]["Average Mass"]),
                    xy=(self._exp_mass_data
                            .iloc[peak_id]["Average Mass"],
                        self._exp_mass_data
                            .iloc[peak_id]["Relative Abundance"]),
                    xytext=(0, 5),
                    textcoords="offset pixels",
                    horizontalalignment="center",
                    bbox=dict(facecolor="white", alpha=.75, linewidth=0))
        self.spectrum_canvas.draw()


    def show_results(self, selected_peaks=None):
        """
        Update the spectrum, the list of masses and the result tree
        after the selection or parameters influencing the selection
        have changed.

        :param selected_peaks: list of selected peaks
        :return: nothing
        """

        if selected_peaks is None:
            selected_peaks = [i.row() for i in self.lwPeaks.selectedIndexes()]

        if self.bgSpectrum.checkedButton() == self.btModeSelection:
            self.highlight_selected_peaks(selected_peaks)
        else:
            selected_peaks = self.highlight_delta_series()

        # update the selection in the mass list
        self.lwPeaks.blockSignals(True)
        for i in range(self.lwPeaks.count()):
            self.lwPeaks.item(i).setSelected(i in selected_peaks)
        self.lwPeaks.blockSignals(False)

        # show the first currently selected mass
        # or the mass of the central peak in the delta series
        # in the mass list
        if self.spectrum_picked_peak is not None:
            self.lwPeaks.scrollToItem(
                self.lwPeaks.item(self.spectrum_picked_peak))
        else:
            self.lwPeaks.scrollToItem(self.lwPeaks.item(selected_peaks[0]))

        # fill the single mass spin box with the currently selected mass
        try:
            self.sbSingleMass.setValue(
                float(self.lwPeaks.currentItem().text()))
        except AttributeError:  # occurs when second peak file is loaded
            pass

        # activate the polymer hit filter if more than one peak is selected
        if len(selected_peaks) > 1:
            self.chFilterStructureHits.blockSignals(True)
            self.chFilterStructureHits.setChecked(True)
            self.chFilterStructureHits.blockSignals(False)

        # update the results tree if available
        if self._monomer_hits is not None:
            self.twResults.clear()
            self.twResults.setUpdatesEnabled(False)

            missing_color = QColor(255, 185, 200)

            if (self.chFilterStructureHits.isChecked()
                    and self._polymer_hits is not None):
                df_hit = self._polymer_hits
            else:
                df_hit = self._monomer_hits

            # set column headers
            header_labels = ["Exp. Mass", "%"]
            header_labels.extend(df_hit.columns)
            self.twResults.setColumnCount(len(header_labels))
            self.twResults.setHeaderLabels(header_labels)

            for mass_index in selected_peaks:
                # generate root item (experimental mass, relative abundance)
                root_item = SortableTreeWidgetItem(self.twResults)
                root_item.setText(
                    0, "{:.2f}".format(self._exp_mass_data
                                       .loc[mass_index]["Average Mass"]))
                root_item.setTextAlignment(0, Qt.AlignRight)
                root_item.setText(
                    1, "{:.1f}".format(self._exp_mass_data
                                       .loc[mass_index]["Relative Abundance"]))
                root_item.setTextAlignment(1, Qt.AlignRight)
                root_item.mass_index = mass_index

                if mass_index not in df_hit.index.levels[0]:
                    root_item.setBackground(0, QBrush(missing_color))
                    root_item.setBackground(1, QBrush(missing_color))
                else:
                    # generate one child item per possible combination
                    for _, hit in df_hit.loc[mass_index].iterrows():
                        child_item = SortableTreeWidgetItem(root_item)
                        monomers = (hit[:df_hit.columns.get_loc("Exp. Mass")]
                                    .index)

                        # monomer counts
                        for j, monomer in enumerate(monomers):
                            child_item.setText(
                                2 + j, "{:.0f}".format(hit[monomer]))
                            child_item.setTextAlignment(2 + j, Qt.AlignHCenter)

                        # hit properties
                        for j, label in enumerate(["Exp. Mass", "Theo. Mass",
                                                   "Da.", "ppm"]):
                            child_item.setText(
                                len(monomers) + 2 + j,
                                "{:.2f}".format(hit[label]))
                            child_item.setTextAlignment(
                                len(monomers) + 2 + j,
                                Qt.AlignRight)

                        if (self.chFilterStructureHits.isChecked()
                                and self._polymer_hits is not None):
                            # polymer composition
                            sites = (hit[df_hit.columns.get_loc("ppm") + 1: -1]
                                     .index)
                            for j, site in enumerate(sites):
                                child_item.setText(
                                    len(monomers) + 6 + j,
                                    "{}".format(hit[site]))
                                child_item.setTextAlignment(
                                    len(monomers) + 6 + j,
                                    Qt.AlignHCenter)

                            # polymer abundance
                            child_item.setText(
                                len(monomers) + 7 + j,
                                "{:.2f}".format(hit["Abundance"]))
                            child_item.setTextAlignment(
                                len(monomers) + 7 + j,
                                Qt.AlignHCenter)

            self.twResults.expandAll()
            self.twResults.header().setSectionResizeMode(
                QHeaderView.ResizeToContents)
            self.twResults.header().setStretchLastSection(False)
            self.twResults.setUpdatesEnabled(True)


    def save_settings(self):
        """
        Export the current settings as XML file.

        :return: nothing
        """
        filename = QFileDialog.getSaveFileName(
            self,
            "Save settings",
            self._path,
            "ModFinder XML settings (*.xml)")[0]
        self._path = os.path.split(filename)[0]
        if filename:
            if not filename.endswith(".xml"):
                filename = filename + ".xml"

            root = ETree.Element("settings")
            settings = [("sequence", self.teSequence.toPlainText()),
                        ("disulfides", self.sbDisulfides.value()),
                        ("pngasef", self.chPngase.isChecked()),
                        ("tolerance-value", self.sbTolerance.value()),
                        ("tolerance-flavor", self.cbTolerance.currentIndex())]
            for child, text in settings:
                ETree.SubElement(root, child).text = str(text)

            dataframes = [("masslist", self._exp_mass_data),
                          ("monomers", self.table_to_df(which="monomers")),
                          ("polymers", self.table_to_df(which="polymers"))]
            for child, df in dataframes:
                ETree.SubElement(root, child).append(dataframe_to_xml(df))

            prettify_xml(root)
            try:
                ETree.ElementTree(root).write(
                    filename,
                    encoding="utf-8",
                    xml_declaration=True)
            except OSError:
                QMessageBox.critical(
                    self,
                    "Error",
                    "Error when writing to " + filename + OSError.args)


    def load_settings(self):
        """
        Load settings form an XML file.

        :return: nothing
        """

        filename = QFileDialog.getOpenFileName(
            self,
            "Load settings",
            self._path,
            "ModFinder XML settings (*.xml)")[0]
        self._path = os.path.split(filename)[0]
        if filename:
            try:
                root = ETree.ElementTree().parse(filename)
            except IOError:
                QMessageBox.critical(
                    self,
                    "Error",
                    "Error loading " + filename + OSError.args)
                return

            self.teSequence.setText(root.find("sequence").text)
            self.sbDisulfides.setValue(bool(root.find("disulfides").text))
            self.chPngase.setChecked(bool(root.find("pngasef").text))
            self.calculate_protein_mass()

            self._monomer_hits = None
            self._polymer_hits = None
            self.chFilterStructureHits.setEnabled(False)
            self.twResults.clear()
            self._exp_mass_data = dataframe_from_xml(root.find("masslist"))
            if self._exp_mass_data is not None:
                self.fill_peak_list(self._exp_mass_data["Average Mass"])
                self.draw_spectrum()
                self.spectrum_picked_peak = 0
                self.show_results()

            self.cbTolerance.setCurrentIndex(
                int(root.find("tolerance-flavor").text))
            self.sbTolerance.setValue(float(root.find("tolerance-value").text))

            self.table_from_df(
                dataframe_from_xml(root.find("monomers")),
                self.tbMonomers,
                _monomer_table_columns)
            self.table_from_df(
                dataframe_from_xml(root.find("polymers")),
                self.tbPolymers,
                _polymer_table_columns)
            self.calculate_mod_mass()


    def save_spectrum(self):
        filename = QFileDialog.getSaveFileName(
            self,
            "Save spectrum",
            self._path,
            "Portable network graphics (*.png)")[0]
        self._path = os.path.split(filename)[0]
        if filename:
            try:
                self.spectrum_fig.savefig(filename, dpi=200)
            except OSError:
                QMessageBox.critical(
                    self,
                    "Error",
                    "Error when writing to " + filename + OSError.args)


def main():  # possible argument: argv=sys.argv
    app = QApplication(sys.argv)
    frame = MainWindow()
    frame.show()
    app.exec_()


if __name__ == "__main__":
    main()
