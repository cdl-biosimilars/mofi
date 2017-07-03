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

import distutils.util
import math
import os
import pickle
import re
import sys
import time

from qtpy.QtWidgets import (QApplication, QMainWindow, QMenu, QActionGroup, QVBoxLayout, QTableWidgetItem, QCheckBox,
                            QMessageBox, QFileDialog, QTreeWidgetItem, QHeaderView, QSpinBox, QDoubleSpinBox,
                            QWidget, QHBoxLayout, QAction, QToolBar, QProgressBar, QLabel, QSizePolicy, QButtonGroup)
from qtpy.QtGui import QColor, QBrush
from qtpy.QtCore import Qt

import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 5000)

import matplotlib
import matplotlib.text
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.widgets import SpanSelector, RectangleSelector
from matplotlib.figure import Figure

import configure
import mass_tools
import modification_search
import sequence_tools
from ModFinder_UI import Ui_ModFinder


class CollapsingRectangleSelector(RectangleSelector):
    """
    Select a rectangular region of an axes.
    The rectangle collapses to a line if a dimension is less than minspanx and minspany.
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
        Overrides method from parent: Calculate the coordinates of the drawn rectangle.

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
    Additionally, the attribute mass_index stores mass index of a top level widget.
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
        return tuple((e if i % 2 == 0 else float(e)) for i, e in enumerate(parts))


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
        self.acSaveAnnotation.triggered.connect(self.save_search_results)
        self.acSaveSettings.triggered.connect(self.save_settings)

        self.btClearMonomers.clicked.connect(lambda: self.table_clear(self.tbMonomers))
        self.btClearPolymers.clicked.connect(lambda: self.table_clear(self.tbPolymers))
        self.btDeleteRowMonomers.clicked.connect(lambda: self.table_delete_row(self.tbMonomers))
        self.btDeleteRowPolymers.clicked.connect(lambda: self.table_delete_row(self.tbPolymers))
        self.btFindModifications.clicked.connect(self.sample_modifications)
        self.btInsertRowAboveMonomers.clicked.connect(lambda: self.table_insert_row(self.tbMonomers, above=True))
        self.btInsertRowAbovePolymers.clicked.connect(lambda: self.table_insert_row(self.tbPolymers, above=True))
        self.btInsertRowBelowMonomers.clicked.connect(lambda: self.table_insert_row(self.tbMonomers, above=False))
        self.btInsertRowBelowPolymers.clicked.connect(lambda: self.table_insert_row(self.tbPolymers, above=False))
        self.btLabelPeaks.clicked.connect(lambda: self.display_selected_peaks())
        self.btLoadMonomers.clicked.connect(self.load_monomers)
        self.btLoadPolymers.clicked.connect(self.load_polymers)
        self.btResetZoom.clicked.connect(self.reset_zoom)
        self.btSaveMonomers.clicked.connect(self.save_monomers)
        self.btSavePolymers.clicked.connect(self.save_polymers)
        self.btUpdateMass.clicked.connect(self.calculate_protein_mass)

        self.cbTolerance.activated.connect(self.choose_tolerance_units)

        self.chDelta1.clicked.connect(self.toggle_delta_series)
        self.chDelta2.clicked.connect(self.toggle_delta_series)
        self.chFilterStructureHits.clicked.connect(self.filter_structure_hits)
        self.chPngase.clicked.connect(self.calculate_protein_mass)

        self.lwPeaks.itemSelectionChanged.connect(self.select_peaks_in_list)

        self.sbDeltaRepetition1.valueChanged.connect(lambda: self.display_selected_peaks())
        self.sbDeltaRepetition2.valueChanged.connect(lambda: self.display_selected_peaks())
        self.sbDeltaTolerance1.valueChanged.connect(lambda: self.display_selected_peaks())
        self.sbDeltaTolerance2.valueChanged.connect(lambda: self.display_selected_peaks())
        self.sbDeltaValue1.valueChanged.connect(lambda: self.display_selected_peaks())
        self.sbDeltaValue2.valueChanged.connect(lambda: self.display_selected_peaks())
        self.sbDisulfides.valueChanged.connect(self.calculate_protein_mass)

        self.teSequence.textChanged.connect(
            lambda: self.teSequence.setStyleSheet("QTextEdit { background-color: rgb(255, 225, 225) }"))

        self.twResults.itemClicked.connect(self.select_peaks_in_tree)

        self.menubar.setVisible(False)

        # add toolbars
        self.tlMassSet = QToolBar(self)
        self.tlMassSet.setObjectName("tlMassSet")
        self.addToolBar(Qt.TopToolBarArea, self.tlMassSet)
        self.tlHelp = QToolBar(self)
        self.tlHelp.setObjectName("tlHelp")
        self.tlHelp.addActions([self.acAbout, self.acHelp])
        self.addToolBar(Qt.TopToolBarArea, self.tlHelp)

        # generate mass set selectors from config file
        self.agSelectMassSet = QActionGroup(self.menuAtomicMasses)
        set_id = 0
        for mass_set in configure.mass_sets:
            ac_select_set = QAction(self)
            ac_select_set.setCheckable(True)
            if set_id == 0:
                ac_select_set.setChecked(True)
            ac_select_set.setObjectName("acSelectMassSet{:d}".format(set_id))
            ac_select_set.setText(mass_set)
            ac_select_set.setToolTip(configure.mass_sets[mass_set].get("description", ""))
            self.menuAtomicMasses.addAction(ac_select_set)
            self.tlMassSet.addAction(ac_select_set)
            self.agSelectMassSet.addAction(ac_select_set)
            set_id += 1
        # noinspection PyUnresolvedReferences
        self.agSelectMassSet.triggered.connect(self.choose_mass_set)

        # init statusbar
        sb_widget = QWidget(self)
        sb_layout = QHBoxLayout(sb_widget)
        sb_layout.setContentsMargins(0, 0, 0, 0)
        sb_widget.setLayout(sb_layout)
        sb_layout.addWidget(QLabel("Progress:"))
        self.pbSearchProgress = QProgressBar(self)
        self.pbSearchProgress.setValue(0)
        self.pbSearchProgress.setMinimumWidth(300)
        self.pbSearchProgress.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Fixed)
        sb_layout.addWidget(self.pbSearchProgress)
        self.statusbar.addPermanentWidget(sb_widget)

        # set up the mass spectrum plot
        layout = QVBoxLayout(self.spectrumView)
        self.spectrum_fig = Figure(dpi=100, frameon=False, tight_layout={"pad": 0}, edgecolor="white")
        self.spectrum_canvas = FigureCanvas(self.spectrum_fig)
        self.spectrum_canvas.setParent(self.spectrumView)
        layout.addWidget(self.spectrum_canvas)
        self.spectrum_axes = None  # single Axes of the figure
        self.spectrum_peak_lines = None  # LineCollection object representing the peaks in the spectrum
        self.spectrum_pick_selector = None  # picker connection to select a single peak
        self.spectrum_span_selector = None  # SpanSelector to select multiple peaks
        self.spectrum_rectangle_selector = None  # CollapsingRectangleSelector to zoom the spectrum
        self.spectrum_x_limits = None  # original limits of the x axis when the plot is drawn
        self.spectrum_picked_peak = None  # the single peak that was picked, or (in a span) a representative singe peak

        # init button group for specrum interaction mode
        self.bgSpectrum = QButtonGroup()
        self.bgSpectrum.addButton(self.btModeDelta)
        self.bgSpectrum.addButton(self.btModeSelection)
        # noinspection PyUnresolvedReferences
        self.bgSpectrum.buttonClicked.connect(self.toggle_spectrum_mode)

        # init the monomer table and associated buttons
        self.tbMonomers.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        for col, width in [(0, 40), (2, 130), (3, 45), (4, 45), (5, 25)]:
            self.tbMonomers.setColumnWidth(col, width)
        self.tbMonomers.verticalHeader().setSectionResizeMode(QHeaderView.Fixed)
        self.tbMonomers.verticalHeader().setDefaultSectionSize(22)

        menu = QMenu()
        for library in configure.default_monomer_libraries:
            menu.addAction(library, self.load_default_monomers)
        self.btDefaultModsMonomers.setMenu(menu)

        # init the polymer table and associated buttons
        self.tbPolymers.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        for col, width in [(0, 40), (2, 130), (3, 80), (4, 50)]:
            self.tbPolymers.setColumnWidth(col, width)
        self.tbPolymers.verticalHeader().setSectionResizeMode(QHeaderView.Fixed)
        self.tbPolymers.verticalHeader().setDefaultSectionSize(22)

        menu = QMenu()
        for library in configure.default_polymer_libraries:
            menu.addAction(library, self.load_default_polymers)
        self.btDefaultModsPolymers.setMenu(menu)

        # initialize private members
        self._monomer_hits = None  # results from the monomer search
        self._delta_lines = [None, None]  # delta lines in the mass spectrum
        self._exp_mass_data = None  # pandas dataframe containing the contents of the mass file
        self._known_mods_mass = 0  # mass of known modification
        self._mass_filename = None  # name of the mass file
        self._path = configure.path  # last path selected in a file chooser dialog
        self._polymer_hits = None  # results from the polymer search
        self._protein = None  # a Protein object representing the input sequence with disulfides and PNGase F digest
        self._protein_mass = 0  # mass of the current Protein object


    def _monomer_table_create_row(self, row_id, active=False, name="", composition="",
                                  min_count=0, max_count=-1, part_of_polymer=False):
        """
        Create a new row in the table of modifications.
        This row describes a modification with the given parameters.

        Private function, to be called by general row insertion functions and
        loaders for modification libraries.

        :param row_id: Row index passed to QTableWidget.insertRow()
        :param active: true if the monomer should be used in the combinatorial search
        :param name: name of the modification (str)
        :param composition: composition or mass of the modification (str)
        :param min_count: minimum number of modifications (int)
        :param max_count: maximum number of modifications (int)
        :param part_of_polymer: true if the monomer should be used in the polymer search
        :return: nothing
        """
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

        monomer_checkbox = QCheckBox()
        monomer_checkbox.setChecked(part_of_polymer)
        monomer_ch_widget = QWidget()
        monomer_ch_layout = QHBoxLayout(monomer_ch_widget)
        monomer_ch_layout.addWidget(monomer_checkbox)
        monomer_ch_layout.setAlignment(Qt.AlignCenter)
        monomer_ch_layout.setContentsMargins(0, 0, 0, 0)
        self.tbMonomers.setCellWidget(row_id, 5, monomer_ch_widget)


    def _polymer_table_create_row(self, row_id, active=True, name="", composition="", sites="", abundance=0.0):
        """
        Create a new row in the table of modifications.
        This row describes a modification with the given parameters.

        Private function, to be called by general row insertion functions and
        loaders for modification libraries.

        :param row_id: Row index passed to QTableWidget.insertRow()
        :param active: true if the polymer should be used in the combinatorial search
        :param name: name of the modification (str)
        :param composition: composition or mass of the modification (str)
        :param sites: identifier for the modification site; empty means any site (str)
        :param abundance: relative abundance (float)
        :return: nothing
        """
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


    def table_insert_row(self, table_widget, above=True):
        """
        Insert a row into the table of modifications.
        The row will be inserted relative to the current selection (if one exists) or to all rows otherwise.

        :param table_widget: the QTableWidget to modify
        :param above: True if the row should be inserted above the current selection
        :return: nothing
        """
        if table_widget.selectionModel().selectedRows():
            if above:
                last_row = table_widget.selectionModel().selectedRows()[0].row()
            else:
                last_row = table_widget.selectionModel().selectedRows()[-1].row() + 1
        else:
            if above:
                last_row = 0
            else:
                last_row = table_widget.rowCount()
        if table_widget == self.tbMonomers:
            self._monomer_table_create_row(last_row)
        else:
            self._polymer_table_create_row(last_row)


    @staticmethod
    def table_clear(table_widget):
        """
        Delete all rows in the table of modifications.

        :param table_widget: the QTableWidget to modify
        :return: nothing
        """
        while table_widget.rowCount() > 0:
            table_widget.removeRow(0)


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


    def show_about(self):
        """
        Show the about dialog.

        :return: nothing
        """
        v = "ModFinder version: {}\n".format(configure.version)
        r = configure.rights
        c = "\nContact: {}".format(configure.contact)
        QMessageBox.about(self, "About ModFinder", v + r + c)


    def show_help(self):
        """
        Show the help dialog.

        :return: nothing
        """
        with open("../docs/help.txt", "r") as helpfile:
            QMessageBox.about(self, "Help", helpfile.read())


    def choose_tolerance_units(self):
        """
        Adjust the settings of the tolerance spin box when PPM or Da are selected.

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
        Opens a FASTA file and displays its contents in the QTextEdit

        Changes:
            self._path to directory of selected file
            text of self.teSequence to FASTA string
            contents and maximum value of self.sbDisulfides to half the number of cysteines in the sequence

        :return: nothing
        """

        filename = QFileDialog.getOpenFileName(self,
                                               "Open FASTA file",
                                               self._path,
                                               "Sequence files (*.fasta *.txt)")[0]
        self._path = os.path.split(filename)[0]
        ext = os.path.splitext(filename)[1]
        if filename:
            if ext in [".fasta", ".txt"]:
                with open(filename) as f:
                    fasta_input = f.read()
                self.teSequence.setText(fasta_input)
                self.calculate_protein_mass()
            else:
                QMessageBox.warning(self, "Error loading sequence", "Not a valid file format: {}".format(ext))


    def load_mass_file(self):
        """
        Opens a mass list as generated by Thermo BioPharma Finder and display its contents.

        :return: nothing
        """

        filename = QFileDialog.getOpenFileName(self,
                                               "Open mass list",
                                               self._path,
                                               "Excel files (*.xlsx *.xls);; CSV files (*.csv *.txt)")[0]
        self._path = os.path.split(filename)[0]
        ext = os.path.splitext(filename)[1]
        if filename:
            self._mass_filename = os.path.split(filename)[1]
            mass_data = mass_tools.read_massfile(filename, sort_by="Average Mass")
            if mass_data is None:
                QMessageBox.warning(self, "Error loading mass file", "Not a valid file format: {}".format(ext))
                return

            self._exp_mass_data = mass_data
            self.twResults.clear()
            self._monomer_hits = None
            self._polymer_hits = None
            self.chFilterStructureHits.setEnabled(False)
            self.lwPeaks.blockSignals(True)
            self.lwPeaks.clear()
            self.lwPeaks.addItems(["{:.2f}".format(i) for i in self._exp_mass_data["Average Mass"]])
            self.lwPeaks.setCurrentRow(0)
            self.lwPeaks.blockSignals(False)
            self.draw_spectrum()


    def save_search_results(self):
        """
        Write the results from a combinatorial search to a CSV file.

        :return: nothing
        """
        outfilename = QFileDialog.getSaveFileName(self,
                                                  "Save results",
                                                  self._path,
                                                  "Comma-separated value (*.csv)")[0]
        self._path = os.path.split(outfilename)[0]
        if outfilename:
            if not outfilename.endswith(".csv"):
                outfilename += ".csv"
            parameters = [("Disulfides", self.sbDisulfides.value()),
                          ("PNGaseF", self.chPngase.isChecked()),
                          ("Tolerance", self.sbTolerance.value()),
                          ("Tol. Type", self.cbTolerance.currentText())]

            for row_id in range(self.tbMonomers.rowCount()):
                ch = self.tbMonomers.cellWidget(row_id, 0).findChild(QCheckBox)
                name = self.tbMonomers.item(row_id, 1).text()
                min_count = self.tbMonomers.cellWidget(row_id, 3).value()
                max_count = self.tbMonomers.cellWidget(row_id, 4).value()
                if ch.isChecked():
                    parameters.append((name + "_min", min_count))
                    parameters.append((name + "_max", max_count))

            mindex = self._monomer_hits.reset_index(level=0)["Massindex"]  # TODO also save polymer hits; optimize code
            ra = mindex.map(self._exp_mass_data["Relative Abundance"]).reset_index()
            self._monomer_hits["Relative Abundance"] = list(ra["Massindex"])

            try:
                parameter_string = ["{}: {}".format(k, v) for k, v in parameters]
                with open(outfilename, "w") as f:
                    f.write("# Combinatorial search results by ModFinder\n")
                    f.write("# Date: " + time.strftime("%c") + "\n")
                    f.write("# Parameters: ")
                    f.write(", ".join(parameter_string))
                    f.write("\n")
                    self._monomer_hits.to_csv(f)
            except IOError:
                QMessageBox.warning(self, "Warning", "Permission denied for {}".format(outfilename))


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

        Changes:
            self._protein to a protein with given sequence, disulfides and PNGase F modifications
            self._protein_mass to the mass of self._protein
            updates the value of self.lbMassProtein, self.lbMassMods and self.lbMassTotal

        :return: True if there was no error, otherwise False
        """

        protein_sequence = self.teSequence.toPlainText()
        chains, sequence = sequence_tools.read_fasta_string(protein_sequence)
        try:
            self._protein = sequence_tools.Protein(sequence,
                                                   chains,
                                                   self.sbDisulfides.value(),
                                                   self.chPngase.isChecked())
        except KeyError as e:
            QMessageBox.critical(self,
                                 "Error",
                                 "Error when parsing sequence: {} is not a valid symbol".format(e.args[0]))
            return False

        self.sbDisulfides.setEnabled(True)
        self.chPngase.setEnabled(True)
        self.sbDisulfides.setMaximum(self._protein.amino_acid_composition["C"] / 2)
        self.teSequence.setStyleSheet("QTextEdit { background-color: rgb(240, 251, 240) }")
        self._protein_mass = self._protein.mass
        self.lbMassProtein.setText("{:,.2f}".format(self._protein_mass))
        self.lbMassMods.setText("{:,.2f}".format(self._known_mods_mass))
        self.lbMassTotal.setText("{:,.2f}".format(self._protein_mass + self._known_mods_mass))
        return True


    def calculate_mod_mass(self, return_values=""):
        """
        Calculate the mass of known modifications.

        Changes:
            self._known_mods_mass to the mass of known modifications
            updates the value of self.lbMassProtein, self.lbMassMods and self.lbMassTotal

        :param return_values: influences the return value (see below)
        :return: if return_values is "all": list of (checked, name, composition, min count,
                                                     max count, include in polymer search) tuples
                 otherwise list of (name, mass, min count, max count, include in polymer search) tuples
        """

        self._known_mods_mass = 0
        result = []

        # add min counts for glycans to the theoretical mass
        for row_id in range(self.tbMonomers.rowCount()):
            ch = self.tbMonomers.cellWidget(row_id, 0).findChild(QCheckBox)
            name = self.tbMonomers.item(row_id, 1).text()
            composition = self.tbMonomers.item(row_id, 2).text()
            min_count = self.tbMonomers.cellWidget(row_id, 3).value()
            max_count = self.tbMonomers.cellWidget(row_id, 4).value()
            is_poly = self.tbMonomers.cellWidget(row_id, 5).findChild(QCheckBox).isChecked()
            if return_values == "all":
                result.append((ch.isChecked(), name, composition, min_count, max_count, is_poly))
            else:
                if ch.isChecked():
                    mass = 0
                    try:  # 'composition' could be a mass
                        mass = float(composition)
                    except ValueError:  # 'composition' could be a formula
                        try:
                            monomer_formula = mass_tools.Formula(composition)
                            mass = monomer_formula.mass
                        except ValueError:  # 'composition' is invalid
                            QMessageBox.critical(self,
                                                 "Error",
                                                 "Invalid formula in row {:d}: {}".format(row_id + 1, composition))
                    self._known_mods_mass += mass * min_count
                    result.append((name, mass, min_count, max_count, is_poly))

        self.lbMassProtein.setText("{:,.2f}".format(self._protein_mass))
        self.lbMassMods.setText("{:,.2f}".format(self._known_mods_mass))
        self.lbMassTotal.setText("{:,.2f}".format(self._protein_mass + self._known_mods_mass))
        return result


    def sample_modifications(self):
        """
        Prepare data for the two-stage search and process its results.

        :return: nothing
        """

        self.calculate_protein_mass()

        # calculate required input if a single mass was entered (i.e., no peak list was loaded)
        if self.rbSingleMass.isChecked():
            self._exp_mass_data = pd.DataFrame({"Average Mass": self.sbSingleMass.value(),
                                                "Relative Abundance": 100.0}, index=[0])
            self._mass_filename = "Input Mass: {:.2f}".format(self.sbSingleMass.value())
            self.lwPeaks.blockSignals(True)
            self.lwPeaks.clear()
            self.lwPeaks.addItems(["{:.2f}".format(i) for i in self._exp_mass_data["Average Mass"]])
            self.lwPeaks.setCurrentRow(0)
            self.lwPeaks.blockSignals(False)
            self.draw_spectrum()

        if self._exp_mass_data is None:
            QMessageBox.critical(self, "Error", "No mass list loaded. Aborting search.")
            return

        monomers = self.calculate_mod_mass()
        modifications = []  # list of modifications to be used in the combinatorial search
        explained_mass = self._protein_mass + self._known_mods_mass
        unknown_masses: pd.DataFrame = self._exp_mass_data["Average Mass"] - explained_mass

        if self.cbTolerance.currentIndex() == 0:  # that is, "Da."
            # calculate largest mass plus tolerance
            max_tol_mass = max(self._exp_mass_data["Average Mass"]) + self.sbTolerance.value()
            mass_tolerance = self.sbTolerance.value()
        else:
            max_tol_mass = max(self._exp_mass_data["Average Mass"]) * (1 + self.sbTolerance.value() / 1_000_000)
            mass_tolerance = []  # calculate a mass tolerance for each peak if we're working with ppm tolerance
            for _, m in self._exp_mass_data["Average Mass"].iteritems():
                mass_tolerance.append(m * self.sbTolerance.value() / 1_000_000)

        # calculate the upper limit of glycans that may appear
        # add all checked single glycans to the list of modifications
        monomers_for_polymer_search = []
        for name, mass, min_count, max_count, is_poly in monomers:
            if max_count == -1:
                max_count = min(int((max_tol_mass - self._protein_mass) / mass), configure.maxmods)
            modifications.append((name, mass, max_count - min_count))
            if is_poly:
                monomers_for_polymer_search.append(name)

        if not modifications:
            QMessageBox.critical(self, "Error", "List of monomers is empty. Aborting search.")
            return

        # prepare list of polymers for search stage 2
        polymers = {}
        for row_id in range(self.tbPolymers.rowCount()):
            if self.tbPolymers.cellWidget(row_id, 0).findChild(QCheckBox).isChecked():
                name = self.tbPolymers.item(row_id, 1).text()
                composition = self.tbPolymers.item(row_id, 2).text()
                sites = self.tbPolymers.item(row_id, 3).text()
                abundance = self.tbPolymers.cellWidget(row_id, 4).value()
                polymers[name] = (composition, sites, abundance)

        if polymers:  # dict contains at least one entry
            df_polymers = pd.DataFrame.from_dict(polymers, orient="index")
            df_polymers.columns = ["Composition", "Sites", "Abundance"]
            monomers_in_library = modification_search.get_monomers_from_library(df_polymers)
        else:  # dict remained empty, i.e., there are no polymers
            monomers_in_library = []
            df_polymers = pd.DataFrame()

        if sorted(monomers_for_polymer_search) != sorted(monomers_in_library):  # TODO make this check automatic!
            monomers_only_checked = set(monomers_for_polymer_search) - set(monomers_in_library)
            monomers_only_in_library = set(monomers_in_library) - set(monomers_for_polymer_search)
            error_message = []
            if monomers_only_checked:
                error_message += ["The following monomers are checked for polymer search ",
                                  "but do not appear in the library: "]
                error_message += " ".join(sorted(monomers_only_checked))
                error_message.append(".\n")
            if monomers_only_in_library:
                error_message += ["The following monomers appear in the library ",
                                  "but are not checked for polymer search: "]
                error_message += " ".join(sorted(monomers_only_in_library))
                error_message.append(".\n")
            QMessageBox.critical(self,
                                 "Error",
                                 "".join(error_message))
            return

        self.statusbar.showMessage("Starting monomer search (stage 1) ...")
        print("Experimental Masses:", self._exp_mass_data["Average Mass"].head(), sep="\n")
        print("Explained mass (protein + known modifications):", explained_mass)
        print("Unknown masses searched:", unknown_masses.head(), sep="\n")
        print("Mass tolerance: {:f} {}".format(self.sbTolerance.value(), self.cbTolerance.currentText()))
        print("Monomers used in search:\nName\tMass\tmax")
        for m in modifications:
            print("{}\t{:.2f}\t{:d}".format(*m))
        print("Monomers checked for polymer search:")
        print(", ".join(monomers_for_polymer_search))
        print("Monomers extracted from library:")
        print(", ".join(monomers_in_library))

        # stage 1: monomer search
        self._monomer_hits = modification_search.find_monomers(
            modifications,
            list(unknown_masses),
            mass_tolerance=mass_tolerance,
            explained_mass=explained_mass,
            progress_bar=self.pbSearchProgress)

        self.statusbar.showMessage("Monomer search done! Starting polymer search (stage 2) ...")
        # the modification search was not successful
        if self._monomer_hits is None:
            QMessageBox.critical(self, "Error", "Combinatorial search was unsuccessful.")
            return

        # add the minimum monomer counts to the result data frame
        for name, _, min_count, _, _ in monomers:
            self._monomer_hits[name] += min_count

        # stage 2: polymer search
        if polymers:
            self._polymer_hits = modification_search.find_polymers(self._monomer_hits,
                                                                   glycan_library=df_polymers,
                                                                   monomers=monomers_for_polymer_search,
                                                                   progress_bar=self.pbSearchProgress)
            self.statusbar.showMessage("Polymer search DONE!", 5000)
        self.chFilterStructureHits.setEnabled(True)
        self.display_selected_peaks()


    def draw_spectrum(self):
        """
        Show a mass spectrum after loading a peak list.

        :return: nothing
        """

        self.spectrum_fig.clear()
        self.spectrum_axes = self.spectrum_fig.add_subplot(111)
        self.spectrum_axes.set_xmargin(.02)
        self.spectrum_peak_lines = self.spectrum_axes.vlines(x=self._exp_mass_data["Average Mass"],
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
        self.spectrum_pick_selector = self.spectrum_fig.canvas.mpl_connect("pick_event", self.select_peaks_by_pick)
        self.spectrum_span_selector = SpanSelector(self.spectrum_axes,
                                                   self.select_peaks_by_span,
                                                   "horizontal",
                                                   minspan=10,  # in order not to interfere with the picker
                                                   rectprops=dict(alpha=.1),
                                                   button=1)  # only the left mouse button spans
        self.spectrum_rectangle_selector = CollapsingRectangleSelector(self.spectrum_axes,
                                                                       self.zoom_spectrum,
                                                                       collapsex=50,
                                                                       collapsey=5,
                                                                       button=3,
                                                                       rectprops=dict(facecolor=(1, 0, 0, .1),
                                                                                      edgecolor="black",
                                                                                      fill=True,
                                                                                      linewidth=1.5))


    def zoom_spectrum(self, start, stop):
        """
        Set the new limits of the x and y axis after drawing the zoom rectangle.

        :param start: MouseEvent that describes the coordinates where the mouse button was pressed.
        :param stop: MouseEvent that describes the coordinates where the mouse button was released.
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
        Enable the delta series spin boxes whose correspondihc checkbox is checked.

        :return: nothing
        """
        for ch, sb_delta, sb_tolerance, sb_repetitions in \
            [(self.chDelta1, self.sbDeltaValue1, self.sbDeltaTolerance1, self.sbDeltaRepetition1),
             (self.chDelta2, self.sbDeltaValue2, self.sbDeltaTolerance2, self.sbDeltaRepetition2)]:
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
        Ensure that only a single mass is displayed in the result tree if this box is unchecked.

        :return: nothing
        """
        if not self.chFilterStructureHits.isChecked():
            i = 0
            for i in range(self.lwPeaks.count()):
                if self.lwPeaks.item(i).isSelected():
                    break
            self.display_selected_peaks([i])
        else:
            self.display_selected_peaks()



    def select_peaks_in_list(self):
        """
        Select one or several peaks in the peak list.

        :return: nothing
        """
        self.spectrum_picked_peak = self.lwPeaks.currentRow()
        self.display_selected_peaks()


    def select_peaks_in_tree(self):
        """
        Select the peak corresponding to an entry in the result tree.

        :return: nothing
        """

        toplevel_item = self.twResults.currentItem()
        while toplevel_item.parent():
            toplevel_item = toplevel_item.parent()
        self.spectrum_picked_peak = toplevel_item.mass_index
        self.display_selected_peaks([toplevel_item.mass_index])


    def select_peaks_by_pick(self, event):
        """
        Select a peak picked by a mouseclick on the spectrum.

        :param event: PickEvent from the canvas
        :return: nothing
        """

        self.spectrum_picked_peak = event.ind[0]
        self.display_selected_peaks(event.ind)


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
        self.display_selected_peaks(peak_indices)


    @staticmethod
    def concat_interval_names(row):
        """
        Concatenate rows from a dataframe indicating interval names.
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

        main_mass = float(self._exp_mass_data.iloc[self.spectrum_picked_peak]["Average Mass"])
        min_mass = float(min(self._exp_mass_data["Average Mass"]))
        max_mass = float(max(self._exp_mass_data["Average Mass"]))

        # dataframe with the mass index as index, the delta series index as columns, and the intervall names als values
        df_delta_distances = pd.DataFrame(index=self._exp_mass_data.index, dtype=int)

        # dataframe with the following values:
        # 0 - peaks not in any delta mass series
        # 1 - peaks in series 1
        # 2 - peaks in series 2
        # 3 - peaks in both series
        # 4 - selected (central) peak
        df_delta_peaks = pd.DataFrame(index=self._exp_mass_data.index, dtype=int)

        # the algorithm is the same for both delta series; hence, use this loop
        for delta_id, ch_delta, sb_value, sb_tolerance, sb_repetitions in \
                [(1, self.chDelta1, self.sbDeltaValue1, self.sbDeltaTolerance1, self.sbDeltaRepetition1),
                 (2, self.chDelta2, self.sbDeltaValue2, self.sbDeltaTolerance2, self.sbDeltaRepetition2)]:
            if not ch_delta.isChecked():
                continue

            if sb_repetitions.value() == -1:
                max_iterations = int(sb_value.value() / sb_tolerance.value() / 2)
            else:
                max_iterations = sb_repetitions.value()
            intervals = {}  # will be a {number of mass differences: (interval start, interval end)} dict

            # calculate possible intervals for peaks separated from the selected peaks
            # by an integer number of mass differences
            # for each added difference, increase the interval size by the original tolerance times two
            current_mass = main_mass
            tolerance = sb_tolerance.value()
            i = 1
            while i <= max_iterations and current_mass > min_mass:
                current_mass -= sb_value.value()
                intervals[-i] = (current_mass - tolerance, current_mass + tolerance)
                tolerance += sb_tolerance.value()
                i += 1

            current_mass = main_mass
            tolerance = sb_tolerance.value()
            i = 1
            while i <= max_iterations and current_mass < max_mass:
                current_mass += sb_value.value()
                intervals[i] = (current_mass - tolerance, current_mass + tolerance)
                tolerance += sb_tolerance.value()
                i += 1

            df_delta_distances[delta_id] = self._exp_mass_data["Average Mass"] \
                .apply(self.find_in_intervals, intervals=intervals)
            df_delta_distances[delta_id][self.spectrum_picked_peak] = "0"

            df_delta_peaks[delta_id] = df_delta_distances[delta_id] \
                .map(lambda x: delta_id if x else 0)

        # calculate "total" columns describing the results from both delta series searches
        df_delta_peaks["total"] = df_delta_peaks.sum(axis="columns")
        df_delta_peaks["total"][self.spectrum_picked_peak] = 4
        df_delta_distances["total"] = df_delta_distances.apply(self.concat_interval_names, axis=1)

        # color the peaks in the delta series and increase their line width
        color_set = np.array(["#b3b3b3",   # light gray
                              "#aa0088",   # violet
                              "#2ca05a",   # greenish
                              "#6b5071",   # mixture of the prevoius two
                              "#ff0000"])  # red

        linewidth_set = np.array([1, 2, 2, 2, 3])
        self.spectrum_peak_lines.set_color(color_set[df_delta_peaks["total"]])
        self.spectrum_peak_lines.set_linewidth(linewidth_set[df_delta_peaks["total"]])

        # annotate the peaks in the delta series
        for annotation in self.spectrum_axes.findobj(matplotlib.text.Annotation):
            annotation.remove()
        for peak_id in np.where(df_delta_peaks["total"] > 0)[0]:
            self.spectrum_axes.annotate(s=df_delta_distances["total"][peak_id],
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
        Simple O(n) algorithm to determine whether a value falls into a set of intervals.
        Example: value=12,  intervals=[(1, 6), (9, 14)] -> True
                 value=8, intervals=[(1, 6), (9, 14)] -> False

        :param value: Value to search
        :param intervals: {interval name: (lower interval boundary, upper interval boundary)} dict
        :return: Name of the interval containing the value; "" if no such interval exists
        """
        for name, (lower, upper) in intervals.items():
            if lower <= value <= upper:
                return name
        return ""


    def highlight_selected_peaks(self, peak_indices):
        """
        Highlight selected peaks in the spectrum.

        :param peak_indices: list of indices of peaks that should be highlighted
        :return: nothing
        """
        polymer_peaks = np.zeros(self.lwPeaks.count(), dtype=int)
        try:
            polymer_peaks[self._polymer_hits.index.levels[0]] = 2
        except AttributeError:
            pass

        selected_peaks = np.zeros(self.lwPeaks.count(), dtype=int)
        selected_peaks[peak_indices] = 1

        # peak colors will be an array with one entry per peak:
        # no polymers: 0 - not selected, 1 - selected
        # polymers:    2 - not selected, 3 - selected
        # alternative colors: orange [1, .49, .16, 1.0], light red [1, .66, .66, 1.0]
        peak_colors = selected_peaks + polymer_peaks
        color_set = np.array(["#000000",   # black for unselected peaks without polymers
                              "#ffcc00",   # yellow for selected peaks without polymers
                              "#00c000",   # green for unselected peaks with polymers
                              "#ff0000"])  # red for selected peaks with polymers
        self.spectrum_peak_lines.set_color(color_set[peak_colors])
        self.spectrum_peak_lines.set_linewidth(1)

        # label the selection by masses
        for annotation in self.spectrum_axes.findobj(matplotlib.text.Annotation):
            annotation.remove()
        if self.btLabelPeaks.isChecked():
            for peak_id in peak_indices:
                self.spectrum_axes.annotate(s="{:.2f}".format(self._exp_mass_data.iloc[peak_id]["Average Mass"]),
                                            xy=(self._exp_mass_data.iloc[peak_id]["Average Mass"],
                                                self._exp_mass_data.iloc[peak_id]["Relative Abundance"]),
                                            xytext=(0, 5),
                                            textcoords="offset pixels",
                                            horizontalalignment="center",
                                            bbox=dict(facecolor="white", alpha=.75, linewidth=0))
        self.spectrum_canvas.draw()


    def display_selected_peaks(self, selected_peaks=None):
        """
        Update the spectrum, the list of masses and the result tree after the selection
        or parameters influencing the selection have changed.

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

        # show the first currently selected mass or the mass of the central peak of the delta series in the mass list
        if self.spectrum_picked_peak is not None:
            self.lwPeaks.scrollToItem(self.lwPeaks.item(self.spectrum_picked_peak))
        else:
            self.lwPeaks.scrollToItem(self.lwPeaks.item(selected_peaks[0]))

        # fill the single mass spin box with the currently selected mass
        try:
            self.sbSingleMass.setValue(float(self.lwPeaks.currentItem().text()))
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
            if self.chFilterStructureHits.isChecked():
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
                root_item.setText(0, "{:.2f}".format(self._exp_mass_data.loc[mass_index]["Average Mass"]))
                root_item.setTextAlignment(0, Qt.AlignRight)
                root_item.setText(1, "{:.1f}".format(self._exp_mass_data.loc[mass_index]["Relative Abundance"]))
                root_item.setTextAlignment(1, Qt.AlignRight)
                root_item.mass_index = mass_index

                if mass_index not in df_hit.index.levels[0]:
                    root_item.setBackground(0, QBrush(missing_color))
                    root_item.setBackground(1, QBrush(missing_color))
                else:
                    # generate child items, one per possible combination of modifications
                    for _, hit in df_hit.loc[mass_index].iterrows():
                        child_item = SortableTreeWidgetItem(root_item)
                        monomers = hit[:df_hit.columns.get_loc("Exp. Mass")].index

                        # monomer counts
                        for j, monomer in enumerate(monomers):
                            child_item.setText(2 + j, "{:.0f}".format(hit[monomer]))
                            child_item.setTextAlignment(2 + j, Qt.AlignHCenter)

                        # hit properties
                        for j, label in enumerate(["Exp. Mass", "Theo. Mass", "Da.", "ppm"]):
                            child_item.setText(len(monomers) + 2 + j, "{:.2f}".format(hit[label]))
                            child_item.setTextAlignment(len(monomers) + 2 + j, Qt.AlignRight)

                        if self.chFilterStructureHits.isChecked():
                            # polymer composition
                            sites = hit[df_hit.columns.get_loc("ppm")+1:-1].index
                            for j, site in enumerate(sites):
                                child_item.setText(len(monomers) + 6 + j, "{}".format(hit[site]))
                                child_item.setTextAlignment(len(monomers) + 6 + j, Qt.AlignHCenter)

                            # polymer abundance
                            child_item.setText(len(monomers) + 7 + j, "{:.2f}".format(hit["Abundance"]))
                            child_item.setTextAlignment(len(monomers) + 7 + j, Qt.AlignHCenter)


            self.twResults.expandAll()
            self.twResults.header().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.twResults.header().setStretchLastSection(False)
            self.twResults.setUpdatesEnabled(True)


    def save_settings(self):
        """
        Dump the current settings via pickle.

        :return: nothing
        """
        settings_filename = QFileDialog.getSaveFileName(self,
                                                        "Save settings",
                                                        self._path,
                                                        "ModFinder settings (*.mofi)")[0]
        self._path = os.path.split(settings_filename)[0]
        if settings_filename:
            if not settings_filename.endswith(".mofi"):
                settings_filename = settings_filename + ".mofi"

            # get monomers from table
            monomers = []
            for row_id in range(self.tbMonomers.rowCount()):
                is_used = self.tbMonomers.cellWidget(row_id, 0).findChild(QCheckBox).isChecked()
                name = self.tbMonomers.item(row_id, 1).text()
                composition = self.tbMonomers.item(row_id, 2).text()
                min_count = self.tbMonomers.cellWidget(row_id, 3).value()
                max_count = self.tbMonomers.cellWidget(row_id, 4).value()
                is_poly = self.tbMonomers.cellWidget(row_id, 5).findChild(QCheckBox).isChecked()
                monomers.append((is_used, name, composition, min_count, max_count, is_poly))

            # get polymers from table
            polymers = []
            for row_id in range(self.tbPolymers.rowCount()):
                is_used = self.tbPolymers.cellWidget(row_id, 0).findChild(QCheckBox).isChecked()
                name = self.tbPolymers.item(row_id, 1).text()
                composition = self.tbPolymers.item(row_id, 2).text()
                sites = self.tbPolymers.item(row_id, 3).text()
                abundance = self.tbPolymers.cellWidget(row_id, 4).value()
                polymers.append((is_used, name, composition, sites, abundance))

            settings = {"sequence": self.teSequence.toPlainText(),
                        "exp mass data": self._exp_mass_data,
                        "mass filename": self._mass_filename,
                        "disulfides": self.sbDisulfides.value(),
                        "pngase f": self.chPngase.isChecked(),
                        "tolerance value": self.sbTolerance.value(),
                        "tolerance flavor": self.cbTolerance.currentIndex(),
                        "monomers": monomers,
                        "polymers": polymers}
            with open(settings_filename, "wb") as f:
                pickle.dump(settings, f)


    def load_settings(self):
        settings_filename = QFileDialog.getOpenFileName(self,
                                                        "Load settings",
                                                        self._path,
                                                        "ModFinder settings (*.mofi)")[0]
        self._path = os.path.split(settings_filename)[0]
        if settings_filename:
            with open(settings_filename, "rb") as f:
                settings = pickle.load(f)

            self.teSequence.setText(settings["sequence"])
            if self.teSequence.toPlainText():
                self.calculate_protein_mass()
            self.sbDisulfides.setValue(settings["disulfides"])
            self.chPngase.setChecked(settings["pngase f"])
            self.calculate_protein_mass()

            self._mass_filename = settings["mass filename"]

            self._monomer_hits = None
            self._polymer_hits = None
            self.chFilterStructureHits.setEnabled(False)
            self.twResults.clear()
            if settings["exp mass data"] is not None:
                self._exp_mass_data = settings["exp mass data"]
                self.lwPeaks.blockSignals(True)
                self.lwPeaks.clear()
                self.lwPeaks.addItems(["{:.2f}".format(i) for i in self._exp_mass_data["Average Mass"]])
                self.lwPeaks.setCurrentRow(0)
                self.lwPeaks.blockSignals(False)
                self.draw_spectrum()
                self.spectrum_picked_peak = 0
                self.display_selected_peaks()

            self.cbTolerance.setCurrentIndex(settings["tolerance flavor"])
            self.sbTolerance.setValue(settings["tolerance value"])

            self.table_clear(self.tbMonomers)
            for row_id, (is_used, name, composition,
                         min_count, max_count, is_poly) in enumerate(settings["monomers"]):
                self._monomer_table_create_row(row_id, is_used, name, composition, min_count, max_count, is_poly)

            self.table_clear(self.tbPolymers)
            for row_id, (is_used, name, composition, sites, abundance) in enumerate(settings["polymers"]):
                self._polymer_table_create_row(row_id, is_used, name, composition, sites, abundance)

            self.calculate_mod_mass()


    def save_monomers(self):
        """
        Export the contents of the monomer table.

        :return: nothing
        """
        filename = QFileDialog.getSaveFileName(self,
                                               "Export monomers",
                                               self._path,
                                               "Comma-separated value (*.csv)")[0]
        self._path = os.path.split(filename)[0]
        if filename:
            if not filename.endswith(".csv"):
                filename += ".csv"
            df_monomers = pd.DataFrame(self.calculate_mod_mass(return_values="all"),
                                       columns=["Checked", "Name", "Composition", "Min", "Max", "Poly?"])
            try:
                df_monomers.to_csv(filename, index=False)
            except IOError:
                QMessageBox.critical(self, "Error", "Error when writing to " + filename + IOError.args)


    def load_monomers(self):
        """
        Import the contents of the monomer table.
        Columns labelled "Name" and "Composition" must exists in the input file.
        The following columns are optional and are filled with the indicated default values if missing:
        "Checked" (False), "Min" (0), "Max" (-1, i.e., inf) and "Poly?" (False).

        :return: nothing
        """
        filename = QFileDialog.getOpenFileName(self,
                                               "Import monomers",
                                               self._path,
                                               "Excel files (*.xlsx *.xls);; CSV files (*.csv *.txt)")[0]
        self._path = os.path.split(filename)[0]
        ext = os.path.splitext(filename)[1]
        if ext in [".xls", ".xlsx"]:
            df_monomers = pd.read_excel(filename)
        elif ext in [".txt", ".csv"]:
            df_monomers = pd.read_csv(filename)
        else:
            return

        if "Name" not in df_monomers.columns:
            QMessageBox.warning(self, "Warning", "Column 'Name' missing in monomer input. No data imported.")
            return
        if "Composition" not in df_monomers.columns:
            QMessageBox.warning(self, "Warning", "Column 'Composition' missing in monomer input. No data imported.")
            return

        if "Checked" not in df_monomers.columns:
            df_monomers["Checked"] = False
        if "Min" not in df_monomers.columns:
            df_monomers["Min"] = 0
        if "Max" not in df_monomers.columns:
            df_monomers["Max"] = -1
        if "Poly?" not in df_monomers.columns:
            df_monomers["Poly?"] = False

        self.table_clear(self.tbMonomers)
        for row_id, data in df_monomers.iterrows():
            self._monomer_table_create_row(row_id, data["Checked"], data["Name"], data["Composition"],
                                           data["Min"], data["Max"], data["Poly?"])
        self.calculate_mod_mass()


    def load_default_monomers(self):
        """
        Load a default monomer library.

        :return: nothing
        """
        self.table_clear(self.tbMonomers)
        library = self.sender().text()
        row_id = 0

        for name, data in configure.default_monomer_libraries[library].items():
            self._monomer_table_create_row(row_id,
                                           name=name,
                                           composition=data["composition"],
                                           min_count=int(data["min"]),
                                           max_count=int(data["max"]),
                                           part_of_polymer=distutils.util.strtobool(data["glycan"]))
            row_id += 1


    def save_polymers(self):
        """
        Export the contents of the polymer table.

        :return: nothing
        """
        filename = QFileDialog.getSaveFileName(self,
                                               "Export polymers",
                                               self._path,
                                               "Comma-separated value (*.csv)")[0]
        self._path = os.path.split(filename)[0]
        if filename:
            if not filename.endswith(".csv"):
                filename += ".csv"

            polymer_data = []
            for row_id in range(self.tbPolymers.rowCount()):
                ch = self.tbPolymers.cellWidget(row_id, 0).findChild(QCheckBox)
                name = self.tbPolymers.item(row_id, 1).text()
                composition = self.tbPolymers.item(row_id, 2).text()
                sites = self.tbPolymers.item(row_id, 3).text()
                abundance = self.tbPolymers.cellWidget(row_id, 4).value()

                polymer_data.append((ch.isChecked(), name, composition, sites, abundance))

            df_polymers = pd.DataFrame(polymer_data,
                                       columns=["Checked", "Name", "Composition", "Sites", "Abundance"])
            try:
                df_polymers.to_csv(filename, index=False)
            except IOError:
                QMessageBox.critical(self, "Error", "Error when writing to " + filename + IOError.args)


    def load_polymers(self):
        """
        Import the contents of the polymer table.
        Columns labelled "Name", "Composition" and "Sites" must exists in the input file.
        The following columns are optional and are filled with the indicated default values if missing:
        "Checked" (False) and "Abundance" (0.0).

        :return: nothing
        """
        filename = QFileDialog.getOpenFileName(self,
                                               "Import polymers",
                                               self._path,
                                               "Excel files (*.xlsx *.xls);; CSV files (*.csv *.txt)")[0]
        self._path = os.path.split(filename)[0]
        ext = os.path.splitext(filename)[1]
        if ext in [".xls", ".xlsx"]:
            df_polymers = pd.read_excel(filename)
        elif ext in [".txt", ".csv"]:
            df_polymers = pd.read_csv(filename)
        else:
            return

        if "Name" not in df_polymers.columns:
            QMessageBox.warning(self, "Warning", "Column 'Name' missing in polymer input. No data imported.")
            return
        if "Composition" not in df_polymers.columns:
            QMessageBox.warning(self, "Warning", "Column 'Composition' missing in polymer input. No data imported.")
            return
        if "Sites" not in df_polymers.columns:
            QMessageBox.warning(self, "Warning", "Column 'Sites' missing in polymer input. No data imported.")
            return

        if "Checked" not in df_polymers.columns:
            df_polymers["Checked"] = False
        if "Abundance" not in df_polymers.columns:
            df_polymers["Abundance"] = 0.0

        self.table_clear(self.tbPolymers)
        for row_id, data in df_polymers.iterrows():
            self._polymer_table_create_row(row_id, data["Checked"], data["Name"], data["Composition"],
                                           data["Sites"], data["Abundance"])

    def load_default_polymers(self):
        """
        Load a default polymer library.

        :return: nothing
        """
        self.table_clear(self.tbPolymers)
        library = self.sender().text()
        row_id = 0

        for name, data in configure.default_polymer_libraries[library].items():
            self._polymer_table_create_row(row_id, name=name, composition=data["composition"], sites=data["sites"])
            row_id += 1


if __name__ == "__main__":
    app = QApplication(sys.argv)
    frame = MainWindow()
    frame.show()
    app.exec_()
