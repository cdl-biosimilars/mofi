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

import sys
import os
import re
import pickle
import time
import qtpy
from qtpy.QtWidgets import (QApplication, QMainWindow, QMenu, QActionGroup, QVBoxLayout, QTableWidgetItem, QCheckBox,
                            QMessageBox, QFileDialog, QTreeWidgetItem, QHeaderView, QSpinBox, QDoubleSpinBox,
                            QWidget, QHBoxLayout, QAction)
from qtpy.QtGui import QColor, QBrush

from qtpy.QtCore import Qt
AlignHCenter = qtpy.QtCore.Qt.AlignHCenter
AlignRight = qtpy.QtCore.Qt.AlignRight

import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 5000)

import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import (FigureCanvasQTAgg as FigureCanvas,
                                                NavigationToolbar2QT as NavigationToolbar)
from matplotlib.backend_bases import PickEvent


from matplotlib.figure import Figure

import configure
import mass_tools
import modification_search
import sequence_tools

from ModFinder_UI import Ui_ModFinder


class SortableTreeWidgetItem(QTreeWidgetItem):
    """
    A QTreeWidget which supports numerical sorting.
    """
    def __init__(self, parent=None):
        super(SortableTreeWidgetItem, self).__init__(parent)

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

        self.acAbout.triggered.connect(self.show_about)
        self.acHelp.triggered.connect(self.show_help)
        self.acLoadSettings.triggered.connect(self.load_settings)
        self.acOpenFasta.triggered.connect(self.read_fasta_file)
        self.acOpenPeaks.triggered.connect(self.read_mass_file)
        self.acQuit.triggered.connect(QApplication.instance().quit)
        self.acSaveAnnotation.triggered.connect(self.save_csv)
        self.acSaveSettings.triggered.connect(self.save_settings)

        self.btCalcCombinations.clicked.connect(self.sample_modifications)
        self.btClearMonomers.clicked.connect(lambda: self.table_clear(self.tbMonomers))
        self.btClearPolymers.clicked.connect(lambda: self.table_clear(self.tbPolymers))
        self.btDeleteRowMonomers.clicked.connect(lambda: self.table_delete_row(self.tbMonomers))
        self.btDeleteRowPolymers.clicked.connect(lambda: self.table_delete_row(self.tbPolymers))
        self.btInsertRowAboveMonomers.clicked.connect(lambda: self.table_insert_row(self.tbMonomers, above=True))
        self.btInsertRowAbovePolymers.clicked.connect(lambda: self.table_insert_row(self.tbPolymers, above=True))
        self.btInsertRowBelowMonomers.clicked.connect(lambda: self.table_insert_row(self.tbMonomers, above=False))
        self.btInsertRowBelowPolymers.clicked.connect(lambda: self.table_insert_row(self.tbPolymers, above=False))
        self.btLoadMonomers.clicked.connect(self.load_monomers)
        self.btLoadPolymers.clicked.connect(self.load_polymers)
        self.btSaveMonomers.clicked.connect(self.save_monomers)
        self.btSavePolymers.clicked.connect(self.save_polymers)
        self.btUpdateMass.clicked.connect(self.calculate_protein_mass)
        self.btShowPolymerHits.clicked.connect(lambda: self.set_result_tree(show_monomers=False))

        self.cbTolerance.activated.connect(self.choose_tolerance_units)

        self.chDelta1.clicked.connect(self.show_deltas1)
        self.chDelta2.clicked.connect(self.show_deltas2)
        self.chOnlyPolymerResults.clicked.connect(lambda: self.set_result_tree(show_monomers=True))
        self.chPngase.clicked.connect(self.calculate_protein_mass)

        self.lwPeaks.currentItemChanged.connect(lambda: self.set_result_tree(show_monomers=True))

        self.sbDelta1.valueChanged.connect(self.show_deltas1)
        self.sbDelta2.valueChanged.connect(self.show_deltas2)
        self.sbDeltaTolerance.valueChanged.connect(self.check_tolerance)
        self.sbDisulfides.valueChanged.connect(self.calculate_protein_mass)

        self.teSequence.textChanged.connect(
            lambda: self.teSequence.setStyleSheet("QTextEdit { background-color: rgb(255, 225, 225) }"))

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
            self.agSelectMassSet.addAction(ac_select_set)
            set_id += 1
        # noinspection PyUnresolvedReferences
        self.agSelectMassSet.triggered.connect(self.choose_mass_set)

        # set up the plot
        layout = QVBoxLayout(self.spectrumView)
        self.fig = Figure(dpi=100, frameon=False, tight_layout={"pad": 0}, edgecolor="white")
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.spectrumView)
        layout.addWidget(self.canvas)
        # layout.addWidget(NavigationToolbar(self.canvas, self.spectrumView))  # TODO: vertical?
        self.peak_lines = None

        # init the monomer table and associated buttons
        self.tbMonomers.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        for col, width in [(0, 25), (2, 130), (3, 45), (4, 45), (5, 40)]:
            self.tbMonomers.setColumnWidth(col, width)
        self.tbMonomers.verticalHeader().setSectionResizeMode(QHeaderView.Fixed)
        self.tbMonomers.verticalHeader().setDefaultSectionSize(22)

        menu = QMenu()
        for library in configure.default_monomer_libraries:
            menu.addAction(library, self.load_default_monomers)
        self.btDefaultModsMonomers.setMenu(menu)

        # init the polymer table and associated buttons
        self.tbPolymers.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        for col, width in [(0, 25), (2, 130), (3, 80), (4, 60)]:
            self.tbPolymers.setColumnWidth(col, width)
        self.tbPolymers.verticalHeader().setSectionResizeMode(QHeaderView.Fixed)
        self.tbPolymers.verticalHeader().setDefaultSectionSize(22)

        menu = QMenu()
        for library in configure.default_polymer_libraries:
            menu.addAction(library, self.load_default_polymers)
        self.btDefaultModsPolymers.setMenu(menu)

        self.choose_tolerance_units()

        # initialize private members
        # these variables completely describe the state of the program
        self._monomer_hits = None  # results from the monomer search
        self._delta_lines_1 = None  # delta 1 lines in the mass spectrum
        self._delta_lines_2 = None  # delta 2 lines in the mass spectrum
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

        score_spinbox = QDoubleSpinBox()
        score_spinbox.setMinimum(0)
        score_spinbox.setMaximum(100)
        score_spinbox.setSingleStep(.1)
        score_spinbox.setFrame(False)
        score_spinbox.setValue(abundance)
        score_spinbox.setStyleSheet(configure.double_spin_box_flat_style)
        self.tbPolymers.setCellWidget(row_id, 4, score_spinbox)


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


    def load_default_monomers(self):
        """
        Load a default monomer library.

        :return: nothing
        """
        self.table_clear(self.tbMonomers)
        library = self.sender().text()
        row_id = 0

        for name, data in configure.default_monomer_libraries[library].items():
            self._monomer_table_create_row(row_id, name=name, composition=data["composition"],
                                           min_count=int(data["min"]), max_count=int(data["max"]), part_of_polymer=True)
            row_id += 1


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


    def check_tolerance(self):
        """
        Update delta lines in the plot if the corresponding tolerance value has been changed.

        :return: nothing
        """
        self.show_deltas1()
        self.show_deltas1()


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


    def read_fasta_file(self):
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


    def read_mass_file(self):
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
            self.lwPeaks.clear()
            self.twResults.clear()
            self._monomer_hits = None
            self._polymer_hits = None
            self.lwPeaks.addItems(["{:.2f}".format(i) for i in self._exp_mass_data["Average Mass"]])
            self.lwPeaks.setCurrentRow(0)
            self.draw_naked_plot()


    def save_csv(self):
        """
        Write the results from a combinatorial search to a CSV file.

        :return: nothing
        """
        outfilename = QFileDialog.getSaveFileName(self,
                                                  "Save Annotation",
                                                  self._path,
                                                  "MoFi Annotation (*.csv)")[0]
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

        # prepare list of polymers for search stage 2
        polymers = {}
        for row_id in range(self.tbPolymers.rowCount()):
            if self.tbPolymers.cellWidget(row_id, 0).findChild(QCheckBox).isChecked():
                name = self.tbPolymers.item(row_id, 1).text()
                composition = self.tbPolymers.item(row_id, 2).text()
                sites = self.tbPolymers.item(row_id, 3).text()
                score = self.tbPolymers.cellWidget(row_id, 4).value()
                polymers[name] = (composition, sites, score)

        if polymers:  # dict contains at least one entry
            df_polymers = pd.DataFrame.from_dict(polymers, orient="index")
            df_polymers.columns = ["Composition", "Sites", "Score"]
            monomers_in_library = modification_search.get_monomers_from_library(df_polymers)
        else:  # dict remained empty, i.e., there are no polymers
            monomers_in_library = []
            df_polymers = pd.DataFrame()

        if sorted(monomers_for_polymer_search) != sorted(monomers_in_library):
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

        print("PERFORMING COMBINATORIAL SEARCH...")
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
            explained_mass=explained_mass)

        print("Monomer search DONE!")
        # the modification search was not successful
        if self._monomer_hits is None:
            QMessageBox.critical(self, "Error", "Combinatorial search was unsuccessful.")

        # add the minimum monomer counts to the result data frame
        for name, _, min_count, _, _ in monomers:
            self._monomer_hits[name] += min_count

        # stage 2: polymer search
        if polymers:
            self._polymer_hits = modification_search.find_polymers(self._monomer_hits,
                                                                   glycan_library=df_polymers,
                                                                   monomers=monomers_for_polymer_search)
        print("Polymer search DONE!")
        # self._monomer_hits.to_csv("csv/monomer_hits.csv")
        # self._polymer_hits.to_csv("csv/polymer_hits.csv")

        if self._monomer_hits is not None:
            self.set_result_tree(show_monomers=True)
            # self.draw_annotated_plot()  TODO


    def pick_peak(self, event):
        """
        If the user clicks on a peak of the spectrum, color this peak in red and select
        the corresponding mass in the list widget.

        :param event: PickEvent from the canvas
        :return: nothing
        """
        peak_index = event.ind[0]
        self.lwPeaks.setCurrentRow(peak_index)

        selected_peaks = np.zeros(len(self._exp_mass_data), dtype=int)
        normal_selected_color = np.array([[0, 0, 0, 1.0], [1, 0, 0, 1.0]])
        selected_peaks[peak_index] = 1
        event.artist.set_color(normal_selected_color[selected_peaks])
        self.canvas.draw()


    def draw_naked_plot(self):
        """
        Show a mass spectrum after loading a peak list (i.e., without annotations).

        :return: nothing
        """

        self.fig.clear()
        axes = self.fig.add_subplot(111)
        self.peak_lines = axes.vlines(x=self._exp_mass_data["Average Mass"],
                                      ymin=0,
                                      ymax=self._exp_mass_data["Relative Abundance"],
                                      linewidth=1,
                                      picker=5)
        axes.set_ylim(0, 115)
        axes.set_xlabel("Mass (Da)")
        axes.set_ylabel("Relative Abundance (%)")
        axes.spines["right"].set_visible(False)
        axes.spines["top"].set_visible(False)
        axes.yaxis.set_ticks_position("left")
        axes.xaxis.set_ticks_position("bottom")
        axes.tick_params(direction="out")
        self.fig.canvas.mpl_connect("pick_event", self.pick_peak)
        self.show_deltas1()
        self.show_deltas2()
        self.canvas.draw()
        print(self.peak_lines)


    def draw_annotated_plot(self):
        """
        Show an annotated mass spectrum after the modification search.

        :return: nothing
        """
        self.fig.clear()
        axes = self.fig.add_subplot(111)
        axes.vlines(self._exp_mass_data["Average Mass"],
                    0,
                    self._exp_mass_data["Relative Abundance"],
                    lw=1)
        axes.axhline(0, c="k", lw=1)
        axes.set_ylim(0, 115)
        axes.set_xlabel("Mass (Da.)")
        axes.set_ylabel("Relative Abundance (%)")
        axes.set_title(self._mass_filename, fontsize=6)
        axes.spines["right"].set_visible(False)
        axes.spines["top"].set_visible(False)
        axes.yaxis.set_ticks_position("left")
        axes.xaxis.set_ticks_position("bottom")
        axes.tick_params(direction="out")

        for i in self._monomer_hits.index.levels[0]:
            x = self._exp_mass_data.loc[i, "Average Mass"]
            y = self._exp_mass_data.loc[i, "Relative Abundance"]
            first_hit = self._monomer_hits.loc[i].iloc[0]
            s = "%.2f" % x
            s += "(%.1f %s)\n" % (first_hit["ppm"], "ppm")
            axes.annotate(s,
                          xy=(x, y),
                          fontsize=8,
                          rotation=0,
                          va="bottom",
                          ha="left",
                          bbox=dict(fc=(0.9, 0.9, 0.9, 1), edgecolor="none", pad=0))
            axes.plot([x], [y], ".", c="orange")

        self.show_deltas1()
        self.show_deltas2()
        self.canvas.draw()


    def set_result_tree(self, show_monomers=True):
        try:  # to set the value of the single mass calculation spinbox
            self.sbSingleMass.setValue(float(self.lwPeaks.currentItem().text()))
        except AttributeError:  # occurs when second peak file is loaded
            pass

        try:  # to highlight the peak that belongs to the selected mass
            # noinspection PyTypeChecker
            self.pick_peak(PickEvent("", None, None, self.peak_lines, ind=[self.lwPeaks.currentRow()]))
        except AttributeError:  # occurs when peak file is loaded
            pass

        if self._monomer_hits is not None:
            self.twResults.clear()
            self.twResults.setUpdatesEnabled(False)

            missing_color = QColor(255, 185, 200)
            if show_monomers:
                indices = [self.lwPeaks.currentRow()]
                df_hit = self._monomer_hits
                self.chOnlyPolymerResults.setEnabled(True)
            else:
                indices = self._polymer_hits.index.levels[0]
                df_hit = self._polymer_hits
                self.chOnlyPolymerResults.setEnabled(False)

            # set column headers
            header_labels = ["Exp. Mass", "%"]
            header_labels.extend(df_hit.columns)
            self.twResults.setColumnCount(len(header_labels))
            self.twResults.setHeaderLabels(header_labels)

            for mass_index in indices:
                # generate root item (experimental mass, relative abundance)
                root_item = SortableTreeWidgetItem(self.twResults)
                root_item.setText(0, "{:.2f}".format(self._exp_mass_data.loc[mass_index]["Average Mass"]))
                root_item.setTextAlignment(0, AlignRight)
                root_item.setText(1, "{:.1f}".format(self._exp_mass_data.loc[mass_index]["Relative Abundance"]))
                root_item.setTextAlignment(1, AlignRight)

                if mass_index not in df_hit.index.levels[0]:
                    root_item.setBackground(0, QBrush(missing_color))
                    root_item.setBackground(1, QBrush(missing_color))
                else:
                    # generate child items, one per possible combination of modifications
                    for _, hit in df_hit.loc[mass_index].iterrows():
                        child_item = SortableTreeWidgetItem(root_item)

                        monomers = hit[:df_hit.columns.get_loc("Exp. Mass")].index
                        for j, monomer in enumerate(monomers):
                            child_item.setText(2 + j, "{:.0f}".format(hit[monomer]))
                            child_item.setTextAlignment(2 + j, AlignHCenter)

                        for j, label in enumerate(["Exp. Mass", "Theo. Mass", "Da.", "ppm"]):
                            child_item.setText(len(monomers) + 2 + j, "{:.2f}".format(hit[label]))
                            child_item.setTextAlignment(len(monomers) + 2 + j, AlignRight)

                        sites = hit[df_hit.columns.get_loc("ppm")+1:].index
                        for j, site in enumerate(sites):
                            child_item.setText(len(monomers) + 6 + j, "{}".format(hit[site]))
                            child_item.setTextAlignment(len(monomers) + 6 + j, AlignHCenter)


            self.twResults.expandAll()
            self.twResults.header().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.twResults.header().setStretchLastSection(False)
            self.twResults.setUpdatesEnabled(True)


    # def selection_plot(self):
    #     # COMMENT: Apparently, whis function was intended to annotate only peaks selected in the tree
    #     self.fig.clear()
    #     axes = self.fig.add_subplot(111)
    #     massdf = self._exp_mass_data
    #     axes.vlines(massdf['Average Mass'], 0,
    #     massdf['Relative Abundance'], lw=1)
    #     axes.axhline(0, c='k', lw=1)
    #     axes.set_ylim(0,115)
    #     axes.set_xlabel('Mass (Da.)')
    #     axes.set_ylabel('Relative Abundance (%)')
    #     axes.set_title(self._mass_filename, fontsize=6)
    #     axes.spines['right'].set_visible(False)
    #     axes.spines['top'].set_visible(False)
    #     axes.yaxis.set_ticks_position('left')
    #     axes.xaxis.set_ticks_position('bottom')
    #     axes.tick_params(direction='out')
    #     selection = self.twResults.selectedItems()
    #     for sel in selection:
    #         mass = round(float(sel.parent().text(0)),4)
    #         r_abundance = round(float(sel.parent().text(1)),4)
    #         x = mass
    #         y = r_abundance
    #         s = '%.2f' % x
    #         seldict = {}
    #         for i in range(2,len(self._hlabels)):
    #             try:
    #                 sel_value = int(sel.text(i))
    #             except:
    #                 sel_value = float(sel.text(i))
    #             seldict[self._hlabels[i]] = sel_value
    #         seldata= pd.Series(seldict)
    #         mods = self._hlabels[2:-4]
    #         s += '(%.1f %s)\n' % (seldata[str(self.cbTolerance.currentText())], self.cbTolerance.currentText())
    #         s += '\n'.join(['%s:%d'%(m,seldata[m]) for m in mods if seldata[m] != 0])
    #         axes.annotate(s, xy=(x,y),fontsize=8, ha='center',
    #                      backgroundcolor='lightgreen')
    #     self.show_deltas1()
    #     self.show_deltas2()
    #     self.canvas.draw()
    #     self.tabWidget.setCurrentIndex(0)


    def show_deltas1(self):  # TODO: much too slow!!!
        """
        Show lines between peaks whose masses differ by the value of the Delta1 spinbox (+- tolerance).

        :return: nothing
        """
        if self._exp_mass_data is not None:
            delta_masses = mass_tools.find_delta_masses(self._exp_mass_data,
                                                        self.sbDelta1.value(),
                                                        self.sbDeltaTolerance.value() / 2)
            delta_plot = self.fig.add_subplot(111)
            y = [mass[5] for mass in delta_masses]
            x_start = [mass[3] for mass in delta_masses]
            x_stop = [mass[4] for mass in delta_masses]
            if self._delta_lines_1 is not None:
                try:
                    self._delta_lines_1.remove()
                except ValueError:
                    pass
            if self.chDelta1.isChecked():
                self._delta_lines_1 = delta_plot.hlines(y, x_start, x_stop, lw=1.5, color="green")
            self.canvas.draw()


    def show_deltas2(self):
        """
            Show lines between peaks whose masses differ by the value of the Delta1 spinbox (+- tolerance).

            :return: nothing
            """
        if self._exp_mass_data is not None:
            delta_masses = mass_tools.find_delta_masses(self._exp_mass_data,
                                                        self.sbDelta2.value(),
                                                        self.sbDeltaTolerance.value() / 2)
            delta_plot = self.fig.add_subplot(111)
            y = [mass[5] * .6 for mass in delta_masses]
            x_start = [mass[3] for mass in delta_masses]
            x_stop = [mass[4] for mass in delta_masses]
            if self._delta_lines_2 is not None:
                try:
                    self._delta_lines_2.remove()
                except ValueError:
                    pass
            if self.chDelta2.isChecked():
                self._delta_lines_2 = delta_plot.hlines(y, x_start, x_stop, lw=1.5, color="purple")
            self.canvas.draw()


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
                score = self.tbPolymers.cellWidget(row_id, 4).value()
                polymers.append((is_used, name, composition, sites, score))

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
            self.twResults.clear()
            if settings["exp mass data"] is not None:
                self._exp_mass_data = settings["exp mass data"]
                self.lwPeaks.clear()
                self.lwPeaks.addItems(["{:.2f}".format(i) for i in self._exp_mass_data["Average Mass"]])
                self.lwPeaks.setCurrentRow(0)
                self.sbSingleMass.setValue(float(self.lwPeaks.currentItem().text()))
                self.draw_naked_plot()

            self.cbTolerance.setCurrentIndex(settings["tolerance flavor"])
            self.sbTolerance.setValue(settings["tolerance value"])

            self.table_clear(self.tbMonomers)
            for row_id, (is_used, name, composition,
                         min_count, max_count, is_poly) in enumerate(settings["monomers"]):
                self._monomer_table_create_row(row_id, is_used, name, composition, min_count, max_count, is_poly)

            self.table_clear(self.tbPolymers)
            for row_id, (is_used, name, composition, sites, score) in enumerate(settings["polymers"]):
                self._polymer_table_create_row(row_id, is_used, name, composition, sites, score)

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
                score = self.tbPolymers.cellWidget(row_id, 4).value()

                polymer_data.append((ch.isChecked(), name, composition, sites, score))

            df_polymers = pd.DataFrame(polymer_data,
                                       columns=["Checked", "Name", "Composition", "Sites", "Score"])
            try:
                df_polymers.to_csv(filename, index=False)
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


    def load_polymers(self):
        """
        Import the contents of the polymer table.
        Columns labelled "Name", "Composition" and "Sites" must exists in the input file.
        The following columns are optional and are filled with the indicated default values if missing:
        "Checked" (False) and "Score" (0.0).

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
        if "Score" not in df_polymers.columns:
            df_polymers["Score"] = 0.0

        self.table_clear(self.tbPolymers)
        for row_id, data in df_polymers.iterrows():
            self._polymer_table_create_row(row_id, data["Checked"], data["Name"], data["Composition"],
                                           data["Sites"], data["Score"])


if __name__ == "__main__":
    app = QApplication(sys.argv)
    frame = MainWindow()
    frame.show()
    app.exec_()
