#! python

# Some commands for UI converting and freezing that need to be documented:

# Converting .ui to .py
# pyside-uic ModFinder_UI.ui -o ModFinder_UI.py

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
import qtpy
from qtpy.QtWidgets import (QApplication, QMainWindow, QMenu, QActionGroup, QVBoxLayout, QTableWidgetItem, QCheckBox,
                            QMessageBox, QFileDialog, QTreeWidgetItem, QHeaderView, QSpinBox, QDoubleSpinBox,
                            QWidget, QHBoxLayout)
from qtpy.QtGui import QColor, QBrush

from qtpy.QtCore import Qt
AlignHCenter = qtpy.QtCore.Qt.AlignHCenter
AlignRight = qtpy.QtCore.Qt.AlignRight

import pandas as pd
pd.set_option('display.max_rows', 5000)

import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import (FigureCanvasQTAgg as FigureCanvas,
                                                NavigationToolbar2QT as NavigationToolbar)


from matplotlib.figure import Figure


import configure
import glyco_tools
import mass_tools
import modification_search
import output_tools
import sequence_tools

from ModFinder_UI import Ui_ModFinder


class SortableTreeWidgetItem(QTreeWidgetItem):
    """ A QTreeWidget which implements numerical sorting """
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

        self.btCalcCombinations.clicked.connect(self.run_modification_search)
        self.btClearMonomers.clicked.connect(lambda: self.table_clear(self.tbMonomers))
        self.btClearPolymers.clicked.connect(lambda: self.table_clear(self.tbPolymers))
        self.btDeleteRowMonomers.clicked.connect(lambda: self.table_delete_row(self.tbMonomers))
        self.btDeleteRowPolymers.clicked.connect(lambda: self.table_delete_row(self.tbPolymers))
        self.btInsertRowAboveMonomers.clicked.connect(lambda: self.table_insert_row(self.tbMonomers, above=True))
        self.btInsertRowAbovePolymers.clicked.connect(lambda: self.table_insert_row(self.tbPolymers, above=True))
        self.btInsertRowBelowMonomers.clicked.connect(lambda: self.table_insert_row(self.tbMonomers, above=False))
        self.btInsertRowBelowPolymers.clicked.connect(lambda: self.table_insert_row(self.tbPolymers, above=False))
        # self.btLoadMods.clicked.connect(self.read_nglycan_file)
        self.btUpdateMass.clicked.connect(self.calculate_disulfides_and_protein_mass)

        self.cbTolerance.activated.connect(self.choose_tolerance_units)

        self.chDelta1.clicked.connect(self.show_deltas1)
        self.chDelta2.clicked.connect(self.show_deltas2)
        self.chPngase.clicked.connect(self.calculate_protein_mass)
        self.chRemoveUnannotated.clicked.connect(self.set_result_tree)

        self.lwPeaks.currentItemChanged.connect(
            lambda: self.sbSingleMass.setValue(float(self.lwPeaks.currentItem().text())))

        self.sbDelta1.valueChanged.connect(self.show_deltas1)
        self.sbDelta2.valueChanged.connect(self.show_deltas2)
        self.sbDeltaTolerance.valueChanged.connect(self.check_tolerance)
        self.sbDisulfides.valueChanged.connect(self.calculate_protein_mass)

        self.teSequence.textChanged.connect(
            lambda: self.teSequence.setStyleSheet("QTextEdit { background-color: rgb(255, 225, 225) }"))

        # iterators for the glycan checkboxes and associated spin boxes
        # self._glycan_checkboxes = [self.chHex, self.chHexnac, self.chNeu5ac, self.chNeu5gc,
        #                            self.chFuc, self.chPent, self.chNcore, self.chOcore]
        # self._glycan_min_spinboxes = [self.sbHexMin, self.sbHexnacMin, self.sbNeu5acMin, self.sbNeu5gcMin,
        #                               self.sbFucMin, self.sbPentMin, self.sbNcoreMin, self.sbOcoreMin]
        # self._glycan_max_spinboxes = [self.sbHexMax, self.sbHexnacMax, self.sbNeu5acMax, self.sbNeu5gcMax,
        #                               self.sbFucMax, self.sbPentMax, self.sbNcoreMax, self.sbOcoreMax]

        # group the mass set selectors
        self.agAverage = QActionGroup(self.menuAtomicMasses)
        self.agAverage.addAction(self.acAverageIupac)
        self.agAverage.addAction(self.acAverageZhang)
        self.agAverage.addAction(self.acMonoisotopic)
        # noinspection PyUnresolvedReferences
        self.agAverage.triggered.connect(self.choose_mass_set)

        # set up the plot
        layout = QVBoxLayout(self.spectrumView)
        self.fig = Figure((5.0, 3.0), dpi=100, frameon=False, tight_layout=True, edgecolor="white")
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.spectrumView)
        layout.addWidget(self.canvas)
        layout.addWidget(NavigationToolbar(self.canvas, self.spectrumView))  # TODO: vertical?

        # init the monomer table and associated buttons
        menu = QMenu()
        menu.addAction("Monosaccharides", lambda: self.load_default_monomers(monomers="monosaccharides"))
        menu.addAction("C-terminal lysines", lambda: self.load_default_monomers(monomers="lysines"))
        self.btDefaultModsMonomers.setMenu(menu)

        self.tbMonomers.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        for col, width in [(0, 25), (2, 130), (3, 45), (4, 45), (5, 40)]:
            self.tbMonomers.setColumnWidth(col, width)
        self.tbMonomers.verticalHeader().setSectionResizeMode(QHeaderView.Fixed)
        self.tbMonomers.verticalHeader().setDefaultSectionSize(22)

        # init the polymer table and associated buttons
        menu = QMenu()
        menu.addAction("Typical mAB glycans", lambda: self.load_default_polymers("mAB glycans"))
        self.btDefaultModsPolymers.setMenu(menu)
        next step: modification search should work again with monomers from table!!!
        self.tbPolymers.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        for col, width in [(0, 25), (2, 130), (3, 80), (4, 60)]:
            self.tbPolymers.setColumnWidth(col, width)
        self.tbPolymers.verticalHeader().setSectionResizeMode(QHeaderView.Fixed)
        self.tbPolymers.verticalHeader().setDefaultSectionSize(22)

        self.choose_tolerance_units()

        # initialize private members
        # these variables completely describe the state of the program
        self._all_hits = None  # results from the modification search
        self._delta_lines_1 = None  # delta 1 lines in the mass spectrum
        self._delta_lines_2 = None  # delta 2 lines in the mass spectrum
        self._exp_mass_data = None  # pandas dataframe containing the contents of the mass file
        self._known_mods_mass = 0  # mass of known modification
        self._mass_filename = None  # name of the mass file
        self._mass_set = configure.default_masses  # currently used atomic masses
        self._nglycans = None  # list of N-glycan modifications
        self._nglycans_data = None  # pandas dataframe containing the N-glycan library
        self._path = configure.path  # last path selected in a file chooser dialog
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
        self.tbMonomers.setCellWidget(row_id, 3, min_spinbox)

        max_spinbox = QSpinBox()
        max_spinbox.setMinimum(-1)
        max_spinbox.setSpecialValueText("inf")
        max_spinbox.setFrame(False)
        max_spinbox.setValue(max_count)
        self.tbMonomers.setCellWidget(row_id, 4, max_spinbox)

        monomer_checkbox = QCheckBox()
        monomer_checkbox.setChecked(part_of_polymer)
        monomer_ch_widget = QWidget()
        monomer_ch_layout = QHBoxLayout(monomer_ch_widget)
        monomer_ch_layout.addWidget(monomer_checkbox)
        monomer_ch_layout.setAlignment(Qt.AlignCenter)
        monomer_ch_layout.setContentsMargins(0, 0, 0, 0)
        self.tbMonomers.setCellWidget(row_id, 5, monomer_ch_widget)


    def _polymer_table_create_row(self, row_id, active=True, name="", composition="", site="", abundance=0.0):
        """
        Create a new row in the table of modifications.
        This row describes a modification with the given parameters.

        Private function, to be called by general row insertion functions and
        loaders for modification libraries.

        :param row_id: Row index passed to QTableWidget.insertRow()
        :param active: true if the polymer should be used in the combinatorial search
        :param name: name of the modification (str)
        :param composition: composition or mass of the modification (str)
        :param site: identifier for the modification site; empty means any site (str)
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

        self.tbPolymers.setItem(row_id, 3, QTableWidgetItem(site))

        abundance_spinbox = QDoubleSpinBox()
        abundance_spinbox.setMinimum(0)
        abundance_spinbox.setMaximum(100)
        abundance_spinbox.setSingleStep(.1)
        abundance_spinbox.setFrame(False)
        abundance_spinbox.setValue(abundance)
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


    def table_clear(self, table_widget):
        """
        Delete all rows in the table of modifications.

        :param table_widget: the QTableWidget to modify
        :return: nothing
        """
        while table_widget.rowCount() > 0:
            table_widget.removeRow(0)


    def table_delete_row(self, table_widget):
        """
        Delete selected rows in the table of modifications.

        :param table_widget: the QTableWidget to modify
        :return: nothing
        """
        if table_widget.selectionModel().selectedRows():
            for i in table_widget.selectionModel().selectedRows()[::-1]:
                table_widget.removeRow(i.row())


    def load_default_monomers(self, monomers="monosaccharides"):
        """
        Create a default modification: 0 to 2 C-terminal lysines

        :param monomers: specifies which set of monomers should be loaded
        :return: nothing
        """
        self.table_clear(self.tbMonomers)

        if monomers == "monosaccharides":
            for row_id, m in enumerate(glyco_tools.monosaccharides):
                self._monomer_table_create_row(row_id, name=m, composition=str(glyco_tools.glycan_formula[m]),
                                               min_count=0, max_count=-1, part_of_polymer=True)
        elif monomers == "lysines":
            lysine = mass_tools.Formula(sequence_tools.amino_acid_compositions["K"])
            self._monomer_table_create_row(0, name="C-terminal Lys", composition=str(lysine),
                                           min_count=0, max_count=2, part_of_polymer=False)
        else:
            pass


    def load_default_polymers(self, polymers="mAB glycans"):
        """
        Create a default modification: Typical mAB glycans

        :param polymers: specifies which set of polymers should be loaded
        :return: nothing
        """
        self.table_clear(self.tbPolymers)

        if polymers == "mAB glycans":
            last_row = 0
            for name, mass, max_count in glyco_tools.glycanlist_generator(glyco_tools.fc_glycans,
                                                                          mono_masses=self.acMonoisotopic.isChecked()):
                self._polymer_table_create_row(last_row, name=name, composition=mass, site="ch_A, ch_B")
                last_row += 1
        else:
            pass


    def show_about(self):
        """
        Show the about dialog.

        :return: nothing
        """
        v = "ModFinder version: %s" % configure.version
        r = configure.rights
        c = "Contact: %s" % configure.contact
        QMessageBox.about(self, "About ModFinder", "%s\n%s\n%s" % (v, r, c))


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
            contents and maximum value of self.disfBox to half the number of cysteines in the sequence

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
                self.calculate_disulfides_and_protein_mass()
            else:
                QMessageBox.warning(self, "Error loading sequence", "Not a valid file format: %s" % ext)


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
                QMessageBox.warning(self, "Error loading mass file", "Not a valid file format: %s" % ext)
                return

            self._exp_mass_data = mass_data
            self.lwPeaks.clear()
            self.twResults.clear()
            self.fig.clear()
            self.lwPeaks.addItems(["{:.2f}".format(i) for i in self._exp_mass_data["Average Mass"]])
            self.lwPeaks.setCurrentRow(0)
            self.sbSingleMass.setValue(float(self.lwPeaks.currentItem().text()))
            self.draw_naked_plot()
            self._all_hits = None


    def read_nglycan_file(self):
        """
        Loads an N-glycan library and stores the data in self._nglycans_data and self._nglycans.

        :return: none
        """
        nglycan_filename = QFileDialog.getOpenFileName(self,
                                                       "Load Nglycans",
                                                       self._path,
                                                       "Nglycan Library (*.ngl)")[0]
        self._path = os.path.split(nglycan_filename)[0]
        if nglycan_filename:
            self._nglycans_data = pd.read_table(nglycan_filename)
            self._nglycans_data["Name"] = self._nglycans_data["Name"].astype(str)
            self._nglycans_data["Nr."] = self._nglycans_data["Nr."].astype(int)
            self._nglycans = glyco_tools.dataframe_to_modlist(self._nglycans_data)
            self.nglycanCheck.show()
            self.nglycanCheck.setChecked(True)


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
            for ch, sp_min, sp_max in zip(self._glycan_checkboxes,
                                          self._glycan_min_spinboxes,
                                          self._glycan_max_spinboxes):
                if ch.isChecked():
                    parameters.append((ch.text() + "_min", sp_min.value()))
                    parameters.append((ch.text() + "_max", sp_max.value()))

            mindex = self._all_hits.reset_index(level=0)["Massindex"]
            ra = mindex.map(self._exp_mass_data["Relative Abundance"]).reset_index()
            self._all_hits["Relative Abundance"] = list(ra["Massindex"])

            try:
                output_tools.write_hits_to_csv(self._all_hits, outfilename, parameters)
            except IOError:
                QMessageBox.warning(self, "Warning", "Permission denied for %s" % outfilename)


    def choose_mass_set(self):
        """
        Chooses the set of atomic weights to use for mass calculations

        :return: nothing
        """
        if self.acAverageIupac.isChecked():
            self._mass_set = "AtomsAverageIUPAC"
            configure.set_average_masses(self._mass_set)
        elif self.acAverageZhang.isChecked():
            self._mass_set = "AtomsAverageZhang"
            configure.set_average_masses(self._mass_set)
        else:
            self._mass_set = "AtomsMonoisotopic"

        if self.teSequence.toPlainText():
            self.calculate_protein_mass()
            self.calculate_mod_mass()


    def calculate_disulfides_and_protein_mass(self):
        """
        Calculates the protein mass and composition, updates the disulfide spinbox and
        recalculates the protein mass (since the disulfides change it).

        :return: nothing
        """
        if not self.calculate_protein_mass() is None:
            self.sbDisulfides.setEnabled(True)
            self.chPngase.setEnabled(True)
            self.sbDisulfides.setMaximum(self._protein.amino_acid_composition["C"] / 2)
            self.sbDisulfides.setValue(self._protein.amino_acid_composition["C"] / 2)
            self.calculate_protein_mass()


    def calculate_protein_mass(self):
        """
        Calculates the mass of the protein and the known modifications from the current data.

        Changes:
            self._protein to a protein with given sequence, disulfides and PNGase F modifications
            self._protein_mass to the mass of self._protein
            updates the value of self.lbMassProtein, self.lbMassMods and self.lbMassTotal

        :return: True if there was no error, otherwise None
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
            return None

        self.teSequence.setStyleSheet("QTextEdit { background-color: rgb(240, 251, 240) }")
        if self._mass_set == "AtomsMonoisotopic":
            self._protein_mass = self._protein.monoisotopic_mass
        else:
            self._protein_mass = self._protein.average_mass
        self.lbMassProtein.setText("{:,.2f}".format(self._protein_mass))
        self.lbMassMods.setText("{:,.2f}".format(self._known_mods_mass))
        self.lbMassTotal.setText("{:,.2f}".format(self._protein_mass + self._known_mods_mass))
        return True


    def calculate_mod_mass(self):
        """
        Calculate the mass of known modifications.

        :return: nothing
        """

        # add min counts for glycans to the theoretical mass
        self._known_mods_mass = 0
        for ch, sp_min, sp_max in zip(self._glycan_checkboxes,
                                      self._glycan_min_spinboxes,
                                      self._glycan_max_spinboxes):
            if ch.isChecked():
                sp_min.setEnabled(True)
                sp_max.setEnabled(True)
                g = ch.text()
                if self._mass_set == "AtomsMonoisotopic":
                    glycan_mass = glyco_tools.glycan_formula[g].monoisotopic_mass
                else:
                    glycan_mass = glyco_tools.glycan_formula[g].average_mass
                self._known_mods_mass += glycan_mass * sp_min.value()
            else:
                sp_min.setEnabled(False)
                sp_max.setEnabled(False)
        self.lbMassProtein.setText("{:,.2f}".format(self._protein_mass))
        self.lbMassMods.setText("{:,.2f}".format(self._known_mods_mass))
        self.lbMassTotal.setText("{:,.2f}".format(self._protein_mass + self._known_mods_mass))


    def sample_modifications(self):
        """
        Main function, which eventually calls the combinatorial search and processes its results.

        :return: nothing
        """

        # calculate required input if a single mass was entered (i.e., no peak list was loaded)
        if self.rbSingleMass.isChecked():
            self._exp_mass_data = pd.DataFrame({"Average Mass": self.sbSingleMass.value(),
                                                "Relative Abundance": 100.0}, index=[0])
            self._mass_filename = "Input Mass: {:.4f}".format(self.sbSingleMass.value())

        self.calculate_mod_mass()
        modifications = []  # list of modifications to be used in the combinatorial search
        unknown_masses: pd.DataFrame = self._exp_mass_data["Average Mass"] - self._protein_mass - self._known_mods_mass

        if self.cbTolerance.currentIndex() == 0:  # == "Da."; max_tol_mass is the largest mass plus tolerance
            max_tol_mass = max(self._exp_mass_data["Average Mass"]) + self.sbTolerance.value()
        else:
            max_tol_mass = max(self._exp_mass_data["Average Mass"]) * (1 + self.sbTolerance.value() / 1000000)

        # calculate the explained mass, i.e., the mass of the protein sequence plus known modifications
        # add all enabled single glycans to the list of modifications
        for ch, sp_min, sp_max in zip(self._glycan_checkboxes,
                                      self._glycan_min_spinboxes,
                                      self._glycan_max_spinboxes):
            if ch.isChecked():
                name = ch.text()
                if self._mass_set == "AtomsMonoisotopic":
                    mass = glyco_tools.glycan_formula[name].monoisotopic_mass
                else:
                    mass = glyco_tools.glycan_formula[name].average_mass

                if sp_max.value() == -1:
                    # determine the upper limit of glycans that may appear
                    max_count = min(int((max_tol_mass - self._protein_mass) / mass),
                                    configure.maxmods)
                else:
                    max_count = sp_max.value()
                modifications.append((name, mass, max_count - sp_min.value(), ""))

        print("PERFORMING COMBINATORIAL SEARCH...")
        print("Experimental Masses:", self._exp_mass_data["Average Mass"].head(), sep="\n")
        print("Explained mass (protein + known modifications):", self._protein_mass + self._known_mods_mass)
        print("Unknown masses searched:", unknown_masses.head(), sep="\n")
        print("Mass tolerance: %f %s" % (self.sbTolerance.value(), self.cbTolerance.currentText()))

        # extend the list of modifications by the table entries
        for row_id in range(self.tbMonomers.rowCount()):
            name = self.tbMonomers.item(row_id, 0).text()
            mass = self.tbMonomers.cellWidget(row_id, 1).value()
            min_count = self.tbMonomers.cellWidget(row_id, 2).value()
            max_count = self.tbMonomers.cellWidget(row_id, 3).value()
            site = self.tbMonomers.item(row_id, 4).text()
            if max_count == -1:
                max_count = min(int((max_tol_mass - self._protein_mass) / mass),
                                configure.maxmods)
            modifications.append((name, mass, max_count - min_count, site))

        print("Modifications List:")
        for m in modifications:
            print("{}\t{:.2f}\t{:d}\t{}".format(*m))

        # the actual combinatorial search
        if self.cbTolerance.currentIndex() == 0:  # that is, Da.
            mass_tolerance = self.sbTolerance.value()
        else:
            # calculate a mass tolerance for each peak if we're working with ppm tolerance
            mass_tolerance = []
            for _, m in self._exp_mass_data["Average Mass"].iteritems():
                mass_tolerance.append(m * self.sbTolerance.value() / 1_000_000)

        self._all_hits = modification_search.fast_find_modifications(
            modifications,
            list(unknown_masses),
            mass_tolerance=mass_tolerance,
            explained_mass=self._protein_mass + self._known_mods_mass)
        # n_glycans=self._nglycans_data,  TODO: remove
        # n_positions=len(self._protein.n_sites)

        # the modification search was not successful
        if self._all_hits is None:
            QMessageBox.warning(self, "Warning", "Combinatorial search was unsuccessful.")

        # add the minimum glycan counts to the result data frame
        for ch, sp_min in zip(self._glycan_checkboxes, self._glycan_min_spinboxes):
            if ch.isChecked():
                self._all_hits[ch.text()] += sp_min.value()


    def draw_naked_plot(self):
        """
        Show a mass spectrum after loading a peak list (i.e., without annotations).

        :return: nothing
        """

        axes = self.fig.add_subplot(111)
        axes.vlines(self._exp_mass_data["Average Mass"], 0,
                    self._exp_mass_data["Relative Abundance"], lw=1, picker=10)
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
        self.show_deltas1()
        self.show_deltas2()
        self.canvas.draw()


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

        for i in self._all_hits.index.levels[0]:
            x = self._exp_mass_data.loc[i, "Average Mass"]
            y = self._exp_mass_data.loc[i, "Relative Abundance"]
            first_hit = self._all_hits.loc[i].iloc[0]
            s = "%.2f" % x
            s += "(%.1f %s)\n" % (first_hit["ppm"], "ppm")
            s += "%s" % "\n".join(first_hit["Modstring"].split("|"))
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


    def set_result_tree(self):
        """
        Displays results from the combinatorial search in the result tree.

        :return: nothing
        """
        number_not_annotated = len(self._exp_mass_data) - len(self._all_hits.index.levels[0])
        self.chRemoveUnannotated.setText("Show unannotated (%d out of %d not annotated)"
                                         % (number_not_annotated, len(self._exp_mass_data)))
        self.twResults.clear()

        # set column headers
        header_labels = ["Exp. Mass", "%"]
        header_labels.extend(self._all_hits.columns)
        self.twResults.setColumnCount(len(header_labels))
        self.twResults.setHeaderLabels(header_labels)

        # fill the tree
        if self.chRemoveUnannotated.isChecked():
            iterindex = self._exp_mass_data.index
        else:
            iterindex = self._all_hits.index.levels[0]
        fcolor = QColor(255, 185, 200)
        for mass_index in iterindex:
            # generate root item (experimental mass, relative abundance)
            root_item = SortableTreeWidgetItem(self.twResults)
            root_item.setText(0, "{:.2f}".format(self._exp_mass_data.loc[mass_index]["Average Mass"]))
            root_item.setTextAlignment(0, AlignRight)
            root_item.setText(1, "{:.1f}".format(self._exp_mass_data.loc[mass_index]["Relative Abundance"]))
            root_item.setTextAlignment(1, AlignRight)

            if mass_index not in self._all_hits.index.levels[0]:
                root_item.setBackground(0, QBrush(fcolor))
                root_item.setBackground(1, QBrush(fcolor))  # i.e., for both columns
            else:
                # generate child items, one per possible combination of modifications
                for _, hit in self._all_hits.loc[mass_index].iterrows():
                    child_item = SortableTreeWidgetItem(root_item)
                    mods = hit[:-5].index
                    for j, mod in enumerate(mods):
                        child_item.setText(j + 2, "{:d}".format(hit[mod]))
                        child_item.setTextAlignment(j + 2, AlignHCenter)
                    child_item.setText(len(mods) + 2, "{:.2f}".format(hit["Exp. Mass"]))
                    child_item.setTextAlignment(len(mods) + 2, AlignRight)
                    child_item.setText(len(mods) + 3, "{:.2f}".format(hit["Theo. Mass"]))
                    child_item.setTextAlignment(len(mods) + 3, AlignRight)
                    child_item.setText(len(mods) + 4, "{:.2f}".format(hit["Da."]))
                    child_item.setTextAlignment(len(mods) + 4, AlignRight)
                    child_item.setText(len(mods) + 5, "{:.1f}".format(hit["ppm"]))
                    child_item.setTextAlignment(len(mods) + 5, AlignRight)
                    child_item.setText(len(mods) + 6, hit["Modstring"])
                    # child_item.setTextAlignment(len(mods) + 6, AlignRight)
        self.twResults.expandAll()
        self.twResults.header().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.twResults.header().setStretchLastSection(False)


    def set_result_text(self, output):
        """
        Displays output in the result text browser.

        :param output: The contents of the text edit.
        :return: nothing
        """
        self.showGlycans.clear()
        self.showGlycans.setText(str(output))


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
                                                        "Save Settings",
                                                        self._path,
                                                        "ModFinder settings (*.mofi)")[0]
        self._path = os.path.split(settings_filename)[0]
        if settings_filename:
            if not settings_filename.endswith(".mofi"):
                settings_filename = settings_filename + ".mofi"
            # get modifications from table
            modifications = []
            for row_id in range(self.tbMonomers.rowCount()):
                name = self.tbMonomers.item(row_id, 0).text()
                mass = self.tbMonomers.cellWidget(row_id, 1).value()
                min_count = self.tbMonomers.cellWidget(row_id, 2).value()
                max_count = self.tbMonomers.cellWidget(row_id, 3).value()
                site = self.tbMonomers.item(row_id, 4).text()
                modifications.append((name, mass, min_count, max_count, site))

            settings = {"sequence": self.teSequence.toPlainText(),
                        "exp mass data": self._exp_mass_data,
                        "mass filename": self._mass_filename,
                        "disulfides": self.sbDisulfides.value(),
                        "pngase f": self.chPngase.isChecked(),
                        "glycan states": [ch.isChecked() for ch in self._glycan_checkboxes],
                        "glycan min counts": [sp_min.value() for sp_min in self._glycan_min_spinboxes],
                        "glycan max counts": [sp_max.value() for sp_max in self._glycan_max_spinboxes],
                        "tolerance value": self.sbTolerance.value(),
                        "tolerance flavor": self.cbTolerance.currentIndex(),
                        "modifications": modifications}
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
                self.sbDisulfides.setEnabled(True)
                self.chPngase.setEnabled(True)
                self.sbDisulfides.setMaximum(self._protein.amino_acid_composition["C"] / 2)
            self.sbDisulfides.setValue(settings["disulfides"])
            self.chPngase.setChecked(settings["pngase f"])
            self.calculate_protein_mass()

            self._mass_filename = settings["mass filename"]

            self._all_hits = None
            self.twResults.clear()
            self.fig.clear()
            if settings["exp mass data"] is not None:
                self._exp_mass_data = settings["exp mass data"]
                self.lwPeaks.clear()
                self.lwPeaks.addItems(["{:.2f}".format(i) for i in self._exp_mass_data["Average Mass"]])
                self.lwPeaks.setCurrentRow(0)
                self.sbSingleMass.setValue(float(self.lwPeaks.currentItem().text()))
                self.draw_naked_plot()

            self.cbTolerance.setCurrentIndex(settings["tolerance flavor"])
            self.sbTolerance.setValue(settings["tolerance value"])

            for ch, is_checked in zip(self._glycan_checkboxes, settings["glycan states"]):
                ch.setChecked(is_checked)
            for sp_min, value in zip(self._glycan_min_spinboxes, settings["glycan min counts"]):
                sp_min.setValue(value)
            for sp_max, value in zip(self._glycan_max_spinboxes, settings["glycan max counts"]):
                sp_max.setValue(value)

            self.modtable_clear()
            for row_id, (name, mass, min_count, max_count, site) in enumerate(settings["modifications"]):
                self._modtable_create_row(row_id, name, mass, min_count, max_count, site)

            self.calculate_mod_mass()


    def run_modification_search(self):
        """
        Prepares everysthing for the main algorithm, runs the combinatorial search,
        and displays the results (annotated peak plot, tree and text browser).

        :return: nothing
        """
        self.calculate_protein_mass()
        self.sample_modifications()
        if self._all_hits is not None:
            self.draw_annotated_plot()
            outstring = ["Protein Mass Assessment:",
                         "Disulfide bonds:\t%s" % self.sbDisulfides.value(),
                         "PNGaseF:\t\t%s" % self.chPngase.isChecked(),
                         "Protein sum formula:\t%s" % self._protein.formula,
                         "Average mass:\t%f Da" % self._protein_mass, "%s\nMASS SEARCH:" % (50 * "-"),
                         self._all_hits.round(2).set_index("Exp. Mass").to_string(max_cols=999)]
            self.set_result_text("\n".join(outstring))
            self.set_result_tree()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    frame = MainWindow()
    frame.show()
    app.exec_()
