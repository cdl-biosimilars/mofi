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
import PySide
from PySide.QtGui import (QApplication, QMainWindow, QWidget, QColor, QActionGroup,
                          QMessageBox, QFileDialog, QTreeWidgetItem, QBrush, QHeaderView)

AlignHCenter = PySide.QtCore.Qt.AlignHCenter
AlignRight = PySide.QtCore.Qt.AlignRight

import pandas as pd
pd.set_option('display.max_rows', 5000)

import matplotlib
matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4'] = 'PySide'
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)

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
        self.runSNOB.clicked.connect(self.run_modification_search)
        self.pngasefBox.clicked.connect(self.calculate_protein_mass)
        self.calculateSequence.clicked.connect(self.calculate_protein_mass)
        self.hexCheck.clicked.connect(self.calculate_mod_mass)
        self.hexnacCheck.clicked.connect(self.calculate_mod_mass)
        self.pentCheck.clicked.connect(self.calculate_mod_mass)
        self.neu5acCheck.clicked.connect(self.calculate_mod_mass)
        self.neu5gcCheck.clicked.connect(self.calculate_mod_mass)
        self.fucCheck.clicked.connect(self.calculate_mod_mass)
        self.ocoreCheck.clicked.connect(self.calculate_mod_mass)
        self.ncoreCheck.clicked.connect(self.calculate_mod_mass)
        self.actionOpenFasta.triggered.connect(self.read_fasta_file)
        self.actionOpenMass.triggered.connect(self.read_mass_file)
        self.actionSaveAnnotationcsv.triggered.connect(self.save_csv)
        self.actionSaveSettings.triggered.connect(self.save_settings)
        self.actionLoadSettings.triggered.connect(self.load_settings)
        self.actionLoadNglycanlibrary.triggered.connect(self.read_nglycan_file)
        self.actionAbout.triggered.connect(self.show_about)
        self.actionHelp.triggered.connect(self.show_help)
        self.toleranceBox.activated.connect(self.choose_tolerance_units)
        self.hexCount.valueChanged.connect(self.calculate_mod_mass)
        self.hexnacCount.valueChanged.connect(self.calculate_mod_mass)
        self.pentCount.valueChanged.connect(self.calculate_mod_mass)
        self.neu5acCount.valueChanged.connect(self.calculate_mod_mass)
        self.neu5gcCount.valueChanged.connect(self.calculate_mod_mass)
        self.fucCount.valueChanged.connect(self.calculate_mod_mass)
        self.ocoreCount.valueChanged.connect(self.calculate_mod_mass)
        self.ncoreCount.valueChanged.connect(self.calculate_mod_mass)
        self.hexCount_max.valueChanged.connect(self.calculate_mod_mass)
        self.hexnacCount_max.valueChanged.connect(self.calculate_mod_mass)
        self.pentCount_max.valueChanged.connect(self.calculate_mod_mass)
        self.neu5acCount_max.valueChanged.connect(self.calculate_mod_mass)
        self.neu5gcCount_max.valueChanged.connect(self.calculate_mod_mass)
        self.fucCount_max.valueChanged.connect(self.calculate_mod_mass)
        self.ocoreCount_max.valueChanged.connect(self.calculate_mod_mass)
        self.ncoreCount_max.valueChanged.connect(self.calculate_mod_mass)
        self.deltaToleranceBox.valueChanged.connect(self.check_tolerance)
        self.disfBox.valueChanged.connect(self.calculate_protein_mass)
        self.delta1Check.clicked.connect(self.show_deltas1)
        self.delta1Box.valueChanged.connect(self.show_deltas1)
        self.delta2Check.clicked.connect(self.show_deltas2)
        self.delta2Box.valueChanged.connect(self.show_deltas2)
        self.treeFilter.clicked.connect(self.set_result_tree)
        self.nglycanCheck.clicked.connect(self.calculate_mod_mass)
        self.nglycanCheck.hide()
        self.actionQuit.triggered.connect(QApplication.instance().quit)

        # group the mass set selectors
        self.agAverage = QActionGroup(self.menuSelect_Ataomic_Masses)
        self.agAverage.addAction(self.actionAverage_IUPAC)
        self.agAverage.addAction(self.actionAverage_Zhang)
        self.agAverage.addAction(self.actionMonoisotopic)
        self.agAverage.triggered.connect(self.choose_mass_set)

        # set up the plot
        self.graph_frame = QWidget()
        self.fig = Figure((5.0, 3.0), dpi=100, frameon=False, tight_layout=True, edgecolor='white')
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.graph_frame)
        self.canvas.setFocus()
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.graph_frame)
        self.choose_tolerance_units()

        # initialize private members
        # these variables completely describe the state of the program
        self._all_hits = None  # results from the modification search
        self._delta_lines_1 = None  # delta 1 lines in the mass spectrum
        self._delta_lines_2 = None  # delta 2 lines in the mass spectrum
        self._exp_mass_data = None  # pandas dataframe containing the contents of the mass file
        self._glycan_max_counts = None  # maximum number of glycans to search
        self._glycan_min_counts = None  # minimum number of glycans to search
        self._glycan_states = None  # indicates which of the glycan checkboxes are active
        self._known_mods_mass = 0  # mass of known modification
        self._mass_filename = None  # name of the mass file
        self._mass_set = configure.default_masses  # currently used atomic masses
        self._nglycans = None  # list of N-glycan modifications
        self._nglycans_data = None  # pandas dataframe containing the N-glycan library
        self._path = configure.path
        self._protein = None  # a Protein object representing the input sequence with disulfides and PNGase F digest
        self._protein_mass = 0  # mass of the current Protein object


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
        if self.toleranceBox.currentIndex() == 0:  # that is, Da.
            self.massRange.setDecimals(2)
            self.massRange.setMinimum(0.0)
            self.massRange.setMaximum(50.0)
            self.massRange.setSingleStep(0.1)
            self.massRange.setValue(configure.default_da)
        else:
            self.massRange.setDecimals(0)
            self.massRange.setMinimum(0)
            self.massRange.setMaximum(150)
            self.massRange.setSingleStep(1)
            self.massRange.setValue(configure.default_ppm)


    def read_fasta_file(self):
        """
        Opens a FASTA file and displays its contents in the QTextEdit

        Changes:
            self._path to directory of selected file
            text of self.sequenceInput to FASTA string
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
                self.sequenceInput.setText(fasta_input)
                self.calculate_protein_mass()
                self.disfBox.setValue(self._protein.amino_acid_composition['C'] / 2)
                self.disfBox.setMaximum(self._protein.amino_acid_composition['C'] / 2)
                self.calculate_protein_mass()
            else:
                QMessageBox.warning(self, "Error loading sequence", "Not a valid file format: %s" % ext)


    def read_mass_file(self):
        """
        Opens a mass list as generated by Thermo BioPharma Finder and display its contents.

        :return: nothing
        """

        self.massList.clear()
        self.resultsTree.clear()
        self.fig.clear()
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
            self._exp_mass_data = mass_data
            self.massList.clear()
            self.massList.insertItems(0, [str(i) for i in self._exp_mass_data["Average Mass"].round(4)])
            self.massList.show()
            self.massBox.setValue(float(self.massList.currentText()))
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
            parameters = [("Disulfides", self.disfBox.value()),
                          ("PNGaseF", self.pngasefBox.isChecked()),
                          ("Lysines", self.lysineCheck.isChecked()),
                          ("Tolerance", self.massRange.value()),
                          ("Tol. Type", self.toleranceBox.currentText())]
            for g in self._glycan_states:
                if self._glycan_states[g]:
                    parameters.append((g + "_min", self._glycan_min_counts[g]))
                    parameters.append((g + "_max", self._glycan_max_counts[g]))
            if self.customCheck.isChecked():
                parameters.append(("custom Name", self.customName.text()))
                parameters.append(("custom Mass", self.customMass.value()))
                parameters.append(("custom Count", self.customCount.value()))

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
        if self.actionAverage_IUPAC.isChecked():
            self._mass_set = "AtomsAverageIUPAC"
            configure.set_average_masses(self._mass_set)
        elif self.actionAverage_Zhang.isChecked():
            self._mass_set = "AtomsAverageZhang"
            configure.set_average_masses(self._mass_set)
        else:
            self._mass_set = "AtomsMonoisotopic"

        if self.sequenceInput.toPlainText():
            self.calculate_protein_mass()
            self.calculate_mod_mass()


    def calculate_protein_mass(self):
        """
        Calculates the mass of the protein and the known modifications from the current data.

        Changes:
            self._protein to a protein with given sequence, disulfides and PNGase F modifications
            self._protein_mass to the mass of self._protein
            updates the value of self.sequenceMass
        """

        protein_sequence = self.sequenceInput.toPlainText()
        chains, sequence = sequence_tools.read_fasta_string(protein_sequence)
        self._protein = sequence_tools.Protein(sequence,
                                               chains,
                                               self.disfBox.value(),
                                               self.pngasefBox.isChecked())

        if self._mass_set == "AtomsMonoisotopic":
            self._protein_mass = self._protein.monoisotopic_mass
        else:
            self._protein_mass = self._protein.average_mass
        self.sequenceMass.display(self._protein_mass + self._known_mods_mass)


    def calculate_mod_mass(self):
        """
        Calculate the mass of known modifications.

        :return: nothing
        """

        self._glycan_states = {"Hex": self.hexCheck.isChecked(),
                               "HexNAc": self.hexnacCheck.isChecked(),
                               "Fuc": self.fucCheck.isChecked(),
                               "Neu5Ac": self.neu5acCheck.isChecked(),
                               "Neu5Gc": self.neu5gcCheck.isChecked(),
                               "Pent": self.pentCheck.isChecked(),
                               "O-core": self.ocoreCheck.isChecked(),
                               "N-core": self.ncoreCheck.isChecked()}

        self._glycan_min_counts = {"Hex": self.hexCount.value(),
                                   "HexNAc": self.hexnacCount.value(),
                                   "Fuc": self.fucCount.value(),
                                   "Neu5Ac": self.neu5acCount.value(),
                                   "Neu5Gc": self.neu5gcCount.value(),
                                   "Pent": self.pentCount.value(),
                                   "O-core": self.ocoreCount.value(),
                                   "N-core": self.ncoreCount.value()}

        self._glycan_max_counts = {"Hex": self.hexCount_max.value(),
                                   "HexNAc": self.hexnacCount_max.value(),
                                   "Fuc": self.fucCount_max.value(),
                                   "Neu5Ac": self.neu5acCount_max.value(),
                                   "Neu5Gc": self.neu5gcCount_max.value(),
                                   "Pent": self.pentCount_max.value(),
                                   "O-core": self.ocoreCount_max.value(),
                                   "N-core": self.ncoreCount_max.value()}

        # add min counts for glycans to the theoretical mass
        self._known_mods_mass = 0
        for g, is_active in self._glycan_states.items():
            if is_active:
                if self._mass_set == "AtomsMonoisotopic":
                    glycan_mass = glyco_tools.glycan_formula[g].monoisotopic_mass
                else:
                    glycan_mass = glyco_tools.glycan_formula[g].average_mass
                self._known_mods_mass += glycan_mass * self._glycan_min_counts[g]
        self.sequenceMass.display(self._protein_mass + self._known_mods_mass)


    def sample_modifications(self):
        """
        Main function, which eventually calls the combinatorial search and processes its results.

        :return: nothing
        """

        # calculate required input if a single mass was entered (i.e., no peak list was loaded)
        if not self._mass_filename and self.massBox.value() > 0:
            self._exp_mass_data = pd.DataFrame(
                {"Average Mass": pd.Series([self.massBox.value()]), "Relative Abundance": 100.0})
            self._mass_filename = "Input Mass: {}".format(self._mass)

        self.calculate_mod_mass()
        modifications = []  # list of modifications to be used in the combinatorial search

        # calculate the explained mass, i.e., the mass of the protein sequence plus known modifications
        # add all enabled single glycans to the list of modifications
        for g, is_active in self._glycan_states.items():
            if is_active:
                if self._mass_set == "AtomsMonoisotopic":
                    glycan_mass = glyco_tools.glycan_formula[g].monoisotopic_mass
                else:
                    glycan_mass = glyco_tools.glycan_formula[g].average_mass

                if self._glycan_max_counts[g] == -1:
                    # determine the upper limit of glycans that may appear
                    if self.toleranceBox.currentIndex() == 0:  # that is, Da.
                        max_tol_mass = max(self._exp_mass_data["Average Mass"]) + self.massRange.value()
                    else:
                        max_tol_mass = max(self._exp_mass_data["Average Mass"]) * (1 + self.massRange.value() / 1000000)

                    glycan_maxcount = int((max_tol_mass - self._protein_mass) / glycan_mass)
                    glycan_maxcount = min(glycan_maxcount, configure.maxmods)
                else:
                    glycan_maxcount = self._glycan_max_counts[g]
                modifications.append((g, glycan_mass, glycan_maxcount - self._glycan_min_counts[g]))

        # add up to two C-terminal lysines to the list of modifications
        if self.lysineCheck.isChecked():
            lysine = mass_tools.Formula(sequence_tools.amino_acid_compositions["K"])
            if self._mass_set == "AtomsMonoisotopic":
                modifications.append(("Lys", lysine.monoisotopic_mass, 2))
            else:
                modifications.append(("Lys", lysine.average_mass, 2))

        # add the single custom modification to the list of modifications
        custom_mod_name = self.customName.text()
        custom_mod_mass = self.customMass.value()
        custom_mod_count = self.customCount.value()
        if self.customCheck.isChecked() and custom_mod_mass != 0 and custom_mod_count > 0:
            modifications.append((custom_mod_name, custom_mod_mass, custom_mod_count))

        # print("PERFORMING COMBINATORIAL SEARCH...")
        # print("Experimental Masses:", self._exp_mass_data["Average Mass"].head(), sep="\n")
        # print("Explained mass (protein + known modifications):", self._protein_mass + self._known_mods_mass)
        # print("Unknown masses searched:", unknown_masses.head(), sep="\n")
        # print("Mass tolerance: %f %s" % (self.massRange.value(), self.toleranceBox.currentText()))

        # extend the list of modifications by mABs or a loaded N-glycan library
        if self.mabCheck.isChecked():
            mabs = glyco_tools.glycanlist_to_modlist(glyco_tools.fc_glycans,
                                                     use_monoisotopic_masses=self.actionMonoisotopic.isChecked())
            modifications.extend(mabs)
            self._nglycans_data = pd.DataFrame(mabs)
            self._nglycans_data.columns = ["Name", "Mass", "Site"]
            self._nglycans_data["Nr."] = self._nglycans_data["Site"]
            del self._nglycans_data["Mass"]
        elif self.nglycanCheck.isChecked():
            modifications.extend(self._nglycans)
        else:
            self._nglycans_data = pd.DataFrame()

        print("Modifications List:")
        for m in modifications:
            print('%s : %.2f : %d' % m)

        # the actual combinatorial search
        unknown_masses = self._exp_mass_data["Average Mass"] - self._protein_mass - self._known_mods_mass
        if self.toleranceBox.currentIndex() == 0:  # that is, Da.
            mass_tolerance = self.massRange.value()
        else:
            # calculate a mass tolernace for each peak if we're working with ppm tolerance
            mass_tolerance = []
            for _, m in self._exp_mass_data["Average Mass"].iteritems():
                mass_tolerance.append(m * self.massRange.value() / 1000000)

        self._all_hits = modification_search.fast_find_modifications(
            modifications,
            list(unknown_masses),
            mass_tolerance=mass_tolerance,
            explained_mass=self._protein_mass + self._known_mods_mass,
            n_glycans=self._nglycans_data,
            n_positions=len(self._protein.n_sites))

        # the modification search was not successful
        if self._all_hits is None:
            QMessageBox.warning(self, "Warning", "Combinatorial search was unsuccessful.")

        # add the minimum glycan counts to the result data frame
        for g in self._glycan_min_counts:
            if self._glycan_states[g]:
                self._all_hits[g] += self._glycan_min_counts[g]


    def draw_naked_plot(self):
        """
        Show a mass spectrum after loading a peak list (i.e., without annotations).

        :return: nothing
        """

        self.spectrumView.addWidget(self.canvas)
        self.spectrumView.addWidget(self.mpl_toolbar)
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
        self.treeFilter.setText("Show unannotated (%d out of %d not annotated)"
                                % (number_not_annotated, len(self._exp_mass_data)))
        self.resultsTree.clear()

        # set column headers
        header_labels = ["Exp. Mass", "%"]
        header_labels.extend(self._all_hits.columns)
        self.resultsTree.setColumnCount(len(header_labels))
        self.resultsTree.setHeaderLabels(header_labels)

        # fill the tree
        if self.treeFilter.isChecked():
            iterindex = self._exp_mass_data.index
        else:
            iterindex = self._all_hits.index.levels[0]
        fcolor = QColor(255, 185, 200)
        for mass_index in iterindex:
            # generate root item (experimental mass, relative abundance)
            root_item = SortableTreeWidgetItem(self.resultsTree)
            root_item.setText(0, "%.1f" % self._exp_mass_data.loc[mass_index]["Average Mass"])
            root_item.setTextAlignment(0, AlignRight)
            root_item.setText(1, "%.1f" % self._exp_mass_data.loc[mass_index]["Relative Abundance"])
            root_item.setTextAlignment(1, AlignRight)

            if mass_index not in self._all_hits.index.levels[0]:
                root_item.setBackground(0, QBrush(fcolor))
                root_item.setBackground(1, QBrush(fcolor))  # i.e., for both columns
            else:
                # generate child items, one per possible combination of modifications
                for _, hit in self._all_hits.loc[mass_index].iterrows():
                    child_item = SortableTreeWidgetItem(root_item)
                    child_item.setText(2, hit["Modstring"])
                    mods = hit[1:-4].index
                    for j, mod in enumerate(mods):
                        child_item.setText(j + 3, "%d" % hit[mod])
                        child_item.setTextAlignment(j + 3, AlignHCenter)
                    child_item.setText(len(mods) + 3, "%.1f" % hit["Exp. Mass"])
                    child_item.setTextAlignment(len(mods) + 3, AlignRight)
                    child_item.setText(len(mods) + 4, "%.1f" % hit["Theo. Mass"])
                    child_item.setTextAlignment(len(mods) + 4, AlignRight)
                    child_item.setText(len(mods) + 5, "%.1f" % hit["Da."])
                    child_item.setTextAlignment(len(mods) + 5, AlignRight)
                    child_item.setText(len(mods) + 6, "%.1f" % hit["ppm"])
                    child_item.setTextAlignment(len(mods) + 6, AlignRight)
        self.resultsTree.expandAll()
        self.resultsTree.header().setResizeMode(QHeaderView.ResizeToContents)
        self.resultsTree.header().setStretchLastSection(False)


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
    #     selection = self.resultsTree.selectedItems()
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
    #         s += '(%.1f %s)\n' % (seldata[str(self.toleranceBox.currentText())], self.toleranceBox.currentText())
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
                                                        self.delta1Box.value(),
                                                        self.deltaToleranceBox.value() / 2)
            delta_plot = self.fig.add_subplot(111)
            y = [mass[5] for mass in delta_masses]
            x_start = [mass[3] for mass in delta_masses]
            x_stop = [mass[4] for mass in delta_masses]
            if self._delta_lines_1 is not None:
                try:
                    self._delta_lines_1.remove()
                except ValueError:
                    pass
            if self.delta1Check.isChecked():
                self._delta_lines_1 = delta_plot.hlines(y, x_start, x_stop, lw=1.5, color="green")
            self.canvas.draw()


    def show_deltas2(self):
        """
            Show lines between peaks whose masses differ by the value of the Delta1 spinbox (+- tolerance).

            :return: nothing
            """
        if self._exp_mass_data is not None:
            delta_masses = mass_tools.find_delta_masses(self._exp_mass_data,
                                                        self.delta2Box.value(),
                                                        self.deltaToleranceBox.value() / 2)
            delta_plot = self.fig.add_subplot(111)
            y = [mass[5] * .6 for mass in delta_masses]
            x_start = [mass[3] for mass in delta_masses]
            x_stop = [mass[4] for mass in delta_masses]
            if self._delta_lines_2 is not None:
                try:
                    self._delta_lines_2.remove()
                except ValueError:
                    pass
            if self.delta2Check.isChecked():
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
            settings = {"sequence": self.sequenceInput.toPlainText(),
                        "exp mass data": self._exp_mass_data,
                        "mass filename": self._mass_filename,
                        "nglycans data": self._nglycans_data,
                        "sequence settings": {"dbonds": self.disfBox.value(),
                                              "pngasef": self.pngasefBox.isChecked()},
                        "glycan_states": self._glycan_states,
                        "glycan min counts": self._glycan_min_counts,
                        "glycan max counts": self._glycan_max_counts,
                        "custom settings": {"custom on": self.customCheck.isChecked(),
                                            "custom name": self.customName.text(),
                                            "custom mass": self.customMass.value(),
                                            "custom count": self.customCount.value()},
                        "lysine sampling": self.lysineCheck.isChecked(),
                        "mabs checked": self.mabCheck.isChecked(),
                        "nglycans checked": self.nglycanCheck.isChecked(),
                        "tolerance settings": {"tolerance value": self.massRange.value(),
                                               "tolerance flavor": self.toleranceBox.currentIndex()}}
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

            self.sequenceInput.setText(settings["sequence"])
            self._mass_filename = settings["mass filename"]

            self.resultsTree.clear()
            self.fig.clear()
            if settings["exp mass data"] is not None:
                self._exp_mass_data = settings["exp mass data"]
                self.massList.clear()
                self.massList.insertItems(0, [str(i) for i in self._exp_mass_data["Average Mass"]])
                self.massList.show()
                self.massBox.setValue(float(self.massList.currentText()))
                self.draw_naked_plot()

            self._all_hits = None
            self._nglycans_data = settings["nglycans data"]
            if self._nglycans_data is not None:
                self._nglycans = glyco_tools.dataframe_to_modlist(self._nglycans_data)
            if settings["nglycans checked"]:
                self.nglycanCheck.show()
                self.nglycanCheck.setChecked(True)

            self.disfBox.setValue(settings["sequence settings"]["dbonds"])
            self.pngasefBox.setChecked(settings["sequence settings"]["pngasef"])
            self.lysineCheck.setChecked(settings["lysine sampling"])
            self.mabCheck.setChecked(settings["mabs checked"])

            self.toleranceBox.setCurrentIndex(settings["tolerance settings"]["tolerance flavor"])
            self.massRange.setValue(settings["tolerance settings"]["tolerance value"])

            self.customCheck.setChecked(settings["custom settings"]["custom on"])
            self.customName.setText(settings["custom settings"]["custom name"])
            self.customMass.setValue(settings["custom settings"]["custom mass"])
            self.customCount.setValue(settings["custom settings"]["custom count"])

            if settings["glycan_states"] is not None:
                self.hexCheck.setChecked(settings["glycan states"]["Hex"])
                self.hexnacCheck.setChecked(settings["glycan states"]["HexNAc"])
                self.fucCheck.setChecked(settings["glycan states"]["Fuc"])
                self.neu5acCheck.setChecked(settings["glycan states"]["Neu5Ac"])
                self.neu5gcCheck.setChecked(settings["glycan states"]["Neu5Gc"])
                self.pentCheck.setChecked(settings["glycan states"]["Pent"])
                self.ocoreCheck.setChecked(settings["glycan states"]["O-core"])
                self.ncoreCheck.setChecked(settings["glycan states"]["N-core"])
                self.hexCount.setValue(settings["glycan min counts"]["Hex"])
                self.hexnacCount.setValue(settings["glycan min counts"]["HexNAc"])
                self.fucCount.setValue(settings["glycan min counts"]["Fuc"])
                self.neu5acCount.setValue(settings["glycan min counts"]["Neu5Ac"])
                self.neu5gcCount.setValue(settings["glycan min counts"]["Neu5Gc"])
                self.pentCount.setValue(settings["glycan min counts"]["Pent"])
                self.ocoreCount.setValue(settings["glycan min counts"]["O-core"])
                self.ncoreCount.setValue(settings["glycan min counts"]["N-core"])
                self.hexCount_max.setValue(settings["glycan max counts"]["Hex"])
                self.hexnacCount_max.setValue(settings["glycan max counts"]["HexNAc"])
                self.fucCount_max.setValue(settings["glycan max counts"]["Fuc"])
                self.neu5acCount_max.setValue(settings["glycan max counts"]["Neu5Ac"])
                self.neu5gcCount_max.setValue(settings["glycan max counts"]["Neu5Gc"])
                self.pentCount_max.setValue(settings["glycan max counts"]["Pent"])
                self.ocoreCount_max.setValue(settings["glycan max counts"]["O-core"])
                self.ncoreCount_max.setValue(settings["glycan max counts"]["N-core"])

            self.calculate_protein_mass()
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
                         "Disulfide bonds:\t%s" % self.disfBox.value(),
                         "PNGaseF:\t\t%s" % self.pngasefBox.isChecked(),
                         "Protein sum formula:\t%s" % self._protein.formula,
                         "Average mass:\t%f Da" % self._protein_mass, "%s\nMASS SEARCH:" % (50 * "-"),
                         self._all_hits.round(2).set_index("Exp. Mass").to_string(max_cols=999)]
            self.set_result_text("\n".join(outstring))
            self.set_result_tree()


if __name__ == "__main__":
    import ctypes
    myappid = u"MoFi-1.0"
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
    
    app = QApplication(sys.argv)
    frame = MainWindow()
    frame.show()
    app.exec_()
