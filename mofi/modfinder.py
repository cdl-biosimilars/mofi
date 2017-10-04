"""
The main module, which includes GUI definitions and major functions.
"""

import math
import os
import re
import sys
import time
import xml.etree.ElementTree as ETree
import webbrowser

from PyQt5.QtWidgets import (QApplication, QMainWindow, QMenu,
                             QTableWidgetItem, QCheckBox, QMessageBox,
                             QHeaderView, QButtonGroup,
                             QSpinBox, QDoubleSpinBox, QWidget, QHBoxLayout,
                             QProgressBar, QLabel, QSizePolicy)
from PyQt5.QtGui import QColor
from PyQt5.QtCore import Qt, QLocale

import numpy as np
import pandas as pd

import matplotlib
import matplotlib.text
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.widgets import SpanSelector
from matplotlib.figure import Figure

from mofi import (configure, mass_tools, io_tools,
                  search_tools, sequence_tools)
from mofi.paths import data_dir, docs_dir
from mofi.modfinder_ui import Ui_ModFinder
from mofi.widgets import (FilterHeader, CollapsingRectangleSelector,
                          SortableTreeWidgetItem, SortableTableWidgetItem,
                          get_filename)

_version_info = """ModFinder v1.0

© 2017 Christian Doppler Laboratory
for Innovative Tools for Biosimilar Characterization

Contact: Wolfgang.Skala@sbg.ac.at

Python version:
{}""".format(sys.version)

# default values for columns in the monomer table
_monomer_table_columns = [
    ("Checked", False),
    ("Name", None),
    ("Composition", None),
    ("Min", 0),
    ("Max", -1)
]

# default values for columns in the polymer table
_polymer_table_columns = [
    ("Checked", True),
    ("Name", None),
    ("Composition", None),
    ("Sites", ""),
    ("Abundance", 0.0)
]

# sorting key for column 1 of the results trees
_default_col_key = {1: lambda x: [int(i) for i in x.split("-")]}


def find_in_intervals(value, intervals):
    """
    Simple :math:`O(n)` algorithm to determine whether a value
    falls into a set of intervals.

    Examples:
    ``value=12, intervals={"a": (1, 6), "b": (9, 14)}`` -> ``"b"``
    ``value=8,  intervals={"a": (1, 6), "b": (9, 14)}`` -> ``""``

    :param float value: Value to search
    :param dict intervals: {interval name: (lower interval boundary,
                                            upper interval boundary)}
    :return: Name of the interval containing the value;
             empty string if no such interval exists
    :rtype: str
    """
    for name, (lower, upper) in intervals.items():
        if lower <= value <= upper:
            return name
    return ""


def table_insert_row(table_widget, above=True):
    """
    Insert a row into a :class:`QTableWidget`.
    The row will be inserted relative to the current selection
    (if one exists) or to all rows otherwise.

    :param QTableWidget table_widget: the table widget to modify
    :param bool above: True if the row should be inserted
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


def create_child_items(df, root_item, column_count, monomers, sites):
    """
    Fill the results table with child items (hit and permutation)
    for each peak.

    :param pd.DataFrame df: data for the child items
    :param SortableTreeWidgetItem root_item: root item for the child
    :param int column_count: total number of columns in results table
    :param list monomers: list of monomer column names
    :param list sites: list of glycan site column names
    :return: nothing
    """

    if sites:  # show stage 2 results
        # add the best annotation to the parent row
        # peak_optimum is a (hit_ID, permutation_ID) tuple
        try:
            peak_optimum = df.loc[df["Permutation score"].idxmax()]
            peak_optimum.name = (peak_optimum.name[1],
                                 peak_optimum.name[2])
            create_hit_columns(root_item, peak_optimum, monomers, sites)
        except ValueError:  # if df is empty, like due to filtering
            pass

        # (1) child items for all hits per peak
        for stage2_id, hit in (df.reset_index("Isobar", drop=True)
                                 .groupby("Stage2_hit")):
            hit_item = SortableTreeWidgetItem(
                root_item,
                default_key=float,
                col_key=_default_col_key)
            hit_optimum = hit.loc[hit["Permutation score"].idxmax()]
            hit_item.setText(
                1, "{}-{}".format(root_item.text(1), stage2_id))
            hit_item.setTextAlignment(1, Qt.AlignLeft)
            hit_item.setFlags(hit_item.flags() | Qt.ItemIsTristate)
            hit_item.setCheckState(0, Qt.Unchecked)
            pos = create_hit_columns(hit_item, hit_optimum,
                                     monomers, sites)

            # background color
            if stage2_id % 2 == 0:
                hit_item.setTotalBackground(
                    QColor(configure.colors["table_hit_even"]),
                    column_count)
            else:
                hit_item.setTotalBackground(
                    QColor(configure.colors["table_hit_odd"]),
                    column_count)

            # (2) child items for all permutations per hit
            # only if there are at least two permutations
            if hit.shape[0] > 1:
                for perm_id, perm in (hit.reset_index("Stage2_hit")
                                         .iterrows()):
                    perm.name = (0, perm.name)
                    perm_item = SortableTreeWidgetItem(
                        hit_item,
                        default_key=float,
                        col_key=_default_col_key)
                    perm_item.setText(
                        1, "{}-{}-{}".format(root_item.text(1),
                                             stage2_id,
                                             perm_id))
                    perm_item.setTextAlignment(1, Qt.AlignLeft)
                    perm_item.setCheckState(0, Qt.Unchecked)
                    create_site_columns(perm_item, pos, perm, sites)

    else:  # stage 1 results
        try:
            create_hit_columns(root_item, df.iloc[0], monomers, [])
        except ValueError:  # if df is empty, like due to filtering
            pass
        # child items for all hits per peak
        stage1_id = 0
        for hit in df.reset_index().itertuples():
            hit_item = SortableTreeWidgetItem(
                root_item,
                default_key=float,
                col_key=_default_col_key)
            hit_item.setText(
                1, "{}-{}".format(root_item.text(1), stage1_id))
            hit_item.setTextAlignment(1, Qt.AlignLeft)
            hit_item.setCheckState(0, Qt.Unchecked)
            create_hit_columns(hit_item, hit, monomers, sites)
            stage1_id += 1

            # background color
            if hit.Index % 2 == 0:
                hit_item.setTotalBackground(
                    QColor(configure.colors["table_hit_even"]),
                    column_count)
            else:
                hit_item.setTotalBackground(
                    QColor(configure.colors["table_hit_odd"]),
                    column_count)


def create_hit_columns(item, hit, monomers, sites):
    """
    Create columns that contain information on a stage 1/2 hit
    in the results table.

    :param SortableTreeWidgetItem item: row to fill
    :param pd.Series/namedtuple hit: hit data
    :param list monomers: list of monosaccharides
    :param list sites: list of glycosylation sites
    :return: the position of the first unused column
    :rtype: int
    """

    pos = 4
    if sites:  # stage 2 results
        # hit index
        item.setText(pos, "{}".format(hit.name[0]))
        item.setTextAlignment(pos, Qt.AlignRight)
        pos += 1

        # hit properties
        for label, form in [("Hit score", "{:.2f}"),
                            ("Permutations", "{}"),
                            ("Theo_Mass", "{:.2f}"),
                            ("Da", "{:.2f}"),
                            ("ppm", "{:.2f}")]:
            item.setText(pos, form.format(hit[label]))
            item.setTextAlignment(pos, Qt.AlignRight)
            pos += 1

        # monomer counts
        for monomer in monomers:
            item.setText(pos, "{}".format(hit[monomer]))
            item.setTextAlignment(pos, Qt.AlignRight)
            pos += 1

        # highest-scoring glycan combination
        create_site_columns(item, pos, hit, sites)

    else:  # stage 1 results
        # hit properties
        for label, form in [("Theo_Mass", "{:.2f}"),
                            ("Da", "{:.2f}"),
                            ("ppm", "{:.2f}")]:
            item.setText(pos, form.format(getattr(hit, label)))
            item.setTextAlignment(pos, Qt.AlignRight)
            pos += 1

        # monomer counts
        for monomer in monomers:
            item.setText(pos, "{:.0f}".format(getattr(hit, monomer)))
            item.setTextAlignment(pos, Qt.AlignRight)
            pos += 1

    return pos


def create_site_columns(item, pos, hit, sites):
    """
    Create columns that contain information on glycan sites
    in the results table.

    :param SortableTreeWidgetItem item: row to fill
    :param int pos: position of the first column
    :param pd.Series hit: hit data
    :param list sites: list of glycosylation sites
    :return: nothing
    """

    # permutation index
    item.setText(pos, "{}".format(hit.name[1]))
    item.setTextAlignment(pos, Qt.AlignRight)
    pos += 1

    # permutation score
    if not np.isnan(hit["Permutation score"]):
        item.setText(pos, "{:.2f}".format(hit["Permutation score"]))
        item.setTextAlignment(pos, Qt.AlignRight)
    pos += 1

    # glycan sites
    for site in sites:
        item.setText(pos, "{}".format(hit[site]))
        item.setTextAlignment(pos, Qt.AlignRight)
        pos += 1


def table_clear(table_widget):
    """
    Delete all rows in a :class:`QTableWidget`.

    :param QTableWidget table_widget: the table widget to modify
    :return: nothing
    """

    table_widget.clearContents()
    table_widget.setRowCount(0)


def table_delete_row(table_widget):
    """
    Delete selected rows in a :class:`QTableWidget`.

    :param QTableWidget table_widget: the table widget to modify
    :return: nothing
    """
    if table_widget.selectionModel().selectedRows():
        for i in table_widget.selectionModel().selectedRows()[::-1]:
            table_widget.removeRow(i.row())


def item_is_shown(query, item, col):
    """
    Determine whether a table widget item is shown
    after applying a filter.

    :param list query: list of query tuples
    :param SortableTreeWidgetItem item: the item to check
    :param int col: the first column containing a monosaccharide
    :return: True if all columns satisfy the filtering condition
             False otherwise
    """

    for lower, upper in query:
        if not (lower <= int(item.text(col)) <= upper):
            return False
        col += 1
    return True


class MainWindow(QMainWindow, Ui_ModFinder):
    def __init__(self, parent=None):

        # initialize the GUI
        super(MainWindow, self).__init__(parent)
        self.setupUi(self)

        # connect signals to slots
        self.acAbout.triggered.connect(self.show_about)
        self.acLoadSettings.triggered.connect(self.load_settings)
        self.acManual.triggered.connect(
            lambda: webbrowser.open("file://"
                                    + os.path.join(docs_dir, "index.html")))
        self.acQuit.triggered.connect(QApplication.instance().quit)
        self.acSaveSettings.triggered.connect(self.save_settings)

        self.btCheckAll.clicked.connect(
            lambda: self.check_results_tree(check=True))
        self.btClearFilters.clicked.connect(self.clear_filters)
        self.btClearSequence.clicked.connect(self.teSequence.clear)
        self.btClearMonomers.clicked.connect(
            lambda: table_clear(self.tbMonomers))
        self.btClearPolymers.clicked.connect(
            lambda: table_clear(self.tbPolymers))
        self.btCollapseAll.clicked.connect(
            lambda: self.expand_results_tree(expand=False))
        self.btDeleteRowMonomers.clicked.connect(
            lambda: table_delete_row(self.tbMonomers))
        self.btDeleteRowPolymers.clicked.connect(
            lambda: table_delete_row(self.tbPolymers))
        self.btExpandParents.clicked.connect(
            lambda: self.expand_results_tree(expand=True, depth=0))
        self.btExpandAll.clicked.connect(
            lambda: self.expand_results_tree(expand=True, depth=1))
        self.btFindModifications.clicked.connect(self.sample_modifications)
        self.btInsertRowAboveMonomers.clicked.connect(
            lambda: table_insert_row(self.tbMonomers, above=True))
        self.btInsertRowAbovePolymers.clicked.connect(
            lambda: table_insert_row(self.tbPolymers, above=True))
        self.btInsertRowBelowMonomers.clicked.connect(
            lambda: table_insert_row(self.tbMonomers, above=False))
        self.btInsertRowBelowPolymers.clicked.connect(
            lambda: table_insert_row(self.tbPolymers, above=False))
        self.btLabelPeaks.clicked.connect(lambda: self.update_selection())
        self.btLoadMonomers.clicked.connect(
            lambda: self.load_table(
                default=False,
                dialog_title="Import modifications",
                extensions=["csv", "xls", "xlsx"],
                table_widget=self.tbMonomers,
                cols=_monomer_table_columns
            )
        )
        self.btLoadPolymers.clicked.connect(
            lambda: self.load_table(
                default=False,
                dialog_title="Import glycans",
                extensions=["csv", "xls", "xlsx"],
                table_widget=self.tbPolymers,
                cols=_polymer_table_columns
            )
        )
        self.btLoadPeaks.clicked.connect(self.load_spectrum)
        self.btLoadSequence.clicked.connect(self.load_sequence)
        self.btResetZoom.clicked.connect(self.reset_zoom)
        self.btSaveMonomers.clicked.connect(
            lambda: self.save_table("Export modifications", "monomers"))
        self.btSavePolymers.clicked.connect(
            lambda: self.save_table("Export glycans", "polymers"))
        self.btSaveSequence.clicked.connect(self.save_sequence)
        self.btSaveSpectrum.clicked.connect(self.save_spectrum)
        self.btUncheckAll.clicked.connect(
            lambda: self.check_results_tree(check=False))
        self.btUpdateMass.clicked.connect(self.calculate_protein_mass)

        self.cbTolerance.activated.connect(self.choose_tolerance_units)

        self.chCombineDelta.clicked.connect(lambda: self.update_selection())
        self.chDelta1.clicked.connect(self.toggle_delta_series)
        self.chDelta2.clicked.connect(self.toggle_delta_series)
        self.chPngase.clicked.connect(self.calculate_protein_mass)

        self.sbDeltaRepetition1.valueChanged.connect(
            lambda: self.update_selection())
        self.sbDeltaRepetition2.valueChanged.connect(
            lambda: self.update_selection())
        self.sbDeltaTolerance1.valueChanged.connect(
            lambda: self.update_selection())
        self.sbDeltaTolerance2.valueChanged.connect(
            lambda: self.update_selection())
        self.sbDeltaValue1.valueChanged.connect(
            lambda: self.update_selection())
        self.sbDeltaValue2.valueChanged.connect(
            lambda: self.update_selection())
        self.sbDisulfides.valueChanged.connect(self.calculate_protein_mass)

        self.taResults.currentChanged.connect(self.set_save_results_menu)

        self.tbMonomers.cellChanged.connect(self.calculate_mod_mass)

        self.teSequence.textChanged.connect(
            lambda: self.teSequence.setStyleSheet(
                "QTextEdit { background-color: rgb(255, 225, 225) }"))

        self.twResults1.itemClicked.connect(
            lambda item: self.update_selection(
                clicked_item=item, clicked_tree=self.twResults1))
        self.twResults2.itemClicked.connect(
            lambda item: self.update_selection(
                clicked_item=item, clicked_tree=self.twResults2))

        # fill the mass sets combobox
        for set_id, mass_set in enumerate(configure.mass_sets):
            self.cbMassSet.addItem(mass_set)
            self.cbMassSet.setItemData(
                set_id,
                configure.mass_sets[mass_set].get("description", ""),
                Qt.ToolTipRole)
        self.cbMassSet.currentTextChanged.connect(self.choose_mass_set)

        # status bar: (1) progress bar
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

        # (2) mass information
        sb_widget = QWidget(self)
        sb_layout = QHBoxLayout(sb_widget)
        sb_layout.setContentsMargins(0, 0, 0, 0)
        sb_widget.setLayout(sb_layout)
        self.lbProteinMass = QLabel()
        self.lbKnownModMass = QLabel()
        self.lbTotalMass = QLabel()
        sb_layout.addWidget(self.lbProteinMass)
        sb_layout.addSpacing(20)
        sb_layout.addWidget(self.lbKnownModMass)
        sb_layout.addSpacing(20)
        sb_layout.addWidget(self.lbTotalMass)
        self.statusbar.addWidget(sb_widget)

        # mass spectrum plot
        self.spectrum_fig = Figure(dpi=100, frameon=False,
                                   tight_layout={"pad": 0}, edgecolor="white")
        self.spectrum_canvas = FigureCanvas(self.spectrum_fig)
        self.vlSpectrumView.addWidget(self.spectrum_canvas)
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

        # button group for spectrum interaction mode
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

        # headers with filters for the results trees
        for tree_widget in (self.twResults1, self.twResults2):
            header = FilterHeader(tree_widget)
            header.setDefaultSectionSize(80)
            header.setHighlightSections(True)
            header.setMinimumSectionSize(0)
            tree_widget.setHeader(header)

        # statistics table
        self.tbStatistics.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeToContents)
        self.tbStatistics.horizontalHeader().setDefaultAlignment(
            Qt.AlignRight)
        self.tbPolymers.verticalHeader().setSectionResizeMode(
            QHeaderView.Fixed)

        # menu for save results button
        self.save_results_menu = {"tree": QMenu(), "table": QMenu()}
        self.save_results_menu["tree"].addAction(
            "Save all entries …",
            lambda: self.save_search_results("all"))
        self.save_results_menu["tree"].addAction(
            "Save checked entries …",
            lambda: self.save_search_results("checked"))
        self.save_results_menu["tree"].addAction(
            "Save checked entries with parents …",
            lambda: self.save_search_results("checked_parent"))
        self.save_results_menu["table"].addAction(
            "Save statistics table (wide format) …",
            lambda: self.save_search_results("stats_wide"))
        self.save_results_menu["table"].addAction(
            "Save statistics table (long format) …",
            lambda: self.save_search_results("stats_long"))
        self.set_save_results_menu(1)

        # private members
        self._current_selection = [0]  # list of currently selected peaks
        self._disulfide_mass = 0  # mass of the current number of disulfides
        self._exp_mass_data = None  # peak list (mass + relative abundance)
        self._known_mods_mass = 0  # mass of known modification
        self._monomer_hits = None  # results from the monomer search
        self._path = configure.path  # last path selected in a file dialog
        self._polymer_hits = None  # results from the polymer search
        self._protein_mass = 0  # mass of the current Protein object
        self._results_tree_headers = [[], []]  # results tables headers
        self._search_statistics = None  # search statistics


    def _monomer_table_create_row(self, row_id, active=False, name="",
                                  composition="", min_count=0, max_count=-1):
        """
        Create a new row in the table of modifications.
        This row describes a modification with the given parameters.

        Private function, to be called by general row insertion functions and
        loaders for modification libraries.

        :param int row_id: Row index passed to :meth:`QTableWidget.insertRow()`
        :param bool active: True if the monomer is to be used in search stage 1
        :param str name: name of the modification
        :param str composition: composition or mass of the modification
        :param int min_count: minimum number of modifications
        :param int max_count: maximum number of modifications
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
        min_spinbox.setStyleSheet(configure.spin_box_flat_style(bg="white"))
        # noinspection PyUnresolvedReferences
        min_spinbox.valueChanged.connect(self.calculate_mod_mass)
        self.tbMonomers.setCellWidget(row_id, 3, min_spinbox)

        max_spinbox = QSpinBox()
        max_spinbox.setMinimum(-1)
        max_spinbox.setSpecialValueText("inf")
        max_spinbox.setFrame(False)
        max_spinbox.setValue(max_count)
        max_spinbox.setStyleSheet(configure.spin_box_flat_style(bg="white"))
        # noinspection PyUnresolvedReferences
        max_spinbox.valueChanged.connect(self.calculate_mod_mass)
        self.tbMonomers.setCellWidget(row_id, 4, max_spinbox)

        self.tbMonomers.blockSignals(False)


    def _polymer_table_create_row(self, row_id, active=True, name="",
                                  composition="", sites="", abundance=0.0):
        """
        Create a new row in the table of modifications.
        This row describes a modification with the given parameters.

        Private function, to be called by general row insertion functions and
        loaders for modification libraries.

        :param int row_id: row index passed to :meth:`QTableWidget.insertRow()`
        :param bool active: True if the polymer is to be used in search stage 2
        :param str name: name of the modification
        :param str composition: composition or mass of the modification
        :param str sites: identifier for the modification site;
                          empty means any site
        :param float abundance: relative abundance
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
        abundance_spinbox.setStyleSheet(
            configure.spin_box_flat_style(double=True))
        self.tbPolymers.setCellWidget(row_id, 4, abundance_spinbox)

        self.tbPolymers.blockSignals(False)


    def table_to_df(self, which="monomers"):
        """
        Create a dataframe from the contents of the monomer or polymer table.

        :param str which: table to convert, "monomers" or "polymers"
        :return: a dataframe
        :rtype: pd.DataFrame
        :raises ValueError: if an unknown string is passed to :arg:`which`
        """

        if which == "monomers":
            columns = ["Checked", "Name", "Composition", "Mass", "Min", "Max"]
            dtypes = [bool, str, str, np.float64, np.int64, np.int64]
            df = pd.DataFrame(self.calculate_mod_mass(), columns=columns)
        elif which == "polymers":
            columns = ["Checked", "Name", "Composition", "Sites", "Abundance"]
            dtypes = [bool, str, str, str, np.float64]
            df = pd.DataFrame(self.get_polymers(), columns=columns)
        else:
            raise ValueError("Invalid value for 'which': {}".format(which))

        return df.astype(dict(zip(columns, dtypes)))


    def save_table(self, dialog_title=None, which="monomers"):
        """
        Export the contents of the monomer or polymer table.

        :param str dialog_title: title of the file dialog
        :param str which: table to export; "monomers" or "polymers"
        :return: nothing
        """

        filename, self._path = get_filename(
            self, "save", dialog_title, self._path, ["csv"])
        if filename is None:
            return

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

        :param df: a :class:`pandas.DataFrame`
        :param table_widget: :class:`QTableWidget` to fill with values
        :param cols: list of (column header, default value) tuples,
                     sorted according to the order of the arguments to
                     _[monomer/polymer]_table_create_row.
                     If the default value is :class:`None`, a column must
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

        table_clear(table_widget)
        for row_id, data in df.iterrows():
            table_widget.create_row(row_id, *[data[c[0]] for c in cols])
        self.calculate_mod_mass()


    def load_table(self, default=False, subdir=None, dialog_title=None,
                   extensions=None, table_widget=None, cols=None):
        """
        Import the contents of the monomer/polymer table.

        :param default: true if a default monomer library should be loaded
                        false if the user should choose a monomer library file
        :param subdir: directory in config containing the default libraries
        :param dialog_title: title of the file dialog
        :param extensions: list of file extensions
                           for :func:`file_extensions()`
        :param table_widget: :class:`QTableWidget` to fill with values
        :param cols: see parameter ``cols``
                     in :meth:`~MainWindow.table_from_df`
        :return: nothing
        """

        if default:
            filename = os.path.join(data_dir, subdir,
                                    self.sender().text() + ".csv")
        else:
            filename, self._path = get_filename(
                self, "open", dialog_title, self._path, extensions)
        if filename is None:
            return

        if filename.endswith("csv"):
            df = pd.read_csv(filename, keep_default_na=False)
        elif filename.endswith("xls") or filename.endswith("xlsx"):
            df = pd.read_excel(filename, keep_default_na=False)
        else:  # TODO detect BPF results
            df = io_tools.read_bpf_library(filename)
        self.table_from_df(df, table_widget, cols)


    def show_about(self):
        """
        Show the about dialog.

        :return: nothing
        """

        QMessageBox.about(self, "About ModFinder", _version_info)


    def set_save_results_menu(self, index):
        """
        Set the appropriate context menu of the save results button.

        :param int index: current index of the results tab widget
        :return: nothing
        """

        if index == 2:
            self.btSaveResults.setMenu(self.save_results_menu["table"])
        else:
            self.btSaveResults.setMenu(self.save_results_menu["tree"])


    def choose_tolerance_units(self):
        """
        Adjust the settings of the tolerance spin box
        when ppm or Da are selected.

        :return: nothing
        """

        if self.cbTolerance.currentText() == "Da":
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


    def load_sequence(self):
        """
        Opens a FASTA file and displays its contents in ``self.teSequence``.

        Changes ``self._path`` to directory of selected file
        and text of ``self.teSequence`` to contents of input file.

        :return: nothing
        """

        filename, self._path = get_filename(
            self, "open", "Open sequence", self._path, ["fasta"])
        if filename is None:
            return

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


    def save_sequence(self):
        """
        Save the contents of ``self.teSequence`` to a FASTA file.

        :return: nothing
        """

        filename, self._path = get_filename(
            self, "save", "Save sequence", self._path, ["fasta"])
        if filename is None:
            return

        try:
            with open(filename, "w") as f:
                f.write(self.teSequence.toPlainText())
        except OSError:
            QMessageBox.critical(
                self,
                "Error",
                "Error while writing to {}. No output saved.".format(filename))


    def load_spectrum(self):
        """
        Open a mass list and display its contents.

        :return: nothing
        """

        filename, self._path = get_filename(
            self, "open", "Open mass list", self._path, ["xls", "xlsx", "csv"])
        if filename is None:
            return

        mass_data = mass_tools.read_massfile(filename)
        if mass_data is None:
            QMessageBox.critical(self,
                                 "Error",
                                 "Could not load mass file.")
            return

        self._exp_mass_data = mass_data
        self.twResults1.clear()
        self.twResults2.clear()
        self._monomer_hits = None
        self._polymer_hits = None
        self.draw_spectrum()


    def save_spectrum(self):
        """
        Save the spectrum as graphics file.

        :return: nothing
        """

        file_types = self.spectrum_canvas.get_supported_filetypes()
        filename, self._path = get_filename(
            self, "save", "Save spectrum", self._path,
            sorted(file_types.keys()), file_types)
        if filename is None:
            return

        try:
            self.spectrum_fig.savefig(filename, dpi=300)
        except OSError:
            QMessageBox.critical(
                self,
                "Error",
                "Error when writing to " + filename + OSError.args)


    def save_search_results(self, mode):
        """
        Write the search results to a CSV file.

        :param str mode: Specifies which results should be exported.
                         Possible choices:

             * ``"checked"``: checked entries in results tree
             * ``"checked_parent"``: checked entries in results tree
                                     with (partially checked) parents
             * ``"all"``: all entries in results tree
             * ``"stats_wide"``: statistics table in wide format (as shown)
             * ``"stats_long"``: statistics table in long (tidy) format
        :return: nothing
        """

        filename, self._path = get_filename(
            self, "save", "Save results", self._path, ["csv"])
        if filename is None:
            return

        # retrieve general information about the search parameters
        out = ["# Combinatorial search results by ModFinder\n",
               "# Date: " + time.strftime("%c") + "\n",
               "#   Sequence: "
               + self.teSequence.toPlainText().replace("\n", r"\n")
               + "\n",
               "#   Masses: "
               + ", ".join(list(self._exp_mass_data["avg_mass"].astype(str)))
               + "\n",
               "#   Tolerance: {:.2f} {}\n".format(
                   self.sbTolerance.value(),
                   self.cbTolerance.currentText()),
               "#   Composition:\n"]

        for (is_checked, name, _, mass,
             min_count, max_count) in self.calculate_mod_mass():
            if is_checked:
                if max_count == -1:
                    max_count = "inf"
                out += ["#     {} ({:.2f} Da), ".format(name, mass),
                        "min {}, ".format(min_count),
                        "max {}\n".format(max_count)]

        out.append("#   Structures:\n")
        for (is_checked, name, composition,
             sites, abundance) in self.get_polymers():
            if is_checked:
                out += ["#     {} ({}); ".format(name, composition),
                        "sites: {}; ".format(sites),
                        "abundance: {:.2f}\n".format(abundance)]

        if mode.startswith("stats"):
            out.append("# Search statistics:\n")

            df = pd.DataFrame(
                [[self.tbStatistics.item(row_id, col_id).text()
                  for col_id in range(self.tbStatistics.columnCount())]
                 for row_id in range(self.tbStatistics.rowCount())],
                columns=[self.tbStatistics.horizontalHeaderItem(c).text()
                         for c in range(self.tbStatistics.columnCount())])
            if mode.endswith("long"):
                df = df.melt(id_vars=["ID",
                                      "Exp. Mass",
                                      "%"],
                             value_vars=["Search space size",
                                         "Stage 1 results",
                                         "Stage 2 results",
                                         "Stage 2 uniques"],
                             value_name="Measure",
                             var_name="Value")
        else:
            # set up the BFS algorithm
            def explore(node):
                yield ([node.checkState(0)]
                       + [node.text(c) for c in range(1, node.columnCount())])
                for i in range(node.childCount()):
                    yield from explore(node.child(i))

            # choose the appropriate data source
            if self.taResults.currentIndex() == 0:
                out.append("# Composition search (stage 1), ")
                results_tree = self.twResults1
            else:
                out.append("# Structure search (stage 2), ")
                results_tree = self.twResults2

            # create dataframe from (selected) entries
            header = results_tree.headerItem()
            df = (
                pd.DataFrame(
                    explore(results_tree.invisibleRootItem()),
                    columns=(["checked"]
                             + [header.text(c)
                                for c in range(1, header.columnCount())]))
                .replace("", np.nan)
                .ffill())
            if not df.empty and df.iloc[0, 1] is None:
                df = df.drop(0)
            df = df.drop("checked", axis=1)

            if mode == "all":
                out.append("all results.\n")
            elif mode == "checked_parent":
                out.append("checked results (with parents).\n")
                df = df[df.checked != Qt.Unchecked]
            else:
                out.append("checked results.\n")
                df = df[df.checked == Qt.Checked]

        # write to file
        try:
            with open(filename, "w") as f:
                f.write("".join(out))
                df.to_csv(f, index=False)
        except OSError:
            QMessageBox.critical(
                self,
                "Error",
                "Error while writing to {}. No output saved.".format(filename))


    def choose_mass_set(self, mass_set):
        """
        Chooses the set of atomic weights to use for mass calculations.

        :param str mass_set: mass set to choose
        :return: nothing
        """

        configure.select_mass_set(mass_set)
        if self.teSequence.toPlainText():
            self.calculate_protein_mass()
            self.calculate_mod_mass()


    def calculate_protein_mass(self):
        """
        Calculates the mass of the protein from the current data.
        Changes ``self._protein_mass`` and ``self._disulfide_mass``.

        :return: nothing
        """

        chains, sequence = sequence_tools.read_fasta_string(
            self.teSequence.toPlainText())
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
        self._protein_mass = protein.mass_without_disulfides
        self._disulfide_mass = (mass_tools.Formula("H-2").mass
                                * self.sbDisulfides.value())
        self.update_mass_in_statusbar()


    def calculate_mod_mass(self):
        """
        Calculate the mass of known modifications.

        Changes ``self._known_mods_mass`` to the mass of known modifications.

        :return: list of (checked, name, composition,
                          mass, min count, max count) tuples
        :rtype: list(tuple(bool, str, str, float, int, int))
        """

        self._known_mods_mass = 0
        result = []

        # add min counts for monomers to the theoretical mass
        for row_id in range(self.tbMonomers.rowCount()):
            # extract the variables
            ch = self.tbMonomers.cellWidget(row_id, 0).findChild(QCheckBox)
            name = self.tbMonomers.item(row_id, 1).text()
            composition = self.tbMonomers.item(row_id, 2).text().strip()
            min_count = self.tbMonomers.cellWidget(row_id, 3).value()
            max_count = self.tbMonomers.cellWidget(row_id, 4).value()
            mass = 0

            # get the mass from the formula cell
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

            # set the mass tooltip and color the cell
            if error_in_formula or mass == 0:
                bg_color = QColor(255, 225, 225)
                tooltip = ""
            else:
                bg_color = QColor(255, 255, 255)
                tooltip = "{:.2f} Da".format(mass)
            self.tbMonomers.item(row_id, 2).setBackground(bg_color)
            self.tbMonomers.item(row_id, 2).setToolTip(tooltip)

            # color the Min cell if its value exceeds the one of Max
            if min_count > max_count != -1:
                style = configure.spin_box_flat_style(bg="red")
            else:
                style = configure.spin_box_flat_style(bg="white")
            self.tbMonomers.cellWidget(row_id, 3).setStyleSheet(style)

            if ch.isChecked():
                self._known_mods_mass += mass * min_count
            result.append((ch.isChecked(), name, composition, mass,
                           min_count, max_count))

        self.update_mass_in_statusbar()
        return result


    def update_mass_in_statusbar(self):
        """
        Updates the contents of the mass labels in the status bar.

        :return: nothing
        """

        protein_mass = "<b>Protein:</b> {:.2f} Da".format(
            self._protein_mass)
        mod_mass = "<b>known modifications:</b> {:.2f} Da".format(
            self._known_mods_mass + self._disulfide_mass)
        total_mass = "<b>total:</b> {:.2f} Da".format(
            self._protein_mass + self._disulfide_mass + self._known_mods_mass)
        self.lbProteinMass.setText(protein_mass)
        self.lbKnownModMass.setText(mod_mass)
        self.lbTotalMass.setText(total_mass)


    def get_polymers(self):
        """
        Read the current glycan library from the table and return as list.

        :return: list of (is_checked, name, compos., sites, abundance) tuples
        :rtype: list(tuple(bool, str, str, str, float))
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
                {"avg_mass": self.sbSingleMass.value(),
                 "rel_abundance": 100.0},
                index=[0])
            self.draw_spectrum()

        if self._exp_mass_data is None:
            QMessageBox.critical(
                self, "Error", "No mass list loaded. Aborting search.")
            return

        monomers = [(m[1], m[3], m[4], m[5])
                    for m in self.calculate_mod_mass()
                    if m[0]]
        modifications = []  # list of modifications for search stage 1
        explained_mass = (self._protein_mass
                          + self._disulfide_mass
                          + self._known_mods_mass)
        unknown_masses = (self._exp_mass_data["avg_mass"]
                          - explained_mass)  # type: pd.DataFrame

        if self.cbTolerance.currentText() == "Da":
            # calculate largest mass plus tolerance
            max_tol_mass = (max(self._exp_mass_data["avg_mass"])
                            + self.sbTolerance.value())
            mass_tolerance = self.sbTolerance.value()
        else:  # ppm
            # calculate a mass tolerance for each peak
            max_tol_mass = (max(self._exp_mass_data["avg_mass"])
                            * (1 + self.sbTolerance.value() / 1000000))
            mass_tolerance = []
            for _, m in self._exp_mass_data["avg_mass"].iteritems():
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
                search_tools.get_monomers_from_library(df_polymers))

            monomers_for_polymer_search = [m for m in available_monomers
                                           if m in monomers_in_library]
            try:
                polymer_combs = search_tools.calc_polymer_combinations(
                    df_polymers,
                    monomers_for_polymer_search,
                    self.pbSearchProgress)
            except ValueError as e:
                QMessageBox.critical(self, "Error", str(e))
                return

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
                        int((max_tol_mass
                             - self._protein_mass
                             - self._disulfide_mass)
                            / mass),
                        configure.maxmods)
            modifications.append((name, mass, max_count - min_count))

        if not modifications:
            QMessageBox.critical(
                self, "Error", "List of modifications is empty. "
                               "Aborting search.")
            return

        # check whether all monomers required in search stage 2 are available
        missing_monomers = monomers_in_library - set(available_monomers)
        if missing_monomers:
            error_message = [
                "The following modifications appear in the library ",
                "but do not appear in the composition list: "]
            error_message += " ".join(sorted(missing_monomers))
            error_message.append(".\n")
            QMessageBox.critical(self,
                                 "Error",
                                 "".join(error_message))
            return

        # stage 1: monomer search
        self.statusbar.showMessage(
            "Performing composition search (stage 1) ...")
        (self._monomer_hits,
         self._search_statistics) = search_tools.find_monomers(
            modifications,
            list(unknown_masses),
            mass_tolerance=mass_tolerance,
            explained_mass=explained_mass,
            progress_bar=self.pbSearchProgress)

        if self._monomer_hits is None:
            # the modification search was not successful
            QMessageBox.critical(
                self, "Error", "Composition search was unsuccessful.")
            return

        # add the minimum monomer counts to the result data frame
        for name, _, min_count, _ in monomers:
            self._monomer_hits[name] += min_count

        # stage 2: polymer search
        if polymers:
            self.statusbar.showMessage(
                "Performing structure search (stage 2) ...")
            (self._polymer_hits,
             stage2_statistics) = search_tools.find_polymers(
                self._monomer_hits,
                polymer_combinations=polymer_combs,
                monomers=monomers_for_polymer_search,
                progress_bar=self.pbSearchProgress)
            self._search_statistics = (
                self._search_statistics
                .join(stage2_statistics)
                .fillna(0)
                .astype(int))
        self.fill_results_tables()
        self.tbStatistics.clearContents()
        self.tbStatistics.setRowCount(0)
        self._search_statistics.apply(self.fill_statistics_table, axis=1)


    def draw_spectrum(self):
        """
        Show a mass spectrum after loading a peak list.

        :return: nothing
        """

        self.spectrum_fig.clear()
        self.spectrum_axes = self.spectrum_fig.add_subplot(111)
        self.spectrum_axes.set_xmargin(.02)
        self.spectrum_peak_lines = self.spectrum_axes.vlines(
            x=self._exp_mass_data["avg_mass"],
            ymin=0,
            ymax=self._exp_mass_data["rel_abundance"],
            linewidth=1,
            color=configure.colors["unselect_no_annotation"],
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

        :param MouseEvent start: describes the coordinates
                                 where the mouse button was pressed
        :param MouseEvent stop: describes the coordinates
                                where the mouse button was released
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

        if self.chDelta1.isChecked() and self.chDelta2.isChecked():
            self.chCombineDelta.setEnabled(True)
        else:
            self.chCombineDelta.setEnabled(False)

        self.update_selection()


    def select_peaks_by_pick(self, event):
        """
        Select a peak picked by a mouseclick on the spectrum.

        :param PickEvent event: pick event from the canvas
        :return: nothing
        """
        if event.mouseevent.button == 1:
            self.update_selection(event.ind)


    def select_peaks_by_span(self, min_mass, max_mass):
        """
        Select all peaks that fall within a selected span.

        :param float min_mass: lower end of the span
        :param float max_mass: upper end of the span
        :return: nothing
        """

        peak_indices = self._exp_mass_data.index[
                (min_mass <= self._exp_mass_data["avg_mass"])
                & (self._exp_mass_data["avg_mass"] <= max_mass)
            ].tolist()

        if peak_indices:
            self.update_selection(peak_indices)


    def update_selection(self, new_selection=None,
                         clicked_item=None, clicked_tree=None):
        """
        Update the spectrum and the single mass spinbox
        after the selection has changed.

        :param list new_selection: list of indices of selected peaks
        :param SortableTreeWidgetItem clicked_item: item that was clicked
        :param QTreeWidget clicked_tree: results tree whose item was clicked
        :return: the DataFrame whose contents are shown in the results table
        :rtype: pd.DataFrame
        """

        if self._exp_mass_data is None:  # there's no spectrum
            return

        # if the selection did not change, use the old one
        if new_selection:
            self._current_selection = new_selection[:]

        # an item of a results tree was clicked
        # for top-level items, proceed with selection
        # for child items, do not change the selection in the clicked table
        scroll_tree1 = True
        scroll_tree2 = True
        if clicked_tree:
            item_index = (clicked_tree.invisibleRootItem()
                          .indexOfChild(clicked_item))
            if item_index == -1:
                item_index = clicked_tree.invisibleRootItem().indexOfChild(
                    clicked_item.getTopParent())
                if clicked_tree == self.twResults1:
                    scroll_tree1 = False
                else:
                    scroll_tree2 = False
            self._current_selection = [item_index]


        # (1) update selection in the spectrum
        central_peak = self._current_selection[0]
        if self.bgSpectrum.checkedButton() == self.btModeSelection:
            self.highlight_selected_peaks(self._current_selection)
        else:
            self._current_selection = self.highlight_delta_series(central_peak)

        # (2) fill the single mass spin box with the currently selected mass
        try:
            self.sbSingleMass.setValue(
                self._exp_mass_data.loc[central_peak, "avg_mass"])
        except AttributeError:  # occurs when second peak file is loaded
            pass

        # (3) scroll to item in results trees
        for tree, scroll in [(self.twResults1, scroll_tree1),
                             (self.twResults2, scroll_tree2)]:
            if scroll:
                try:
                    for i in range(tree.invisibleRootItem().childCount()):
                        tree.invisibleRootItem().child(i).setSelected(False)
                    item = tree.invisibleRootItem().child(central_peak)
                    tree.scrollToItem(item)
                    item.setSelected(True)
                except AttributeError:
                    pass


    def find_delta_peaks(self, query_peak, delta, tolerance, iterations):
        """
        Find peaks separated from a given peak
        by a multiple of a given mass difference.

        :param int query_peak: index of the peak at the center of the series
        :param float delta: mass difference
        :param float tolerance: tolerance for finding peaks in the series
        :param int iterations: maximum number of peaks to find at each
                               side of the main peak
        :return: a Series with information where a peak was found
        :rtype: pd.Series
        """

        if iterations == -1:
            iterations = int(delta / tolerance / 2)
        intervals = {}  # a {number of differences: (start, end)} dict

        main_mass = float(self._exp_mass_data.iloc[query_peak]["avg_mass"])
        min_mass = float(min(self._exp_mass_data["avg_mass"]))
        max_mass = float(max(self._exp_mass_data["avg_mass"]))

        # calculate putative intervals
        # increase the interval size by 2 * tolerance in each iteration
        current_mass = main_mass
        current_tolerance = tolerance
        i = 1
        while i <= iterations and current_mass > min_mass:
            current_mass -= delta
            intervals[str(-i)] = (current_mass - current_tolerance,
                                  current_mass + current_tolerance)
            current_tolerance += tolerance
            i += 1

        current_mass = main_mass
        current_tolerance = tolerance
        i = 1
        while i <= iterations and current_mass < max_mass:
            current_mass += delta
            intervals[str(i)] = (current_mass - current_tolerance,
                                 current_mass + current_tolerance)
            current_tolerance += tolerance
            i += 1

        interval_per_peak = (self._exp_mass_data["avg_mass"]
                             .apply(find_in_intervals,
                                    intervals=intervals))
        interval_per_peak[query_peak] = "0"
        return interval_per_peak


    def highlight_delta_series(self, central_peak):
        """
        Highlights a series of peaks that differ by a given mass.

        :param int central_peak: index of the central peak
        :return: list of peak indices in the delta series
        :rtype: list(int)
        """

        if self._exp_mass_data is None:  # there's no spectrum
            return []

        # a dataframe that shares its index with self._exp_mass_data
        # columns 1 and 2 indicate how many delta masses each peak is away
        # from the selected peak in a series
        df_counts = pd.DataFrame(
            index=self._exp_mass_data.index,
            dtype=str)

        # column "1" is straight forward
        if self.chDelta1.isChecked():
            df_counts["1"] = self.find_delta_peaks(
                central_peak,
                self.sbDeltaValue1.value(),
                self.sbDeltaTolerance1.value(),
                self.sbDeltaRepetition1.value())
        else:
            df_counts["1"] = ""

        # column "2" has to take into account whether the second delta series
        # should start at each peak of the first one
        if self.chDelta2.isChecked():
            if self.chDelta1.isChecked() and self.chCombineDelta.isChecked():
                delta_series = []
                for peak_id in np.flatnonzero(df_counts["1"]):
                    subseries = self.find_delta_peaks(
                        peak_id,
                        self.sbDeltaValue2.value(),
                        self.sbDeltaTolerance2.value(),
                        self.sbDeltaRepetition2.value())
                    df_counts.loc[np.where(subseries)[0], "1"] = (
                        df_counts.loc[peak_id, "1"])  # correct labels
                    delta_series.append(subseries)
                df_counts["2"] = (pd.concat(delta_series, axis=1)
                                  .apply(lambda x: ", ".join(filter(None, x)),
                                         axis=1))
            else:
                df_counts["2"] = self.find_delta_peaks(
                    central_peak,
                    self.sbDeltaValue2.value(),
                    self.sbDeltaTolerance2.value(),
                    self.sbDeltaRepetition2.value())
        else:
            df_counts["2"] = ""

        # calculate "label" and "color" columns; values of the latter:
        # 0 - peaks not in any delta mass series
        # 1 - peaks in series 1
        # 2 - peaks in series 2
        # 3 - peaks in both series
        # 4 - selected (central) peak
        df_counts["label"] = df_counts.apply(
            lambda x: "/".join(filter(None, x)), axis=1)
        if (self.chDelta1.isChecked()
                and self.chDelta2.isChecked()
                and self.chCombineDelta.isChecked()):
            # calculate colors for main and sub delta series
            df_counts["color"] = (
                (
                    (df_counts["1"] != "")
                    & (df_counts["2"] == "0")
                ) * 1
                + (
                    (df_counts["1"] != "")
                    & (df_counts["2"] != "0")
                ) * 2
            )
        else:
            # calculate colors for both primary delat series
            df_counts["color"] = ((df_counts["1"] != "") * 1
                                  + (df_counts["2"] != "") * 2)
        df_counts.loc[central_peak, "color"] = 4

        # color the peaks in the delta series and increase their line width
        color_set = np.array([configure.colors["delta_other"],
                              configure.colors["delta_1"],
                              configure.colors["delta_2"],
                              configure.colors["delta_both"],
                              configure.colors["delta_main"]])

        lw_set = np.array([1, 2, 2, 2, 3])
        self.spectrum_peak_lines.set_color(color_set[df_counts["color"]])
        self.spectrum_peak_lines.set_linewidth(lw_set[df_counts["color"]])

        # annotate the peaks in the delta series
        for annotation in self.spectrum_axes.findobj(
                matplotlib.text.Annotation):
            annotation.remove()

        peaks_in_series = list(np.flatnonzero(df_counts["color"]))
        for peak_id in peaks_in_series:
            label = df_counts["label"][peak_id]
            if self.btLabelPeaks.isChecked():
                label += " ({:.2f})".format(
                    self._exp_mass_data["avg_mass"][peak_id])
            self.spectrum_axes.annotate(
                s=label,
                xy=(self._exp_mass_data.iloc[peak_id]["avg_mass"],
                    self._exp_mass_data.iloc[peak_id]["rel_abundance"]),
                xytext=(0, 5),
                textcoords="offset pixels",
                horizontalalignment="center",
                bbox=dict(facecolor="white", alpha=.75,
                          linewidth=0, pad=0))
        self.spectrum_canvas.draw()

        return peaks_in_series


    def highlight_selected_peaks(self, peak_indices):
        """
        Highlight selected peaks in the spectrum.

        :param list(int) peak_indices: list of indices of peaks
                                       that should be highlighted
        :return: nothing
        """

        peaks_with_result = np.full(self._exp_mass_data.shape[0], 2, dtype=int)
        try:
            peaks_with_result[self._polymer_hits
                                  .swaplevel(0, 1)
                                  .loc[-1]
                                  .index.labels[0]] = 0
        except KeyError:
            pass  # results for all peaks found
        except AttributeError:
            try:
                peaks_with_result[self._monomer_hits
                                      .swaplevel(0, 1)
                                      .loc[-1]
                                      .index.labels[0]] = 0
            except KeyError:
                pass  # results for all peaks found
            except AttributeError:
                peaks_with_result = np.zeros(
                    self._exp_mass_data.shape[0], dtype=int)

        selected_peaks = np.zeros(self._exp_mass_data.shape[0], dtype=int)
        selected_peaks[peak_indices] = 1

        # peak colors will be an array with one entry per peak:
        # no polymers: 0 - not selected, 1 - selected
        # polymers:    2 - not selected, 3 - selected
        peak_colors = selected_peaks + peaks_with_result
        color_set = np.array([configure.colors["unselect_no_annotation"],
                              configure.colors["select_no_annotation"],
                              configure.colors["unselect_annotation"],
                              configure.colors["select_annotation"]])
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
                                      .iloc[peak_id]["avg_mass"]),
                    xy=(self._exp_mass_data
                            .iloc[peak_id]["avg_mass"],
                        self._exp_mass_data
                            .iloc[peak_id]["rel_abundance"]),
                    xytext=(0, 5),
                    textcoords="offset pixels",
                    horizontalalignment="center",
                    bbox=dict(facecolor="white", alpha=.75,
                              linewidth=0, pad=0))
        self.spectrum_canvas.draw()


    def _fill_results_table(self, stage, tree_widget, cols, df_hit):
        """
        Fills a results table with rows.

        :param int stage: search stage (0 = stage 1, 1 = stage 2)
        :param QTreeWidget tree_widget: results table
        :param list(list(str), list(str), list(str)) cols: column headers
        :param pd.DataFrame df_hit: data for rows
        :return: nothing
        """

        tree_widget.clear()
        tree_widget.setUpdatesEnabled(False)

        # set column headers
        mono_columns = list(
            df_hit.columns[:df_hit.columns.get_loc("Exp_Mass")])
        self._results_tree_headers[stage] = (cols[0]
                                             + mono_columns
                                             + cols[1]
                                             + cols[2])
        column_count = len(self._results_tree_headers[stage])
        tree_widget.setColumnCount(column_count)
        tree_widget.setHeaderLabels(self._results_tree_headers[stage])
        tree_widget.header().setDefaultAlignment(Qt.AlignRight)

        self.pbSearchProgress.setValue(0)
        index_length = len(self._exp_mass_data.index)
        for index_count, mass_index in enumerate(self._exp_mass_data.index):
            # generate root item (experimental mass, relative abundance)
            root_item = SortableTreeWidgetItem(
                tree_widget,
                default_key=float,
                col_key=_default_col_key)
            root_item.setFlags(root_item.flags() | Qt.ItemIsTristate)
            root_item.setCheckState(0, Qt.Unchecked)
            root_item.setText(1, "{}".format(mass_index))
            root_item.setTextAlignment(1, Qt.AlignLeft)
            root_item.setText(
                2, "{:.2f}".format(self._exp_mass_data
                                   .loc[mass_index, "avg_mass"]))
            root_item.setTextAlignment(2, Qt.AlignRight)
            root_item.setText(
                3, "{:.1f}".format(self._exp_mass_data
                                   .loc[mass_index, "rel_abundance"]))
            root_item.setTextAlignment(3, Qt.AlignRight)

            # color root item
            if df_hit.loc[mass_index].index.values[0][0] == -1:
                if mass_index % 2 == 0:
                    root_item.setTotalBackground(
                        QColor(configure.colors["table_root_missing_even"]),
                        column_count)
                else:
                    root_item.setTotalBackground(
                        QColor(configure.colors["table_root_missing_odd"]),
                        column_count)
            else:
                if mass_index % 2 == 0:
                    root_item.setTotalBackground(
                        QColor(configure.colors["table_root_annotated_even"]),
                        column_count)
                else:
                    root_item.setTotalBackground(
                        QColor(configure.colors["table_root_annotated_odd"]),
                        column_count)

                # create child items
                df_hit.loc[mass_index].pipe(create_child_items,
                                            root_item,
                                            column_count,
                                            mono_columns,
                                            cols[2])
            self.pbSearchProgress.setValue(int((index_count + 1)
                                               / index_length * 100))

        tree_widget.header().setSectionResizeMode(
            QHeaderView.ResizeToContents)
        tree_widget.header().setSectionResizeMode(0, QHeaderView.Fixed)
        tree_widget.header().setStretchLastSection(False)
        tree_widget.header().setStyleSheet(
            """
            QHeaderView::section {
                padding-top: 2px;
                padding-bottom: 2px;
                padding-left: 12px;
            }
            """)
        tree_widget.setUpdatesEnabled(True)

        # create filters
        self.taResults.setCurrentIndex(2)  # filters not drawn on current tab
        filter_index = [self._results_tree_headers[stage].index("ppm") + i + 1
                        for i in range(df_hit.columns.get_loc("Exp_Mass"))]
        tree_widget.header().setFilterBoxes(filter_index)
        tree_widget.header().filterChanged.connect(
            lambda: self.filter_results_table(stage))
        self.taResults.setCurrentIndex(1)


    def fill_results_tables(self):
        """
        Populate both results tables with rows.
        Prepare column labels and then call the actual function that
        adds rows to each results table.

        :return: nothing
        """

        self.update_selection()

        if self._monomer_hits is not None:
            self.statusbar.showMessage(
                "Preparing results table for stage 1 ...")
            cols = [
                ["", "ID", "Exp. Mass", "%", "Theo. Mass", "Da", "ppm"],  # hit
                [],  # perm
                []  # site
            ]
            self._fill_results_table(
                stage=0,
                tree_widget=self.twResults1,
                cols=cols,
                df_hit=self._monomer_hits)

        if self._polymer_hits is not None:
            self.statusbar.showMessage(
                "Preparing results table for stage 2 ...")
            cols = [
                ["", "ID", "Exp. Mass", "%", "Hit", "Hit Score", "# Perms",
                 "Theo. Mass", "Da", "ppm"],
                ["Perm", "Perm Score"],
                list(
                    self._polymer_hits.columns[
                        self._polymer_hits.columns.get_loc("ppm") + 1:
                        self._polymer_hits.columns.get_loc(
                            "Permutation score")]
                    )
            ]
            self._fill_results_table(
                stage=1,
                tree_widget=self.twResults2,
                cols=cols,
                df_hit=self._polymer_hits)

        self.statusbar.showMessage("Results tables ready!", 5000)


    def expand_results_tree(self, expand=True, depth=0):
        """
        Expand or collapse the results trees.

        :param bool expand: expand (True) or collapse (False) the tree
        :param int depth: expand to this depth
        :return: nothing
        """

        if self.taResults.currentIndex() == 0:
            tree_widget = self.twResults1
        elif self.taResults.currentIndex() == 1:
            tree_widget = self.twResults2
        else:
            return

        if expand:
            tree_widget.expandToDepth(depth)
        else:
            tree_widget.collapseAll()


    def check_results_tree(self, check=True):
        """
        Check or uncheck all entries of a results tree.

        :param bool check: whether to check or uncheck the entries
        :return: nothing
        """
        if self.taResults.currentIndex() == 0:
            tree_widget = self.twResults1
        elif self.taResults.currentIndex() == 1:
            tree_widget = self.twResults2
        else:
            return

        for i in range(tree_widget.invisibleRootItem().childCount()):
                tree_widget.invisibleRootItem().child(i).setCheckState(
                    0, Qt.Checked if check else Qt.Unchecked)


    def filter_results_table(self, stage):
        """
        Filter the results table.

        :param int stage: search stage (0 = stage 1, 1 = stage 2)
        :return: nothing
        """

        if stage == 0:
            tree_widget = self.twResults1
        else:
            tree_widget = self.twResults2
        start_col = self._results_tree_headers[stage].index("ppm") + 1

        # Generate a list of query conditions for filtering the results trees.
        # Each condition is a (lower, upper) tuple, indicating the lower
        # and upper limit for the count of a given monosaccharide.
        re_filter = re.compile("(\d*)(-?)(\d*)")
        query = []
        for _, child in tree_widget.header().allFilters():
            f = re_filter.match(child.text()).groups()
            if "".join(f):
                if f[0] and not f[1] and not f[2]:
                    query.append((int(f[0]), int(f[0])))  # x == value
                elif f[0] and f[1] and not f[2]:
                    query.append((int(f[0]), math.inf))  # x >= lower
                elif f[0] and f[1] and f[2]:
                    query.append((int(f[0]), int(f[2])))  # lower <= x <= upper
                else:
                    query.append((0, int(f[2])))  # x <= upper
            else:
                query.append((0, math.inf))  # accept any value

        # check filter conditions for all level 2 items in the tree widget
        # if all of its level 2 items are hidden, also hide the parent
        root = tree_widget.invisibleRootItem()
        for i in range(root.childCount()):
            only_hidden_children = True
            for j in range(root.child(i).childCount()):
                child = root.child(i).child(j)
                if item_is_shown(query, child, start_col):
                    child.setHidden(False)
                    only_hidden_children = False
                else:
                    child.setHidden(True)

            if only_hidden_children:
                root.child(i).setHidden(True)
            else:
                root.child(i).setHidden(False)


    def fill_statistics_table(self, row):
        """
        Fill the statistics table with data from
        :var:`self._search_statistics`. This method is called
        by :meth:`pd.DataFrame.apply()`.

        :param pd.Series row: row of the statistics dataframe
        :return: nothing
        """

        row_id = self.tbStatistics.rowCount()
        self.tbStatistics.insertRow(row_id)

        # (1) running counter
        item = SortableTableWidgetItem(str(row_id))
        item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
        item.setFlags(Qt.ItemIsEnabled)
        self.tbStatistics.setItem(row_id, 0, item)

        # (2) peak information
        for col_id, label, form in [(1, "avg_mass", "{:.2f}"),
                                    (2, "rel_abundance", "{:.1f}")]:
            item = SortableTableWidgetItem(
                    form.format(self._exp_mass_data.loc[row_id, label]))
            item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            item.setFlags(Qt.ItemIsEnabled)
            self.tbStatistics.setItem(row_id, col_id, item)

        # (3) peak statistics
        for col_id, label in [(3, "search_space_size"),
                              (4, "stage1_results"),
                              (5, "stage2_results"),
                              (6, "stage2_uniques")]:
            try:
                item = SortableTableWidgetItem(str(row[label]))
                item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                item.setFlags(Qt.ItemIsEnabled)
                self.tbStatistics.setItem(row_id, col_id, item)
            except KeyError:
                pass



    def clear_filters(self):
        """
        Clear the contents of the filter line edits.

        :return: nothing
        """

        if self.taResults.currentIndex() == 0:
            self.twResults1.header().clearFilters()
        elif self.taResults.currentIndex() == 1:
            self.twResults2.header().clearFilters()


    def save_settings(self):
        """
        Export the current settings as XML file.

        :return: nothing
        """

        filename, self._path = get_filename(
            self, "save", "Save settings", self._path, ["xml"])
        if filename is None:
            return

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
            ETree.SubElement(root, child).append(
                io_tools.dataframe_to_xml(df))

        io_tools.prettify_xml(root)
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

        filename, self._path = get_filename(
            self, "open", "Load settings", self._path, ["xml"])
        if filename is None:
            return

        try:
            root = ETree.ElementTree().parse(filename)
        except IOError:
            QMessageBox.critical(
                self,
                "Error",
                "Error loading " + filename + OSError.args)
            return

        self.teSequence.setText(root.find("sequence").text)
        disulfides = int(root.find("disulfides").text)
        self.sbDisulfides.setMaximum(disulfides)
        self.sbDisulfides.setValue(disulfides)
        self.chPngase.setChecked(root.find("pngasef").text == "True")
        self.calculate_protein_mass()

        self._monomer_hits = None
        self._polymer_hits = None
        self.twResults1.clear()
        self.twResults2.clear()
        self._exp_mass_data = io_tools.dataframe_from_xml(
            root.find("masslist"))
        if not self._exp_mass_data.empty:
            self.draw_spectrum()

        self.cbTolerance.setCurrentIndex(
            int(root.find("tolerance-flavor").text))
        self.sbTolerance.setValue(float(root.find("tolerance-value").text))

        self.table_from_df(
            io_tools.dataframe_from_xml(root.find("monomers")),
            self.tbMonomers,
            _monomer_table_columns)
        self.table_from_df(
            io_tools.dataframe_from_xml(root.find("polymers")),
            self.tbPolymers,
            _polymer_table_columns)
        self.calculate_mod_mass()


def main():
    """
    Execute the main application loop.

    :return: nothing
    """
    app = QApplication(sys.argv)
    QLocale.setDefault(QLocale.c())
    frame = MainWindow()
    frame.show()
    app.exec_()


if __name__ == "__main__":
    main()
