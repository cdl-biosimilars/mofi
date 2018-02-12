"""
The main module, which includes GUI definitions and major functions.
"""

import math
import os
import re
import sys
import time
import xml.etree.ElementTree as ETree
from urllib.request import pathname2url
import webbrowser

from PyQt5.QtWidgets import (QApplication, QMainWindow, QMenu,
                             QTableWidgetItem, QCheckBox, QMessageBox,
                             QHeaderView, QButtonGroup,
                             QSpinBox, QDoubleSpinBox, QWidget, QHBoxLayout,
                             QProgressBar, QLabel, QSizePolicy)
from PyQt5.QtGui import QColor, QCursor, QIcon, QPixmap
from PyQt5.QtCore import Qt, QLocale, QEvent

import numpy as np
import pandas as pd

import matplotlib
import matplotlib.collections
import matplotlib.markers
import matplotlib.text
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure

from mofi import (configure, io_tools, mass_tools,
                  search_tools, sequence_tools)
from mofi.paths import data_dir, docs_dir
from mofi.mofi_ui import Ui_MoFi
from mofi.widgets import (FilterHeader, CollapsingRectangleSelector,
                          SortableTreeWidgetItem, SortableTableWidgetItem,
                          ImportTabDataDialog, FileTypes, get_filename)

_version_info = """MoFi v1.0

© 2017 Christian Doppler Laboratory
for Innovative Tools for Biosimilar Characterization

Contact: Wolfgang.Skala@sbg.ac.at

Python version:
{}""".format(sys.version)

# default values for columns in the monomer table
_monomer_table_columns = [
    ("Checked", False, "Use?"),
    ("Name", None, "Name"),
    ("Composition", None, "Formula/Mass"),
    ("Min", 0, "Min"),
    ("Max", -1, "Max")
]

# default values for columns in the polymer table
_polymer_table_columns = [
    ("Checked", True, "Use?"),
    ("Name", None, "Name"),
    ("Composition", None, "Composition"),
    ("Sites", "", "Sites"),
    ("Abundance", 0.0, "Abundance")
]

# default values for columns in the polymer table
_pepmap_columns = [
    ("Checked", True, "Use?"),
    ("Modification", None, "Modification"),
    ("Abundance", 0.0, "Abundance")
]

# sorting key for column 1 of the results trees
_default_col_key = {1: lambda x: [int(i) for i in x.split("-")]}

# matplotlib markers used in the spectrum
_delta_symbols = pd.DataFrame.from_records([
    ("off", ""),
    ("Point", "."),
    ("Circle", "o"),
    ("TriangleDown", "v"),
    ("TriangleUp", "^"),
    ("TriangleLeft", "<"),
    ("TriangleRight", ">"),
    ("TriDown", "1"),
    ("TriUp", "2"),
    ("TriLeft", "3"),
    ("TriRight", "4"),
    ("Octagon", "8"),
    ("Square", "s"),
    ("Pentagon", "p"),
    ("PlusFilled", "P"),
    ("Star", "*"),
    ("Hexagon1", "h"),
    ("Hexagon2", "H"),
    ("Plus", "+"),
    ("Cross", "x"),
    ("CrossFilled", "X"),
    ("Diamond", "D"),
    ("ThinDiamond", "d"),
    ("VLine", "|"),
    ("HLine", "_")
], columns=["name", "marker"])

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
            peak_optimum.name = (peak_optimum.name[0],
                                 peak_optimum.name[1])
            create_hit_columns(root_item, peak_optimum, monomers, sites)
        except ValueError:  # if df is empty, like due to filtering
            pass

        # (1) child items for all hits per peak
        for stage2_id, hit in (df.groupby("Stage2_hit")):
            hit_item = SortableTreeWidgetItem(
                root_item,
                default_key=float,
                col_key=_default_col_key)
            hit_optimum = hit.loc[hit["Permutation score"].idxmax()]
            hit_item.setText(
                1, "{}-{}-{}".format(root_item.text(1),
                                     hit_optimum["Stage1_ID"],
                                     hit_optimum["Stage2_ID"]))
            hit_item.setTextAlignment(1, Qt.AlignLeft)
            hit_item.setFlags(hit_item.flags() | Qt.ItemIsTristate)
            hit_item.setCheckState(0, Qt.Unchecked)
            pos = create_hit_columns(hit_item, hit_optimum,
                                     monomers, sites)

            # background color
            hit_item.setTotalBackground(
                even_color=configure.colors["table"]["child_even"],
                odd_color=configure.colors["table"]["child_odd"],
                column_count=column_count,
                alternating=stage2_id)

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
                        1, "{}-{}".format(hit_item.text(1), perm_id))
                    perm_item.setTextAlignment(1, Qt.AlignLeft)
                    perm_item.setCheckState(0, Qt.Unchecked)
                    create_site_columns(perm_item, pos, perm, sites)

                    perm_item.setTotalBackground(
                        even_color=configure.colors["table"][
                            "grandchild_even"],
                        odd_color=configure.colors["table"]["grandchild_odd"],
                        column_count=column_count,
                        alternating=perm_id)

    else:  # stage 1 results
        try:
            create_hit_columns(root_item, df.iloc[0], monomers, [])
        except ValueError:  # if df is empty, like due to filtering
            pass
        # child items for all hits per peak
        stage1_id = 0
        for hit_id, hit in df.reset_index().iterrows():
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
            hit_item.setTotalBackground(
                even_color=configure.colors["table"]["child_even"],
                odd_color=configure.colors["table"]["child_odd"],
                column_count=column_count,
                alternating=hit_id)


def create_hit_columns(item, hit, monomers, sites):
    """
    Create columns that contain information on a stage 1/2 hit
    in the results table.

    :param SortableTreeWidgetItem item: row to fill
    :param pd.Series hit: hit data
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
                            ("Theo_Mass", configure.dec_places()),
                            ("Da", configure.dec_places()),
                            ("ppm", "{:.0f}")]:
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
        for label, form in [("Theo_Mass", configure.dec_places()),
                            ("Da", configure.dec_places()),
                            ("ppm", "{:.0f}")]:
            item.setText(pos, form.format(hit[label]))
            item.setTextAlignment(pos, Qt.AlignRight)
            pos += 1

        # monomer counts
        for monomer in monomers:
            item.setText(pos, "{:.0f}".format(hit[monomer]))
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
    :return: True if all columns satisfy the filtering condition;
             False otherwise
    """

    for lower, upper in query:
        if not (lower <= int(item.text(col)) <= upper):
            return False
        col += 1
    return True


class MainWindow(QMainWindow, Ui_MoFi):
    """
    The main window.

    .. automethod:: __init__
    """

    def __init__(self, parent=None):
        """
        Initialize the main window.

        :param QWidget parent: parent widget
        """

        # initialize the GUI
        super().__init__(parent)
        self.setupUi(self)
        QApplication.instance().installEventFilter(self)

        # initialize private members
        self._exp_mass_data = None  # peak list (mass + relative abundance)
        self._in_context_help_mode = False  # True if in context help mode
        self._known_mods_formula = None  # formula of known modifications
        self._known_mods_mass = 0  # mass of known modification
        self._monomer_hits = None  # results from the monomer search
        self._osa_checked = [False, False]  # check status of btOSA
        self._path = configure.defaults["path"]  # last selected path
        self._polymer_hits = None  # results from the polymer search
        self._protein_formula = mass_tools.Formula()  # f. of the sequence
        self._results_tree_headers = [[], []]  # results tables headers
        self._search_statistics = None  # search statistics
        self._selected_peak = 0  # currently selected peak

        # connect signals to slots
        self.acAbout.triggered.connect(self.show_about)
        self.acLoadSettings.triggered.connect(self.load_settings)
        self.acManual.triggered.connect(self.show_context_help)
        self.acQuit.triggered.connect(QApplication.instance().quit)
        self.acSaveSettings.triggered.connect(self.save_settings)
        self.acWhatIs.triggered.connect(self.enter_context_help_mode)

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
                extensions=["csv", "xls", "xlsx", "xls-bpf"],
                table_widget=self.tbPolymers,
                cols=_polymer_table_columns
            )
        )
        self.btLoadPeaks.clicked.connect(self.load_spectrum)
        self.btLoadSequence.clicked.connect(self.load_sequence)
        self.btOnlyShowUnannotated.clicked.connect(self.only_show_unannotated)
        self.btResetZoom.clicked.connect(self.reset_zoom)
        self.btRestoreSortOrder.clicked.connect(
            self.restore_original_sort_order)
        self.btSaveMonomers.clicked.connect(
            lambda: self.save_table("Export modifications", "monomers"))
        self.btSavePolymers.clicked.connect(
            lambda: self.save_table("Export glycans", "polymers"))
        self.btSaveSequence.clicked.connect(self.save_sequence)
        self.btSaveSpectrum.clicked.connect(self.save_spectrum)
        self.btUncheckAll.clicked.connect(
            lambda: self.check_results_tree(check=False))
        self.btUpdateMass.clicked.connect(self.calculate_protein_mass)

        self.cbDeltaSymbol1.currentIndexChanged.connect(
            self.toggle_delta_series)
        self.cbDeltaSymbol2.currentIndexChanged.connect(
            self.toggle_delta_series)
        self.cbTolerance.activated.connect(self.choose_tolerance_units)

        self.chCombineDelta.clicked.connect(lambda: self.update_selection())
        self.chIndexDelta.clicked.connect(lambda: self.update_selection())
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
        self.sbDisulfides.valueChanged.connect(self.calculate_mod_mass)

        self.taResults.currentChanged.connect(self.set_save_results_menu)

        self.tbMonomers.cellChanged.connect(self.calculate_mod_mass)
        self.tbStatistics.itemClicked.connect(
            lambda item: self.update_selection(clicked_table_item=item))

        self.teSequence.textChanged.connect(
            lambda: self.teSequence.setStyleSheet(
                "QTextEdit {{ background-color: {} }}"
                .format(configure.colors["widgets"]["bg_error"])))

        self.twResults1.itemClicked.connect(
            lambda item: self.update_selection(
                clicked_tree_item=item, clicked_tree=self.twResults1))
        self.twResults2.itemClicked.connect(
            lambda item: self.update_selection(
                clicked_tree_item=item, clicked_tree=self.twResults2))

        # fill the mass sets combobox
        for set_id, mass_set in enumerate(configure.mass_sets):
            self.cbMassSet.addItem(mass_set)
            self.cbMassSet.setItemData(
                set_id,
                configure.mass_sets[mass_set]["tooltip"],
                Qt.ToolTipRole)
        self.cbMassSet.currentTextChanged.connect(self.choose_mass_set)

        # spinbox with upper limit for each modification
        self.sbMaxMods.setValue(configure.defaults["maxmods"])

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

        # delta series symbols
        for cb, start_id in [(self.cbDeltaSymbol1, 3),
                             (self.cbDeltaSymbol2, 0)]:
            cb.blockSignals(True)
            cb.addItem("off")
            cb.setItemData(0, "series disabled", Qt.ToolTipRole)
            for i, symb in _delta_symbols["name"][1:].iteritems():
                resource = ":/mofi resource/images/Symbol-{}.png".format(symb)
                cb.addItem(QIcon(QPixmap(resource)), "")
                cb.setItemData(i, symb, Qt.ToolTipRole)
            cb.setCurrentIndex(start_id)
            cb.blockSignals(False)

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
        self.btDefaultMonomers.setMenu(menu)

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
        self.btDefaultPolymers.setMenu(menu)

        # headers with filters for the results trees
        for tree_widget in (self.twResults1, self.twResults2):
            header = FilterHeader(tree_widget)
            header.setDefaultSectionSize(80)
            header.setHighlightSections(True)
            header.setMinimumSectionSize(0)
            if (QApplication.style().metaObject()
                            .className() == "QWindowsVistaStyle"):
                header.setFixedHeight(50)
                header._vertical_padding = 8
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
            "Save in wide format …",
            lambda: self.save_search_results("stats_wide"))
        self.save_results_menu["table"].addAction(
            "Save in long format …",
            lambda: self.save_search_results("stats_long"))
        self.set_save_results_menu(1)


    # noinspection PyUnusedLocal,PyPep8Naming
    def eventFilter(self, obj, event):
        """
        Filter any mouseclick if in context help mode
        and report the widget that was clicked.

        :param QObject obj: watched object
        :param QEvent event: filtered event
        :return: True if the event was filtered, False otherwise
        :rtype: bool
        """

        if event.type() == QEvent.MouseButtonRelease:
            if self._in_context_help_mode:
                widget = QApplication.instance().widgetAt(QCursor.pos())
                self._in_context_help_mode = False
                QApplication.instance().restoreOverrideCursor()
                self.show_context_help(widget)
                return True
        return False


    def keyPressEvent(self, event):
        """
        Get the widget that currently has focus if the F1 key is pressed.

        :param QKeyEvent event: catched event
        :return: nothing
        """
        if event.key() == Qt.Key_F1:
            widget = QApplication.instance().focusWidget()
            self.show_context_help(widget)



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
        min_spinbox.setStyleSheet(configure.spin_box_flat_style())
        # noinspection PyUnresolvedReferences
        min_spinbox.valueChanged.connect(self.calculate_mod_mass)
        self.tbMonomers.setCellWidget(row_id, 3, min_spinbox)

        max_spinbox = QSpinBox()
        max_spinbox.setMinimum(-1)
        max_spinbox.setSpecialValueText("max")
        max_spinbox.setFrame(False)
        max_spinbox.setValue(max_count)
        max_spinbox.setStyleSheet(configure.spin_box_flat_style())
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

        :param str which: table to convert, ``"monomers"`` or ``"polymers"``
        :return: a dataframe
        :rtype: pd.DataFrame
        :raises ValueError: if an unknown string is passed to ``which``
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
        :param str which: table to export; ``"monomers"`` or ``"polymers"``
        :return: nothing
        """

        filename, _, self._path = get_filename(
            self, "save", dialog_title, self._path, FileTypes(["csv"]))
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

        :param pd.DataFrame df: the dataframe to be converted to a table
        :param QTableWidget table_widget: widget to fill with values
        :param list(tuple) cols: list of (column header, default value,
                     label in import dialog) tuples,
                     sorted according to the order of the arguments to
                     :meth:`~MainWindow._monomer_table_create_row()` and
                     :meth:`~MainWindow._polymer_table_create_row()`.
                     If the default value is ``None``, a column must
                     exist in the input file; otherwise, it will be filled
                     with the default value if missing in the input file.
        :return: nothing
        """

        for label, default, _ in cols:
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


    def load_table(self, default=False, subdir=None, dialog_title=None,
                   extensions=None, table_widget=None, cols=None):
        """
        Import the contents of the monomer/polymer table.

        :param bool default: true: load a default monomer library;
                             false: the user should choose a library file
        :param str subdir: directory containing the default libraries
        :param str dialog_title: title of the file dialog
        :param list(str) extensions: list of file extensions
        :param QTableWidget table_widget: table widget to fill with values
        :param cols: see parameter ``cols``
                     in :meth:`~MainWindow.table_from_df`
        :return: nothing
        """

        if default:
            filename = os.path.join(data_dir, subdir,
                                    self.sender().text() + ".csv")
            filetype = "DEFAULT"
        else:
            filename, filetype, self._path = get_filename(
                self, "open", dialog_title, self._path, FileTypes(extensions))
        if filename is None:
            return

        if filetype == "DEFAULT":
            df = ImportTabDataDialog.get_data(self, filename, cols, "csv",
                                              force_accept=True)
        elif filename.endswith("csv"):
            df = ImportTabDataDialog.get_data(self, filename, cols, "csv")
        elif filetype == "BioPharma Finder PepMap results":
            df = ImportTabDataDialog.get_data(self, filename,
                                              _pepmap_columns, "xls")
            df = io_tools.read_bpf_library(df)
        else:
            df = ImportTabDataDialog.get_data(self, filename, cols, "xls")
        if df is not None:
            self.table_from_df(df, table_widget, cols)
            self.calculate_mod_mass()


    def enter_context_help_mode(self):
        """
        Enter context help mode.

        :return: nothing
        """
        self._in_context_help_mode = True
        QApplication.instance().setOverrideCursor(QCursor(Qt.WhatsThisCursor))


    def show_context_help(self, widget):
        """
        Open the manual in a browser and jump to the appropriate anchor.

        :param QWidget widget: widget about which help was requested
        :return: nothing
        """

        if widget in [self.btLoadSequence, self.btSaveSequence,
                      self.btClearSequence, self.btUpdateMass,
                      self.teSequence, self.sbDisulfides, self.chPngase]:
            site = "workflow"
            anchor = "load-seq"
        elif widget == self.cbMassSet:
            site = "workflow"
            anchor = "mass-sets"
        elif widget in [self.statusbar, self.lbProteinMass,
                        self.lbKnownModMass, self.lbTotalMass]:
            site = "workflow"
            anchor = "status-bar"
        elif widget in [self.btLoadMonomers, self.btSaveMonomers,
                        self.btDefaultMonomers,
                        self.btInsertRowAboveMonomers,
                        self.btInsertRowBelowMonomers,
                        self.btClearMonomers, self.btDeleteRowMonomers,
                        self.tbMonomers]:
            site = "workflow"
            anchor = "mod-list"
        elif widget in [self.btLoadPolymers, self.btSavePolymers,
                        self.btDefaultPolymers,
                        self.btInsertRowAbovePolymers,
                        self.btInsertRowBelowPolymers,
                        self.btClearPolymers, self.btDeleteRowPolymers,
                        self.tbPolymers]:
            site = "workflow"
            anchor = "glycan-library"
        elif widget in [self.btLoadPeaks, self.btSaveSpectrum,
                        self.btResetZoom, self.btLabelPeaks,
                        self.spectrum_canvas]:
            site = "workflow"
            anchor = "spectrum"
        elif widget in [self.btModeSelection, self.btModeDelta,
                        self.cbDeltaSymbol1, self.sbDeltaValue1,
                        self.sbDeltaTolerance1, self.sbDeltaRepetition1,
                        self.cbDeltaSymbol2, self.sbDeltaValue2,
                        self.sbDeltaTolerance2, self.sbDeltaRepetition2,
                        self.chCombineDelta, self.chIndexDelta]:
            site = "workflow"
            anchor = "delta-series"
        elif widget in [self.btFindModifications, self.rbAllPeaks,
                        self.rbSingleMass, self.sbSingleMass,
                        self.sbTolerance, self.cbTolerance]:
            site = "workflow"
            anchor = "perform-search"
        elif widget in [self.twResults1, self.stage1Tab]:
            site = "results"
            anchor = "stage-1-results"
        elif widget in [self.twResults1, self.stage2Tab]:
            site = "results"
            anchor = "stage-2-results"
        elif widget in [self.tbStatistics, self.statisticsTab]:
            site = "results"
            anchor = "statistics"
        elif widget == self.btRestoreSortOrder:
            site = "results"
            anchor = "sort-results"
        elif widget in [self.btClearFilters, self.btOnlyShowUnannotated]:
            site = "results"
            anchor = "filter-results"
        elif widget in [self.btCollapseAll, self.btExpandParents,
                        self.btExpandAll]:
            site = "results"
            anchor = "expand-results"
        elif widget in [self.btSaveResults, self.btCheckAll,
                        self.btUncheckAll]:
            site = "results"
            anchor = "save-results"
        else:
            site = "index"
            anchor = ""

        # create the URL of the manual including the anchor
        site += ".html"
        path = os.path.abspath(os.path.join(docs_dir, "html", site))
        url = "#".join([pathname2url(path), anchor])
        webbrowser.open("file:{}".format(url))


    def show_about(self):
        """
        Show the about dialog.

        :return: nothing
        """

        QMessageBox.about(self, "About MoFi", _version_info)


    def set_save_results_menu(self, index):
        """
        Set the appropriate context menu of the save results button
        and initiate recoloring of the peaks.
        Also enable the 'Only show unannotated' button accordingly
        and set its check state.

        :param int index: current index of the results tab widget
        :return: nothing
        """

        if index == 2:
            self.btSaveResults.setMenu(self.save_results_menu["table"])
            self.btOnlyShowUnannotated.setEnabled(False)
            self.btOnlyShowUnannotated.setChecked(False)
        else:
            self.btSaveResults.setMenu(self.save_results_menu["tree"])
            self.btOnlyShowUnannotated.setEnabled(True)
            self.btOnlyShowUnannotated.setChecked(self._osa_checked[index])
        self.update_selection()


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
            self.sbTolerance.setValue(configure.defaults["da"])
        else:
            self.sbTolerance.setDecimals(0)
            self.sbTolerance.setMinimum(0)
            self.sbTolerance.setMaximum(150)
            self.sbTolerance.setSingleStep(1)
            self.sbTolerance.setValue(configure.defaults["ppm"])


    def load_sequence(self):
        """
        Opens a FASTA file and displays its contents
        in :attr:`~MainWindow.teSequence`.

        Changes :attr:`~MainWindow._path` to directory of selected file
        and text of :attr:`~MainWindow.teSequence` to contents of input file.

        :return: nothing
        """

        filename, _, self._path = get_filename(
            self, "open", "Open sequence", self._path, FileTypes(["fasta"]))
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
        Save the contents of :attr:`~MainWindow.teSequence` to a FASTA file.

        :return: nothing
        """

        filename, _, self._path = get_filename(
            self, "save", "Save sequence", self._path, FileTypes(["fasta"]))
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
        Open a mass list.

        :return: nothing
        """

        filename, _, self._path = get_filename(
            self, "open", "Open mass list", self._path,
            FileTypes(["csv", "xlsx", "xls"]))
        if filename is None:
            return

        cols = [("avg_mass", None, "Average mass"),
                ("rel_abundance", 100, "Relative abundance")]
        if filename.endswith("csv"):
            df = ImportTabDataDialog.get_data(self, filename, cols, "csv")
        else:
            df = ImportTabDataDialog.get_data(self, filename, cols, "xls")

        if df is not None:
            self._exp_mass_data = df
            self.twResults1.clear()
            self.twResults2.clear()
            self._monomer_hits = None
            self._polymer_hits = None
            self.draw_spectrum()


    def save_spectrum(self):
        """
        Save the spectrum as CSV or graphics file.

        :return: nothing
        """

        if self._exp_mass_data is None:
            return

        types = self.spectrum_canvas.get_supported_filetypes()
        types.update({"csv": "Comma-separated value"})
        filename, filetype, self._path = get_filename(
            self, "save", "Save spectrum", self._path, FileTypes(types))
        if filename is None:
            return

        try:
            if filetype == "Comma-separated value":
                self._exp_mass_data.to_csv(filename, index=False)
            else:
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

        filename, _, self._path = get_filename(
            self, "save", "Save results", self._path,
            FileTypes(["csv", "xlsx"]))
        if filename is None:
            return

        # retrieve general information about the search parameters
        out_general = [
            ("Combinatorial search results by MoFi", ),
            ("Date", time.strftime("%c")),
            ("Sequence",
             self.teSequence.toPlainText().replace("\n", r"\n")),
            ("Masses",
             ", ".join(list(self._exp_mass_data["avg_mass"].astype(str)))),
            ("Tolerance",
             "{:.1f} {}".format(self.sbTolerance.value(),
                                self.cbTolerance.currentText()))]

        out_mass_set = [("atom", "atomic weight")]
        for (atom, weight) in configure.current_mass_set.items():
            if atom not in ["description", "tooltip"]:
                out_mass_set.append((atom, str(weight)))

        out_composition = [("name", "formula", "mass",
                            "min_count", "max_count")]
        for (is_checked, name, formula, mass,
             min_count, max_count) in self.calculate_mod_mass():
            if is_checked:
                if max_count == -1:
                    max_count = "max"
                out_composition.append((name,
                                        formula,
                                        configure.dec_places().format(mass),
                                        str(min_count),
                                        str(max_count)))

        out_structures = [("name", "composition", "sites", "abundance")]
        for (is_checked, name, composition,
             sites, abundance) in self.get_polymers():
            if is_checked:
                out_structures.append((name,
                                       composition,
                                       sites,
                                       "{:.2f}".format(abundance)))

        if mode.startswith("stats"):
            out_contents = ["Search statistics"]

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
                                         "Stage 2 permutations",
                                         "Stage 2 hits"],
                             value_name="Value",
                             var_name="Measure")
        else:
            # set up the depth-first search algorithm
            def explore(node):
                yield ([node.checkState(0)]
                       + [node.text(c) for c in range(1, node.columnCount())])
                for i in range(node.childCount()):
                    yield from explore(node.child(i))

            # choose the appropriate data source
            if self.taResults.currentIndex() == 0:
                out_contents = ["Composition search (stage 1)"]
                results_tree = self.twResults1
            else:
                out_contents = ["Structure search (stage 2)"]
                results_tree = self.twResults2

            # create dataframe from (selected) entries
            header = results_tree.headerItem()
            df = (
                pd.DataFrame(
                    explore(results_tree.invisibleRootItem()),
                    columns=(["checked"]
                             + [header.text(c)
                                for c in range(1, header.columnCount())]))
                .replace("", np.nan))
            df[["Exp. Mass", "%"]] = df[["Exp. Mass", "%"]].ffill()
            if not df.empty and df.iloc[0, 1] is None:
                df = df.drop(0)

            if mode == "all":
                out_contents.append("all results")
            elif mode == "checked_parent":
                out_contents.append("checked results (with parents)")
                df = df[df.checked != Qt.Unchecked]
            else:
                out_contents.append("checked results")
                df = df[df.checked == Qt.Checked]
            df = df.drop("checked", axis=1)
        out_general.insert(1, (", ".join(out_contents), ))

        # write to csv or xlsx file
        try:
            if filename.endswith("csv"):
                with open(filename, "w") as f:
                    out_general = [": ".join(i) for i in out_general]
                    f.write("# ")
                    f.write("\n# ".join(out_general))

                    f.write("\n# Mass set: ")
                    f.write(self.cbMassSet.currentText())
                    f.write("\n#   ")
                    out_mass_set = ["; ".join(i) for i in out_mass_set]
                    f.write("\n#   ".join(out_mass_set))

                    f.write("\n# Composition:\n#   ")
                    out_composition = ["; ".join(i) for i in out_composition]
                    f.write("\n#   ".join(out_composition))

                    f.write("\n# Structure:\n#   ")
                    out_structures = ["; ".join(i) for i in out_structures]
                    f.write("\n#   ".join(out_structures))
                    f.write("\n")

                    df.to_csv(f, index=False, na_rep="NA")
            else:
                with pd.ExcelWriter(filename, engine="xlsxwriter") as f:
                    df.to_excel(f, sheet_name="MoFi results", index=False)

                    for sheet, source in [("Parameters", out_general),
                                          ("Mass set", out_mass_set),
                                          ("Composition", out_composition),
                                          ("Structures", out_structures)]:
                        worksheet = f.book.add_worksheet(sheet)
                        for row_id, row in enumerate(source):
                            for col_id, value in enumerate(row):
                                worksheet.write(row_id, col_id, value)
                    f.save()

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
        Changes :attr:`~MainWindow._protein_formula`.

        :return: nothing
        """

        chains, sequence = sequence_tools.read_fasta_string(
            self.teSequence.toPlainText())
        try:
            protein = sequence_tools.Protein(
                sequence,
                chains=chains,
                pngasef=self.chPngase.isChecked())
        except KeyError as e:
            QMessageBox.critical(
                self,
                "Error",
                "Error when parsing sequence: "
                + "{} is not a valid symbol".format(e.args[0]))
            return

        self.sbDisulfides.setEnabled(True)
        self.chPngase.setEnabled(True)
        self.sbDisulfides.setMaximum(protein.amino_acid_composition["C"] / 2)
        self.teSequence.setStyleSheet(
            "QTextEdit {{ background-color: {} }}"
            .format(configure.colors["widgets"]["bg_ok"]))
        self._protein_formula = protein.formula
        self.update_mass_in_statusbar()


    def calculate_mod_mass(self):
        """
        Calculate the mass of known modifications.

        Changes :attr:`~MainWindow._known_mods_mass`
        and, if possible, :attr:`~MainWindow._known_mods_formula`.

        :return: list of (checked, name, composition,
                          mass, min count, max count) tuples
        :rtype: list(tuple(bool, str, str, float, int, int))
        """

        total_mass = 0
        total_formula = mass_tools.Formula()
        only_formulas = True
        result = []

        # add min counts for monomers to the theoretical mass
        for row_id in range(self.tbMonomers.rowCount()):
            # extract the variables
            ch = self.tbMonomers.cellWidget(row_id, 0).findChild(QCheckBox)
            name = self.tbMonomers.item(row_id, 1).text()
            composition = self.tbMonomers.item(row_id, 2).text().strip()
            min_count = self.tbMonomers.cellWidget(row_id, 3).value()
            max_count = self.tbMonomers.cellWidget(row_id, 4).value()
            formula = mass_tools.Formula()
            mass = 0
            is_formula = True

            # get the mass from the formula cell
            error_in_formula = False
            try:  # composition could be a formula ...
                formula = mass_tools.Formula(composition)
                mass = formula.mass
            except ValueError:
                try:  # ... or a mass in Da
                    mass = float(composition)
                    is_formula = False
                except ValueError:
                    error_in_formula = True

            # set the mass tooltip and color the cell
            if error_in_formula or mass == 0:
                bg_color = QColor(configure.colors["widgets"]["bg_error"])
                tooltip = ""
            else:
                bg_color = QColor(configure.colors["widgets"]["bg_ok"])
                tooltip = configure.dec_places().format(mass) + " Da"
            self.tbMonomers.item(row_id, 2).setBackground(bg_color)
            self.tbMonomers.item(row_id, 2).setToolTip(tooltip)

            # color the Min cell if its value exceeds the one of Max
            if min_count > max_count != -1:
                style = configure.spin_box_flat_style(bg="bg_error")
            else:
                style = configure.spin_box_flat_style()
            self.tbMonomers.cellWidget(row_id, 3).setStyleSheet(style)

            if ch.isChecked():
                total_mass += mass * min_count
                total_formula += formula * min_count
                if not is_formula:
                    only_formulas = False
            result.append((ch.isChecked(), name, composition, mass,
                           min_count, max_count))

        # include disulfides
        ss_hydrogens = mass_tools.Formula("H2") * self.sbDisulfides.value()
        total_mass -= ss_hydrogens.mass
        total_formula -= ss_hydrogens

        self._known_mods_mass = total_mass
        if only_formulas:
            self._known_mods_formula = total_formula
        else:
            self._known_mods_formula = mass_tools.Formula()
        self.update_mass_in_statusbar()
        return result


    def update_mass_in_statusbar(self):
        """
        Updates the contents of the mass labels in the status bar.

        :return: nothing
        """

        # (1) protein mass and formula
        if self._protein_formula:
            protein_formula = self._protein_formula
        else:
            protein_formula = "N/A"
        protein_mass = ("<b>Protein:</b> "
                        + configure.dec_places()
                        .format(self._protein_formula.mass)
                        + " Da <font color='#808080'>({})</font>"
                        .format(protein_formula))
        self.lbProteinMass.setText(protein_mass)

        # (2) known modifications mass and (if available) formula
        if self._known_mods_formula:
            known_mods_formula = self._known_mods_formula
            total_formula = self._protein_formula + self._known_mods_formula
        else:
            known_mods_formula = "N/A"
            total_formula = "N/A"
        mod_mass = ("<b>known modifications:</b> "
                    + configure.dec_places()
                    .format(self._known_mods_mass)
                    + " Da <font color='#808080'>({})</font>"
                    .format(known_mods_formula))
        self.lbKnownModMass.setText(mod_mass)

        # (3) total mass and (if available) formula
        total_mass = ("<b>total:</b> "
                      + configure.dec_places()
                      .format(self._protein_formula.mass
                              + self._known_mods_mass)
                      + " Da <font color='#808080'>({})</font>"
                      .format(total_formula))
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
        explained_mass = (self._protein_formula.mass
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
                        int((max_tol_mass - self._protein_formula.mass)
                            / mass),
                        self.sbMaxMods.value())
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
            linewidth=1.5,
            color=configure.colors["spectrum"]["unselected"],
            picker=5)
        self.spectrum_axes.set_ylim(0, 115)
        self.spectrum_axes.set_xlabel("Mass (Da)")
        self.spectrum_axes.set_ylabel("Relative Abundance (%)")
        self.spectrum_axes.yaxis.set_ticks_position("left")
        self.spectrum_axes.xaxis.set_ticks_position("bottom")
        self.spectrum_axes.tick_params(direction="out")
        self.spectrum_axes.autoscale(enable=False)
        self.spectrum_x_limits = self.spectrum_axes.get_xlim()
        self.spectrum_canvas.draw()

        # set up the pick and span selectors.
        self.spectrum_fig.canvas.mpl_connect(
            "pick_event", self.select_peak)
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
        Toggle between single peak and delta series selection.

        :return: nothing
        """

        try:
            if self.bgSpectrum.checkedButton() == self.btModeSelection:
                self.gbDeltaSeries.setEnabled(False)
            else:
                self.gbDeltaSeries.setEnabled(True)
                self.toggle_delta_series()
        except AttributeError:  # i.e., there is currently no spectrum
            pass


    def toggle_delta_series(self):
        """
        Enable the delta series spin boxes
        whose corresponding checkbox is checked.

        :return: nothing
        """

        for cb, sb_delta, sb_tolerance, sb_repetitions in [
                (self.cbDeltaSymbol1, self.sbDeltaValue1,
                 self.sbDeltaTolerance1, self.sbDeltaRepetition1),
                (self.cbDeltaSymbol2, self.sbDeltaValue2,
                 self.sbDeltaTolerance2, self.sbDeltaRepetition2)]:
            if cb.currentIndex() != 0:
                sb_delta.setEnabled(True)
                sb_tolerance.setEnabled(True)
                sb_repetitions.setEnabled(True)
            else:
                sb_delta.setEnabled(False)
                sb_tolerance.setEnabled(False)
                sb_repetitions.setEnabled(False)

        if (self.cbDeltaSymbol1.currentIndex() != 0
                and self.cbDeltaSymbol2.currentIndex() != 0):
            self.chCombineDelta.setEnabled(True)
        else:
            self.chCombineDelta.setEnabled(False)

        self.update_selection()


    def select_peak(self, event):
        """
        Select a peak picked by a mouseclick on the spectrum.

        :param PickEvent event: pick event from the canvas
        :return: nothing
        """
        if event.mouseevent.button == 1:
            self.update_selection(clicked_peak=event.ind[len(event.ind)//2])


    def update_selection(self, clicked_peak=None, clicked_tree_item=None,
                         clicked_tree=None, clicked_table_item=None):
        """
        Update the spectrum and the single mass spinbox
        after the selection has changed.

        :param int clicked_peak: peak in the spectrum that was clicked
        :param SortableTreeWidgetItem clicked_tree_item: item that was clicked
        :param QTreeWidget clicked_tree: results tree whose item was clicked
        :param SortableTableWidgetItem clicked_table_item: item of the
               statistics table that was clicked
        :return: nothing
        """

        if self._exp_mass_data is None:  # there's no spectrum
            return

        # (A) Get the selected peak:
        # (1) a peak was picked in the spectrum
        if clicked_peak is not None:
            self._selected_peak = clicked_peak

        # (2) an item of a results tree was clicked
        #     -> extract the peak index from the annotation index
        scroll_tree1 = True
        scroll_tree2 = True
        if clicked_tree:
            if clicked_tree == self.twResults1:
                scroll_tree1 = False
            else:
                scroll_tree2 = False
            self._selected_peak = int(clicked_tree_item.text(1).split("-")[0])

        # (3) a row of the statistics table was clicked
        #     -> extract the peak index from the leftmost cell
        scroll_table = True
        if clicked_table_item is not None:
            self._selected_peak = int(
                self.tbStatistics.item(clicked_table_item.row(), 0).text())
            scroll_table = False

        # (B) Update the GUI:
        # (1) highlight selected peaks in the spectrum
        if self.bgSpectrum.checkedButton() == self.btModeSelection:
            self.highlight_selected_peak(self._selected_peak)
        else:
            self.highlight_delta_series(self._selected_peak)

        # (2) fill the single mass spin box with the currently selected mass
        try:
            self.sbSingleMass.setValue(
                self._exp_mass_data.loc[self._selected_peak, "avg_mass"])
        except AttributeError:  # occurs when second peak file is loaded
            pass

        # (3) scroll to item in results trees
        for tree, scroll in [(self.twResults1, scroll_tree1),
                             (self.twResults2, scroll_tree2)]:
            if scroll:
                try:
                    for i in range(tree.invisibleRootItem().childCount()):
                        tree.invisibleRootItem().child(i).setSelected(False)
                    for i in range(tree.invisibleRootItem().childCount()):
                        item = tree.invisibleRootItem().child(i)
                        peak_id = int(item.text(1).split("-")[0])
                        if peak_id == self._selected_peak:
                            tree.scrollToItem(item)
                            item.setSelected(True)
                            break
                except AttributeError:
                    pass

        # (4) scroll to row in statistics table
        if scroll_table:
            for row in range(self.tbStatistics.rowCount()):
                peak_id = int(self.tbStatistics.item(row, 0).text())
                if peak_id == self._selected_peak:
                    self.tbStatistics.selectRow(row)
                    break


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
        :return: nothing
        """

        if self._exp_mass_data is None:  # there's no spectrum
            return

        # a dataframe that shares its index with self._exp_mass_data
        # columns 1 and 2 indicate how many delta masses each peak is away
        # from the selected peak in a series
        df_counts = pd.DataFrame(
            index=self._exp_mass_data.index,
            dtype=str)

        # column "1" is straight forward
        if self.cbDeltaSymbol1.currentIndex() != 0:
            df_counts["1"] = self.find_delta_peaks(
                central_peak,
                self.sbDeltaValue1.value(),
                self.sbDeltaTolerance1.value(),
                self.sbDeltaRepetition1.value())
        else:
            df_counts["1"] = ""

        # column "2" has to take into account whether the second delta series
        # should start at each peak of the first one
        if self.cbDeltaSymbol2.currentIndex() != 0:
            if (self.cbDeltaSymbol1.currentIndex() != 0
                    and self.chCombineDelta.isChecked()):
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
        if (self.cbDeltaSymbol1.currentIndex() != 0
                and self.cbDeltaSymbol2.currentIndex() != 0
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
        color_set = np.array([configure.colors["delta"]["other"],
                              configure.colors["delta"]["series_1"],
                              configure.colors["delta"]["series_2"],
                              configure.colors["delta"]["both"],
                              configure.colors["delta"]["main"]])

        lw_set = np.array([1, 2, 2, 2, 3])
        self.spectrum_peak_lines.set_color(color_set[df_counts["color"]])
        self.spectrum_peak_lines.set_linewidth(lw_set[df_counts["color"]])

        # annotate the peaks in the delta series:
        # remove previous annotations and markers
        for obj in self.spectrum_axes.findobj(
            lambda o:
                isinstance(o, matplotlib.text.Annotation) or
                isinstance(o, matplotlib.collections.PathCollection)):
            obj.remove()

        # show labels (if "Index" is checked)
        if self.chIndexDelta.isChecked():
            for peak_id in list(np.flatnonzero(df_counts["color"])):
                label = df_counts["label"][peak_id]
                if self.btLabelPeaks.isChecked():
                    label += " ({})".format(configure.dec_places()).format(
                        self._exp_mass_data["avg_mass"][peak_id])
                self.spectrum_axes.annotate(
                    s=label,
                    xy=(self._exp_mass_data.iloc[peak_id]["avg_mass"],
                        self._exp_mass_data.iloc[peak_id]["rel_abundance"]),
                    xytext=(0, 17),
                    textcoords="offset pixels",
                    horizontalalignment="center",
                    bbox=dict(facecolor="white", alpha=.75,
                              linewidth=0, pad=0))

        # show markers
        for peak_type, cb, y_offset in [(1, self.cbDeltaSymbol1, 4),
                                        (2, self.cbDeltaSymbol2, 4),
                                        (3, self.cbDeltaSymbol1, 10),
                                        (3, self.cbDeltaSymbol2, 4),
                                        (4, self.cbDeltaSymbol1, 4)]:
            mask = df_counts["color"] == peak_type
            marker = _delta_symbols["marker"][cb.currentIndex()]
            if marker in matplotlib.markers.MarkerStyle.filled_markers:
                facecolors = "none"
            else:
                facecolors = color_set[peak_type]
            self.spectrum_axes.scatter(
                x=self._exp_mass_data["avg_mass"][mask],
                y=self._exp_mass_data["rel_abundance"][mask] + y_offset,
                marker=marker,
                facecolors=facecolors,
                edgecolors=color_set[peak_type])
        self.spectrum_canvas.draw()


    def highlight_selected_peak(self, peak_index):
        """
        Highlight the selected peak in the spectrum.

        :param int peak_index: index of the peak to be highlighted
        :return: nothing
        """

        peaks_with_result = np.full(self._exp_mass_data.shape[0], 4, dtype=int)
        if self.taResults.currentIndex() == 0:
            try:
                peaks_with_result[self._monomer_hits
                                      .swaplevel(0, 1)
                                      .loc[-1]
                                      .index.labels[0]] = 2
            except KeyError:
                pass  # results for all peaks found
            except AttributeError:
                peaks_with_result = np.zeros(
                    self._exp_mass_data.shape[0], dtype=int)
        else:
            try:
                peaks_with_result[self._polymer_hits
                                      .swaplevel(0, 1)
                                      .loc[-1]
                                      .index.labels[0]] = 2
            except KeyError:
                pass  # results for all peaks found
            except AttributeError:
                peaks_with_result = np.zeros(
                    self._exp_mass_data.shape[0], dtype=int)

        selected_peaks = np.zeros(self._exp_mass_data.shape[0], dtype=int)
        selected_peaks[peak_index] = 1

        # peak colors will be an array with one entry per peak:
        # before search: 0 - not selected, 1 - selected
        # no annotation: 2 - not selected, 3 - selected
        # annotation:    4 - not selected, 5 - selected
        peak_colors = selected_peaks + peaks_with_result
        color_set = np.array(
            [configure.colors["spectrum"]["unselected"],
             configure.colors["spectrum"]["selected"],
             configure.colors["spectrum"]["unselected_no_hit"],
             configure.colors["spectrum"]["selected_no_hit"],
             configure.colors["spectrum"]["unselected_hit"],
             configure.colors["spectrum"]["selected_hit"]])
        lw_set = np.array([1.5, 1.5, 1.5, 1.5, 1.5, 1.5])
        self.spectrum_peak_lines.set_color(color_set[peak_colors])
        self.spectrum_peak_lines.set_linewidth(lw_set[peak_colors])

        # remove previous annotations and markers
        for obj in self.spectrum_axes.findobj(
            lambda o:
                isinstance(o, matplotlib.text.Annotation) or
                isinstance(o, matplotlib.collections.PathCollection)):
            obj.remove()

        # label selected peaks with masses
        if self.btLabelPeaks.isChecked():
            self.spectrum_axes.annotate(
                s=configure.dec_places().format(
                    self._exp_mass_data.iloc[peak_index]["avg_mass"]),
                xy=(self._exp_mass_data.iloc[peak_index]["avg_mass"],
                    self._exp_mass_data.iloc[peak_index]["rel_abundance"]),
                xytext=(0, 5),
                textcoords="offset pixels",
                horizontalalignment="center",
                bbox=dict(facecolor="white", alpha=.75, linewidth=0, pad=0))
        self.spectrum_canvas.draw()


    def _fill_results_table(self, stage, tree_widget, cols, df_hit):
        """
        Fills a results table with rows.

        :param int stage: search stage (0 = stage 1, 1 = stage 2)
        :param QTreeWidget tree_widget: results table
        :param list(list)) cols: column headers
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
                2, configure.dec_places().format(
                    self._exp_mass_data.loc[mass_index, "avg_mass"]))
            root_item.setTextAlignment(2, Qt.AlignRight)
            root_item.setText(
                3, "{:.1f}".format(self._exp_mass_data
                                   .loc[mass_index, "rel_abundance"]))
            root_item.setTextAlignment(3, Qt.AlignRight)

            # color root item
            if df_hit.loc[mass_index].index.values[0][0] == -1:
                root_item.setTotalBackground(
                    even_color=configure.colors["table"]["parent_no_hit_even"],
                    odd_color=configure.colors["table"]["parent_no_hit_odd"],
                    column_count=column_count,
                    alternating=mass_index)
            else:
                root_item.setTotalBackground(
                    even_color=configure.colors["table"]["parent_hit_even"],
                    odd_color=configure.colors["table"]["parent_hit_odd"],
                    column_count=column_count,
                    alternating=mass_index)

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
        Prepare column labels and then call the actual function that adds rows
        to each results table (:meth:`~MainWindow._fill_results_table()`).

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
                        self._polymer_hits.columns.get_loc("Stage1_ID") + 1:
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


    def only_show_unannotated(self):
        """
        Trigger the results table filter
        to only show items without annotations.

        :return: nothing
        """
        stage = self.taResults.currentIndex()
        if stage != 2:
            self._osa_checked[stage] = self.btOnlyShowUnannotated.isChecked()
            self.filter_results_table(stage)


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
        empty_filter = True
        for _, child in tree_widget.header().allFilters():
            f = re_filter.match(child.text()).groups()
            if "".join(f):
                empty_filter = False
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
            if self.btOnlyShowUnannotated.isChecked():
                # only show items without annotation
                root.child(i).setHidden(root.child(i).childCount() != 0)
            else:  # otherwise use the table filters
                if empty_filter:
                    # simply show all items
                    for j in range(root.child(i).childCount()):
                        root.child(i).child(j).setHidden(False)
                    root.child(i).setHidden(False)
                else:
                    # each item has to be tested
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
        :attr:`~MainWindow._search_statistics`. This method is called
        by :meth:`pd.DataFrame.apply()`.

        :param pd.Series row: row of the statistics dataframe
        :return: nothing
        """

        row_id = self.tbStatistics.rowCount()
        self.tbStatistics.insertRow(row_id)

        # (1) running counter
        item = SortableTableWidgetItem(str(row_id))
        item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
        item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)
        self.tbStatistics.setItem(row_id, 0, item)

        # (2) peak information
        for col_id, label, form in [(1, "avg_mass", configure.dec_places()),
                                    (2, "rel_abundance", "{:.1f}")]:
            item = SortableTableWidgetItem(
                    form.format(self._exp_mass_data.loc[row_id, label]))
            item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)
            self.tbStatistics.setItem(row_id, col_id, item)

        # (3) peak statistics
        for col_id, label in [(3, "search_space_size"),
                              (4, "stage1_results"),
                              (5, "stage2_results"),
                              (6, "stage2_uniques")]:
            try:
                item = SortableTableWidgetItem(str(row[label]))
            except KeyError:
                item = SortableTableWidgetItem("")
            item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)
            self.tbStatistics.setItem(row_id, col_id, item)


    def clear_filters(self):
        """
        Clear the contents of the filter line edits.

        :return: nothing
        """

        if self.taResults.currentIndex() == 0:
            self.twResults1.header().clearFilters()
        elif self.taResults.currentIndex() == 1:
            self.twResults2.header().clearFilters()


    def restore_original_sort_order(self):
        """
        Restore the original order of rows in the results trees.

        :return: nothing
        """

        if self.taResults.currentIndex() == 0:
            # simply sort ascending by ID column
            self.twResults1.sortByColumn(1, Qt.AscendingOrder)
            self.twResults1.header().setSortIndicator(-1, 0)

        elif self.taResults.currentIndex() == 1:
            # sort root items ascending by ID,
            # and child items descending by hit score
            self.twResults2.header().setSortIndicator(-1, 0)
            root = self.twResults2.invisibleRootItem()
            root.sortChildren(1, Qt.AscendingOrder)
            for i in range(root.childCount()):
                root.child(i).sortChildren(5, Qt.DescendingOrder)


    def save_settings(self):
        """
        Export the current settings as XML file.

        :return: nothing
        """

        filename, _, self._path = get_filename(
            self, "save", "Save settings", self._path, FileTypes(["xml"]))
        if filename is None:
            return

        root = ETree.Element("settings")
        settings = [("sequence", self.teSequence.toPlainText()),
                    ("disulfides", self.sbDisulfides.value()),
                    ("pngasef", self.chPngase.isChecked()),
                    ("mass-set", self.cbMassSet.currentText()),
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

        filename, _, self._path = get_filename(
            self, "open", "Load settings", self._path, FileTypes(["xml"]))
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
        self.chPngase.setChecked(root.find("pngasef").text == "True")
        self.cbMassSet.setCurrentText(root.find("mass-set").text)
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
        disulfides = int(root.find("disulfides").text)
        self.sbDisulfides.setMaximum(disulfides)
        self.sbDisulfides.setValue(disulfides)
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
