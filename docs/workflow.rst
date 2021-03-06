.. |bt_clear| image:: ../images/Clear.png
                      :scale: 50 %
.. |bt_clear_table| image:: ../images/ClearTable.png
                            :scale: 50 %
.. |bt_delete_row| image:: ../images/DeleteRow.png
                           :scale: 50 %
.. |bt_insert_row_above| image:: ../images/InsertRowAbove.png
                                 :scale: 50 %
.. |bt_insert_row_below| image:: ../images/InsertRowBelow.png
                                 :scale: 50 %
.. |bt_label_peaks| image:: ../images/Label.png
                            :scale: 50 %
.. |bt_load| image:: ../images/Open.png
                     :scale: 50 %
.. |bt_tools| image:: ../images/Tools.png
                             :scale: 50 %
.. |bt_reset_zoom| image:: ../images/ResetZoom.png
                           :scale: 50 %
.. |bt_save| image:: ../images/Save.png
                     :scale: 50 %
.. |bt_select_delta_series| image:: ../images/DeltaMassMode.png
                                    :scale: 50 %
.. |bt_select_peaks| image:: ../images/SelectMode.png
                             :scale: 50 %
.. |bt_update| image:: ../images/UpdateMass.png
                       :scale: 50 %


********
Workflow
********

The layout of MoFi's main window encourages a workflow that comprises six steps:

1. Load a protein *sequence*.
2. Specify a list of *modifications* for the composition search.
3. Provide a *glycan library* for the structure search (optional).
4. Load or enter *mass values*.
5. Perform the *search(es)*
6. Evaluate and save the results.

.. image:: images/workflow.png
           :alt: Workflow
           :align: center



.. _load-seq:

===========================
(1) Load a protein sequence
===========================

Click |bt_load| *Load sequence*  to load the sequence (in FASTA format) of the protein to be analyzed. Alternatively, paste the sequence into the editor. If the protein has multiple chains, each chain must be specified with its own header line, i.e., a line starting with ``>``. (The mass of one oxygen atom and two hydrogens will be added to the protein mass per chain.) If there is no header line, MoFi assumes that the protein comprises a single chain.

Enter the number of disulfide bonds. For each disulfide bond, the mass of two hydrogen atoms will be removed from the protein.

Specify whether the protein was treated with PNGase F. The asparagine residue in each N-glycosylation site (consensus sequence ``[ST]N[^P][ST]``) will be mutated to aspartate, which corresponds to a mass change of N\ :sub:`–1`\ H\ :sub:`–1`\ O.

.. image:: images/sequence.png
           :alt: Sequence input
           :align: center

.. _mass-sets:

You can select different mass sets in the combobox below the protein sequence. MoFi ships with three predefined sets of atomic masses:

* Average masses as defined by IUPAC in 2013 (Meija, J. *et al.* Atomic weights of the elements 2013 (IUPAC Technical Report). *Pure Appl. Chem.* **88**\ (3), 265–291 (2016))
* Average masses estimated from the isotopic abundance in organic materials (Zhang, Z. *et al.* Mass spectrometry for structural characterization of therapeutic proteins. *Mass Spec. Rev.* **28**, 147–176 (2009))
* Monoisotopic masses from the 2012 atomic mass evaluation (Wang, M. *et al.* The Ame2012 atomic mass evaluation. (II). Tables, graphs and references. *Chin. Phys. C* **36**\ (12), 1603–2014 (2012))

Additional mass sets may be provided in the configuration file ``config/mass_sets.ini``. The atoms and weights in a mass set appear as tooltip.

.. image:: images/mass_sets.png
           :alt: Mass sets
           :align: center

.. _status-bar:

MoFi uses the selected atomic masses for calculating the mass of the protein, the mass of known modifications, and the total mass. These masses and the respective formulas always appear in the status bar.

.. image:: images/statusbar.png
           :alt: Masses in status bar
           :align: center

“N/A” denotes formulas which MoFi could not calculate. This happens, for instance, if you enter the mass of a modification instead of its formula.

Click |bt_save| *Save sequence* to save the contents of the sequence editor, |bt_clear| *Clear sequence* to clear its contents, and |bt_update| *Update masses* to update the masses in the status bar.


.. admonition:: Example
   :class: note
   
   Load the sequence of kadcyla in ``sample data/1_kadcyla.fasta`` by clicking on |bt_load| *Load sequence* next to the sequence editor.

   The antibody comprises two light and two heady chains, which are labeled by FASTA headers (``>Light 1``, ``>Heavy 1``, ``>Light 2`` and ``>Heavy 2``). Enter the correct number of disulfides (16) in the spinbox.

   Using IUPAC average masses, the sequence has a molecular mass of 145199.4360 Da. Including the mass of known modifications (i.e., the disulfide briges), which is –32.2560 Da, the total mass of the protein is 145167.1800 Da.

.. _mod-list:

==============================================================
(2) Specify a list of modifications for the composition search
==============================================================

Load modifications from a CSV or Excel file via |bt_load| *Load modifications*. A modification file must contain at least two columns labeled "Name" and "Composition". The remaining columns are optional; if missing, MoFi assumes the following default values: "Checked" (False), "Min" (0), and "Max" (inf). Choose the columns to use in the :ref:`import dialog <import-dialog>`.

Alternatively, click |bt_tools| *Default modifications and tools → Monosaccharides and frequent modifications* to load a set of default modifications.

Save the current list of modifications to a CSV file via |bt_save| *Save modifications*.

.. image:: images/modification_table.png
           :alt: Table of modifications
           :align: center

The table of modifications contains the following columns:

  :Use?: Check the box for each modification that you want to include in the composition search.
  :Name: Modification names may include any Unicode character.
  :Formula/Mass: accepts either molecular formulas (as shown for Hex) or mass values in Da (as shown for DM1 and MCC).
    
    A molecular formula consists of space-separated ``symbol[count]`` pairs. ``symbol`` is any one- or two-letter atomic symbol whose mass is specified in the current mass set. The optional ``[count]`` is a positive or negative integer. A symbol without count is counted once.

    If you enter a formula and move the mouse cursor over the cell, a tooltip containing the mass of this formula appears (if the syntax of the formula is correct). If there is an error in the formula, the tooltip displays an error message.

  :Min: the minimum …
  :Max: … and maximum number of occurrences, respectively. If the maximum count for a modification is *max*, MoFi calculates it from the glycan library, the mass of the molecule, or the value of *Upper limit for each modification* (see :ref:`below <perform-search>`).

Manipulate the table via the buttons next to it:

  * |bt_insert_row_above| *Insert row above*
  * |bt_insert_row_below| *Insert row below*
  * |bt_delete_row| *Delete row*
  * |bt_clear_table| *Clear table*


.. admonition:: Example
   :class: note
   
   Load modifications from ``sample data/2_modifications.csv`` by clicking on |bt_load| *Load modifications* next to the table of modifications. In the case of kadcyla, the combinatorial search requires the following modifications:

   :Hex, HexNAc, Neu5Ac, Fuc: Monosaccharides that form the N-glycans of the antibody moiety.
   :DM1-MCC: The drug emtansine, coupled to the antibody via a linker. We expect that kadcyla contains up to ten drug molecules.
   :MCC: The linker maleimidylmethyl cyclohexane-1-carboxylate. We expect that kadcyla also contains 'dead' linkers without any attached drug molecule.


.. _glycan-library:

================================================================
(3) Provide a glycan library for the structure search (optional)
================================================================

Load a glycan library from a CSV or Excel file via |bt_load| *Load glycans* and the subsequent :ref:`import dialog <import-dialog>`. A glycan library must contain at least one column labeled "Name". The remaining columns are optional; if missing, MoFi assumes the following default values: "Checked" (True), "Composition" (derive from name; see below), "Sites" (derive from name), and "Abundance" (0.0).

While MoFi supports arbitrary glycan names, it is able to extract information on the monosaccharide composition from abbreviations conforming to the Zhang nomenclature (see Zhang, Z. Large-scale identification and quantification of covalent modifications in therapeutic proteins. *Anal. Chem.* **81**\ (20), 8354–8364 (2009)). For example, the abbreviation "A2G0F" is converted to the composition "3 Hex, 4 HexNAc, 1 Fuc". Moreover, the program recognizes modifications like "N300+A2G0F", which appear in peptide mapping results from *Thermo Fisher BioPharma Finder*. From such an abbrevation, MoFi additionally extracts the glycosylation site (here, "N300").

Alternatively, click |bt_tools| *Default glycans and tools → Default mAb glycans* to load a default glycan library.

Save the current list of glycans to a CSV file via |bt_save| *Save glycans*.

.. image:: images/glycan_table.png
           :alt: Table of glycans
           :align: center

The table of glycans contains the following columns:

  :Use?: Check the box for each glycan that you want to include in the structure search.
  :Name: contains the name of the glycan. Press :kbd:`Shift+Return` after entering a glycan abbreviation following the Zhang or BPF nomenclature to automatically extract the monosaccharide composition and site, if applicable.
  :Composition: accepts a comma-separated list of modifications, all of which must appear in the table of modifications.

    If you enter a composition and move the mouse cursor over the cell, a tooltip containing its mass appears (if the syntax of composition is correct). If there is an error in the composition, the tooltip displays an error message.
  :Sites: accepts a comma-separated list of glycosylation sites.
  :Abundance: may contain relative abundances as determined, e.g., by peptide mapping. MoFi calculates the score of a glycan combination from these values.

Manipulate the table via the buttons next to it:

  * |bt_insert_row_above| *Insert row above*
  * |bt_insert_row_below| *Insert row below*
  * |bt_delete_row| *Delete row*
  * |bt_clear_table| *Clear table*


.. admonition:: Example
   :class: note
   
   Load the glycan library from ``sample data/3_glycan_library.csv`` by clicking on |bt_load| *Load glycans* next to the table of glycans. Note that MoFi also accepts unglycosylated sites (here, the structure 'no_glycan'). We arbitrarily named the glycosylation sites 'ch_A' and 'ch_B', but any other name will also work.

   Alternatively, load the glycan library in ``sample data/3_glycan_library_BPF.xls``. This file contains the results of a peptide mapping analysis in Thermo BioPharma Finder and was directly exported from this program. MoFi automatically extracts the name of the glycosylation site (here, 'N300') and the glycan composition from the column 'Modification' in the XLS file. (For instance, the abbreviation 'A2S1G1F' denotes a glycan comprising 5 Hex, 4 HexNAc, 1 Neu5Ac and 1 Fuc.)

   NB: Since each heavy chain harbors a glycosylation site at N300, you have to change the values in column 'Site' of the table of glycans to 'ch_A, ch_B' or similar.


.. _spectrum:

=============================
(4) Load or enter mass values
=============================

Click |bt_load| *Load mass list* to load a peak list (in CSV or Excel format) that represents a mass spectrum. The file must contain at least one column labeled "Average Mass" or "Average Mass (mean)". If a column labeled "Relative Abundance" is present, MoFi will interpret its values as peak heights. Again, choose columns to be used in the :ref:`import dialog <import-dialog>`.

Click |bt_save| *Save spectrum* to save the spectrum in CSV format or in one of several image file formats (e.g., jpg, pdf, png, …).

.. image:: images/spectrum.png
           :alt: Spectrum
           :align: center

|bt_label_peaks| *Label peaks* turns labels (peak masses) on or off.

If |bt_select_peaks| *Select peaks* is active, you may interact with the spectrum by

* Clicking onto a single peak with the left mouse button, which highlights that peak.
* Dragging a line or rectangle with the right mouse button, which zooms into the selected region of the spectrum. |bt_reset_zoom| *Reset zoom* shows the entire spectrum.

.. image:: images/selection.png
           :alt: Interaction with the spectrum
           :align: center

.. _delta-series:

|bt_select_delta_series| *Select delta series* enters delta series selection mode: Select a single peak to mark it as the main peak (highlighted in red). All peaks that are separated from the main peak by equal distances are highlighted in blue.

.. image:: images/delta_series.png
           :alt: Delta series
           :align: center
 
You can display a second delta series by selecting a marker symbol to the right of the spectrum. The peaks in the second series are highlighted in yellow.

.. image:: images/delta_series_2.png
           :alt: Second delta series
           :align: center

For each series, you may set the following parameters:

* Mass differences between neighboring peaks
* Tolerance of the mass differences
* Maximum repetitions (i.e., the maximum number of labeled peaks on each side of the main peak)

.. image:: images/delta_series_parameters.png
           :alt: Delta series parameters
           :align: center

Select *Index* to display the delta series index above each peak. The main peak is numbered 0, the other peaks in the series are consecutively numbered 1, 2, … (increasing masses) and –1, –2, … (decreasing masses). The selection status of |bt_label_peaks| *Label peaks* determines whether MoFi also displays masses next to the indices (compare the left half of the figure below to its right half).

.. image:: images/delta_series_labels.png
           :alt: Delta series labels
           :align: center

It is also possible to combine the delta series (check button *Combine*). In this case, the second delta series will start at each peak in the first delta series.

.. image:: images/delta_series_combined.png
           :alt: Combining delta series
           :align: center

The following table summarizes the color scheme for delta series:

.. image:: images/colortable_delta.png
           :alt: Delta series color scheme
           :align: center


.. admonition:: Example
   :class: note
   
   Load the mass spectrum of kadcyla from ``sample data/4_spectrum.csv`` or ``sample date/4_spectrum.xls`` by clicking on |bt_load| *Load mass list* next to the delta series parameters.

   Apparently, the spectrum contains group of peaks whose largest peaks are separated by equal masses. Highlight those peaks by clicking |bt_select_delta_series| *Select delta series* and choosing the following parameters for series 1: Mass difference, 957.5 Da (i.e., one DM1-MCC molecule); tolerance: 5.0 Da; maximum repetitions: auto.

   Within each group, the major peaks also differ by equal masses. Highlight those peaks by activating the second delta series, entering a mass difference of 162.1 Da (i.e., one hexose) and two maximum repetitions, and finally checking *Combine*.


.. _perform-search:

==========================
(5) Perform the search(es)
==========================

.. image:: images/search_parameters.png
           :alt: Search parameters
           :align: center

Click onto *Find modifications* to start the composition search, possibly followed by the structure search if you specified a list of glycans in step 3. You may

* analyze either all peaks in the spectrum or a single mass.
* set the tolerance for acceptable annotations in Da or ppm.
* specify an upper limit for each modification to be used in the absence of a glycan library.

.. admonition:: Example
   :class: note
   
   Search for modifications in kadcyla by clicking *Find modifications*.


.. _import-dialog:

===========================
The "Import CSV/XLS" dialog
===========================

.. image:: images/import_dialog.png
           :alt: Import dialog
           :align: center

Since the CSV file format is not standardized, MoFi allows you to set CSV parameters upon importing such a file. Some of these parameters also apply to Excel files, which is why a similar dialog appears if you import data from a spreadsheet.

Importing tabular data comprises four steps:


------------------------------------
1. View the contents of the raw file
------------------------------------

This step is only possible for CSV files. It may help you to choose correct parameters in the second step.


---------------------------
2. Choose import parameters
---------------------------

  :Delimiter: Single character (e.g., ``,``) or regular expression (e.g., ``\t``), separates each field in a line (CSV files only).
  :Comment char: Single character, denotes lines that should be ignored (CSV files only).
  :Quote char: Single character, encloses fields that contain the delimiter as normal character (CSV files only).
  :Header: Indicates whether the first line should be interpreted as column header (available for CSV and XLS files).
  :Decimal point: Single character, separates the integer part from the fractional part of a decimal number (CSV files only).
  :Thousands separator: Single character, used in digit grouping (CSV files only).
  :Skip rows: Number of rows to skip at the top of the file (available for CSV and XLS files).
  :Sheet name: The sheet to be imported (XLS files only).

Click *Apply* to apply the current settings for the import parameters.


--------------------------
3. Preview the parsed file
--------------------------

This table allows you to check whether the import parameters are correct.


-------------------------------
4. Select the columns to import
-------------------------------

Depending on the file contents, there will be a different set of mandatory and optional columns to be filled with values. These columns are described in the corresponding sections of the manual (see :ref:`list of modifications <mod-list>`, :ref:`glycan library <glycan-library>`, and :ref:`spectrum <spectrum>`).

The selection in the image above prompts MoFi to

* use default values for columns *Use?* and *Min*
* derive values for column *Formula/Mass* (*Name*, *Max*) in the table of modifications from column *Composition* (*Name*, *Max*) in the CSV file.

Click *OK* once you have chosen the correct columns.


.. _mutation-dialog:

==================================
The "Create point mutation" dialog
==================================

Click |bt_tools| *Default modifications and tools → Create point mutation …* to open the *Create point mutation* dialog.

.. image:: images/mutation_dialog.png
           :alt: Create point mutation dialog
           :align: center

Use this dialog to quickly create a modification corresponding to the mass change associated with a single amino acid exchange. Enter the original residue and the new one, using their single-letter abbreviations. MoFi will immediately display the formula describing this mutation and its mass, which is calculated from the current mass set in the main window.

Click *OK* to generate a modification. For the Ser→Val point mutation, this modification will be

.. image:: images/mutation_modification.png
           :alt: Modification corresponding to the Ser→Val point mutation
           :align: center


.. _truncation-dialog:

=======================================
The "Create terminal truncation" dialog
=======================================

Click |bt_tools| *Default modifications and tools → Create terminal truncation …* or |bt_tools| *Default glycans and tools → Create terminal truncation …* to open the *Create terminal truncation* dialog.

.. image:: images/truncation_dialog.png
           :alt: Create terminal truncation
           :align: center

Use this dialog to quickly create terminal truncations to be used as modifications for the composition search or structures for the second search stage:

1. Enter the sequence of the terminus. Alternatively, import the sequence of a single chain or all chains from the main window (button *From parameters*).
2. Choose whether the sequence should be truncated from the *N*- or *C*-terminus.
3. Select cleavage positions. The sequence may either be cleaved after every *n*-th residue or after residues given as a list (e.g., ``1, 3-5, 8, 10-13``).

Click *Preview* to display a list of the cleavage products that will be generated, including their sequence, formula and mass.

If you opened this dialog next to the table of modifications, MoFi will create one modification for each possible truncation. Cleaving the sequence SLSPG after every residue, starting at the *N*-terminus, yields

.. image:: images/truncation_modifications.png
           :alt: Truncation modifications
           :align: center

By contrast, if you opened this dialog next to the table of glycans, MoFi will create one structure for each possible truncation and also generate the appropriate monomers in the table of modifications:

.. image:: images/truncation_structures.png
           :alt: Truncation structures
           :align: center



========
Settings
========

.. image:: images/menu_file.png
           :alt: File menu
           :align: center

* *Save settings …* (Ctrl+S) saves the current settings (i.e., all parameters which you specified in steps (1) to (5) above) as an XML file.
* *Load settings …* (Ctrl+O) loads settings from a previously generated XML file.
* *Quit* (Ctrl+Q) closes MoFi.
