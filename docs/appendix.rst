**************************************
Appendix: Configuration and data files
**************************************

=================
config/config.ini
=================

This file contains default parameters and colors for MoFi. Default values are given in parentheses. All colors must be specified as HTML hex strings.


.. rubric:: defaults

General parameters.

  ``da``
    default tolerance value in Dalton (5.0)

  ``ppm``
    default tolerance value in parts per million (40)

  ``maxmods``
    maximum number of modifications to consider in the composition search (100)

  ``path``
    default path for file dialogs (``docs/sample data``)


.. rubric:: colors.widgets

Widget colors.

  ``bg_ok``
    background color for widgets whose contents are correct (``#ffffff``)

  ``bg_error``
    background color for widgets whose contents are incorrect (``#fde0ef``)


.. rubric:: colors.spectrum

The default color scheme for the peaks corresponds to ColorBrewer's 4-class *PiYG*, which is color blind friendly. See http://colorbrewer2.org/ for details.

  ``unselected``
    unselected peaks before any search (``#000000``)

  ``unselected_no_hit``
    unselected peaks for which MoFi found no annotations (``#d01c8b``)

  ``unselected_hit``
    unselected peaks for which MoFi found at least one annotation (``#4dac26``)

  ``selected``
    selected peaks before any search (``#b3b3b3``)

  ``selected_no_hit``
    selected peaks for which MoFi found no annotations (``#f1b6da``)

  ``selected_hit``
    selected peaks for which MoFi found at least one annotation (``#b8e186``)


.. rubric:: colors.delta

The default color scheme for delta series is 4-class *RdYlBu*.

  ``other``
    peaks that do not belong to any delta series (``#b3b3b3``)

  ``series_1``
    peaks that belong to the first delta series (``#2c7bb6``)

  ``series_2``
    peaks that belong to the second delta series (``#fdae61``)

  ``both``
    peaks that belong to both delta series (``#abd9e9``)

  ``main``
    main peak in the delta series (``#d7191c``)


.. rubric:: colors.table

The default color scheme for the background of the result tables derives from 9-class *PiYG*. Alternating colors for rows with even and odd line numbers are possible.

  ``parent_no_hit_even`` and  ``parent_no_hit_odd``
    parent rows for which MoFi found no annotations (``#f1b6da``)

  ``parent_hit_even`` and ``parent_hit_odd``
    parent rows for which MoFi found at least one annotation (``#b8e186``)

  ``child_even`` and ``child_odd``
    hit rows (``#e6f5d0``)

  ``grandchild_even`` and ``grandchild_odd``
    permutation rows (``#f7f7f7``)



====================
config/mass_sets.ini
====================

This file contains all mass sets that are available in the mass sets combobox. Each mass set starts by its name in brackets and should at least contain masses for C, H, N, O, P and S. The key ``description`` is optional and will be displayed as tooltip of the combobox.

.. literalinclude:: ../mofi/config/mass_sets.ini
                   :language: ini


=========================
data/modifications/\*.csv
=========================

Each CSV file present in this directory will be available as a set of default modifications. The name of the respective menu item derives from the file name.

Example (``Monosaccharides and frequent modifications.csv``):

.. literalinclude:: ../mofi/data/modifications/Monosaccharides and frequent modifications.csv


===================
data/glycans/\*.csv
===================

Each CSV file present in this directory will be available as a set of default glycans. The name of the respective menu item derives from the file name.

Example (``Default mAB glycans.csv``):

.. literalinclude:: ../mofi/data/glycans/Default mAB glycans.csv
