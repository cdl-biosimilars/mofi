**************************************
Appendix: Configuration and data files
**************************************

=================
config/config.ini
=================

This file contains general parameters for MoFi.

* Section **Defaults**:

  ``da``
    default tolerance value in Dalton

  ``ppm``
    default tolerance value in parts per million

  ``maxmods``
    maximum number of modifications to consider in the composition search

  ``path``
    default path for file dialogs

* Section **Colors**:

  ``unselect_no_annotation``
    Color for peaks that are currently not selected and for which MoFi found no annotations.
    All colors in this section must be specified as HTML hex strings.
    Default: black (``#000000``).

  ``unselect_annotation``
    Color for peaks that are currently not selected and for which MoFi found at least one annotation.
    Default: green (``#00c000``).

  ``select_no_annotation``
    Color for peaks that are currently selected and for which MoFi found no annotations.
    Default: yellow (``#ffcc00``).

  ``select_annotation``
    Color for peaks that are currently selected and for which MoFi found at least one annotation.
    Default: red (``#ff0000``).

  ``delta_other``
    Color for peaks that do not belong to any delta series.
    Default: gray (``#b3b3b3``).

  ``delta_1``
    Color for peaks that belong to the first delta series.
    Default: violet (``#aa0088``).

  ``delta_2``
    Color for peaks that belong to the second delta series.
    Default: dark green (``#2ca05a``).

  ``delta_both``
    Color for peaks that belong to both delta series.
    Default: mixture between violet and green (``#6b5071``).

  ``delta_main``
    Color for the main peak in the delta series.
    Default: red (``#ff0000``).


====================
config/mass_sets.ini
====================

This file contains all mass sets that are available via *Options → Atomic masses*. Each mass set starts by its name in brackets and should at least contain masses for C, H, N, O, P and S. The key ``description`` is optional and will be displayed as tooltip of the menu.

Example:

.. code-block:: ini

   [Average (IUPAC)]
   description: from Table 2/3 in IUPACs Atomic weights of the elements (2013)
   C:  12.011
   H:   1.008
   N:  14.007
   O:  15.999
   P:  30.974
   S:  32.06
   Cl: 35.45


=========================
data/modifications/\*.csv
=========================

Each CSV file present in this directory will be available as a set of default modifications (button *Load default modifications* to the left of the table of modifications). The name of the respective menu item derives from the file name.

Example (``Monosaccharides and frequent modifications.csv``)::

    Checked,Name,Composition,Mass,Min,Max
    False,Hex,C6 H10 O5,162.14084,0,-1
    False,HexNAc,C8 H13 O5 N1,203.19288,0,-1
    False,Neu5Ac,C11 H17 O8 N1,291.25506,0,-1
    False,Neu5Gc,C11 H17 O9 N1,307.25446,0,-1
    False,Fuc,C6 H10 O4,146.14144,0,-1
    False,Pent,C5 H8 O4,132.11482,0,-1
    False,NCore,C34 H56 O25 N2,892.80828,0,-1
    False,OCore,C14 H23 O10 N1,365.33372,0,-1
    False,C_terminal_Lys,C6 H12 N2 O,128.1726,0,2
    False,Oxidation,O,15.9994,0,0
    False,Phosphorylation,P O3 H,79.98,0,0
    False,Succinimide,N-1 H-3,-17.03,0,0


===================
data/glycans/\*.csv
===================

Each CSV file present in this directory will be available as a set of default glycans (button *Load default glycans* to the right of the table of glycans). The name of the respective menu item derives from the file name.

Example (``Default mAB glycans.csv``)::

    Name,Composition,Sites
    Man5,"2 Hex, 1 NCore","ch_A, ch_B"
    G0,"2 HexNAc, 1 NCore","ch_A, ch_B"
    G1,"1 Hex, 2 HexNAc, 1 NCore","ch_A, ch_B"
    G2,"2 Hex, 2 HexNAc, 1 NCore","ch_A, ch_B"
    G0F,"2 HexNAc, 1 Fuc, 1 NCore","ch_A, ch_B"
    G1F,"1 Hex, 2 HexNAc, 1 Fuc, 1 NCore","ch_A, ch_B"
    G2F,"2 Hex, 2 HexNAc, 1 Fuc, 1 NCore","ch_A, ch_B"
    G2FSA1,"2 Hex, 2 HexNAc, 1 Fuc, 1 Neu5Ac, 1 NCore","ch_A, ch_B"
    G2FSA2,"2 Hex, 2 HexNAc, 1 Fuc, 2 Neu5Ac, 1 NCore","ch_A, ch_B"
