.. |bt_check_all| image:: ../images/CheckAll.png
                          :scale: 50 %
.. |bt_clear_filters| image:: ../images/ClearFilters.png
                              :scale: 50 %
.. |bt_collapse_all| image:: ../images/CollapseAll.png
                             :scale: 50 %
.. |bt_expand_parents| image:: ../images/Expand1.png
                               :scale: 50 %
.. |bt_expand_all| image:: ../images/Expand2.png
                           :scale: 50 %
.. |bt_save| image:: ../images/Save.png
                     :scale: 50 %
.. |bt_uncheck_all| image:: ../images/UncheckAll.png
                            :scale: 50 %


*******
Results
*******

========
Spectrum
========

After MoFi has finished the search for modifications, the spectrum gives an overview of the search results:

.. image:: images/results_spectrum.png
           :alt: Spectrum after combinatorial search
           :align: center


The color of each peak indicates whether this peak is selected (lightness) and if MoFi was able to annotate it (hue). Hue derives from the search stage whose results are currently visible. Hence, click on the “Structure” (“Composition”) tab to quickly identify peaks for which MoFi found at least one annotation during structure (composition) search.

.. image:: images/colortable_spectrum.png
           :alt: Spectrum color scheme
           :align: center


==============
Results tables
==============

The hierarchical tables at the bottom of the main window display results of both search stages and overall search statistics. Switch between the tables by clicking the tabs to their right. The results tables share the following color scheme:

.. image:: images/colortable_results.png
           :alt: Results tables color scheme
           :align: center



.. _stage-1-results:

---------------
Stage 1 results
---------------

.. image:: images/results_table_1.png
           :alt: Results table for search stage 1
           :align: center

For each peak, the parent row (dark green) shows its experimental mass (i.e., the mass detected in the spectrum), relative abundance (%), and data of the annotation with the least (absolute) unexplained mass. The format of the row index is ``[peak ID]``.

A child row (light green) with the following columns appears for each possible annotation:

* theoretical mass (i.e., the mass calculated for the current annotation)
* residual unexplained mass in *Da* and *ppm*
* counts for each modification (here: columns *Hex* to *MCC*)

For child rows, the format of the row index is ``[peak ID]-[hit ID]``.



.. _stage-2-results:

---------------
Stage 2 results
---------------

.. image:: images/results_table_2.png
           :alt: Results table for search stage 2
           :align: center

For each peak, the parent row (dark green) shows its experimental mass, relative abundance (%), and data of the annotation with the highest hit score. The format of the row index is ``[peak ID]``.

A child row (light green) with the following columns appears for each possible annotation:

* hit index
* hit score (in percent)
* number of permutations
* theoretical mass
* residual unexplained mass in *Da* and *ppm*
* counts for each modification (here: columns *Hex* to *MCC*)
* data of the permutation with the highest permutation score

For child rows, the format of the row index is ``[peak ID]-[hit ID]``.

Each stage 2 annotation ('hit') may comprise several permutations of glycan assignments. *Permutations* are isobaric annotations that comprise an equal set of glycans, but assign those glycans to different glycosylation sites. For instance, the annotations "G0F at site A, G1F at site B" and "G1F at site A, G0F at site B" are permutations of the general glycan annotation (G0F, G1F).

A grandchild row (white) with the following columns appears for each possible permutation:

* permutation index
* permutation score (in percent)
* one column per glycosylation site (here: *ch_A* and *ch_B*)

For grandchild rows, the format of the row index is ``[peak ID]-[hit ID]-[permutation ID]``.

The *permutation score* is proportional to the normalized probability of a glycan combination (normalization is done peakwise). Assume that MoFi found :math:`K` permutations for a peak, and that the :math:`k`-th permutation (:math:`k = 1, \dots, K`) comprises a set of glycans :math:`g_k \subset \{1, \dots, n\}` with known abundances :math:`p_1, \dots, p_n` (:math:`n`: total number of glycans). Then, the permutation score :math:`S^\mathrm{perm}_k` of the :math:`k`-th permutation is

.. math::

   S^\mathrm{perm}_k = \frac{\prod_{i \in g_k} p_i}{\sum_{k=1}^K \prod_{i \in g_k} p_i } \ .

The *hit score* is the sum of the permutation scores of all possible glycan assignments that constitute an annotation. Assume that MoFi found :math:`H` hits for a peak, and that the :math:`h`-th hit (:math:`h = 1, \dots, H`) comprises a set of permutations :math:`v_h \subset \{1, \dots, K\}`. Then, the hit score :math:`S^\mathrm{hit}_h` of the :math:`h`-th hit is

.. math::

   S^\mathrm{hit}_h = \sum_{i \in v_h} S^\mathrm{perm}_i \ .

Note that both the permutation score and hit score correspond to the contribution of each permutation or hit to the peak height. From the definition of the permutation score, it is evident that

.. math::

   \sum_{k=1}^K S^\mathrm{perm}_k = 1 \ .

Likewise, since each set of permutations :math:`v_h` contains a different subset of all permutations associated with a peak (:math:`\forall h, \eta \in 1, \dots, H: v_h \cap v_\eta = \emptyset`), it also follows that

.. math::

   \sum_{h=1}^H S^\mathrm{hit}_h = 1 \ .


.. _statistics:

----------
Statistics
----------

.. image:: images/statistics_table.png
           :alt: Search statistics
           :align: center

The statistics table lists the following measures for each peak:

* search space size
* number of annotations found in the composition search
* number of permutations found in the structure search
* number of hits found in the structure search




==================================
Manipulation of the results tables
==================================

--------
Sorting
--------

Click on a column header to sort the table by this column. Sorting by ID restores the original order of the rows.

.. _filter-results:

--------
Fitering
--------

You may filter the stage 1/2 results tables by entering a *constraint* for one or several modifications into the filters beneath the table header and then pressing *Enter*. The button |bt_clear_filters| *Clear filters* removes all constraints.

Constraints must have one of the following forms:

* ``N`` selects rows with exactly ``N`` modifications.
* ``L-`` selects rows with at least ``L`` modifications.
* ``L-U`` selects rows with at least ``L`` and at most ``U`` modifications.
* ``-U`` selects rows with at most ``U`` modifications.

.. image:: images/filter.png
           :alt: Results filter
           :align: center

.. _expand-results:

---------
Expanding
---------

Click on the triangle to the left of any row with children to expand this row. Alternatively, use one of the following buttons to the right of the table:

* |bt_collapse_all| *Collapse all* collapses all rows.
* |bt_expand_parents| *Expand parent items* expands all peak rows, leaving the hit rows collapsed.
* |bt_expand_all| *Expand all* expands all peak and hit rows.



.. _save-results:

==============
Saving results
==============

Click |bt_save| *Save results* to save the currently visible table in CSV or XLSX format.

The stage 1/2 results tables support the following saving options:

* *Save all entries* saves all rows in the table.
* *Save checked entries* saves rows whose checkbox is fully checked.
* *Save checked entries with parents* saves rows whose checkbox is at least partially checked.

Click |bt_check_all| *Check all* or |bt_uncheck_all| *Uncheck all* to check or uncheck all rows in the currently visible table.

.. admonition:: Example
   :class: note

   Check the box of permutation row 78-0-1. Its parent rows (78-0 and 78) automatically become partially checked:

   .. image:: images/results_table_checkboxes.png
              :alt: Results table checkboxes
              :align: center

   Hence, *Save checked entries* yields a CSV/Excel file that contains one row (78-0-1), while *Save checked entries with parents* yields a file that contains three rows (78, 78-0 and 78-0-1).

The statistics table supports the following saving options:

* *Save in wide format* saves the statistics table as shown in the main window.
* *Save in long format* saves the statistics table in long (tidy) format. In the exported table, the columns *Search space size* to *Stage 2 hits* will be gathered, which yields two columns *Measure* and *Value*.

.. admonition:: Example
   :class: note

   Saving the statistics table in long format facilitates its analysis by tools that require tidy data. For instance, the R script below employs packages from the tidyverse to plot the search statistics.

   .. literalinclude:: sample data/search_statistics.R
                       :language: R

   .. image:: images/search_statistics.png
              :alt: Search statistics
              :align: center
   
