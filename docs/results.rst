*******
Results
*******

.. image:: images/results_spectrum.png
           :alt: Spectrum after combinatorial search
           :align: center

The peaks in the spectrum are colored according to the search results and selection status:

* No combination of modifications was found: black – peak not selected, yellow – peak selected
* At least one combination of modifications was found: green – peak not selected, red – peak selected

Search results for selected peaks are shown in the table at the bottom.

.. image:: images/results_table.png
           :alt: Results table
           :align: center

For each peak, the parent row shows its experimental mass and relative abundance. A child row appears for each possible composition of modifications, containing

* the counts of each modification (here: columns Hex to MCC)
* the experimental mass, theoretical mass (i.e., the mass calculated for the current combination) and residual unexplained mass (in Da and ppm)
* the glycan combination (one column per site; here: A and B), as well as the overall abundance calculated for this combination from the relative abundances of the glycans

If you select a single peak, you can also display the results of the first search stage by unchecking *Filter structure hits*. The spectrum and the table of results can be saved in PNG/SVG and CSV format, respectively (buttons *Save spectrum …* and *Save results …*).

.. image:: images/save_buttons.png
           :alt: Save buttons
           :align: center

You may filter the results by entering a *constraint* for one or several modifications above the table and then pressing *Enter* or clicking *Apply*. The button *Clear* removes all constraints.

Constraints may have one of the following forms:

* ``N`` selects rows with exactly ``N`` modifications.
* ``L-`` selects rows with at least ``L`` modifications.
* ``L-U`` selects rows with at least ``L`` and at most ``U`` modifications.
* ``-U`` selects rows with at most ``U`` modifications.

.. image:: images/filter.png
           :alt: Results filter
           :align: center
