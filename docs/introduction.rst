************
Introduction
************

MoFi allows you to find post-translational modifications (PTMs) in deconvoluted mass spectra of intact proteins. The program annotates such spectra by performing a two-stage search.

The first search stage, termed the *composition search*, assigns a PTM composition (including monosaccharides and other modifications) to each mass in the spectrum. An exemplary result is

  The residual mass of 3054.58 Da can be explained by the mass of seven hexoses, eight *N*-acetyl hexosamines, and two fucoses, leaving an unexplained mass of 1.78 Da.

The second search stage, termed the *structure search*, reduces the potentially large number of stage 1 annotations by integrating bottom-up data represented by a glycan library. An exemplary result is

  The combination above can be explained by a glycoform carrying the glycans A2G0F on chain A and A2G1F on chain B.

.. image:: images/main_window.png


============
Getting help
============

.. image:: images/menu_help.png
           :alt: Help menu
           :align: center

MoFi provides context-sensitive help: You may open the relevant section of this manual by

* pressing F1 while you enter values via the keyboard
* entering context help mode (*Help → What is …* or Shift+F1) and then clicking on any widget
* Clicking *Help* in any dialog

.. admonition:: Example
   :class: note
   
   This help file also includes a step-by-step tutorial denoted by *Example* boxes like this ones.
   In this tutorial, you will annotate a mass spectrum of kadcyla (ado-trastuzumab emtansine),
   an antibody-drug conjugate used in treatment of HER2-positive metastatic breast cancer.
   The dataset is taken from: Bailey, A. *et al.* Complete characterization of a lysine-linked antibody
   drug conjugate by native LC/MS intact mass analysis and peptide mapping.
   *Thermo Fisher Application Note* 72511 (2017).
   The folder ``sample data`` contains the required data files:
   
   :``0_settings.xml``: MoFi settings     
   :``1_kadcyla.fasta``: sequence of trastuzumab in FASTA format
   :``2_modifications.csv``: PTMs found in the drug
   :``3_glycan_library.csv``: N-glycans determined by peptide mapping
   :``3_glycan_library_BPF.csv``: the same N-glycan library as exported from Thermo BioPharma Finder
   :``4_spectrum.xls``: mass spectrum as exported from Thermo BioPharma Finder
   :``4_spectrum.csv``: the same mass spectrum in comma-separated value format
