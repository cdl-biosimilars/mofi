**************************************
Appendix: Configuration and data files
**************************************

``config/config.ini``
  This file contains general parameters for MoFi:
 
  ``da``
    default tolerance value in Dalton
  
  ``ppm``
    default tolerance value in parts per million
    
  ``maxmods``
    maximum number of modifications to consider in the composition search
 
  ``path``
    default path for file dialogs

``config/mass_sets.ini``
  This file contains all mass sets that are available via *Options -> Atomic masses*.
  Each mass set starts by its name in brackets and should at least contain masses for C, H, N, O, P and S.
  The key description is optional and will be displayed as tooltip of the menu.

``data/modifications/*.csv``
  Each CSV file present in this directory will be available as a set of default modifications (button Load default modifications to the left of the table of modifications).
  The name of the respective menu item derives from the file name.

``data/glycans/*.csv``
  Each CSV file present in this directory will be available as a set of default glycans (button Load default glycans to the right of the table of glycans).
  The name of the respective menu item derives from the file name.
