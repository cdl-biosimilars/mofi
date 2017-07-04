ModFinder - A tool to find molecular modifications (mostly sugars)

Requirements:
    *Python 3.4 (not higher, because of PySide)
    *numpy      |   
    *pandas     |   These packages are included in the anaconda distribution
    *matplotlib |   Anaconda-3-2.3.0 should work. Update packages if necessary
           https://repo.continuum.io/archive/Anaconda3-2.3.0-Windows-x86_64.exe
    *PySide (pip install pyside)
    
    optional:
    *For orphaned features, Thermo MSFileReader bindings are used. Find them
     at: https://github.com/frallain/MSFileReader-Python-bindings. The spectrum
     parser is rudimentary and so far mostly used to extract the last scan of
     a direct infusion experiment (normalized in the mass analyzer).
    *Some other orphaned features use scipy

ModFinder is separated into the following parts:
*The program logic
    -a C++ extension (findmods) that recursively searches through a 
     given list of modifications against a single mass w. given tolerance
    
    -a wrapper for that extension (modification_search.py) that specializes
     certain features for N-glycan search and searches against a list of masses
    
    -helper functions for mass calculations, sequence calculations, glycan
     calculations as well as config tools.
   
   -a graphical interface, realizes using PySide/Qt. PySide at this point 
     only works with Python up to version 3.4. The Interface layout is 
     described in ModFinder_UI.ui and then translated with pyside-uic to a 
     python class (ModFinder_UI.py).The Interface layout was designed in 
     Qt Designer (included in PySide). 
     ModFinder.py then builds the UI and 
     implements its logic by inhereting from the class in ModFinder_UI.

*Building and installing the C++ extension findmods:
    If you're lucky
        python setup.py build
        python setup.py install
    works.
    You'll need the correct version of VisualStudio.
    VisualStudio 2010 worked. There should also be other ways to compile that
    extension.     

*Orphaned features:
    -spectrum_parser.py: Uses A MSFileReader wrapper to read .raw files.
     certain functions were programmed for specific use-cases, such as
     extracting XICs from a LC/MS experiment or getting the last scan from
     a direct infusion experiment.