# ModFinder

A tool to find molecular modifications (mostly sugars)!

## Requirements

* Python 3
* numpy
* pandas
* matplotlib
* PyQt5
* xlrd

For building the C++ extension:

* Unix: gcc / clang
* Windows: Visual C++ 2015 Build Tools with Windows SDK 8.1

## Installing

    python3 setup.py install

## Running in source

Build the extension using `python3 setup.py build` and copy the library from
`build/lib.linux-x86_64-3.5/findmods.cpython-35m-x86_64-linux-gnu.so` (adjust
the path as necessary) to `mofi/`. Then run `run.py`. 

## Components

ModFinder is separated into the following parts:

### Program logic

* A C++ extension (findmods) that recursively searches through a 
  given list of modifications against a single mass w. given tolerance
    
* A Wrapper for that extension (modification_search.py) that specializes
  certain features for N-glycan search and searches against a list of masses

* Helper functions for mass calculations, sequence calculations, glycan
  calculations as well as config tools.

* A graphical interface, implemented using PyQt5. The Interface layout is 
  described in ModFinder_UI.ui and then translated with pyuic5 to a 
  python class (modfinder_ui.py). The Interface layout was designed in 
  Qt Designer. ModFinder.py loads the UI and implements its logic by inhereting
  from the class in ModFinder_UI.
