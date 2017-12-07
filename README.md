# MoFi - A software tool for annotating glycoprotein mass spectra by integrating hybrid data from the intact protein and glycopeptide level

## Requirements

General:

* Python 3.5
* numpy
* pandas
* matplotlib
* PyQt5
* xlrd
* xlsxwriter
* Sphinx (for creating the documentation)
* sphinx_rtd_theme (for creating the documentation)

For building the C++ extension:

* Linux: gcc / clang
* Windows: Visual C++ 2015 Build Tools with Windows SDK 8.1

For freezing:

* pyinstaller



## Installation

To run MoFi directly in source:

* Execute `python setup.py build_ext --inplace` to create the `findmods` library.
* In the `docs` directory, execute `make html` (Linux) or `make.bat html` (Windows) to create the documentation.
* Start MoFi via executing `python run.py`.

To run MoFi from a frozen distribution:

* Extract the contents of the archive and start MoFi by executing `mofi`.



## Creating frozen distributions

* Create the library and documentation as described above.
* Execute `pyinstaller -i images\mofi.ico mofi.spec`.
* The folder `dist/mofi` is now a self-contained MoFi installation.
