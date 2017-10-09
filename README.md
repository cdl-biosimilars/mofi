# ModFinder

A tool to find molecular modifications (like glycans) in mass spectra of intact proteins.



## Requirements

* Python 3.5
* numpy
* pandas
* matplotlib
* PyQt5
* xlrd
* Sphinx (for creating the documentation)
* sphinx_rtd_theme (for creating the documentation)

For building the C++ extension:

* Unix: gcc / clang
* Windows: Visual C++ 2015 Build Tools with Windows SDK 8.1

For freezing:

* pyinstaller (Unix)
* cx_freeze (Windows)



## Installation

You may either run ModFinder directly from source or execute the frozen version from a stand-alone distribution.


### Running in source (Unix/Windows)

* Run `python3 setup.py build` to create the `findmods` library
* Copy the library (`findmods[...].so` on Unix, `findmods[...].pyd` on Windows) from `build/lib[...]/mofi/` to `mofi/`
* Run `make html` (Unix) or `make.bat html` (Windows) in the `docs/` directory to create the documentation
* Start ModFinder via `run.py`




### Running from frozen distribution

* Extract the contents of the archive
* Start ModFinder by executing `ModFinder` (on Windows) or `run` (on Unix).



## Creating frozen distributions and Debian packages

### Freezing

#### General

* Create and copy the library as described above


#### Windows)

* Run `./package.sh`
* The folder `build/mofi-windows/` is now a self-contained ModFinder installation, which can be compressed (if desired) and distributed.
* Start the program by double-clicking `ModFinder`.


#### Unix

* Run `pyinstaller run.spec`
* The folder `dist/run/` is now a self-contained ModFinder installation, which can be compressed (if desired) and distributed.
* Start the program by double-clicking `run`.


### Debian package

* Add new message in `debian/changelog`. Ensure that a key with the appropriate user ID exists.
* Tag that version using `git tag -a v1.0`
* Commit the changes
* Run `gpb buildpackage`
