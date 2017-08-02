# ModFinder

A tool to find molecular modifications (like glycans)
in mass spectra of intact proteins.

## Requirements

* Python 3.5
* numpy
* pandas
* matplotlib
* PyQt5
* xlrd

For building the C++ extension:

* Unix: gcc / clang
* Windows: Visual C++ 2015 Build Tools with Windows SDK 8.1

## Installing

    pip3 install .

## Running in source

Build the extension using `python3 setup.py build` and copy the library from
`build/lib.linux-x86_64-3.5/findmods.cpython-35m-x86_64-linux-gnu.so` (adjust
the path as necessary) to `mofi/`. Then run `run.py`.
