from setuptools import setup, Extension
import glob

require = [
    "PyQt5",
    "numpy",
    "pandas",
    "matplotlib",
    "xlrd",
    "xlsxwriter"
]

c_ext = Extension(
    name="mofi.findmods",
    sources=["mofi/findmods/modsmodule.cpp", "mofi/findmods/modifications.cpp"],
    language="c++"
)

setup(
    name='mofi',
    version='1.0',
    description='A software tool for annotating glycoprotein mass spectra by integrating hybrid data from the intact protein and glycopeptide level',
    packages=['mofi'],
    install_requires=require,
    ext_modules=[c_ext],
    entry_points = {
        'gui_scripts': [
            "mofi = mofi.mofi:main"
        ]
    },
    include_package_data=True,
    package_data={
        "": ["config/*", "data/glycans/*", "data/modifications/*"]
    }
)
