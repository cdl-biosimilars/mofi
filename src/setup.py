from cx_Freeze import setup, Executable


buildOptions = dict(
    includes=[
    #     "math",
    #     "os",
    #     "pickle",
    #     "re",
    #     "sys",
    #     "time",
    #     "webbrowser",
    #     "qtpy.QtWidgets",
    #     "qtpy.QtGui",
    #     "qtpy.QtCore",
    #     "numpy",
    #     "pandas",
    #     "matplotlib",
        "numpy.core._methods",
    ],
    include_files=[
        "config/",
        "findmods_source/"
    ],
    excludes=[
        "pkg_resources",
        "setuptools"
    ]
)

import sys
base = 'Win32GUI' if sys.platform == 'win32' else None

executables = [Executable('ModFinder.py', base=base)]

setup(name='mofi',
      version='1.0',
      description='finds modifications',
      options=dict(build_exe=buildOptions),
      executables=executables)
