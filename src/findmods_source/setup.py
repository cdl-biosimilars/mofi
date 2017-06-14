# The command
#     run python setup.py build_ext --inplace
# creates findmods.pyd which may be imported as a module

from setuptools import setup, Extension

c_ext = Extension(
    name="findmods",
    sources=["modsmodule.cpp", "modifications.cpp"],
    language="c++"
)

setup(
    name="findmods",
    version="1.0",
    description="This package implements a fast combinatorial search",
    ext_modules=[c_ext]
)
