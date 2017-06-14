from setuptools import setup, Extension

c_ext = Extension(
    name="findmods",
    sources=["modsmodule.cpp", "modifications.cpp"],
    language="c++"
)

setup(
    name = "findmods",
    version = "1.0",
    ext_modules = [c_ext],
)
