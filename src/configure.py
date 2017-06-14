"""
    configure.py

    Helper functions/classes for global configurations.

    Authors: Stefan Senn, Wolfgang Skala

    (c) 2017 Christian Doppler Laboratory for Biosimilar Characterization
"""

import configparser

config = configparser.ConfigParser()
config.read("E:/Projects/ModFinder/config/config.ini")
# Config.read('../config/config.ini')


# extract information from config.ini
version = config.get("Version", "Version")
rights = config.get("Version", "Rights")
contact = config.get("Version", "Contact")
default_masses = config.get("Defaults", "DefaultMasses")
default_da = float(config.get("Defaults", "da"))
default_ppm = int(config.get("Defaults", "ppm"))
maxmods = int(config.get("Defaults", "maxmods"))
path = config.get("Defaults", "path")


def _read_atomic_masses(mass_set):
    """
    Read atomic masses from the config file.

    :param mass_set: mass set from config.ini
    :return: dict of atomic masses
    """
    return {a: float(config[mass_set][a]) for a in "CHNOPS"}


def set_average_masses(mass_set):
    """
    Change average masses to be used for mass calculations.

    :param mass_set: mass set from config.ini
    :return: nothing, but sets the global variable atomic_masses
    """
    global average_masses
    average_masses = _read_atomic_masses(mass_set)


# The following dictionaries of atom name=weight pairs can be accessed by other modules.
# atomic_masses stores the current set of atomic masses (e.g., AtomsIUPAC) as defined in config.ini
# monoisotopic_masses stores monoisotopic masses.
average_masses = {}
set_average_masses(default_masses)
monoisotopic_masses = _read_atomic_masses("AtomsMonoisotopic")
