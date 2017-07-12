"""
configure.py

Helper functions/classes for global configurations.

Authors: Stefan Senn, Wolfgang Skala

(c) 2017 Christian Doppler Laboratory for Biosimilar Characterization
"""

from collections import OrderedDict
from configparser import ConfigParser
import os

# (1) general parameters
config = ConfigParser()
config.read("config/config.ini")
version = config.get("Version", "Version")
rights = config.get("Version", "Rights")
contact = config.get("Version", "Contact")
default_da = float(config.get("Defaults", "da"))
default_ppm = int(config.get("Defaults", "ppm"))
maxmods = int(config.get("Defaults", "maxmods"))
path = config.get("Defaults", "path")


# (2) mass sets
mass_set_parser = ConfigParser()
mass_set_parser.optionxform = str  # do not convert keys to lowercase
mass_set_parser.read("config/mass_sets.ini")

mass_sets = OrderedDict()
for set_name in mass_set_parser.sections():
    mass_sets[set_name] = {}
    for (atom, weight) in mass_set_parser.items(set_name):
        if atom != "description":
            mass_sets[set_name][atom] = float(weight)


def select_mass_set(name):
    """
    Change mass set (e.g., average, monoisotopic)
    to be used for mass calculations.

    :param name: section name in config/mass_sets.ini
    :return: nothing, but sets the global variable current_mass_set
    """
    global current_mass_set
    current_mass_set = mass_sets[name]


# An {atom name: weight} dict that can be accessed by other modules
current_mass_set = {}
select_mass_set(mass_set_parser.sections()[0])


# (3) default monomer libraries
def read_default_libraries(subdir):
    """
    Generate the default monomer or polymer libraries.
    Each file in ./config/monomers (or ./config/polymers, respectively)
    contains the description of a library.
    The file name is used as the name of a library.
    Sections indicate names, key/values pairs indicate data
    for the monomer/polymer table.

    :param subdir: either "monomers" or "polymers"
    :return: a dict containing all found libraries
    """
    result = OrderedDict()
    for filename in os.listdir("config/" + subdir):
        library_parser = ConfigParser()
        library_parser.read("config/" + subdir + "/" + filename)
        name = os.path.splitext(filename)[0]
        result[name] = OrderedDict()
        for section in library_parser.sections():
            result[name][section] = OrderedDict()
            for (k, v) in library_parser.items(section):
                result[name][section][k] = v
    return result

default_monomer_libraries = read_default_libraries("monomers")
default_polymer_libraries = read_default_libraries("polymers")

double_spin_box_flat_style = """
    QDoubleSpinBox::up-button {
        background: white
    }
    QDoubleSpinBox::up-arrow {
        image: url(:/mofi resource/images/UpArrow.png);
        width: 7px;
        height: 9px
    }
    QDoubleSpinBox::down-button {
        background: white
    }
    QDoubleSpinBox::down-arrow {
        image: url(:/mofi resource/images/DownArrow.png);
        width: 7px;
        height: 7px
    }
    """

spin_box_flat_style = """
    QSpinBox::up-button {
        background: white
    }
    QSpinBox::up-arrow {
        image: url(:/mofi resource/images/UpArrow.png);
        width: 7px;
        height: 9px
    }
    QSpinBox::down-button {
        background: white
    }
    QSpinBox::down-arrow {
        image: url(:/mofi resource/images/DownArrow.png);
        width: 7px;
        height: 7px
    }
    """