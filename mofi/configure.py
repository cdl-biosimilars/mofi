"""
configure.py

Helper functions/classes for global configurations.

Authors: Stefan Senn, Wolfgang Skala

(c) 2017 Christian Doppler Laboratory for Biosimilar Characterization
"""

from collections import OrderedDict
from configparser import ConfigParser
import os.path
from mofi.paths import config_dir

# (1) general parameters
config = ConfigParser()
with open(os.path.join(config_dir, "config.ini")) as f:
    config.read_file(f)
default_da = float(config.get("Defaults", "da"))
default_ppm = int(config.get("Defaults", "ppm"))
maxmods = int(config.get("Defaults", "maxmods"))
path = config.get("Defaults", "path")


# (2) mass sets
mass_set_parser = ConfigParser()
mass_set_parser.optionxform = str  # do not convert keys to lowercase
with open(os.path.join(config_dir, "mass_sets.ini")) as f:
    mass_set_parser.read_file(f)

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
