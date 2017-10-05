"""
Manage the global configuration.

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
colors = dict(config.items("Colors"))


# (2) mass sets
mass_set_parser = ConfigParser()
mass_set_parser.optionxform = str  # do not convert keys to lowercase
with open(os.path.join(config_dir, "mass_sets.ini")) as f:
    mass_set_parser.read_file(f)

mass_sets = OrderedDict()
for set_name in mass_set_parser.sections():
    mass_sets[set_name] = {}
    for k, v in mass_set_parser.items(set_name):
        try:
            mass_sets[set_name][k] = float(v)
        except ValueError:
            mass_sets[set_name][k] = v


def select_mass_set(name):
    """
    Change mass set (e.g., average, monoisotopic)
    to be used for mass calculations.

    :param str name: section name in ``config/mass_sets.ini``
    :return: nothing, sets the global variable :data:`current_mass_set`
    """
    global current_mass_set
    current_mass_set = mass_sets[name]


# An {atom name: weight} dict that can be accessed by other modules
current_mass_set = {}
select_mass_set(mass_set_parser.sections()[0])


def spin_box_flat_style(bg="bg_ok", double=False):
    """
    Returns a style sheet for a :class:`QSpinBox` or :class:`QDoubleSpinBox`
    in flat style, i.e., without borders.

    :param str bg: background color as defined in ``config/config.ini``
    :param bool double: True: style sheet for a :class:`QDoubleSpinBox`.
                        False: style sheet for a :class:`QSpinBox`.
    :return: style sheet
    :rtype: str
    """

    if double:
        widget = "QDoubleSpinBox"
    else:
        widget = "QSpinBox"

    return """
    {widget} {{
        background: {bg}
    }}
    {widget}::up-button {{
        background: {bg}
    }}
    {widget}::up-arrow {{
        image: url(:/mofi resource/images/UpArrow.png);
        width: 7px;
        height: 9px
    }}
    {widget}::down-button {{
        background: {bg}
    }}
    {widget}::down-arrow {{
        image: url(:/mofi resource/images/DownArrow.png);
        width: 7px;
        height: 7px
    }}
    """.format(bg=colors[bg], widget=widget)
