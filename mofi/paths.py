"""
Set the correct system path for frozen versions.
"""

import sys
import os.path

if getattr(sys, "frozen", False):
    module_dir = os.path.dirname(sys.executable)
else:
    module_dir = os.path.dirname(os.path.realpath(__file__))

config_dir = os.path.join(module_dir, "config")
data_dir = os.path.join(module_dir, "data")
docs_dir = os.path.join(module_dir, "docs")
