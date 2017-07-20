import sys
import os.path

if getattr(sys, "frozen", False):
    module_dir = os.path.dirname(sys.executable)
else:
    module_dir = os.path.dirname(os.path.realpath(__file__))

config_dir = os.path.join(module_dir, "config")
data_dir = os.path.join(module_dir, "data")