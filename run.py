#!/usr/bin/env python3

import sys
import os.path


launch_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, launch_path)

from mofi.modfinder import main
main()
