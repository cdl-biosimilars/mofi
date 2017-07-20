#!/usr/bin/env python3

import sys
import os.path

try:
    from mofi.modfinder import main
except ImportError as e:
    if not getattr(sys, "frozen", False):
        launch_path = os.path.dirname(os.path.realpath(__file__))
        sys.path.insert(0, launch_path)
    else:
        raise e
    from mofi.modfinder import main

main()
