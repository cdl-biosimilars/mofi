import sys
from cx_Freeze import setup, Executable

build_exe_options = {
    "includes": [
        "numpy.core._methods" # Fix import error
    ],
    "include_files": [
        ("mofi/config", "config"),
        ("mofi/data", "data"),
        ("mofi/docs", "docs"),
    ]
}

if sys.platform == "win32":
    base = "Win32GUI"
else:
    base = None

setup(
    name="mofi",
    version="1.0",
    description="finds modifications",
    options={"build_exe": build_exe_options},
    executables=[Executable("run.py", base=base)]
)
