"""
Find the LibXC shared object and imports it as core.
"""

import ctypes
import ctypes.util
import numpy as np
import os

# Attempt to load the compiled C code
__libxc_path = None
core = None

# First check the local folder
__abs_path = os.path.abspath(os.path.dirname(__file__))
for module_file in os.listdir(__abs_path):
    if module_file.startswith("libxc."):
        __libxc_path = os.path.join(__abs_path, module_file)

# If no libxc is local, check LD_LIBRARY_PATHS's
if __libxc_path is None:
    __libxc_path = ctypes.util.find_library("libxc")

# If we still havent found it, give up and throw an error
if __libxc_path is None:
    raise ImportError(
        "LibXC Shared object not found, searched Python module local directory and library paths"
    )

# Load the C object
core = ctypes.CDLL(__libxc_path)


def get_core_path():
    """
    Returns the path of the loaded LibXC shared object.
    """

    return __libxc_path