"""
Binds the LibXC utility functions.
"""

import ctypes

from .core import core

# Set required
core.xc_version.restype = None
core.xc_version.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int)
]

def xc_version():
    """
    Returns the current LibXC version.
    """
    a = ctypes.c_int()
    b = ctypes.c_int()
    c = ctypes.c_int()
    core.xc_version(a, b, c)
    return (int(a.value), int(b.value), int(c.value))

