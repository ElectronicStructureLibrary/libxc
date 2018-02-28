"""
Binds the LibXC utility functions.
"""

import ctypes
import numpy as np

from .core import core
from . import flags

### Set required ctypes bindings

core.xc_version.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int))
core.xc_version.restype = None

core.xc_version_string.restype = ctypes.c_char_p

core.xc_functional_get_number.argtypes = (ctypes.c_char_p, )
core.xc_functional_get_number.restype = ctypes.c_int

core.xc_functional_get_name.argtypes = (ctypes.c_int, )
core.xc_functional_get_name.restype = ctypes.c_char_p

core.xc_family_from_id.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int))
core.xc_family_from_id.restype = ctypes.c_int

core.xc_available_functional_numbers.argtypes = (np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags=("W", "C",
                                                                                                       "A")), )

core.xc_available_functional_names.argtypes = (ctypes.POINTER(ctypes.c_char_p), )

### Build wrapper functions


def xc_version():
    """
    Returns the current LibXC version as semantic versioning tuple.

    Returns
    -------
    version : tuple
        The (major, minor, patch) version of the linked LibXC shared object.

    Examples
    --------
    >>> pylibxc.util.xc_version()
    (4, 0, 1)

    """
    major = ctypes.c_int()
    minor = ctypes.c_int()
    patch = ctypes.c_int()
    core.xc_version(major, minor, patch)
    return (major.value, minor.value, patch.value)


def xc_version_string():
    """
    Returns the current LibXC version as a string.

    Returns
    -------
    version : str
        The string representation of the current LibXC version.

    Examples
    --------
    >>> pylibxc.util.xc_version_string()
    "4.0.1"

    """
    return core.xc_version_string().decode("UTF-8")


def xc_functional_get_number(name):
    """
    Returns the functional ID from a given string.

    Parameters
    ----------
    name : str
        The functional name to find the ID of.

    Returns
    -------
    id : int
        The ID of the requested functional.

    Examples
    --------
    >>> pylibxc.util.xc_functional_get_number("XC_GGA_X_GAM")
    32
    """

    return core.xc_functional_get_number(ctypes.c_char_p(name.encode()))


def xc_functional_get_name(func_id):
    """
    Returns the functional name from a ID.

    Parameters
    ----------
    func_id : int
        The LibXC functional ID

    Returns
    -------
    functional_name : str
        The functional_name of the requested functional.

    Examples
    --------
    >>> pylibxc.util.xc_functional_get_name(32)
    "gga_x_gam"
    """

    ret = core.xc_functional_get_name(func_id)
    if ret is None:
        return ret
    else:
        return ret.decode("UTF-8")


def xc_family_from_id(func_id):
    """
    Returns the family class and family number (?).

    Parameters
    ----------
    func_id : int
        The LibXC functional ID

    Returns
    -------
    functional_family : int
        The family ID.
    functional_number : int
        The functional number within a family.

    Examples
    --------
    >>> pylibxc.util.xc_family_from_id(72)
    (4, 3)

    """
    family = ctypes.c_int()
    number = ctypes.c_int()
    core.xc_family_from_id(func_id, ctypes.pointer(family), ctypes.pointer(number))

    return (family.value, number.value)


def xc_number_of_functionals():
    """
    Returns the totaly number of XC functionals available in LibXC.

    Returns
    -------
    number_of_functinals : int
        The total number of functionals available in LibXC.

    Examples
    --------
    >>> pylibxc.util.xc_family_from_id(72)
    (4, 3)
    """

    return core.xc_number_of_functionals()


def xc_available_functional_numbers():
    """
    Returns a list of all available XC functional IDs


    Returns
    -------
    functional_ids : list of ints
        All available LibXC functional IDs.

    Examples
    --------
    >>> xc_func_list = pylibxc.util.xc_available_functional_numbers()
    np.array([1, 2, ..., 568, 569])
    """

    nfunc = xc_number_of_functionals()

    ret = np.zeros(nfunc, dtype=np.int32)
    core.xc_available_functional_numbers(ret)
    return ret


def xc_available_functional_names():
    """
    Returns a list of all available XC functional names

    Returns
    -------
    functional_names : list of strs
        All available LibXC functional names.

    Examples
    --------
    >>> xc_func_list = pylibxc.util.xc_available_functional_names()
    [1, 2, ..., 568, 569]
    """

    # I give up trying to get char** working, someone else can pick it up.

    func_ids = xc_available_functional_numbers()
    return [xc_functional_get_name(x) for x in func_ids]

