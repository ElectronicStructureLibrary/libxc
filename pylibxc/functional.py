"""
Binds a LibXC Functional struct to a Python object
"""

import ctypes
import numpy as np

from .core import core
from . import flags
from . import util
from . import structs

### Bind required ctypes

# Allocation wrappers
core.xc_func_alloc.restype = ctypes.POINTER(structs.xc_func_type)

core.xc_func_init.argtype = (ctypes.POINTER(structs.xc_func_type), ctypes.c_int, ctypes.c_int)
core.xc_func_init.restype = ctypes.c_int

core.xc_func_end.argtype = (ctypes.POINTER(structs.xc_func_type))

core.xc_func_free.argtype = (ctypes.POINTER(structs.xc_func_type))

# Info wrappers
core.xc_func_get_info.argtype = (ctypes.POINTER(structs.xc_func_type))
core.xc_func_get_info.restype = ctypes.POINTER(structs.xc_func_info_type)

core.xc_func_get_info.argtype = (ctypes.POINTER(structs.xc_func_info_type))

core.xc_func_info_get_kind.argtype = (ctypes.POINTER(structs.xc_func_info_type))

core.xc_func_info_get_name.argtype = (ctypes.POINTER(structs.xc_func_info_type))
core.xc_func_info_get_name.restype = ctypes.c_char_p

core.xc_func_info_get_family.argtype = (ctypes.POINTER(structs.xc_func_info_type))

core.xc_func_info_get_flags.argtype = (ctypes.POINTER(structs.xc_func_info_type))

core.xc_func_info_get_references.argtype = (ctypes.POINTER(structs.xc_func_info_type), ctypes.c_int)
core.xc_func_info_get_references.restype = ctypes.POINTER(structs.func_reference_type)

# Setters
core.xc_func_info_get_n_ext_params.argtype = (ctypes.POINTER(structs.xc_func_info_type))

core.xc_func_info_get_ext_params_description.argtype = (ctypes.POINTER(structs.xc_func_info_type), ctypes.c_int)
core.xc_func_info_get_ext_params_description.restype = ctypes.c_char_p

core.xc_func_info_get_ext_params_default_value.argtype = (ctypes.POINTER(structs.xc_func_info_type), ctypes.c_int)
core.xc_func_info_get_ext_params_default_value.restype = ctypes.c_double

core.xc_func_set_ext_params.argtype = (ctypes.POINTER(structs.xc_func_type), ctypes.POINTER(ctypes.c_double))

core.xc_func_set_dens_threshold.argtype = (ctypes.POINTER(structs.xc_func_type), ctypes.c_double)

### Build LibXCFunctional class


class LibXCFunctional(object):
    def __init__(self, func_name, spin):
        """
        The primary LibXCFunctional class used to build and compute DFT exchange-correlation quantities.

        Parameters
        ----------
        func_name : int or str
            Either the functional name or ID used to create the LibXCFunctional.
        spin : int or str
            The spin of the requested functional either "unpolarized" (1) or polarized" (2).


        Returns
        -------
        func : LibXCFunctional
            A constructed LibXCFunctional.

        Examples
        --------
        >>> pylibxc.util.xc_family_from_id(72)
        (4, 3)


        """
        self.xc_func = None
        self._xc_func_init = False

        # Handle func_name
        if isinstance(func_name, str):
            func_id = util.xc_functional_get_number(func_name)
            if func_id == -1:
                raise KeyError("LibXC Functional name '%s' not found." % func_name)
        elif isinstance(func_name, int):
            func_id = func_name
            if util.xc_functional_get_name(func_name) is None:
                raise KeyError("LibXC Functional ID '%d' not found." % func_name)
        else:
            raise TypeError("func_name must either be a string or int.")

        # Handle spin
        if isinstance(spin, str):
            spin = spin.lower()
            if spin == "polarized":
                spin = 2
            elif spin == "unpolarized":
                spin = 1
            else:
                raise KeyError("Spin must either be 'polarized' or 'unpolarized' if represented by a string.")

        if spin not in [1, 2]:
            raise KeyError("Spin must either be 0 or 1 if represented by a integer.")

        # Build the LibXC functional
        self.xc_func = core.xc_func_alloc()
        ret = core.xc_func_init(self.xc_func, func_id, spin)
        if ret != 0:
            raise ValueError("LibXC Functional construction did not complete. Error code %d" % ret)
        self._xc_func_init = True

        # Unpack functional info
        self.xc_func_info = core.xc_func_get_info(self.xc_func)
        self.__number = core.xc_func_info_get_number(self.xc_func_info)
        self.__kind = core.xc_func_info_get_kind(self.xc_func_info)
        self.__name = core.xc_func_info_get_name(self.xc_func_info).decode("UTF-8")
        self.__family = core.xc_func_info_get_family(self.xc_func_info)
        self.__flags = core.xc_func_info_get_flags(self.xc_func_info)

        # Pull out references
        self.__refs = []
        self.__bibtexs = []
        self.__dois = []

        for pos in range(flags.XC_MAX_REFERENCES):
            ref = core.xc_func_info_get_references(self.xc_func_info, pos)
            if not ref: break

            self.__refs.append(ref.contents.ref.decode("UTF-8"))
            self.__bibtexs.append(ref.contents.bibtex.decode("UTF-8"))
            self.__dois.append(ref.contents.doi.decode("UTF-8"))

    def __del__(self):
        """
        Cleans up the LibXC C struct on deletion
        """
        if self.xc_func is None:
            return

        if self._xc_func_init:
            core.xc_func_end(self.xc_func)

        core.xc_func_free(self.xc_func)

    ### Getters

    def get_number(self):
        """
        Returns the LibXCFunctional ID.
        """

        return self.__number

    def get_kind(self):
        """
        Returns the LibXCFunctional kind.
        """

        return self.__kind

    def get_name(self):
        """
        Returns the LibXCFunctional name.
        """

        return self.__name

    def get_family(self):
        """
        Returns the LibXCFunctional family.
        """

        return self.__family

    def get_flags(self):
        """
        Returns the LibXCFunctional flags.
        """

        return self.__flags

    def get_references(self):
        """
        Returns the LibXCFunctional references.
        """

        return self.__refs

    def get_bibtex(self):
        """
        Returns the LibXCFunctional bibtex references.
        """

        return self.__bibtexs

    def get_doi(self):
        """
        Returns the LibXCFunctional reference DOIs.
        """

        return self.__dois

    ### Setters

    def get_ext_param_descriptions(self):
        """
        Gets the description of all external parameters
        """
        num_param = core.xc_func_info_get_n_ext_params(self.xc_func_info)

        ret = []
        for p in range(num_param):
            tmp = core.xc_func_info_get_ext_params_description(self.xc_func_info, p)
            ret.append(tmp.decode("UTF-8"))

        return ret

    def get_ext_param_default_values(self):
        """
        Gets the default value of all external parameters.
        """
        num_param = core.xc_func_info_get_n_ext_params(self.xc_func_info)

        ret = []
        for p in range(num_param):
            tmp = core.xc_func_info_get_ext_params_default_value(self.xc_func_info, p)
            ret.append(tmp)

        return ret

    def set_ext_params(self, ext_params):
        """
        Sets all external parameters.
        """
        num_param = core.xc_func_info_get_n_ext_params(self.xc_func_info)
        if num_param == 0:
            raise ValueError("The LibXCFunctional '%s' has no extermal parameters to set." % self.get_name())

        if len(ext_params) != num_param:
            raise ValueError("The length of the input external parameters (%d) does not match the length of the Functionals external parameters (%d)." % (len(ext_params), num_param))

        arr = np.array(ext_params, dtype=np.double)
        core.xc_func_set_ext_params(self.xc_func, arr.ctypes.data)

    def set_dens_threshold(self, dens_threshold):
        """
        Sets the density threshold in which densities will not longer be computer.
        """

        if dens_threshold < 0:
            raise ValueError("The density threshold cannot be smaller than 0.")

        core.xc_func_set_dens_threshold(self.xc_func, ctypes.c_double(dens_threshold))

# void xc_func_set_dens_threshold(xc_func_type *p, double dens_threshold);
