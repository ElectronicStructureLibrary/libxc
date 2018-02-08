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

# Build out a few common tmps
__ndptr = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags=("C", "A"))

__ndptr_w = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags=("W", "C", "A"))

__xc_func_p = ctypes.POINTER(structs.xc_func_type)
__xc_func_info_p = ctypes.POINTER(structs.xc_func_info_type)

# Allocation wrappers
core.xc_func_alloc.restype = __xc_func_p

core.xc_func_init.argtypes = (__xc_func_p, ctypes.c_int, ctypes.c_int)
core.xc_func_init.restype = ctypes.c_int

core.xc_func_end.argtypes = (__xc_func_p, )

core.xc_func_free.argtypes = (__xc_func_p, )

# Info wrappers
core.xc_func_get_info.argtypes = (__xc_func_p, )
core.xc_func_get_info.restype = __xc_func_info_p

core.xc_func_info_get_kind.argtypes = (__xc_func_info_p, )

core.xc_func_info_get_name.argtypes = (__xc_func_info_p, )
core.xc_func_info_get_name.restype = ctypes.c_char_p

core.xc_func_info_get_family.argtypes = (__xc_func_info_p, )

core.xc_func_info_get_flags.argtypes = (__xc_func_info_p, )

core.xc_func_info_get_references.argtypes = (__xc_func_info_p, ctypes.c_int)
core.xc_func_info_get_references.restype = ctypes.POINTER(structs.func_reference_type)

# Setters
core.xc_func_info_get_n_ext_params.argtypes = (__xc_func_info_p, )

core.xc_func_info_get_ext_params_description.argtypes = (__xc_func_info_p, ctypes.c_int)
core.xc_func_info_get_ext_params_description.restype = ctypes.c_char_p

core.xc_func_info_get_ext_params_default_value.argtypes = (__xc_func_info_p, ctypes.c_int)
core.xc_func_info_get_ext_params_default_value.restype = ctypes.c_double

core.xc_func_set_ext_params.argtypes = (__xc_func_p, __ndptr)

core.xc_func_set_dens_threshold.argtypes = (__xc_func_p, ctypes.c_double)

# LDA computers

core.xc_lda_exc.argtypes = (__xc_func_p, ctypes.c_int, __ndptr, __ndptr_w)

# core.xc_lda_exc.argtypes = [__xc_func_p,
#                            ctypes.c_int,
#                            np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags=("C", "A"), ),
#                            np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags=("W", "C", "A"))
#                           ]
# void xc_lda        (const xc_func_type *p, int np, const double *rho, double *zk, double *vrho, double *v2rho2, double *v3rho3);
# void xc_lda_exc    (const xc_func_type *p, int np, const double *rho, double *zk);
# void xc_lda_exc_vxc(const xc_func_type *p, int np, const double *rho, double *zk, double *vrho);
# void xc_lda_vxc    (const xc_func_type *p, int np, const double *rho, double *vrho);
# void xc_lda_fxc    (const xc_func_type *p, int np, const double *rho, double *v2rho2);
# void xc_lda_kxc    (const xc_func_type *p, int np, const double *rho, double *v3rho3);

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
        self.xc_func_size_names = [x for x in dir(self.xc_func.contents) if "n_" in x]

        # Set all int attributes to zero
        for attr in self.xc_func_size_names:
            setattr(self.xc_func, attr, 0)

        ret = core.xc_func_init(self.xc_func, func_id, spin)
        if ret != 0:
            raise ValueError("LibXC Functional construction did not complete. Error code %d" % ret)
        self._xc_func_init = True

        # Pull out all sizes
        self.xc_func_sizes = {}
        for attr in self.xc_func_size_names:
            self.xc_func_sizes[attr] = getattr(self.xc_func, attr)

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
            raise ValueError(
                "The length of the input external parameters (%d) does not match the length of the Functionals external parameters (%d)."
                % (len(ext_params), num_param))

        core.xc_func_set_ext_params(self.xc_func, np.asarray(ext_params, dtype=np.double))

    def set_dens_threshold(self, dens_threshold):
        """
        Sets the density threshold in which densities will not longer be computer.
        """

        if dens_threshold < 0:
            raise ValueError("The density threshold cannot be smaller than 0.")

        core.xc_func_set_dens_threshold(self.xc_func, ctypes.c_double(dens_threshold))

    def compute(self, inp, output):

        core.xc_lda_exc(self.xc_func, inp.shape[0], inp, output)


# void xc_func_set_dens_threshold(xc_func_type *p, double dens_threshold);
