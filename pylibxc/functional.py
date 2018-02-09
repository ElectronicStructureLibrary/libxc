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
__ndptr = np.ctypeslib.ndpointer(dtype=np.double, flags=("C", "A"))
__ndptr_w = np.ctypeslib.ndpointer(dtype=np.double, flags=("W", "C", "A"))  # Writable

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
core.xc_lda.argtypes = (__xc_func_p, ctypes.c_int, __ndptr, __ndptr_w, __ndptr_w, __ndptr_w, __ndptr_w)
core.xc_lda_exc_vxc.argtypes = (__xc_func_p, ctypes.c_int, __ndptr, __ndptr_w, __ndptr_w)
core.xc_lda_exc.argtypes = (__xc_func_p, ctypes.c_int, __ndptr, __ndptr_w)
core.xc_lda_vxc.argtypes = (__xc_func_p, ctypes.c_int, __ndptr, __ndptr_w)
core.xc_lda_fxc.argtypes = (__xc_func_p, ctypes.c_int, __ndptr, __ndptr_w)
core.xc_lda_kxc.argtypes = (__xc_func_p, ctypes.c_int, __ndptr, __ndptr_w)

# GGA computers
core.xc_gga.argtypes = (__xc_func_p, ctypes.c_int, __ndptr, __ndptr, *([__ndptr_w] * 10))
core.xc_gga_exc_vxc.argtypes = (__xc_func_p, ctypes.c_int, __ndptr, __ndptr, *([__ndptr_w] * 3))
core.xc_gga_exc.argtypes = (__xc_func_p, ctypes.c_int, __ndptr, __ndptr, *([__ndptr_w] * 1))
core.xc_gga_vxc.argtypes = (__xc_func_p, ctypes.c_int, __ndptr, __ndptr, *([__ndptr_w] * 2))
core.xc_gga_fxc.argtypes = (__xc_func_p, ctypes.c_int, __ndptr, __ndptr, *([__ndptr_w] * 3))
core.xc_gga_kxc.argtypes = (__xc_func_p, ctypes.c_int, __ndptr, __ndptr, *([__ndptr_w] * 4))

# MGGA computers
core.xc_mgga.argtypes = (__xc_func_p, ctypes.c_int, *([__ndptr] * 4), *([__ndptr_w] * 15))
core.xc_mgga_exc_vxc.argtypes = (__xc_func_p, ctypes.c_int, *([__ndptr] * 4), *([__ndptr_w] * 5))
core.xc_mgga_exc.argtypes = (__xc_func_p, ctypes.c_int, *([__ndptr] * 4), *([__ndptr_w] * 1))
core.xc_mgga_vxc.argtypes = (__xc_func_p, ctypes.c_int, *([__ndptr] * 4), *([__ndptr_w] * 4))
core.xc_mgga_fxc.argtypes = (__xc_func_p, ctypes.c_int, *([__ndptr] * 4), *([__ndptr_w] * 10))

### Build LibXCFunctional class


def _check_arrays(current_arrays, required, sizes, factor):
    """
    A specialized function built to construct and check the sizes of arrays given to the LibXCFunctional class.
    """

    # Nothing supplied so we build it out
    if current_arrays is None:
        current_arrays = {}
        for label in required:
            size = sizes["n_" + label]
            current_arrays[label] = np.zeros((size, factor))

    # Supplied arrays, check sizes
    else:
        missing = set(required) - set(current_arrays)
        if len(missing):
            raise KeyError("Missing the following output arrays: %s" % ", ".join(missing))

        for label in required:
            size = sizes["n_" + label] * factor
            if size != current_arrays[label].size:
                raise ValueError("Supplied output array '%s' does not have the correct shape number of points by %d" %
                                 (label, size))

    ret = [current_arrays[x] for x in required]
    return ret


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
                self.__spin = 2
            elif spin == "unpolarized":
                self.__spin = 1
            else:
                raise KeyError("Spin must either be 'polarized' or 'unpolarized' if represented by a string.")
        else:
            self.__spin = spin

        if self.__spin not in [1, 2]:
            raise KeyError("Spin must either be 1 or 2 if represented by a integer.")

        # Build the LibXC functional
        self.xc_func = core.xc_func_alloc()
        self.xc_func_size_names = [x for x in dir(self.xc_func.contents) if "n_" in x]

        # Set all int attributes to zero
        for attr in self.xc_func_size_names:
            setattr(self.xc_func.contents, attr, 0)

        ret = core.xc_func_init(self.xc_func, func_id, self.__spin)
        if ret != 0:
            raise ValueError("LibXC Functional construction did not complete. Error code %d" % ret)
        self._xc_func_init = True

        # Pull out all sizes
        self.xc_func_sizes = {}
        for attr in self.xc_func_size_names:
            self.xc_func_sizes[attr] = getattr(self.xc_func.contents, attr)

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

    def compute(self, inp, output=None, do_exc=True, do_vxc=True, do_fxc=False, do_kxc=False):

        # Parse input
        if isinstance(inp, np.ndarray):
            inp = {"rho": np.asarray(inp, dtype=np.double)}
        elif isinstance(inp, dict):
            inp = {k: np.asarray(v, dtype=np.double) for k, v in inp.items()}
        else:
            raise KeyError("Input must have a 'rho' variable or a single array.")

        # How long are we?
        npoints = int(inp["rho"].size / self.__spin)
        if (inp["rho"].size % self.__spin):
            raise ValueError("Rho input has an invalid shape, must be divisible by %d" % self.__spin)

        # Find the right compute function
        args = [self.xc_func, ctypes.c_int(npoints)]
        if self.get_family() == flags.XC_FAMILY_LDA:

            # Build input args
            required_input = ["rho"]
            args.extend(_check_arrays(inp, required_input, self.xc_func_sizes, npoints))
            input_num_args = len(args)

            # Hybrid computers
            if do_exc and do_vxc and do_fxc and do_kxc:
                required_fields = ["zk", "vrho", "v2rho2", "v3rho3"]
                args.extend(_check_arrays(output, required_fields, self.xc_func_sizes, npoints))
                core.xc_lda(*args)
                do_exc = do_vxc = do_fxc = do_kxc = False

            if do_exc and do_vxc:
                required_fields = ["zk", "vrho"]
                args.extend(_check_arrays(output, required_fields, self.xc_func_sizes, npoints))
                core.xc_lda_exc_vxc(*args)
                do_exc = do_vxc = False

            # Individual computers
            if do_exc:
                required_fields = ["zk"]
                args.extend(_check_arrays(output, required_fields, self.xc_func_sizes, npoints))
                core.xc_lda_exc(*args)
            if do_vxc:
                required_fields = ["vrho"]
                args.extend(_check_arrays(output, required_fields, self.xc_func_sizes, npoints))
                core.xc_lda_vxc(*args)
            if do_fxc:
                required_fields = ["v2rho2"]
                args.extend(_check_arrays(output, required_fields, self.xc_func_sizes, npoints))
                core.xc_lda_fxc(*args)
            if do_kxc:
                required_fields = ["v3rho3"]
                args.extend(_check_arrays(output, required_fields, self.xc_func_sizes, npoints))
                core.xc_lda_kxc(*args)

        elif self.get_family() in [flags.XC_FAMILY_GGA, flags.XC_FAMILY_HYB_GGA]:

            # Build input args
            required_input = ["rho", "sigma"]
            args.extend(_check_arrays(inp, required_input, self.xc_func_sizes, npoints))
            input_num_args = len(args)

            # Hybrid computers
            if do_exc and do_vxc and do_fxc and do_kxc:
                required_fields = [
                    "zk", "vrho", "vsigma", "v2rho2", "v2rhosigma", "v2sigma2", "v3rho3", "v3rho2sigma", "v3rhosigma2",
                    "v3sigma3"
                ]
                args.extend(_check_arrays(output, required_fields, self.xc_func_sizes, npoints))
                core.xc_gga(*args)
                do_exc = do_vxc = do_fxc = do_kxc = False

            if do_exc and do_vxc:
                required_fields = ["zk", "vrho", "vsigma"]
                args.extend(_check_arrays(output, required_fields, self.xc_func_sizes, npoints))
                core.xc_gga_exc_vxc(*args)
                do_exc = do_vxc = False

            # Individual computers
            if do_exc:
                required_fields = ["zk"]
                args.extend(_check_arrays(output, required_fields, self.xc_func_sizes, npoints))
                core.xc_gga_exc(*args)
            if do_vxc:
                required_fields = ["vrho", "vsigma"]
                args.extend(_check_arrays(output, required_fields, self.xc_func_sizes, npoints))
                core.xc_gga_vxc(*args)
            if do_fxc:
                required_fields = ["v2rho2", "v2rhosigma", "v2sigma2"]
                args.extend(_check_arrays(output, required_fields, self.xc_func_sizes, npoints))
                core.xc_gga_fxc(*args)
            if do_kxc:
                required_fields = ["v3rho3", "v3rho2sigma", "v3rhosigma2", "v3sigma3"]
                args.extend(_check_arrays(output, required_fields, self.xc_func_sizes, npoints))
                core.xc_gga_kxc(*args)

        elif self.get_family() in [flags.XC_FAMILY_MGGA, flags.XC_FAMILY_HYB_MGGA]:
            # Build input args
            required_input = ["rho", "sigma", "lapl", "tau"]
            args.extend(_check_arrays(inp, required_input, self.xc_func_sizes, npoints))
            input_num_args = len(args)

            # Hybrid computers

            # Wait until FXC and KXC are available
            # if do_exc and do_vxc and do_fxc and do_kxc:
            #     required_fields = [
            #         "zk", "vrho", "vsigma", "vlapl", "vtau", "v2rho2", "v2sigma2", "v2lapl2", "v2tau2", "v2rhosigma",
            #         "v2rholapl", "v2rhotau", "v2sigmalapl", "v2sigmatau", "v2lapltau"
            #     ]
            #     args.extend(_check_arrays(output, required_fields, self.xc_func_sizes, npoints))
            #     core.xc_mgga(*args)
            #     do_exc = do_vxc = do_fxc = do_kxc = False

            if do_exc and do_vxc:
                required_fields = ["zk", "vrho", "vsigma", "vlapl", "vtau"]
                args.extend(_check_arrays(output, required_fields, self.xc_func_sizes, npoints))
                core.xc_mgga_exc_vxc(*args)
                do_exc = do_vxc = False

            # Individual computers
            if do_exc:
                required_fields = ["zk"]
                args.extend(_check_arrays(output, required_fields, self.xc_func_sizes, npoints))
                core.xc_mgga_exc(*args)
            if do_vxc:
                required_fields = ["vrho", "vsigma", "vlapl", "vtau"]
                args.extend(_check_arrays(output, required_fields, self.xc_func_sizes, npoints))
                core.xc_mgga_vxc(*args)
            if do_fxc:
                raise KeyError("FXC quantities (2rd derivitives) are not defined for MGGA's! (%d)")
                # required_fields = [
                #     "v2rho2", "v2sigma2", "v2lapl2", "v2tau2", "v2rhosigma", "v2rholapl", "v2rhotau", "v2sigmalapl",
                #     "v2sigmatau", "v2lapltau"
                # ]
                # args.extend(_check_arrays(output, required_fields, self.xc_func_sizes, npoints))
                # core.xc_mgga_fxc(*args)
            if do_kxc:
                raise KeyError("KXC quantities (3rd derivitives) are not defined for MGGA's! (%d)")
                # required_fields = ["v3rho3", "v3rho2sigma", "v3rhosigma2", "v3sigma3"]
                # args.extend(_check_arrays(output, required_fields, self.xc_func_sizes, npoints))
                # core.xc_gga_kxc(*args)
        else:
            raise KeyError("Functional kind not recognized! (%d)" % self.get_kind())

        # Return a dictionary
        return {k: v for k, v in zip(required_fields, args[input_num_args:])}
