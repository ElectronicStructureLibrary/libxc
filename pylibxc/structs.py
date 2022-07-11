"""
Wrappers to the LibXC C structs
"""

import ctypes
import numpy as np

xc_inputs = ["rho sigma lapl tau exx"]
xc_outputs = [
    # order 0
    "zk",
    # order 1
    "vrho vsigma vlapl vtau vexx",
    # order 2
    "v2rho2 v2rhosigma v2rholapl v2rhotau v2rhoexx " +
    "v2sigma2 v2sigmalapl v2sigmatau v2sigmaexx " +
    "v2lapl2 v2lapltau v2laplexx " +
    "v2tau2 v2tauexx " +
    "v2exx2",
    # order 3
    "v3rho3 v3rho2sigma v3rho2lapl v3rho2tau v3rho2exx " +
    "v3rhosigma2 v3rhosigmalapl v3rhosigmatau v3rhosigmaexx " +
    "v3rholapl2 v3rholapltau v3rholaplexx " +
    "v3rhotau2 v3rhotauexx " +
    "v3rhoexx2 " +
    "v3sigma3 v3sigma2lapl v3sigma2tau v3sigma2exx " +
    "v3sigmalapl2 v3sigmalapltau v3sigmalaplexx " +
    "v3sigmatau2 v3sigmatauexx " +
    "v3sigmaexx2 " +
    "v3lapl3 v3lapl2tau v3lapl2exx " +
    "v3lapltau2 v3lapltauexx " +
    "v3laplexx2 " +
    "v3tau3 v3tau2exx v3tauexx2 v3exx3",
    # order 4
    "v4rho4 v4rho3sigma v4rho3lapl v4rho3tau v4rho3exx " +
    "v4rho2sigma2 v4rho2sigmalapl v4rho2sigmatau v4rho2sigmaexx " +
    "v4rho2lapl2 v4rho2lapltau v4rho2laplexx " +
    "v4rho2tau2 v4rho2tauexx " +
    "v4rho2exx2 " +
    "v4rhosigma3 v4rhosigma2lapl v4rhosigma2tau v4rhosigma2exx " +
    "v4rhosigmalapl2 v4rhosigmalapltau v4rhosigmalaplexx " +
    "v4rhosigmatau2 v4rhosigmatauexx " +
    "v4rhosigmaexx2 " +
    "v4rholapl3, v4rholapl2tau, v4rholapl2exx " +
    "v4rholapltau2, v4rholapltauexx " +
    "v4rholaplexx2 " +
    "v4rhotau3 v4rhotau2exx v4rhoexx3 " +
    "v4sigma4 v4sigma3lapl v4sigma3tau v4sigma3exx " +
    "v4sigma2lapl2 v4sigma2lapltau v4sigma2laplexx " +
    "v4sigma2tau2 v4sigma2tauexx " +
    "v4sigma2exx2 " +
    "v4sigmalapl3 v4sigmalapl2tau v4sigmalapl2exx " +
    "v4sigmalapltau2 v4sigmalapltauexx " +
    "v4sigmalaplexx2 " +
    "v4sigmatau3 v4sigmatau2exx v4sigmatauexx2 v4sigmaexx3 " +
    "v4lapl4 v4lapl3tau v4lapl3exx " +
    "v4lapl2tau2 v4lapl2tauexx v4lapl2exx2 " +
    "v4lapltau3 v4lapltau2exx v4lapltauexx2 v4laplexx3 " +
    "v4tau4 v4tau3exx v4tauexx3 v4exx4"
]


class func_reference_type(ctypes.Structure):
    """
    Holds reference data for the LibXC functional
    """
    _fields_ = [("ref", ctypes.c_char_p), ("doi", ctypes.c_char_p),
                ("bibtex", ctypes.c_char_p)]


class func_params_type(ctypes.Structure):
    """
    Holds user defined parameters and their description.
    """
    _fields_ = [("n", ctypes.c_int), ("names", ctypes.POINTER(ctypes.c_char_p)),
        ("descriptions", ctypes.POINTER(ctypes.c_char_p)), ("values", ctypes.POINTER(ctypes.c_double)),
        ("set", ctypes.c_void_p)
    ]


class xc_functionals_work_variants(ctypes.Structure):
    _fields_ = [("unpol", ctypes.c_void_p * 5), ("pol", ctypes.c_void_p * 5)]


class xc_func_info_type(ctypes.Structure):
    """
    Holds LibXC information about a single functional primitive.
    """
    _fields_ = [
        ("number", ctypes.c_int),
        ("kind", ctypes.c_int),
        ("name", ctypes.c_char_p),
        ("family", ctypes.c_int),
        ("refs", ctypes.POINTER(func_reference_type) * 5),
        ("flags", ctypes.c_int),
        ("dens_threshold", ctypes.c_double),
        ("ext_params", func_params_type),
        ("init", ctypes.c_void_p),
        ("end", ctypes.c_void_p),
        ("work", xc_functionals_work_variants),
    ]


class xc_dimensions_expl(ctypes.Structure):  
    """
    Named fields
    """
    _fields_ = [(i, ctypes.c_int) for i in " ".join(xc_inputs + xc_outputs).split()]
   
class xc_dimensions(ctypes.Union):
    """
    Holds dimensions of the several arrays.
    """

    _fields_ = [("named", xc_dimensions_expl), ("vec",  ctypes.POINTER(ctypes.c_int))]

class xc_func_type(ctypes.Structure):
    """
    The primary xc_func_type used to hold all data pertaining to a given
    LibXC functional
    """
    _fields_ = [
        ("info", ctypes.POINTER(xc_func_info_type)),  # const xc_func_info_type *info;
        ("nspin", ctypes.c_int),
        ("n_func_aux", ctypes.c_int),
        ("xc_func_type", ctypes.c_void_p),
        ("mix_coef", ctypes.POINTER(ctypes.c_double)),

        # Hybrids
        ("hyb_number_terms", ctypes.c_int),
        ("hyb_type", ctypes.POINTER(ctypes.c_int)),
        ("hyb_coeff", ctypes.POINTER(ctypes.c_double)),
        ("hyb_omega", ctypes.POINTER(ctypes.c_double)),

        # VV10
        ("nlc_b", ctypes.c_double),
        ("nlc_C", ctypes.c_double),

        ("dim", ctypes.POINTER(xc_dimensions)),
        
        # parameters
        ("ext_params", ctypes.POINTER(ctypes.c_double)),
        ("params", ctypes.c_void_p),  # void *params;
        ("dens_threshold", ctypes.c_double),
        ("zeta_threshold", ctypes.c_double),
        ("sigma_threshold", ctypes.c_double),
        ("tau_threshold", ctypes.c_double),
    ]
