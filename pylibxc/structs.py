"""
Wrappers to the LibXC C structs
"""

import ctypes
import numpy as np


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
    _fields_ = [("value", ctypes.c_double), ("description", ctypes.c_char_p)]


class xc_func_info_type(ctypes.Structure):
    """
    Holds LibXC information about a single functional primitive.
    """
    _fields_ = [
        ("number", ctypes.c_int),
        ("kind", ctypes.c_int),
        ("name", ctypes.c_char_p),
        ("family", ctypes.c_int),
        ("refs", ctypes.POINTER(func_reference_type)),
        ("flags", ctypes.c_int),
        ("dens_threshold", ctypes.c_double),
        ("n_ext_params", ctypes.c_int),
        ("ext_params", ctypes.POINTER(func_params_type)),
        ("set_ext_params", ctypes.c_void_p),
        ("init", ctypes.c_void_p),
        ("end", ctypes.c_void_p),
        ("lda", ctypes.c_void_p),
        ("gga", ctypes.c_void_p),
        ("mgga", ctypes.c_void_p),
    ]


class xc_dimensions(ctypes.Structure):
     """
    Holds dimensions of the several arrays.
    """
     _fields_ = [
         ("rho", ctypes.c_size_t),
         ("sigma", ctypes.c_size_t),
         ("lapl", ctypes.c_size_t),
         ("tau", ctypes.c_size_t),

         ("zk", ctypes.c_size_t),

         ("vrho", ctypes.c_size_t),
         ("vsigma", ctypes.c_size_t),
         ("vlapl", ctypes.c_size_t),
         ("vtau", ctypes.c_size_t), 

         ("v2rho2", ctypes.c_size_t),
         ("v2rhosigma", ctypes.c_size_t),
         ("v2rholapl", ctypes.c_size_t),
         ("v2rhotau", ctypes.c_size_t),
         ("v2sigma2", ctypes.c_size_t),
         ("v2sigmalapl", ctypes.c_size_t),
         ("v2sigmatau", ctypes.c_size_t),
         ("v2lapl2", ctypes.c_size_t),
         ("v2lapltau", ctypes.c_size_t),
         ("v2tau2", ctypes.c_size_t),

         ("v3rho3", ctypes.c_size_t),
         ("v3rho2sigma", ctypes.c_size_t),
         ("v3rho2lapl", ctypes.c_size_t),
         ("v3rho2tau", ctypes.c_size_t),
         ("v3rhosigma2", ctypes.c_size_t),
         ("v3rhosigmalapl", ctypes.c_size_t),
         ("v3rhosigmatau", ctypes.c_size_t),
         ("v3rholapl2", ctypes.c_size_t),
         ("v3rholapltau", ctypes.c_size_t),
         ("v3rhotau2", ctypes.c_size_t),
         ("v3sigma3", ctypes.c_size_t),
         ("v3sigma2lapl", ctypes.c_size_t),
         ("v3sigma2tau", ctypes.c_size_t),
         ("v3sigmalapl2", ctypes.c_size_t),
         ("v3sigmalapltau", ctypes.c_size_t),
         ("v3sigmatau2", ctypes.c_size_t),
         ("v3lapl3", ctypes.c_size_t),
         ("v3lapl2tau", ctypes.c_size_t),
         ("v3lapltau2", ctypes.c_size_t),
         ("v3tau3", ctypes.c_size_t),

         ("v4rho4", ctypes.c_size_t),
         ("v4rho3sigma", ctypes.c_size_t),
         ("v4rho3lapl", ctypes.c_size_t),
         ("v4rho3tau", ctypes.c_size_t),
         ("v4rho2sigma2", ctypes.c_size_t),
         ("v4rho2sigmalapl", ctypes.c_size_t),
         ("v4rho2sigmatau", ctypes.c_size_t),
         ("v4rho2lapl2", ctypes.c_size_t),
         ("v4rho2lapltau", ctypes.c_size_t),
         ("v4rho2tau2", ctypes.c_size_t),
         ("v4rhosigma3", ctypes.c_size_t),
         ("v4rhosigma2lapl", ctypes.c_size_t),
         ("v4rhosigma2tau", ctypes.c_size_t),
         ("v4rhosigmalapl2", ctypes.c_size_t),
         ("v4rhosigmalapltau", ctypes.c_size_t),
         ("v4rhosigmatau2", ctypes.c_size_t),
         ("v4rholapl3", ctypes.c_size_t),
         ("v4rholapl2tau", ctypes.c_size_t),
         ("v4rholapltau2", ctypes.c_size_t),
         ("v4rhotau3", ctypes.c_size_t),
         ("v4sigma4", ctypes.c_size_t),
         ("v4sigma3lapl", ctypes.c_size_t),
         ("v4sigma3tau", ctypes.c_size_t),
         ("v4sigma2lapl2", ctypes.c_size_t),
         ("v4sigma2lapltau", ctypes.c_size_t),
         ("v4sigma2tau2", ctypes.c_size_t),
         ("v4sigmalapl3", ctypes.c_size_t),
         ("v4sigmalapl2tau", ctypes.c_size_t),
         ("v4sigmalapltau2", ctypes.c_size_t),
         ("v4sigmatau3", ctypes.c_size_t),
         ("v4lapl4", ctypes.c_size_t),
         ("v4lapl3tau", ctypes.c_size_t),
         ("v4lapl2tau2", ctypes.c_size_t),
         ("v4lapltau3", ctypes.c_size_t),
         ("v4tau4", ctypes.c_size_t),
     ]
   
class xc_hybrid_params_type(ctypes.Structure):
    """
    Holds user defined parameters and their description.
    """
    # The data is stored in an union but it's easier to just use it as raw
    _fields_ = [("raw0", ctypes.c_double), ("raw1", ctypes.c_double), ("raw2", ctypes.c_double)]


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
        ("hyb_params", ctypes.POINTER(xc_hybrid_params_type)),

        # VV10
        ("nlc_b", ctypes.c_double),
        ("nlc_C", ctypes.c_double),

        ("dim", xc_dimensions),
        
        # parameters
        ("params", ctypes.c_void_p),  # void *params;
        ("dens_threshold", ctypes.c_double),
    ]
