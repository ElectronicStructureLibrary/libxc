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
        ("name", ctypes.c_char_p),
        ("kind", ctypes.c_int),
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
         ("rho", ctypes.c_int),
         ("sigma", ctypes.c_int),
         ("lapl", ctypes.c_int),
         ("tau", ctypes.c_int),

         ("zk", ctypes.c_int),

         ("vrho", ctypes.c_int),
         ("vsigma", ctypes.c_int),
         ("vlapl", ctypes.c_int),
         ("vtau", ctypes.c_int), 

         ("v2rho2", ctypes.c_int),
         ("v2rhosigma", ctypes.c_int),
         ("v2rholapl", ctypes.c_int),
         ("v2rhotau", ctypes.c_int),
         ("v2sigma2", ctypes.c_int),
         ("v2sigmalapl", ctypes.c_int),
         ("v2sigmatau", ctypes.c_int),
         ("v2lapl2", ctypes.c_int),
         ("v2lapltau", ctypes.c_int),
         ("v2tau2", ctypes.c_int),

         ("v3rho3", ctypes.c_int),
         ("v3rho2sigma", ctypes.c_int),
         ("v3rho2lapl", ctypes.c_int),
         ("v3rho2tau", ctypes.c_int),
         ("v3rhosigma2", ctypes.c_int),
         ("v3rhosigmalapl", ctypes.c_int),
         ("v3rhosigmatau", ctypes.c_int),
         ("v3rholapl2", ctypes.c_int),
         ("v3rholapltau", ctypes.c_int),
         ("v3rhotau2", ctypes.c_int),
         ("v3sigma3", ctypes.c_int),
         ("v3sigma2lapl", ctypes.c_int),
         ("v3sigma2tau", ctypes.c_int),
         ("v3sigmalapl2", ctypes.c_int),
         ("v3sigmalapltau", ctypes.c_int),
         ("v3sigmatau2", ctypes.c_int),
         ("v3lapl3", ctypes.c_int),
         ("v3lapl2tau", ctypes.c_int),
         ("v3lapltau2", ctypes.c_int),
         ("v3tau3", ctypes.c_int),
     ]
   

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

        # CAM
        ("cam_omega", ctypes.c_double),
        ("cam_alpha", ctypes.c_double),
        ("cam_beta", ctypes.c_double),

        # VV10
        ("nlc_b", ctypes.c_double),
        ("nlc_C", ctypes.c_double),

        ("dim", xc_dimensions),
        
        # parameters
        ("params", ctypes.c_void_p),  # void *params;
        ("dens_threshold", ctypes.c_double),
    ]
