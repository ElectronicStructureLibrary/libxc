"""
Tests the LibXCFunctional class.
"""

import pytest
import numpy as np

import pylibxc

compute_test_dim = 5

def test_libxc_functional_build():

    pylibxc.LibXCFunctional(1, 1)
    pylibxc.LibXCFunctional(1, 2)

    pylibxc.LibXCFunctional("XC_LDA_C_VWN", "polarized")
    pylibxc.LibXCFunctional("lda_c_vwn", "unpolarized")

    # Check functional edge cases
    with pytest.raises(KeyError):
        pylibxc.LibXCFunctional("something", 1)

    with pytest.raises(KeyError):
        pylibxc.LibXCFunctional(5000, 1)

    # Check spin edge cases
    with pytest.raises(KeyError):
        pylibxc.LibXCFunctional("lda_c_vwn", 10)

    with pytest.raises(KeyError):
        pylibxc.LibXCFunctional("lda_c_vwn", "something")


def test_libxc_functional_info():

    func = pylibxc.LibXCFunctional(1, 1)
    assert func.get_number() == 1
    assert func.get_kind() == 0
    assert func.get_name() == "Slater exchange"
    assert func.get_family() == 1
    assert func.get_flags() == 143
    assert len(func.get_bibtex()) == 2
    assert len(func.get_references()) == 2
    assert len(func.get_doi()) == 2

    func = pylibxc.LibXCFunctional("XC_HYB_MGGA_XC_WB97M_V", 1)
    assert func.get_number() == 531
    assert func.get_kind() == 2
    assert func.get_name() == "wB97M-V exchange-correlation functional"
    assert func.get_family() == 64
    assert func.get_flags() == 1411
    assert len(func.get_bibtex()) == 1
    assert len(func.get_references()) == 1
    assert len(func.get_doi()) == 1


def test_ext_params():

    func = pylibxc.LibXCFunctional(1, 1)
    assert 0 == len(func.get_ext_param_descriptions())
    assert 0 == len(func.get_ext_param_default_values())
    func.set_dens_threshold(1.e-16)
    func.set_dens_threshold(5)
    with pytest.raises(ValueError):
        func.set_ext_params([])

    func = pylibxc.LibXCFunctional("XC_HYB_GGA_XC_HSE06", 1)
    assert 3 == len(func.get_ext_param_descriptions())
    assert 3 == len(func.get_ext_param_default_values())
    assert all("param" in x for x in func.get_ext_param_descriptions())
    func.set_dens_threshold(1.e-16)
    func.set_dens_threshold(5)

    # Segfaults, need to check it out
    func.set_ext_params([5, 3, 3])

    with pytest.raises(ValueError):
        func.set_ext_params([5, 3])

    with pytest.raises(ValueError):
        func.set_dens_threshold(-1)

def test_lda_compute():
    inp = np.random.random((compute_test_dim)) 
    out = np.zeros((compute_test_dim)) 

    func = pylibxc.LibXCFunctional("lda_c_vwn", "unpolarized")
    func.compute(inp, out)
    
