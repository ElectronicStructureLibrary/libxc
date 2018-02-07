"""
Tests the LibXCFunctional class.
"""

import pytest

import pylibxc


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

    func = pylibxc.LibXCFunctional("XC_HYB_MGGA_XC_WB97M_V", 1)
    assert func.get_number() == 531
    assert func.get_kind() == 2
    assert func.get_name() == "wB97M-V exchange-correlation functional"
    assert func.get_family() == 64
    assert func.get_flags() == 1411

