"""
Tests the LibXCFunctional class.
"""

import pytest
import numpy as np

import pylibxc

compute_test_dim = 5

np.random.seed(0)

def _dict_array_comp(test, ref, keys):
    for key in keys:
        tmp = np.testing.assert_allclose(test[key], ref[key])
        # print(key, np.allclose(test[key], ref[key]), np.linalg.norm(test[key] - ref[key]))
    return True


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

    # Test polarized
    for polar, ndim in [("unpolarized", 1), ("polarized", 2)]:
        inp = {}
        inp["rho"] = np.random.random((compute_test_dim * ndim))

        func = pylibxc.LibXCFunctional("lda_c_vwn", polar)

        ret_full = func.compute(inp, do_exc=True, do_vxc=True, do_fxc=True, do_kxc=True)
        ret_ev = func.compute(inp, do_exc=True, do_vxc=True, do_fxc=False, do_kxc=False)
        ret_e = func.compute(inp, do_exc=True, do_vxc=False, do_fxc=False, do_kxc=False)
        ret_v = func.compute(inp, do_exc=False, do_vxc=True, do_fxc=False, do_kxc=False)
        ret_f = func.compute(inp, do_exc=False, do_vxc=False, do_fxc=True, do_kxc=False)
        ret_k = func.compute(inp, do_exc=False, do_vxc=False, do_fxc=False, do_kxc=True)

        assert ret_full["zk"].size == compute_test_dim
        assert ret_full["vrho"].size == compute_test_dim * ndim

        assert np.allclose(ret_full["zk"], ret_ev["zk"])
        assert np.allclose(ret_full["vrho"], ret_ev["vrho"])

        assert np.allclose(ret_full["zk"], ret_e["zk"])
        assert np.allclose(ret_full["vrho"], ret_v["vrho"])
        assert np.allclose(ret_full["v2rho2"], ret_f["v2rho2"])
        assert np.allclose(ret_full["v3rho3"], ret_k["v3rho3"])


def test_gga_compute():

    # Test polarized
    for polar, ndim in [("unpolarized", (1, 1)), ("polarized", (2, 3))]:
        inp = {}
        inp["rho"] = np.random.random((compute_test_dim * ndim[0]))
        inp["sigma"] = np.random.random((compute_test_dim * ndim[1]))

        func = pylibxc.LibXCFunctional("gga_c_pbe", polar)

        ret_full = func.compute(inp, do_exc=True, do_vxc=True, do_fxc=True, do_kxc=True)
        ret_ev = func.compute(inp, do_exc=True, do_vxc=True, do_fxc=False, do_kxc=False)
        ret_e = func.compute(inp, do_exc=True, do_vxc=False, do_fxc=False, do_kxc=False)
        ret_v = func.compute(inp, do_exc=False, do_vxc=True, do_fxc=False, do_kxc=False)
        ret_f = func.compute(inp, do_exc=False, do_vxc=False, do_fxc=True, do_kxc=False)
        ret_k = func.compute(inp, do_exc=False, do_vxc=False, do_fxc=False, do_kxc=True)

        assert ret_full["zk"].size == compute_test_dim
        assert ret_full["vrho"].size == compute_test_dim * ndim[0]
        assert ret_full["vsigma"].size == compute_test_dim * ndim[1]

        assert _dict_array_comp(ret_full, ret_ev, ["zk", "vrho", "vsigma"])

        assert _dict_array_comp(ret_full, ret_e, ["zk"])
        assert _dict_array_comp(ret_full, ret_v, ["vrho", "vsigma"])
        assert _dict_array_comp(ret_full, ret_f, ["v2rho2", "v2rhosigma", "v2sigma2"])
        assert _dict_array_comp(ret_full, ret_k, ["v3rho3", "v3rho2sigma", "v3rhosigma2", "v3sigma3"])
