"""
Tests the utility LibXC functions.
"""

import pylibxc

def test_xc_version():

    assert (4, 0, 1) == pylibxc.util.xc_version()

