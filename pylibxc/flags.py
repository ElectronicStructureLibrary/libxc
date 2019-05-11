"""
A port of the LibXC internal flags to Python. Header defines must be physically copied.
"""

XC_UNPOLARIZED = 1
XC_POLARIZED = 2

XC_NON_RELATIVISTIC = 0
XC_RELATIVISTIC = 1

XC_EXCHANGE = 0
XC_CORRELATION = 1
XC_EXCHANGE_CORRELATION = 2
XC_KINETIC = 3

XC_FAMILY_UNKNOWN = -1
XC_FAMILY_LDA = 1
XC_FAMILY_GGA = 2
XC_FAMILY_MGGA = 4
XC_FAMILY_LCA = 8
XC_FAMILY_OEP = 16
XC_FAMILY_HYB_GGA = 32
XC_FAMILY_HYB_MGGA = 64
XC_FAMILY_HYB_LDA = 128

XC_FLAGS_HAVE_EXC = (1 << 0)  # = 1
XC_FLAGS_HAVE_VXC = (1 << 1)  # = 2
XC_FLAGS_HAVE_FXC = (1 << 2)  # = 4
XC_FLAGS_HAVE_KXC = (1 << 3)  # = 8
XC_FLAGS_HAVE_LXC = (1 << 4)  # = 16
XC_FLAGS_1D = (1 << 5)  # = 32
XC_FLAGS_2D = (1 << 6)  # = 64
XC_FLAGS_3D = (1 << 7)  # = 128
XC_FLAGS_HYB_CAM = (1 << 8)  # = 256
XC_FLAGS_HYB_CAMY = (1 << 9)  # = 512
XC_FLAGS_VV10 = (1 << 10)  #  1024
XC_FLAGS_HYB_LC = (1 << 11)  #  2048
XC_FLAGS_HYB_LCY = (1 << 12)  #  4096
XC_FLAGS_STABLE = (1 << 13)  #  8192
XC_FLAGS_DEVELOPMENT = (1 << 14)  # 16384
XC_FLAGS_NEEDS_LAPLACIAN = (1 << 15)  # 32768

XC_TAU_EXPLICIT = 0
XC_TAU_EXPANSION = 1

XC_MAX_REFERENCES = 5

