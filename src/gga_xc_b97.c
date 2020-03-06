/*
 Copyright (C) 2006-2007 M.A.L. Marques
               2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_GGA_XC_HCTH_93     161 /* HCTH functional fitted to  93 molecules  */
#define XC_GGA_XC_HCTH_120    162 /* HCTH functional fitted to 120 molecules  */
#define XC_GGA_XC_HCTH_147    163 /* HCTH functional fitted to 147 molecules  */
#define XC_GGA_XC_HCTH_407    164 /* HCTH functional fitted to 407 molecules  */
#define XC_HYB_GGA_XC_B97     407 /* Becke 97                                 */
#define XC_HYB_GGA_XC_B97_1   408 /* Becke 97-1                               */
#define XC_HYB_GGA_XC_B97_2   410 /* Becke 97-2                               */
#define XC_GGA_XC_B97_D       170 /* Grimme functional to be used with C6 vdW term */
#define XC_HYB_GGA_XC_B97_K   413 /* Boese-Martin for Kinetics                */
#define XC_HYB_GGA_XC_B97_3   414 /* Becke 97-3                               */
#define XC_HYB_GGA_XC_SB98_1A 420 /* Schmider-Becke 98 parameterization 1a    */
#define XC_HYB_GGA_XC_SB98_1B 421 /* Schmider-Becke 98 parameterization 1b    */
#define XC_HYB_GGA_XC_SB98_1C 422 /* Schmider-Becke 98 parameterization 1c    */
#define XC_HYB_GGA_XC_SB98_2A 423 /* Schmider-Becke 98 parameterization 2a    */
#define XC_HYB_GGA_XC_SB98_2B 424 /* Schmider-Becke 98 parameterization 2b    */
#define XC_HYB_GGA_XC_SB98_2C 425 /* Schmider-Becke 98 parameterization 2c    */
#define XC_GGA_XC_B97_GGA1     96 /* Becke 97 GGA-1                           */
#define XC_GGA_XC_HCTH_P14     95 /* HCTH p=1/4                               */
#define XC_GGA_XC_HCTH_P76     94 /* HCTH p=7/6                               */
#define XC_GGA_XC_HCTH_407P    93 /* HCTH/407+                                */
#define XC_HYB_GGA_XC_B97_1P  266 /* version of B97 by Cohen and Handy        */
#define XC_GGA_XC_HLE16       545 /* high local exchange 2016                 */

typedef struct {
  double c_x[5], c_ss[5], c_ab[5];
} gga_xc_b97_params;

static const func_params_type ext_params_par_hcth_93[] = {
  { "_cx0",     1.0932, "u^0 coefficient for exchange"},
  { "_cx1",  -0.744056, "u^1 coefficient for exchange"},
  { "_cx2",     5.5992, "u^2 coefficient for exchange"},
  { "_cx3",   -6.78549, "u^3 coefficient for exchange"},
  { "_cx4",    4.49357, "u^4 coefficient for exchange"},
  {"_css0",   0.222601, "u^0 coefficient for same-spin correlation"},
  {"_css1", -0.0338622, "u^1 coefficient for same-spin correlation"},
  {"_css2",  -0.012517, "u^2 coefficient for same-spin correlation"},
  {"_css3",  -0.802496, "u^3 coefficient for same-spin correlation"},
  {"_css4",    1.55396, "u^4 coefficient for same-spin correlation"},
  {"_cos0",   0.729974, "u^0 coefficient for opposite-spin correlation"},
  {"_cos1",    3.35287, "u^1 coefficient for opposite-spin correlation"},
  {"_cos2",    -11.543, "u^2 coefficient for opposite-spin correlation"},
  {"_cos3",    8.08564, "u^3 coefficient for opposite-spin correlation"},
  {"_cos4",   -4.47857, "u^4 coefficient for opposite-spin correlation"},
};

static const func_params_type ext_params_par_hcth_120[] = {
  { "_cx0",    1.09163, "u^0 coefficient for exchange"},
  { "_cx1",  -0.747215, "u^1 coefficient for exchange"},
  { "_cx2",    5.07833, "u^2 coefficient for exchange"},
  { "_cx3",   -4.10746, "u^3 coefficient for exchange"},
  { "_cx4",    1.17173, "u^4 coefficient for exchange"},
  {"_css0",   0.489508, "u^0 coefficient for same-spin correlation"},
  {"_css1",  -0.260699, "u^1 coefficient for same-spin correlation"},
  {"_css2",   0.432917, "u^2 coefficient for same-spin correlation"},
  {"_css3",   -1.99247, "u^3 coefficient for same-spin correlation"},
  {"_css4",    2.48531, "u^4 coefficient for same-spin correlation"},
  {"_cos0",    0.51473, "u^0 coefficient for opposite-spin correlation"},
  {"_cos1",    6.92982, "u^1 coefficient for opposite-spin correlation"},
  {"_cos2",   -24.7073, "u^2 coefficient for opposite-spin correlation"},
  {"_cos3",    23.1098, "u^3 coefficient for opposite-spin correlation"},
  {"_cos4",   -11.3234, "u^4 coefficient for opposite-spin correlation"},
};

static const func_params_type ext_params_par_hcth_147[] = {
  { "_cx0",    1.09025, "u^0 coefficient for exchange"},
  { "_cx1",  -0.799194, "u^1 coefficient for exchange"},
  { "_cx2",    5.57212, "u^2 coefficient for exchange"},
  { "_cx3",    -5.8676, "u^3 coefficient for exchange"},
  { "_cx4",    3.04544, "u^4 coefficient for exchange"},
  {"_css0",   0.562576, "u^0 coefficient for same-spin correlation"},
  {"_css1",  0.0171436, "u^1 coefficient for same-spin correlation"},
  {"_css2",   -1.30636, "u^2 coefficient for same-spin correlation"},
  {"_css3",    1.05747, "u^3 coefficient for same-spin correlation"},
  {"_css4",   0.885429, "u^4 coefficient for same-spin correlation"},
  {"_cos0",   0.542352, "u^0 coefficient for opposite-spin correlation"},
  {"_cos1",    7.01464, "u^1 coefficient for opposite-spin correlation"},
  {"_cos2",   -28.3822, "u^2 coefficient for opposite-spin correlation"},
  {"_cos3",    35.0329, "u^3 coefficient for opposite-spin correlation"},
  {"_cos4",   -20.4284, "u^4 coefficient for opposite-spin correlation"},
};

static const func_params_type ext_params_par_hcth_407[] = {
  { "_cx0",    1.08184, "u^0 coefficient for exchange"},
  { "_cx1",  -0.518339, "u^1 coefficient for exchange"},
  { "_cx2",    3.42562, "u^2 coefficient for exchange"},
  { "_cx3",   -2.62901, "u^3 coefficient for exchange"},
  { "_cx4",    2.28855, "u^4 coefficient for exchange"},
  {"_css0",    1.18777, "u^0 coefficient for same-spin correlation"},
  {"_css1",   -2.40292, "u^1 coefficient for same-spin correlation"},
  {"_css2",    5.61741, "u^2 coefficient for same-spin correlation"},
  {"_css3",   -9.17923, "u^3 coefficient for same-spin correlation"},
  {"_css4",    6.24798, "u^4 coefficient for same-spin correlation"},
  {"_cos0",   0.589076, "u^0 coefficient for opposite-spin correlation"},
  {"_cos1",    4.42374, "u^1 coefficient for opposite-spin correlation"},
  {"_cos2",   -19.2218, "u^2 coefficient for opposite-spin correlation"},
  {"_cos3",    42.5721, "u^3 coefficient for opposite-spin correlation"},
  {"_cos4",   -42.0052, "u^4 coefficient for opposite-spin correlation"},
};

static const func_params_type ext_params_par_b97[] = {
  { "_cx0",     0.8094, "u^0 coefficient for exchange"},
  { "_cx1",     0.5073, "u^1 coefficient for exchange"},
  { "_cx2",     0.7481, "u^2 coefficient for exchange"},
  { "_cx3",        0.0, "u^3 coefficient for exchange"},
  { "_cx4",        0.0, "u^4 coefficient for exchange"},
  {"_css0",     0.1737, "u^0 coefficient for same-spin correlation"},
  {"_css1",     2.3487, "u^1 coefficient for same-spin correlation"},
  {"_css2",    -2.4868, "u^2 coefficient for same-spin correlation"},
  {"_css3",        0.0, "u^3 coefficient for same-spin correlation"},
  {"_css4",        0.0, "u^4 coefficient for same-spin correlation"},
  {"_cos0",     0.9454, "u^0 coefficient for opposite-spin correlation"},
  {"_cos1",     0.7471, "u^1 coefficient for opposite-spin correlation"},
  {"_cos2",    -4.5961, "u^2 coefficient for opposite-spin correlation"},
  {"_cos3",        0.0, "u^3 coefficient for opposite-spin correlation"},
  {"_cos4",        0.0, "u^4 coefficient for opposite-spin correlation"},
  {"_cxx",      0.1943, "coefficient for exact exchange"}
};

static const func_params_type ext_params_par_b97_1[] = {
  { "_cx0",   0.789518, "u^0 coefficient for exchange"},
  { "_cx1",   0.573805, "u^1 coefficient for exchange"},
  { "_cx2",   0.660975, "u^2 coefficient for exchange"},
  { "_cx3",        0.0, "u^3 coefficient for exchange"},
  { "_cx4",        0.0, "u^4 coefficient for exchange"},
  {"_css0",  0.0820011, "u^0 coefficient for same-spin correlation"},
  {"_css1",    2.71681, "u^1 coefficient for same-spin correlation"},
  {"_css2",   -2.87103, "u^2 coefficient for same-spin correlation"},
  {"_css3",        0.0, "u^3 coefficient for same-spin correlation"},
  {"_css4",        0.0, "u^4 coefficient for same-spin correlation"},
  {"_cos0",   0.955689, "u^0 coefficient for opposite-spin correlation"},
  {"_cos1",   0.788552, "u^1 coefficient for opposite-spin correlation"},
  {"_cos2",   -5.47869, "u^2 coefficient for opposite-spin correlation"},
  {"_cos3",        0.0, "u^3 coefficient for opposite-spin correlation"},
  {"_cos4",        0.0, "u^4 coefficient for opposite-spin correlation"},
  {"_cxx",        0.21, "coefficient for exact exchange"}
};

static const func_params_type ext_params_par_b97_2[] = {
  { "_cx0",   0.827642, "u^0 coefficient for exchange"},
  { "_cx1",    0.04784, "u^1 coefficient for exchange"},
  { "_cx2",    1.76125, "u^2 coefficient for exchange"},
  { "_cx3",        0.0, "u^3 coefficient for exchange"},
  { "_cx4",        0.0, "u^4 coefficient for exchange"},
  {"_css0",   0.585808, "u^0 coefficient for same-spin correlation"},
  {"_css1",  -0.691682, "u^1 coefficient for same-spin correlation"},
  {"_css2",   0.394796, "u^2 coefficient for same-spin correlation"},
  {"_css3",        0.0, "u^3 coefficient for same-spin correlation"},
  {"_css4",        0.0, "u^4 coefficient for same-spin correlation"},
  {"_cos0",   0.999849, "u^0 coefficient for opposite-spin correlation"},
  {"_cos1",    1.40626, "u^1 coefficient for opposite-spin correlation"},
  {"_cos2",    -7.4406, "u^2 coefficient for opposite-spin correlation"},
  {"_cos3",        0.0, "u^3 coefficient for opposite-spin correlation"},
  {"_cos4",        0.0, "u^4 coefficient for opposite-spin correlation"},
  {"_cxx",        0.21, "coefficient for exact exchange"}
};

static const func_params_type ext_params_par_b97_d[] = {
  { "_cx0",    1.08662, "u^0 coefficient for exchange"},
  { "_cx1",   -0.52127, "u^1 coefficient for exchange"},
  { "_cx2",    3.25429, "u^2 coefficient for exchange"},
  { "_cx3",        0.0, "u^3 coefficient for exchange"},
  { "_cx4",        0.0, "u^4 coefficient for exchange"},
  {"_css0",     0.2234, "u^0 coefficient for same-spin correlation"},
  {"_css1",   -1.56208, "u^1 coefficient for same-spin correlation"},
  {"_css2",    1.94293, "u^2 coefficient for same-spin correlation"},
  {"_css3",        0.0, "u^3 coefficient for same-spin correlation"},
  {"_css4",        0.0, "u^4 coefficient for same-spin correlation"},
  {"_cos0",    0.69041, "u^0 coefficient for opposite-spin correlation"},
  {"_cos1",     6.3027, "u^1 coefficient for opposite-spin correlation"},
  {"_cos2",   -14.9712, "u^2 coefficient for opposite-spin correlation"},
  {"_cos3",        0.0, "u^3 coefficient for opposite-spin correlation"},
  {"_cos4",        0.0, "u^4 coefficient for opposite-spin correlation"},
};

static const func_params_type ext_params_par_b97_k[] = {
  { "_cx0",   0.507863, "u^0 coefficient for exchange"},
  { "_cx1",    1.46873, "u^1 coefficient for exchange"},
  { "_cx2",   -1.51301, "u^2 coefficient for exchange"},
  { "_cx3",        0.0, "u^3 coefficient for exchange"},
  { "_cx4",        0.0, "u^4 coefficient for exchange"},
  {"_css0",    0.12355, "u^0 coefficient for same-spin correlation"},
  {"_css1",    2.65399, "u^1 coefficient for same-spin correlation"},
  {"_css2",   -3.20694, "u^2 coefficient for same-spin correlation"},
  {"_css3",        0.0, "u^3 coefficient for same-spin correlation"},
  {"_css4",        0.0, "u^4 coefficient for same-spin correlation"},
  {"_cos0",    1.58613, "u^0 coefficient for opposite-spin correlation"},
  {"_cos1",   -6.20977, "u^1 coefficient for opposite-spin correlation"},
  {"_cos2",    6.46106, "u^2 coefficient for opposite-spin correlation"},
  {"_cos3",        0.0, "u^3 coefficient for opposite-spin correlation"},
  {"_cos4",        0.0, "u^4 coefficient for opposite-spin correlation"},
  {"_cxx",        0.42, "coefficient for exact exchange"}
};

static const func_params_type ext_params_par_b97_3[] = {
  { "_cx0",  0.7334648, "u^0 coefficient for exchange"},
  { "_cx1",   0.292527, "u^1 coefficient for exchange"},
  { "_cx2",   3.338789, "u^2 coefficient for exchange"},
  { "_cx3",  -10.51158, "u^3 coefficient for exchange"},
  { "_cx4",   10.60907, "u^4 coefficient for exchange"},
  {"_css0",  0.5623649, "u^0 coefficient for same-spin correlation"},
  {"_css1",   -1.32298, "u^1 coefficient for same-spin correlation"},
  {"_css2",   6.359191, "u^2 coefficient for same-spin correlation"},
  {"_css3",  -7.464002, "u^3 coefficient for same-spin correlation"},
  {"_css4",   1.827082, "u^4 coefficient for same-spin correlation"},
  {"_cos0",    1.13383, "u^0 coefficient for opposite-spin correlation"},
  {"_cos1",  -2.811967, "u^1 coefficient for opposite-spin correlation"},
  {"_cos2",   7.431302, "u^2 coefficient for opposite-spin correlation"},
  {"_cos3",  -1.969342, "u^3 coefficient for opposite-spin correlation"},
  {"_cos4",  -11.74423, "u^4 coefficient for opposite-spin correlation"},
  {"_cxx", 2.692880E-01, "coefficient for exact exchange"}
};

static const func_params_type ext_params_par_sb98_1a[] = {
  { "_cx0",   0.845975, "u^0 coefficient for exchange"},
  { "_cx1",   0.228183, "u^1 coefficient for exchange"},
  { "_cx2",   0.749949, "u^2 coefficient for exchange"},
  { "_cx3",        0.0, "u^3 coefficient for exchange"},
  { "_cx4",        0.0, "u^4 coefficient for exchange"},
  {"_css0",  -0.817637, "u^0 coefficient for same-spin correlation"},
  {"_css1",  -0.054676, "u^1 coefficient for same-spin correlation"},
  {"_css2",   0.592163, "u^2 coefficient for same-spin correlation"},
  {"_css3",        0.0, "u^3 coefficient for same-spin correlation"},
  {"_css4",        0.0, "u^4 coefficient for same-spin correlation"},
  {"_cos0",   0.975483, "u^0 coefficient for opposite-spin correlation"},
  {"_cos1",   0.398379, "u^1 coefficient for opposite-spin correlation"},
  {"_cos2",    -3.7354, "u^2 coefficient for opposite-spin correlation"},
  {"_cos3",        0.0, "u^3 coefficient for opposite-spin correlation"},
  {"_cos4",        0.0, "u^4 coefficient for opposite-spin correlation"},
  {"_cxx",    0.229015, "coefficient for exact exchange"}
};

static const func_params_type ext_params_par_sb98_1b[] = {
  { "_cx0",   0.800103, "u^0 coefficient for exchange"},
  { "_cx1",  -0.084192, "u^1 coefficient for exchange"},
  { "_cx2",    1.47742, "u^2 coefficient for exchange"},
  { "_cx3",        0.0, "u^3 coefficient for exchange"},
  { "_cx4",        0.0, "u^4 coefficient for exchange"},
  {"_css0",    1.44946, "u^0 coefficient for same-spin correlation"},
  {"_css1",   -2.37073, "u^1 coefficient for same-spin correlation"},
  {"_css2",    2.13564, "u^2 coefficient for same-spin correlation"},
  {"_css3",        0.0, "u^3 coefficient for same-spin correlation"},
  {"_css4",        0.0, "u^4 coefficient for same-spin correlation"},
  {"_cos0",   0.977621, "u^0 coefficient for opposite-spin correlation"},
  {"_cos1",   0.931199, "u^1 coefficient for opposite-spin correlation"},
  {"_cos2",   -4.76973, "u^2 coefficient for opposite-spin correlation"},
  {"_cos3",        0.0, "u^3 coefficient for opposite-spin correlation"},
  {"_cos4",        0.0, "u^4 coefficient for opposite-spin correlation"},
  {"_cxx",    0.199352, "coefficient for exact exchange"}
};

static const func_params_type ext_params_par_sb98_1c[] = {
  { "_cx0",   0.810936, "u^0 coefficient for exchange"},
  { "_cx1",    0.49609, "u^1 coefficient for exchange"},
  { "_cx2",   0.772385, "u^2 coefficient for exchange"},
  { "_cx3",        0.0, "u^3 coefficient for exchange"},
  { "_cx4",        0.0, "u^4 coefficient for exchange"},
  {"_css0",   0.262077, "u^0 coefficient for same-spin correlation"},
  {"_css1",    2.12576, "u^1 coefficient for same-spin correlation"},
  {"_css2",   -2.30465, "u^2 coefficient for same-spin correlation"},
  {"_css3",        0.0, "u^3 coefficient for same-spin correlation"},
  {"_css4",        0.0, "u^4 coefficient for same-spin correlation"},
  {"_cos0",   0.939269, "u^0 coefficient for opposite-spin correlation"},
  {"_cos1",   0.898121, "u^1 coefficient for opposite-spin correlation"},
  {"_cos2",   -4.91276, "u^2 coefficient for opposite-spin correlation"},
  {"_cos3",        0.0, "u^3 coefficient for opposite-spin correlation"},
  {"_cos4",        0.0, "u^4 coefficient for opposite-spin correlation"},
  {"_cxx",    0.192416, "coefficient for exact exchange"}
};

static const func_params_type ext_params_par_sb98_2a[] = {
  { "_cx0",     0.7492, "u^0 coefficient for exchange"},
  { "_cx1",   0.402322, "u^1 coefficient for exchange"},
  { "_cx2",   0.620779, "u^2 coefficient for exchange"},
  { "_cx3",        0.0, "u^3 coefficient for exchange"},
  { "_cx4",        0.0, "u^4 coefficient for exchange"},
  {"_css0",    1.26686, "u^0 coefficient for same-spin correlation"},
  {"_css1",    1.67146, "u^1 coefficient for same-spin correlation"},
  {"_css2",   -1.22565, "u^2 coefficient for same-spin correlation"},
  {"_css3",        0.0, "u^3 coefficient for same-spin correlation"},
  {"_css4",        0.0, "u^4 coefficient for same-spin correlation"},
  {"_cos0",   0.964641, "u^0 coefficient for opposite-spin correlation"},
  {"_cos1",   0.050527, "u^1 coefficient for opposite-spin correlation"},
  {"_cos2",   -3.01966, "u^2 coefficient for opposite-spin correlation"},
  {"_cos3",        0.0, "u^3 coefficient for opposite-spin correlation"},
  {"_cos4",        0.0, "u^4 coefficient for opposite-spin correlation"},
  {"_cxx",    0.232055, "coefficient for exact exchange"}
};

static const func_params_type ext_params_par_sb98_2b[] = {
  { "_cx0",   0.770587, "u^0 coefficient for exchange"},
  { "_cx1",   0.180767, "u^1 coefficient for exchange"},
  { "_cx2",   0.955246, "u^2 coefficient for exchange"},
  { "_cx3",        0.0, "u^3 coefficient for exchange"},
  { "_cx4",        0.0, "u^4 coefficient for exchange"},
  {"_css0",   0.170473, "u^0 coefficient for same-spin correlation"},
  {"_css1",    1.24051, "u^1 coefficient for same-spin correlation"},
  {"_css2",  -0.862711, "u^2 coefficient for same-spin correlation"},
  {"_css3",        0.0, "u^3 coefficient for same-spin correlation"},
  {"_css4",        0.0, "u^4 coefficient for same-spin correlation"},
  {"_cos0",   0.965362, "u^0 coefficient for opposite-spin correlation"},
  {"_cos1",     0.8633, "u^1 coefficient for opposite-spin correlation"},
  {"_cos2",   -4.61778, "u^2 coefficient for opposite-spin correlation"},
  {"_cos3",        0.0, "u^3 coefficient for opposite-spin correlation"},
  {"_cos4",        0.0, "u^4 coefficient for opposite-spin correlation"},
  {"_cxx",    0.237978, "coefficient for exact exchange"}
};

static const func_params_type ext_params_par_sb98_2c[] = {
  { "_cx0",   0.790194, "u^0 coefficient for exchange"},
  { "_cx1",   0.400271, "u^1 coefficient for exchange"},
  { "_cx2",   0.832857, "u^2 coefficient for exchange"},
  { "_cx3",        0.0, "u^3 coefficient for exchange"},
  { "_cx4",        0.0, "u^4 coefficient for exchange"},
  {"_css0",  -0.120163, "u^0 coefficient for same-spin correlation"},
  {"_css1",    2.82332, "u^1 coefficient for same-spin correlation"},
  {"_css2",   -2.59412, "u^2 coefficient for same-spin correlation"},
  {"_css3",        0.0, "u^3 coefficient for same-spin correlation"},
  {"_css4",        0.0, "u^4 coefficient for same-spin correlation"},
  {"_cos0",   0.934715, "u^0 coefficient for opposite-spin correlation"},
  {"_cos1",    1.14105, "u^1 coefficient for opposite-spin correlation"},
  {"_cos2",   -5.33398, "u^2 coefficient for opposite-spin correlation"},
  {"_cos3",        0.0, "u^3 coefficient for opposite-spin correlation"},
  {"_cos4",        0.0, "u^4 coefficient for opposite-spin correlation"},
  {"_cxx",    0.219847, "coefficient for exact exchange"}
};

static const func_params_type ext_params_par_b97_gga1[] = {
  { "_cx0",     1.1068, "u^0 coefficient for exchange"},
  { "_cx1",    -0.8765, "u^1 coefficient for exchange"},
  { "_cx2",     4.2639, "u^2 coefficient for exchange"},
  { "_cx3",        0.0, "u^3 coefficient for exchange"},
  { "_cx4",        0.0, "u^4 coefficient for exchange"},
  {"_css0",     0.4883, "u^0 coefficient for same-spin correlation"},
  {"_css1",     -2.117, "u^1 coefficient for same-spin correlation"},
  {"_css2",     2.3235, "u^2 coefficient for same-spin correlation"},
  {"_css3",        0.0, "u^3 coefficient for same-spin correlation"},
  {"_css4",        0.0, "u^4 coefficient for same-spin correlation"},
  {"_cos0",     0.7961, "u^0 coefficient for opposite-spin correlation"},
  {"_cos1",      5.706, "u^1 coefficient for opposite-spin correlation"},
  {"_cos2",    -14.982, "u^2 coefficient for opposite-spin correlation"},
  {"_cos3",        0.0, "u^3 coefficient for opposite-spin correlation"},
  {"_cos4",        0.0, "u^4 coefficient for opposite-spin correlation"},
};

static const func_params_type ext_params_par_hcth_p14[] = {
  { "_cx0",    1.03161, "u^0 coefficient for exchange"},
  { "_cx1",  -0.360781, "u^1 coefficient for exchange"},
  { "_cx2",    3.51994, "u^2 coefficient for exchange"},
  { "_cx3",   -4.95944, "u^3 coefficient for exchange"},
  { "_cx4",    2.41165, "u^4 coefficient for exchange"},
  {"_css0",    2.82414, "u^0 coefficient for same-spin correlation"},
  {"_css1",  0.0318843, "u^1 coefficient for same-spin correlation"},
  {"_css2",   -1.78512, "u^2 coefficient for same-spin correlation"},
  {"_css3",    2.39795, "u^3 coefficient for same-spin correlation"},
  {"_css4",  -0.876909, "u^4 coefficient for same-spin correlation"},
  {"_cos0",  0.0821827, "u^0 coefficient for opposite-spin correlation"},
  {"_cos1",    4.56466, "u^1 coefficient for opposite-spin correlation"},
  {"_cos2",   -13.5529, "u^2 coefficient for opposite-spin correlation"},
  {"_cos3",     13.382, "u^3 coefficient for opposite-spin correlation"},
  {"_cos4",   -3.17493, "u^4 coefficient for opposite-spin correlation"},
};

static const func_params_type ext_params_par_hcth_p76[] = {
  { "_cx0",    1.16525, "u^0 coefficient for exchange"},
  { "_cx1",  -0.583033, "u^1 coefficient for exchange"},
  { "_cx2",    2.51769, "u^2 coefficient for exchange"},
  { "_cx3",    3.81278, "u^3 coefficient for exchange"},
  { "_cx4",   -5.45906, "u^4 coefficient for exchange"},
  {"_css0",   -3.92143, "u^0 coefficient for same-spin correlation"},
  {"_css1",   -1.10098, "u^1 coefficient for same-spin correlation"},
  {"_css2",  -0.091405, "u^2 coefficient for same-spin correlation"},
  {"_css3",  -0.859723, "u^3 coefficient for same-spin correlation"},
  {"_css4",    2.07184, "u^4 coefficient for same-spin correlation"},
  {"_cos0",   0.192949, "u^0 coefficient for opposite-spin correlation"},
  {"_cos1",   -5.73335, "u^1 coefficient for opposite-spin correlation"},
  {"_cos2",    50.8757, "u^2 coefficient for opposite-spin correlation"},
  {"_cos3",    135.475, "u^3 coefficient for opposite-spin correlation"},
  {"_cos4",    101.268, "u^4 coefficient for opposite-spin correlation"},
};

static const func_params_type ext_params_par_hcth_407p[] = {
  { "_cx0",    1.08018, "u^0 coefficient for exchange"},
  { "_cx1",    -0.4117, "u^1 coefficient for exchange"},
  { "_cx2",     2.4368, "u^2 coefficient for exchange"},
  { "_cx3",      1.389, "u^3 coefficient for exchange"},
  { "_cx4",    -1.3529, "u^4 coefficient for exchange"},
  {"_css0",    0.80302, "u^0 coefficient for same-spin correlation"},
  {"_css1",    -1.0479, "u^1 coefficient for same-spin correlation"},
  {"_css2",     4.9807, "u^2 coefficient for same-spin correlation"},
  {"_css3",     -12.89, "u^3 coefficient for same-spin correlation"},
  {"_css4",     9.6446, "u^4 coefficient for same-spin correlation"},
  {"_cos0",    0.73604, "u^0 coefficient for opposite-spin correlation"},
  {"_cos1",      3.027, "u^1 coefficient for opposite-spin correlation"},
  {"_cos2",    -10.075, "u^2 coefficient for opposite-spin correlation"},
  {"_cos3",     20.611, "u^3 coefficient for opposite-spin correlation"},
  {"_cos4",    -29.418, "u^4 coefficient for opposite-spin correlation"},
};

static const func_params_type ext_params_par_b97_1p[] = {
  { "_cx0",     0.8773, "u^0 coefficient for exchange"},
  { "_cx1",     0.2149, "u^1 coefficient for exchange"},
  { "_cx2",     1.5204, "u^2 coefficient for exchange"},
  { "_cx3",        0.0, "u^3 coefficient for exchange"},
  { "_cx4",        0.0, "u^4 coefficient for exchange"},
  {"_css0",     0.2228, "u^0 coefficient for same-spin correlation"},
  {"_css1",     1.3678, "u^1 coefficient for same-spin correlation"},
  {"_css2",    -1.5068, "u^2 coefficient for same-spin correlation"},
  {"_css3",        0.0, "u^3 coefficient for same-spin correlation"},
  {"_css4",        0.0, "u^4 coefficient for same-spin correlation"},
  {"_cos0",     0.9253, "u^0 coefficient for opposite-spin correlation"},
  {"_cos1",      2.027, "u^1 coefficient for opposite-spin correlation"},
  {"_cos2",    -7.3431, "u^2 coefficient for opposite-spin correlation"},
  {"_cos3",        0.0, "u^3 coefficient for opposite-spin correlation"},
  {"_cos4",        0.0, "u^4 coefficient for opposite-spin correlation"},
  {"_cxx",        0.15, "coefficient for exact exchange"}
};

static const func_params_type ext_params_par_hle16[] = {
  { "_cx0",     1.3523, "u^0 coefficient for exchange"},
  { "_cx1", -0.64792375, "u^1 coefficient for exchange"},
  { "_cx2",   4.282025, "u^2 coefficient for exchange"},
  { "_cx3", -3.2862625, "u^3 coefficient for exchange"},
  { "_cx4",  2.8606875, "u^4 coefficient for exchange"},
  {"_css0",   0.593885, "u^0 coefficient for same-spin correlation"},
  {"_css1",   -1.20146, "u^1 coefficient for same-spin correlation"},
  {"_css2",   2.808705, "u^2 coefficient for same-spin correlation"},
  {"_css3",  -4.589615, "u^3 coefficient for same-spin correlation"},
  {"_css4",    3.12399, "u^4 coefficient for same-spin correlation"},
  {"_cos0",   0.294538, "u^0 coefficient for opposite-spin correlation"},
  {"_cos1",    2.21187, "u^1 coefficient for opposite-spin correlation"},
  {"_cos2",    -9.6109, "u^2 coefficient for opposite-spin correlation"},
  {"_cos3",   21.28605, "u^3 coefficient for opposite-spin correlation"},
  {"_cos4",   -21.0026, "u^4 coefficient for opposite-spin correlation"},
};


static void
gga_xc_b97_init(xc_func_type *p)
{
  gga_xc_b97_params *params;

  assert(p->params == NULL);
  p->params = libxc_malloc(sizeof(gga_xc_b97_params));
  params = (gga_xc_b97_params *)(p->params);
}

static void
set_ext_params_pure(xc_func_type *p, const double *ext_params)
{
  gga_xc_b97_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_xc_b97_params *) (p->params);

  params->c_x[0] = get_ext_param(p->info->ext_params, ext_params, 0);
  params->c_x[1] = get_ext_param(p->info->ext_params, ext_params, 1);
  params->c_x[2] = get_ext_param(p->info->ext_params, ext_params, 2);
  params->c_x[3] = get_ext_param(p->info->ext_params, ext_params, 3);
  params->c_x[4] = get_ext_param(p->info->ext_params, ext_params, 4);
  params->c_ss[0] = get_ext_param(p->info->ext_params, ext_params, 5);
  params->c_ss[1] = get_ext_param(p->info->ext_params, ext_params, 6);
  params->c_ss[2] = get_ext_param(p->info->ext_params, ext_params, 7);
  params->c_ss[3] = get_ext_param(p->info->ext_params, ext_params, 8);
  params->c_ss[4] = get_ext_param(p->info->ext_params, ext_params, 9);
  params->c_ab[0] = get_ext_param(p->info->ext_params, ext_params, 10);
  params->c_ab[1] = get_ext_param(p->info->ext_params, ext_params, 11);
  params->c_ab[2] = get_ext_param(p->info->ext_params, ext_params, 12);
  params->c_ab[3] = get_ext_param(p->info->ext_params, ext_params, 13);
  params->c_ab[4] = get_ext_param(p->info->ext_params, ext_params, 14);
}

static void
set_ext_params_hybrid(xc_func_type *p, const double *ext_params)
{
  gga_xc_b97_params *params;

  assert(p != NULL && p->params != NULL);
  params = (gga_xc_b97_params *) (p->params);

  params->c_x[0] = get_ext_param(p->info->ext_params, ext_params, 0);
  params->c_x[1] = get_ext_param(p->info->ext_params, ext_params, 1);
  params->c_x[2] = get_ext_param(p->info->ext_params, ext_params, 2);
  params->c_x[3] = get_ext_param(p->info->ext_params, ext_params, 3);
  params->c_x[4] = get_ext_param(p->info->ext_params, ext_params, 4);
  params->c_ss[0] = get_ext_param(p->info->ext_params, ext_params, 5);
  params->c_ss[1] = get_ext_param(p->info->ext_params, ext_params, 6);
  params->c_ss[2] = get_ext_param(p->info->ext_params, ext_params, 7);
  params->c_ss[3] = get_ext_param(p->info->ext_params, ext_params, 8);
  params->c_ss[4] = get_ext_param(p->info->ext_params, ext_params, 9);
  params->c_ab[0] = get_ext_param(p->info->ext_params, ext_params, 10);
  params->c_ab[1] = get_ext_param(p->info->ext_params, ext_params, 11);
  params->c_ab[2] = get_ext_param(p->info->ext_params, ext_params, 12);
  params->c_ab[3] = get_ext_param(p->info->ext_params, ext_params, 13);
  params->c_ab[4] = get_ext_param(p->info->ext_params, ext_params, 14);
  p->cam_alpha = get_ext_param(p->info->ext_params, ext_params, 15);
}

#include "decl_gga.h"
#include "maple2c/gga_exc/gga_xc_b97.c"
#include "work_gga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_b97 = {
  XC_HYB_GGA_XC_B97,
  XC_EXCHANGE_CORRELATION,
  "Becke 97",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Becke1997_8554, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-21,
  16, ext_params_par_b97, set_ext_params_hybrid,
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_b97_1 = {
  XC_HYB_GGA_XC_B97_1,
  XC_EXCHANGE_CORRELATION,
  "Becke 97-1",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Hamprecht1998_6264, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-21,
  16, ext_params_par_b97_1, set_ext_params_hybrid,
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_b97_2 = {
  XC_HYB_GGA_XC_B97_2,
  XC_EXCHANGE_CORRELATION,
  "Becke 97-2",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Wilson2001_9233, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-21,
  16, ext_params_par_b97_2, set_ext_params_hybrid,
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_b97_d = {
  XC_GGA_XC_B97_D,
  XC_EXCHANGE_CORRELATION,
  "Becke 97-D",
  XC_FAMILY_GGA,
  {&xc_ref_Grimme2006_1787, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-21,
  15, ext_params_par_b97_d, set_ext_params_pure,
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_b97_k = {
  XC_HYB_GGA_XC_B97_K,
  XC_EXCHANGE_CORRELATION,
  "Boese-Martin for Kinetics",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Boese2004_3405, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-21,
  16, ext_params_par_b97_k, set_ext_params_hybrid,
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_b97_3 = {
  XC_HYB_GGA_XC_B97_3,
  XC_EXCHANGE_CORRELATION,
  "Becke 97-3",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Keal2005_121103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-21,
  16, ext_params_par_b97_3, set_ext_params_hybrid,
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_hcth_93 = {
  XC_GGA_XC_HCTH_93,
  XC_EXCHANGE_CORRELATION,
  "HCTH/93",
  XC_FAMILY_GGA,
  {&xc_ref_Hamprecht1998_6264, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-21,
  15, ext_params_par_hcth_93, set_ext_params_pure,
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_hcth_120 = {
  XC_GGA_XC_HCTH_120,
  XC_EXCHANGE_CORRELATION,
  "HCTH/120",
  XC_FAMILY_GGA,
  {&xc_ref_Boese2000_1670, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-21,
  15, ext_params_par_hcth_120, set_ext_params_pure,
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_hcth_147 = {
  XC_GGA_XC_HCTH_147,
  XC_EXCHANGE_CORRELATION,
  "HCTH/147",
  XC_FAMILY_GGA,
  {&xc_ref_Boese2000_1670, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-21,
  15, ext_params_par_hcth_147, set_ext_params_pure,
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_hcth_407 = {
  XC_GGA_XC_HCTH_407,
  XC_EXCHANGE_CORRELATION,
  "HCTH/407",
  XC_FAMILY_GGA,
  {&xc_ref_Boese2001_5497, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-21,
  15, ext_params_par_hcth_407, set_ext_params_pure,
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_sb98_1a = {
  XC_HYB_GGA_XC_SB98_1A,
  XC_EXCHANGE_CORRELATION,
  "SB98 (1a)",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-21,
  16, ext_params_par_sb98_1a, set_ext_params_hybrid,
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_sb98_1b = {
  XC_HYB_GGA_XC_SB98_1B,
  XC_EXCHANGE_CORRELATION,
  "SB98 (1b)",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-21,
  16, ext_params_par_sb98_1b, set_ext_params_hybrid,
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_sb98_1c = {
  XC_HYB_GGA_XC_SB98_1C,
  XC_EXCHANGE_CORRELATION,
  "SB98 (1c)",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-21,
  16, ext_params_par_sb98_1c, set_ext_params_hybrid,
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_sb98_2a = {
  XC_HYB_GGA_XC_SB98_2A,
  XC_EXCHANGE_CORRELATION,
  "SB98 (2a)",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-21,
  16, ext_params_par_sb98_2a, set_ext_params_hybrid,
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_sb98_2b = {
  XC_HYB_GGA_XC_SB98_2B,
  XC_EXCHANGE_CORRELATION,
  "SB98 (2b)",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-21,
  16, ext_params_par_sb98_2b, set_ext_params_hybrid,
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_sb98_2c = {
  XC_HYB_GGA_XC_SB98_2C,
  XC_EXCHANGE_CORRELATION,
  "SB98 (2c)",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Schmider1998_9624, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-21,
  16, ext_params_par_sb98_2c, set_ext_params_hybrid,
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_b97_gga1 = {
  XC_GGA_XC_B97_GGA1,
  XC_EXCHANGE_CORRELATION,
  "Becke 97 GGA-1",
  XC_FAMILY_GGA,
  {&xc_ref_Cohen2000_160, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-21,
  15, ext_params_par_b97_gga1, set_ext_params_pure,
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_hcth_p14 = {
  XC_GGA_XC_HCTH_P14,
  XC_EXCHANGE_CORRELATION,
  "HCTH p=1/4",
  XC_FAMILY_GGA,
  {&xc_ref_Menconi2001_3958, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-21,
  15, ext_params_par_hcth_p14, set_ext_params_pure,
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_hcth_p76 = {
  XC_GGA_XC_HCTH_P76,
  XC_EXCHANGE_CORRELATION,
  "HCTH p=7/6",
  XC_FAMILY_GGA,
  {&xc_ref_Menconi2001_3958, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-21,
  15, ext_params_par_hcth_p76, set_ext_params_pure,
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_hcth_407p = {
  XC_GGA_XC_HCTH_407P,
  XC_EXCHANGE_CORRELATION,
  "HCTH/407+",
  XC_FAMILY_GGA,
  {&xc_ref_Boese2003_5965, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-21,
  15, ext_params_par_hcth_407p, set_ext_params_pure,
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_gga_xc_b97_1p = {
  XC_HYB_GGA_XC_B97_1P,
  XC_EXCHANGE_CORRELATION,
  "version of B97 by Cohen and Handy",
  XC_FAMILY_HYB_GGA,
  {&xc_ref_Cohen2000_160, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-21,
  16, ext_params_par_b97_1p, set_ext_params_hybrid,
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_gga_xc_hle16 = {
  XC_GGA_XC_HLE16,
  XC_EXCHANGE_CORRELATION,
  "high local exchange 2016",
  XC_FAMILY_GGA,
  {&xc_ref_Verma2017_380, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | MAPLE2C_FLAGS,
  1e-21,
  15, ext_params_par_hle16, set_ext_params_pure,
  gga_xc_b97_init, NULL,
  NULL, work_gga, NULL
};
