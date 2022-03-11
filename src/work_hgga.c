/*
 Copyright (C) 2006-2018 M.A.L. Marques
 Copyright (C) 2019 X. Andrade

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/**
 * @file work_hgga.c
 * @brief This file is to be included in hGGA functionals.
 */

/* define auxiliary functions to NULL in case they are not available */
#if defined(XC_DONT_COMPILE_EXC) || maple2c_order < 0 || defined(XC_NO_EXC)
#define work_hgga_exc_unpol NULL
#define work_hgga_exc_pol NULL
#else
#define ORDER_TXT exc
#define SPIN_TXT  unpol
#include "work_hgga_inc.c"
#undef SPIN_TXT
#define SPIN_TXT  pol
#include "work_hgga_inc.c"
#undef SPIN_TXT
#undef ORDER_TXT
#endif

#if defined(XC_DONT_COMPILE_VXC) || maple2c_order < 1
#define work_hgga_vxc_unpol NULL
#define work_hgga_vxc_pol NULL
#else
#define ORDER_TXT vxc
#define SPIN_TXT  unpol
#include "work_hgga_inc.c"
#undef SPIN_TXT
#define SPIN_TXT  pol
#include "work_hgga_inc.c"
#undef SPIN_TXT
#undef ORDER_TXT
#endif

#if defined(XC_DONT_COMPILE_FXC) || maple2c_order < 2
#define work_hgga_fxc_unpol NULL
#define work_hgga_fxc_pol NULL
#else
#define ORDER_TXT fxc
#define SPIN_TXT  unpol
#include "work_hgga_inc.c"
#undef SPIN_TXT
#define SPIN_TXT  pol
#include "work_hgga_inc.c"
#undef SPIN_TXT
#undef ORDER_TXT
#endif

#if defined(XC_DONT_COMPILE_KXC) || maple2c_order < 3
#define work_hgga_kxc_unpol NULL
#define work_hgga_kxc_pol NULL
#else
#define ORDER_TXT kxc
#define SPIN_TXT  unpol
#include "work_hgga_inc.c"
#undef SPIN_TXT
#define SPIN_TXT  pol
#include "work_hgga_inc.c"
#undef SPIN_TXT
#undef ORDER_TXT
#endif

#if defined(XC_DONT_COMPILE_LXC) || maple2c_order < 4
#define work_hgga_lxc_unpol NULL
#define work_hgga_lxc_pol NULL
#else
#define ORDER_TXT lxc
#define SPIN_TXT  unpol
#include "work_hgga_inc.c"
#undef SPIN_TXT
#define SPIN_TXT  pol
#include "work_hgga_inc.c"
#undef SPIN_TXT
#undef ORDER_TXT
#endif

/* we construct a structure containing all variants */
static xc_hgga_funcs_variants work_hgga =
  {
   {work_hgga_exc_unpol, work_hgga_vxc_unpol, work_hgga_fxc_unpol, work_hgga_kxc_unpol, work_hgga_lxc_unpol},
   {work_hgga_exc_pol,   work_hgga_vxc_pol,   work_hgga_fxc_pol,   work_hgga_kxc_pol,   work_hgga_lxc_pol}
  };