/*
  This file was generated automatically with scripts/maple2c.py.
  Do not edit this file directly as it can be overwritten!!

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.

  Maple version     : Maple 2020 (X86 64 LINUX)
  Maple source      : ./maple/lda_exc/lda_c_gk72.mpl
  Type of functional: lda_exc
*/

#define maple2c_order 4
#define MAPLE2C_FLAGS (XC_FLAGS_I_HAVE_EXC | XC_FLAGS_I_HAVE_VXC | XC_FLAGS_I_HAVE_FXC | XC_FLAGS_I_HAVE_KXC | XC_FLAGS_I_HAVE_LXC)


#ifndef XC_DONT_COMPILE_EXC
GPU_DEVICE_FUNCTION static inline void
func_exc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t1, t2, t3, t4, t5, t6, t7, t9;
  double t10, t11, t12, t13, t20, t23, t25, t29;
  double t30, t35, t37, t38, t42, t44, t48, tzk0;


  t1 = M_CBRT3;
  t2 = 0.1e1 / M_PI;
  t3 = POW_1_3(t2);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = POW_1_3(rho[0]);
  t9 = t6 / t7;
  t10 = t4 * t9;
  t11 = t10 / 0.4e1;
  t12 = t11 < 0.7e0;
  t13 = log(t11);
  t20 = t11 < 0.1e2;
  t23 = t1 * t1;
  t25 = t23 / t3;
  t29 = sqrt(0.4e1);
  t30 = sqrt(t10);
  t35 = t3 * t3;
  t37 = t1 / t35;
  t38 = t7 * t7;
  t42 = t23 * t35;
  t44 = t5 / t38;
  t48 = 0.1e1 / t30 / t42 / t44 / 0.4e1;
  tzk0 = my_piecewise5(t12, 0.311e-1 * t13 - 0.48e-1 + 0.225e-2 * t4 * t9 * t13 - 0.425e-2 * t10, t20, -0.6156e-1 + 0.1898e-1 * t13, 0.146e0 * t25 * t5 * t7 + 0.53e1 * t29 / t30 / t10 - 0.49e0 * t37 * t6 * t38 - 0.64e1 * t29 * t48);

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->out_dim->zk + 0] += tzk0;

}

#endif


#ifndef XC_DONT_COMPILE_VXC
GPU_DEVICE_FUNCTION static inline void
func_vxc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t1, t2, t3, t4, t5, t6, t7, t9;
  double t10, t11, t12, t13, t20, t23, t25, t29;
  double t30, t35, t37, t38, t42, t44, t48, tzk0;

  double t52, t55, t56, t66, t67, t68, t76, t77;
  double t81, tvrho0;


  t1 = M_CBRT3;
  t2 = 0.1e1 / M_PI;
  t3 = POW_1_3(t2);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = POW_1_3(rho[0]);
  t9 = t6 / t7;
  t10 = t4 * t9;
  t11 = t10 / 0.4e1;
  t12 = t11 < 0.7e0;
  t13 = log(t11);
  t20 = t11 < 0.1e2;
  t23 = t1 * t1;
  t25 = t23 / t3;
  t29 = sqrt(0.4e1);
  t30 = sqrt(t10);
  t35 = t3 * t3;
  t37 = t1 / t35;
  t38 = t7 * t7;
  t42 = t23 * t35;
  t44 = t5 / t38;
  t48 = 0.1e1 / t30 / t42 / t44 / 0.4e1;
  tzk0 = my_piecewise5(t12, 0.311e-1 * t13 - 0.48e-1 + 0.225e-2 * t4 * t9 * t13 - 0.425e-2 * t10, t20, -0.6156e-1 + 0.1898e-1 * t13, 0.146e0 * t25 * t5 * t7 + 0.53e1 * t29 / t30 / t10 - 0.49e0 * t37 * t6 * t38 - 0.64e1 * t29 * t48);

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->out_dim->zk + 0] += tzk0;

  t52 = 0.1e1 / rho[0];
  t55 = 0.1e1 / t7 / rho[0];
  t56 = t6 * t55;
  t66 = pow(0.4e1, 0.1e1 / 0.6e1);
  t67 = t66 * t48;
  t68 = t4 * t55;
  t76 = 0.1e1 / t30 / t2 / t52 / 0.48e2;
  t77 = t66 * t76;
  t81 = my_piecewise5(t12, -0.10366666666666666667e-1 * t52 - 0.75e-3 * t4 * t56 * t13 + 0.6666666666666666667e-3 * t4 * t56, t20, -0.63266666666666666667e-2 * t52, 0.48666666666666666667e-1 * t25 * t44 + 0.106e2 * t67 * t68 - 0.32666666666666666667e0 * t37 * t9 - 0.21333333333333333333e2 * t77 * t68);
  tvrho0 = rho[0] * t81 + tzk0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->out_dim->vrho + 0] += tvrho0;

}

#endif


#ifndef XC_DONT_COMPILE_FXC
GPU_DEVICE_FUNCTION static inline void
func_fxc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t1, t2, t3, t4, t5, t6, t7, t9;
  double t10, t11, t12, t13, t20, t23, t25, t29;
  double t30, t35, t37, t38, t42, t44, t48, tzk0;

  double t52, t55, t56, t66, t67, t68, t76, t77;
  double t81, tvrho0;

  double t84, t85, t88, t89, t99, t102, t103, t104;
  double t105, t107, t108, t111, t121, t122, t128, tv2rho20;


  t1 = M_CBRT3;
  t2 = 0.1e1 / M_PI;
  t3 = POW_1_3(t2);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = POW_1_3(rho[0]);
  t9 = t6 / t7;
  t10 = t4 * t9;
  t11 = t10 / 0.4e1;
  t12 = t11 < 0.7e0;
  t13 = log(t11);
  t20 = t11 < 0.1e2;
  t23 = t1 * t1;
  t25 = t23 / t3;
  t29 = sqrt(0.4e1);
  t30 = sqrt(t10);
  t35 = t3 * t3;
  t37 = t1 / t35;
  t38 = t7 * t7;
  t42 = t23 * t35;
  t44 = t5 / t38;
  t48 = 0.1e1 / t30 / t42 / t44 / 0.4e1;
  tzk0 = my_piecewise5(t12, 0.311e-1 * t13 - 0.48e-1 + 0.225e-2 * t4 * t9 * t13 - 0.425e-2 * t10, t20, -0.6156e-1 + 0.1898e-1 * t13, 0.146e0 * t25 * t5 * t7 + 0.53e1 * t29 / t30 / t10 - 0.49e0 * t37 * t6 * t38 - 0.64e1 * t29 * t48);

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->out_dim->zk + 0] += tzk0;

  t52 = 0.1e1 / rho[0];
  t55 = 0.1e1 / t7 / rho[0];
  t56 = t6 * t55;
  t66 = pow(0.4e1, 0.1e1 / 0.6e1);
  t67 = t66 * t48;
  t68 = t4 * t55;
  t76 = 0.1e1 / t30 / t2 / t52 / 0.48e2;
  t77 = t66 * t76;
  t81 = my_piecewise5(t12, -0.10366666666666666667e-1 * t52 - 0.75e-3 * t4 * t56 * t13 + 0.6666666666666666667e-3 * t4 * t56, t20, -0.63266666666666666667e-2 * t52, 0.48666666666666666667e-1 * t25 * t44 + 0.106e2 * t67 * t68 - 0.32666666666666666667e0 * t37 * t9 - 0.21333333333333333333e2 * t77 * t68);
  tvrho0 = rho[0] * t81 + tzk0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->out_dim->vrho + 0] += tvrho0;

  t84 = rho[0] * rho[0];
  t85 = 0.1e1 / t84;
  t88 = 0.1e1 / t7 / t84;
  t89 = t6 * t88;
  t99 = t5 / t38 / rho[0];
  t102 = t66 * t66;
  t103 = t102 * t102;
  t104 = t103 * t66;
  t105 = t104 * t76;
  t107 = 0.1e1 / t38 / t84;
  t108 = t42 * t107;
  t111 = t4 * t88;
  t121 = 0.1e1 / t30 / t1 / t3 / t2 / t56 / 0.48e2;
  t122 = t104 * t121;
  t128 = my_piecewise5(t12, 0.10366666666666666667e-1 * t85 + 0.1e-2 * t4 * t89 * t13 - 0.63888888888888888893e-3 * t4 * t89, t20, 0.63266666666666666667e-2 * t85, -0.32444444444444444445e-1 * t25 * t99 + 0.88333333333333333333e1 * t105 * t108 - 0.14133333333333333333e2 * t67 * t111 + 0.10888888888888888889e0 * t37 * t56 - 0.24888888888888888889e2 * t122 * t108 + 0.28444444444444444444e2 * t77 * t111);
  tv2rho20 = rho[0] * t128 + 0.2e1 * t81;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->out_dim->v2rho2 + 0] += tv2rho20;

}

#endif


#ifndef XC_DONT_COMPILE_KXC
GPU_DEVICE_FUNCTION static inline void
func_kxc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t1, t2, t3, t4, t5, t6, t7, t9;
  double t10, t11, t12, t13, t20, t23, t25, t29;
  double t30, t35, t37, t38, t42, t44, t48, tzk0;

  double t52, t55, t56, t66, t67, t68, t76, t77;
  double t81, tvrho0;

  double t84, t85, t88, t89, t99, t102, t103, t104;
  double t105, t107, t108, t111, t121, t122, t128, tv2rho20;

  double t131, t132, t135, t136, t147, t148, t149, t150;
  double t154, t155, t158, t168, t169, t177, tv3rho30;


  t1 = M_CBRT3;
  t2 = 0.1e1 / M_PI;
  t3 = POW_1_3(t2);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = POW_1_3(rho[0]);
  t9 = t6 / t7;
  t10 = t4 * t9;
  t11 = t10 / 0.4e1;
  t12 = t11 < 0.7e0;
  t13 = log(t11);
  t20 = t11 < 0.1e2;
  t23 = t1 * t1;
  t25 = t23 / t3;
  t29 = sqrt(0.4e1);
  t30 = sqrt(t10);
  t35 = t3 * t3;
  t37 = t1 / t35;
  t38 = t7 * t7;
  t42 = t23 * t35;
  t44 = t5 / t38;
  t48 = 0.1e1 / t30 / t42 / t44 / 0.4e1;
  tzk0 = my_piecewise5(t12, 0.311e-1 * t13 - 0.48e-1 + 0.225e-2 * t4 * t9 * t13 - 0.425e-2 * t10, t20, -0.6156e-1 + 0.1898e-1 * t13, 0.146e0 * t25 * t5 * t7 + 0.53e1 * t29 / t30 / t10 - 0.49e0 * t37 * t6 * t38 - 0.64e1 * t29 * t48);

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->out_dim->zk + 0] += tzk0;

  t52 = 0.1e1 / rho[0];
  t55 = 0.1e1 / t7 / rho[0];
  t56 = t6 * t55;
  t66 = pow(0.4e1, 0.1e1 / 0.6e1);
  t67 = t66 * t48;
  t68 = t4 * t55;
  t76 = 0.1e1 / t30 / t2 / t52 / 0.48e2;
  t77 = t66 * t76;
  t81 = my_piecewise5(t12, -0.10366666666666666667e-1 * t52 - 0.75e-3 * t4 * t56 * t13 + 0.6666666666666666667e-3 * t4 * t56, t20, -0.63266666666666666667e-2 * t52, 0.48666666666666666667e-1 * t25 * t44 + 0.106e2 * t67 * t68 - 0.32666666666666666667e0 * t37 * t9 - 0.21333333333333333333e2 * t77 * t68);
  tvrho0 = rho[0] * t81 + tzk0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->out_dim->vrho + 0] += tvrho0;

  t84 = rho[0] * rho[0];
  t85 = 0.1e1 / t84;
  t88 = 0.1e1 / t7 / t84;
  t89 = t6 * t88;
  t99 = t5 / t38 / rho[0];
  t102 = t66 * t66;
  t103 = t102 * t102;
  t104 = t103 * t66;
  t105 = t104 * t76;
  t107 = 0.1e1 / t38 / t84;
  t108 = t42 * t107;
  t111 = t4 * t88;
  t121 = 0.1e1 / t30 / t1 / t3 / t2 / t56 / 0.48e2;
  t122 = t104 * t121;
  t128 = my_piecewise5(t12, 0.10366666666666666667e-1 * t85 + 0.1e-2 * t4 * t89 * t13 - 0.63888888888888888893e-3 * t4 * t89, t20, 0.63266666666666666667e-2 * t85, -0.32444444444444444445e-1 * t25 * t99 + 0.88333333333333333333e1 * t105 * t108 - 0.14133333333333333333e2 * t67 * t111 + 0.10888888888888888889e0 * t37 * t56 - 0.24888888888888888889e2 * t122 * t108 + 0.28444444444444444444e2 * t77 * t111);
  tv2rho20 = rho[0] * t128 + 0.2e1 * t81;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->out_dim->v2rho2 + 0] += tv2rho20;

  t131 = t84 * rho[0];
  t132 = 0.1e1 / t131;
  t135 = 0.1e1 / t7 / t131;
  t136 = t6 * t135;
  t147 = t29 * t121;
  t148 = t84 * t84;
  t149 = 0.1e1 / t148;
  t150 = t2 * t149;
  t154 = 0.1e1 / t38 / t131;
  t155 = t42 * t154;
  t158 = t4 * t135;
  t168 = 0.1e1 / t30 / t23 / t35 / t2 / t99 / 0.192e3;
  t169 = t29 * t168;
  t177 = my_piecewise5(t12, -0.20733333333333333334e-1 * t132 - 0.23333333333333333333e-2 * t4 * t136 * t13 + 0.11574074074074074075e-2 * t4 * t136, t20, -0.12653333333333333333e-1 * t132, 0.54074074074074074075e-1 * t25 * t5 * t107 + 0.12366666666666666667e3 * t147 * t150 - 0.35333333333333333332e2 * t105 * t155 + 0.32977777777777777777e2 * t67 * t158 - 0.14518518518518518519e0 * t37 * t89 - 0.448e3 * t169 * t150 + 0.99555555555555555556e2 * t122 * t155 - 0.66370370370370370369e2 * t77 * t158);
  tv3rho30 = rho[0] * t177 + 0.3e1 * t128;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->out_dim->v3rho3 + 0] += tv3rho30;

}

#endif


#ifndef XC_DONT_COMPILE_LXC
GPU_DEVICE_FUNCTION static inline void
func_lxc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t1, t2, t3, t4, t5, t6, t7, t9;
  double t10, t11, t12, t13, t20, t23, t25, t29;
  double t30, t35, t37, t38, t42, t44, t48, tzk0;

  double t52, t55, t56, t66, t67, t68, t76, t77;
  double t81, tvrho0;

  double t84, t85, t88, t89, t99, t102, t103, t104;
  double t105, t107, t108, t111, t121, t122, t128, tv2rho20;

  double t131, t132, t135, t136, t147, t148, t149, t150;
  double t154, t155, t158, t168, t169, t177, tv3rho30;

  double t182, t183, t196, t200, t204, t209, t212, t217;
  double t234, tv4rho40;


  t1 = M_CBRT3;
  t2 = 0.1e1 / M_PI;
  t3 = POW_1_3(t2);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = POW_1_3(rho[0]);
  t9 = t6 / t7;
  t10 = t4 * t9;
  t11 = t10 / 0.4e1;
  t12 = t11 < 0.7e0;
  t13 = log(t11);
  t20 = t11 < 0.1e2;
  t23 = t1 * t1;
  t25 = t23 / t3;
  t29 = sqrt(0.4e1);
  t30 = sqrt(t10);
  t35 = t3 * t3;
  t37 = t1 / t35;
  t38 = t7 * t7;
  t42 = t23 * t35;
  t44 = t5 / t38;
  t48 = 0.1e1 / t30 / t42 / t44 / 0.4e1;
  tzk0 = my_piecewise5(t12, 0.311e-1 * t13 - 0.48e-1 + 0.225e-2 * t4 * t9 * t13 - 0.425e-2 * t10, t20, -0.6156e-1 + 0.1898e-1 * t13, 0.146e0 * t25 * t5 * t7 + 0.53e1 * t29 / t30 / t10 - 0.49e0 * t37 * t6 * t38 - 0.64e1 * t29 * t48);

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->out_dim->zk + 0] += tzk0;

  t52 = 0.1e1 / rho[0];
  t55 = 0.1e1 / t7 / rho[0];
  t56 = t6 * t55;
  t66 = pow(0.4e1, 0.1e1 / 0.6e1);
  t67 = t66 * t48;
  t68 = t4 * t55;
  t76 = 0.1e1 / t30 / t2 / t52 / 0.48e2;
  t77 = t66 * t76;
  t81 = my_piecewise5(t12, -0.10366666666666666667e-1 * t52 - 0.75e-3 * t4 * t56 * t13 + 0.6666666666666666667e-3 * t4 * t56, t20, -0.63266666666666666667e-2 * t52, 0.48666666666666666667e-1 * t25 * t44 + 0.106e2 * t67 * t68 - 0.32666666666666666667e0 * t37 * t9 - 0.21333333333333333333e2 * t77 * t68);
  tvrho0 = rho[0] * t81 + tzk0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->out_dim->vrho + 0] += tvrho0;

  t84 = rho[0] * rho[0];
  t85 = 0.1e1 / t84;
  t88 = 0.1e1 / t7 / t84;
  t89 = t6 * t88;
  t99 = t5 / t38 / rho[0];
  t102 = t66 * t66;
  t103 = t102 * t102;
  t104 = t103 * t66;
  t105 = t104 * t76;
  t107 = 0.1e1 / t38 / t84;
  t108 = t42 * t107;
  t111 = t4 * t88;
  t121 = 0.1e1 / t30 / t1 / t3 / t2 / t56 / 0.48e2;
  t122 = t104 * t121;
  t128 = my_piecewise5(t12, 0.10366666666666666667e-1 * t85 + 0.1e-2 * t4 * t89 * t13 - 0.63888888888888888893e-3 * t4 * t89, t20, 0.63266666666666666667e-2 * t85, -0.32444444444444444445e-1 * t25 * t99 + 0.88333333333333333333e1 * t105 * t108 - 0.14133333333333333333e2 * t67 * t111 + 0.10888888888888888889e0 * t37 * t56 - 0.24888888888888888889e2 * t122 * t108 + 0.28444444444444444444e2 * t77 * t111);
  tv2rho20 = rho[0] * t128 + 0.2e1 * t81;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->out_dim->v2rho2 + 0] += tv2rho20;

  t131 = t84 * rho[0];
  t132 = 0.1e1 / t131;
  t135 = 0.1e1 / t7 / t131;
  t136 = t6 * t135;
  t147 = t29 * t121;
  t148 = t84 * t84;
  t149 = 0.1e1 / t148;
  t150 = t2 * t149;
  t154 = 0.1e1 / t38 / t131;
  t155 = t42 * t154;
  t158 = t4 * t135;
  t168 = 0.1e1 / t30 / t23 / t35 / t2 / t99 / 0.192e3;
  t169 = t29 * t168;
  t177 = my_piecewise5(t12, -0.20733333333333333334e-1 * t132 - 0.23333333333333333333e-2 * t4 * t136 * t13 + 0.11574074074074074075e-2 * t4 * t136, t20, -0.12653333333333333333e-1 * t132, 0.54074074074074074075e-1 * t25 * t5 * t107 + 0.12366666666666666667e3 * t147 * t150 - 0.35333333333333333332e2 * t105 * t155 + 0.32977777777777777777e2 * t67 * t158 - 0.14518518518518518519e0 * t37 * t89 - 0.448e3 * t169 * t150 + 0.99555555555555555556e2 * t122 * t155 - 0.66370370370370370369e2 * t77 * t158);
  tv3rho30 = rho[0] * t177 + 0.3e1 * t128;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->out_dim->v3rho3 + 0] += tv3rho30;

  t182 = 0.1e1 / t7 / t148;
  t183 = t6 * t182;
  t196 = t148 * rho[0];
  t200 = 0.1e1 / t7 / t196 * t1 * t3;
  t204 = t2 / t196;
  t209 = t42 / t38 / t148;
  t212 = t4 * t182;
  t217 = M_PI * M_PI;
  t234 = my_piecewise5(t12, 0.62200000000000000002e-1 * t149 + 0.77777777777777777777e-2 * t4 * t183 * t13 - 0.30802469135802469139e-2 * t4 * t183, t20, 0.37959999999999999999e-1 * t149, -0.14419753086419753087e0 * t25 * t5 * t154 + 0.74200000000000000002e3 * t66 * t168 * t2 * t200 - 0.98933333333333333333e3 * t147 * t204 + 0.15703703703703703703e3 * t105 * t209 - 0.10992592592592592592e3 * t67 * t212 + 0.33876543209876543211e0 * t37 * t136 - 0.14259259259259259259e1 * t66 / t30 * t217 / t85 * t2 * t200 + 0.3584e4 * t169 * t204 - 0.4424691358024691358e3 * t122 * t209 + 0.2212345679012345679e3 * t77 * t212);
  tv4rho40 = rho[0] * t234 + 0.4e1 * t177;

  if(out->v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    out->v4rho4[ip*p->out_dim->v4rho4 + 0] += tv4rho40;

}

#endif


#ifndef XC_DONT_COMPILE_EXC
GPU_DEVICE_FUNCTION static inline void
func_exc_pol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t1, t2, t3, t4, t5, t6, t7, t8;
  double t10, t11, t12, t13, t14, t21, t24, t26;
  double t30, t31, t36, t38, t39, t43, t45, t49;
  double tzk0;


  t1 = M_CBRT3;
  t2 = 0.1e1 / M_PI;
  t3 = POW_1_3(t2);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = rho[0] + rho[1];
  t8 = POW_1_3(t7);
  t10 = t6 / t8;
  t11 = t4 * t10;
  t12 = t11 / 0.4e1;
  t13 = t12 < 0.7e0;
  t14 = log(t12);
  t21 = t12 < 0.1e2;
  t24 = t1 * t1;
  t26 = t24 / t3;
  t30 = sqrt(0.4e1);
  t31 = sqrt(t11);
  t36 = t3 * t3;
  t38 = t1 / t36;
  t39 = t8 * t8;
  t43 = t24 * t36;
  t45 = t5 / t39;
  t49 = 0.1e1 / t31 / t43 / t45 / 0.4e1;
  tzk0 = my_piecewise5(t13, 0.311e-1 * t14 - 0.48e-1 + 0.225e-2 * t4 * t10 * t14 - 0.425e-2 * t11, t21, -0.6156e-1 + 0.1898e-1 * t14, 0.146e0 * t26 * t5 * t8 + 0.53e1 * t30 / t31 / t11 - 0.49e0 * t38 * t6 * t39 - 0.64e1 * t30 * t49);

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->out_dim->zk + 0] += tzk0;

}

#endif


#ifndef XC_DONT_COMPILE_VXC
GPU_DEVICE_FUNCTION static inline void
func_vxc_pol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t1, t2, t3, t4, t5, t6, t7, t8;
  double t10, t11, t12, t13, t14, t21, t24, t26;
  double t30, t31, t36, t38, t39, t43, t45, t49;
  double tzk0;

  double t53, t56, t57, t67, t68, t69, t77, t78;
  double t82, tvrho0, tvrho1;


  t1 = M_CBRT3;
  t2 = 0.1e1 / M_PI;
  t3 = POW_1_3(t2);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = rho[0] + rho[1];
  t8 = POW_1_3(t7);
  t10 = t6 / t8;
  t11 = t4 * t10;
  t12 = t11 / 0.4e1;
  t13 = t12 < 0.7e0;
  t14 = log(t12);
  t21 = t12 < 0.1e2;
  t24 = t1 * t1;
  t26 = t24 / t3;
  t30 = sqrt(0.4e1);
  t31 = sqrt(t11);
  t36 = t3 * t3;
  t38 = t1 / t36;
  t39 = t8 * t8;
  t43 = t24 * t36;
  t45 = t5 / t39;
  t49 = 0.1e1 / t31 / t43 / t45 / 0.4e1;
  tzk0 = my_piecewise5(t13, 0.311e-1 * t14 - 0.48e-1 + 0.225e-2 * t4 * t10 * t14 - 0.425e-2 * t11, t21, -0.6156e-1 + 0.1898e-1 * t14, 0.146e0 * t26 * t5 * t8 + 0.53e1 * t30 / t31 / t11 - 0.49e0 * t38 * t6 * t39 - 0.64e1 * t30 * t49);

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->out_dim->zk + 0] += tzk0;

  t53 = 0.1e1 / t7;
  t56 = 0.1e1 / t8 / t7;
  t57 = t6 * t56;
  t67 = pow(0.4e1, 0.1e1 / 0.6e1);
  t68 = t67 * t49;
  t69 = t4 * t56;
  t77 = 0.1e1 / t31 / t2 / t53 / 0.48e2;
  t78 = t67 * t77;
  t82 = my_piecewise5(t13, -0.10366666666666666667e-1 * t53 - 0.75e-3 * t4 * t57 * t14 + 0.6666666666666666667e-3 * t4 * t57, t21, -0.63266666666666666667e-2 * t53, 0.48666666666666666667e-1 * t26 * t45 + 0.106e2 * t68 * t69 - 0.32666666666666666667e0 * t38 * t10 - 0.21333333333333333333e2 * t78 * t69);
  tvrho0 = t7 * t82 + tzk0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->out_dim->vrho + 0] += tvrho0;

  tvrho1 = tvrho0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->out_dim->vrho + 1] += tvrho1;

}

#endif


#ifndef XC_DONT_COMPILE_FXC
GPU_DEVICE_FUNCTION static inline void
func_fxc_pol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t1, t2, t3, t4, t5, t6, t7, t8;
  double t10, t11, t12, t13, t14, t21, t24, t26;
  double t30, t31, t36, t38, t39, t43, t45, t49;
  double tzk0;

  double t53, t56, t57, t67, t68, t69, t77, t78;
  double t82, tvrho0, tvrho1;

  double t85, t86, t89, t90, t100, t103, t104, t105;
  double t106, t108, t109, t112, t122, t123, t129, tv2rho20;
  double tv2rho21, tv2rho22;


  t1 = M_CBRT3;
  t2 = 0.1e1 / M_PI;
  t3 = POW_1_3(t2);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = rho[0] + rho[1];
  t8 = POW_1_3(t7);
  t10 = t6 / t8;
  t11 = t4 * t10;
  t12 = t11 / 0.4e1;
  t13 = t12 < 0.7e0;
  t14 = log(t12);
  t21 = t12 < 0.1e2;
  t24 = t1 * t1;
  t26 = t24 / t3;
  t30 = sqrt(0.4e1);
  t31 = sqrt(t11);
  t36 = t3 * t3;
  t38 = t1 / t36;
  t39 = t8 * t8;
  t43 = t24 * t36;
  t45 = t5 / t39;
  t49 = 0.1e1 / t31 / t43 / t45 / 0.4e1;
  tzk0 = my_piecewise5(t13, 0.311e-1 * t14 - 0.48e-1 + 0.225e-2 * t4 * t10 * t14 - 0.425e-2 * t11, t21, -0.6156e-1 + 0.1898e-1 * t14, 0.146e0 * t26 * t5 * t8 + 0.53e1 * t30 / t31 / t11 - 0.49e0 * t38 * t6 * t39 - 0.64e1 * t30 * t49);

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->out_dim->zk + 0] += tzk0;

  t53 = 0.1e1 / t7;
  t56 = 0.1e1 / t8 / t7;
  t57 = t6 * t56;
  t67 = pow(0.4e1, 0.1e1 / 0.6e1);
  t68 = t67 * t49;
  t69 = t4 * t56;
  t77 = 0.1e1 / t31 / t2 / t53 / 0.48e2;
  t78 = t67 * t77;
  t82 = my_piecewise5(t13, -0.10366666666666666667e-1 * t53 - 0.75e-3 * t4 * t57 * t14 + 0.6666666666666666667e-3 * t4 * t57, t21, -0.63266666666666666667e-2 * t53, 0.48666666666666666667e-1 * t26 * t45 + 0.106e2 * t68 * t69 - 0.32666666666666666667e0 * t38 * t10 - 0.21333333333333333333e2 * t78 * t69);
  tvrho0 = t7 * t82 + tzk0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->out_dim->vrho + 0] += tvrho0;

  tvrho1 = tvrho0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->out_dim->vrho + 1] += tvrho1;

  t85 = t7 * t7;
  t86 = 0.1e1 / t85;
  t89 = 0.1e1 / t8 / t85;
  t90 = t6 * t89;
  t100 = t5 / t39 / t7;
  t103 = t67 * t67;
  t104 = t103 * t103;
  t105 = t104 * t67;
  t106 = t105 * t77;
  t108 = 0.1e1 / t39 / t85;
  t109 = t43 * t108;
  t112 = t4 * t89;
  t122 = 0.1e1 / t31 / t1 / t3 / t2 / t57 / 0.48e2;
  t123 = t105 * t122;
  t129 = my_piecewise5(t13, 0.10366666666666666667e-1 * t86 + 0.1e-2 * t4 * t90 * t14 - 0.63888888888888888893e-3 * t4 * t90, t21, 0.63266666666666666667e-2 * t86, -0.32444444444444444445e-1 * t26 * t100 + 0.88333333333333333333e1 * t106 * t109 - 0.14133333333333333333e2 * t68 * t112 + 0.10888888888888888889e0 * t38 * t57 - 0.24888888888888888889e2 * t123 * t109 + 0.28444444444444444444e2 * t78 * t112);
  tv2rho20 = t7 * t129 + 0.2e1 * t82;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->out_dim->v2rho2 + 0] += tv2rho20;

  tv2rho21 = tv2rho20;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->out_dim->v2rho2 + 1] += tv2rho21;

  tv2rho22 = tv2rho21;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->out_dim->v2rho2 + 2] += tv2rho22;

}

#endif


#ifndef XC_DONT_COMPILE_KXC
GPU_DEVICE_FUNCTION static inline void
func_kxc_pol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t1, t2, t3, t4, t5, t6, t7, t8;
  double t10, t11, t12, t13, t14, t21, t24, t26;
  double t30, t31, t36, t38, t39, t43, t45, t49;
  double tzk0;

  double t53, t56, t57, t67, t68, t69, t77, t78;
  double t82, tvrho0, tvrho1;

  double t85, t86, t89, t90, t100, t103, t104, t105;
  double t106, t108, t109, t112, t122, t123, t129, tv2rho20;
  double tv2rho21, tv2rho22;

  double t132, t133, t136, t137, t148, t149, t150, t151;
  double t155, t156, t159, t169, t170, t178, tv3rho30, tv3rho31;
  double tv3rho32, tv3rho33;


  t1 = M_CBRT3;
  t2 = 0.1e1 / M_PI;
  t3 = POW_1_3(t2);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = rho[0] + rho[1];
  t8 = POW_1_3(t7);
  t10 = t6 / t8;
  t11 = t4 * t10;
  t12 = t11 / 0.4e1;
  t13 = t12 < 0.7e0;
  t14 = log(t12);
  t21 = t12 < 0.1e2;
  t24 = t1 * t1;
  t26 = t24 / t3;
  t30 = sqrt(0.4e1);
  t31 = sqrt(t11);
  t36 = t3 * t3;
  t38 = t1 / t36;
  t39 = t8 * t8;
  t43 = t24 * t36;
  t45 = t5 / t39;
  t49 = 0.1e1 / t31 / t43 / t45 / 0.4e1;
  tzk0 = my_piecewise5(t13, 0.311e-1 * t14 - 0.48e-1 + 0.225e-2 * t4 * t10 * t14 - 0.425e-2 * t11, t21, -0.6156e-1 + 0.1898e-1 * t14, 0.146e0 * t26 * t5 * t8 + 0.53e1 * t30 / t31 / t11 - 0.49e0 * t38 * t6 * t39 - 0.64e1 * t30 * t49);

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->out_dim->zk + 0] += tzk0;

  t53 = 0.1e1 / t7;
  t56 = 0.1e1 / t8 / t7;
  t57 = t6 * t56;
  t67 = pow(0.4e1, 0.1e1 / 0.6e1);
  t68 = t67 * t49;
  t69 = t4 * t56;
  t77 = 0.1e1 / t31 / t2 / t53 / 0.48e2;
  t78 = t67 * t77;
  t82 = my_piecewise5(t13, -0.10366666666666666667e-1 * t53 - 0.75e-3 * t4 * t57 * t14 + 0.6666666666666666667e-3 * t4 * t57, t21, -0.63266666666666666667e-2 * t53, 0.48666666666666666667e-1 * t26 * t45 + 0.106e2 * t68 * t69 - 0.32666666666666666667e0 * t38 * t10 - 0.21333333333333333333e2 * t78 * t69);
  tvrho0 = t7 * t82 + tzk0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->out_dim->vrho + 0] += tvrho0;

  tvrho1 = tvrho0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->out_dim->vrho + 1] += tvrho1;

  t85 = t7 * t7;
  t86 = 0.1e1 / t85;
  t89 = 0.1e1 / t8 / t85;
  t90 = t6 * t89;
  t100 = t5 / t39 / t7;
  t103 = t67 * t67;
  t104 = t103 * t103;
  t105 = t104 * t67;
  t106 = t105 * t77;
  t108 = 0.1e1 / t39 / t85;
  t109 = t43 * t108;
  t112 = t4 * t89;
  t122 = 0.1e1 / t31 / t1 / t3 / t2 / t57 / 0.48e2;
  t123 = t105 * t122;
  t129 = my_piecewise5(t13, 0.10366666666666666667e-1 * t86 + 0.1e-2 * t4 * t90 * t14 - 0.63888888888888888893e-3 * t4 * t90, t21, 0.63266666666666666667e-2 * t86, -0.32444444444444444445e-1 * t26 * t100 + 0.88333333333333333333e1 * t106 * t109 - 0.14133333333333333333e2 * t68 * t112 + 0.10888888888888888889e0 * t38 * t57 - 0.24888888888888888889e2 * t123 * t109 + 0.28444444444444444444e2 * t78 * t112);
  tv2rho20 = t7 * t129 + 0.2e1 * t82;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->out_dim->v2rho2 + 0] += tv2rho20;

  tv2rho21 = tv2rho20;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->out_dim->v2rho2 + 1] += tv2rho21;

  tv2rho22 = tv2rho21;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->out_dim->v2rho2 + 2] += tv2rho22;

  t132 = t85 * t7;
  t133 = 0.1e1 / t132;
  t136 = 0.1e1 / t8 / t132;
  t137 = t6 * t136;
  t148 = t30 * t122;
  t149 = t85 * t85;
  t150 = 0.1e1 / t149;
  t151 = t2 * t150;
  t155 = 0.1e1 / t39 / t132;
  t156 = t43 * t155;
  t159 = t4 * t136;
  t169 = 0.1e1 / t31 / t24 / t36 / t2 / t100 / 0.192e3;
  t170 = t30 * t169;
  t178 = my_piecewise5(t13, -0.20733333333333333334e-1 * t133 - 0.23333333333333333333e-2 * t4 * t137 * t14 + 0.11574074074074074075e-2 * t4 * t137, t21, -0.12653333333333333333e-1 * t133, 0.54074074074074074075e-1 * t26 * t5 * t108 + 0.12366666666666666667e3 * t148 * t151 - 0.35333333333333333332e2 * t106 * t156 + 0.32977777777777777777e2 * t68 * t159 - 0.14518518518518518519e0 * t38 * t90 - 0.448e3 * t170 * t151 + 0.99555555555555555556e2 * t123 * t156 - 0.66370370370370370369e2 * t78 * t159);
  tv3rho30 = t7 * t178 + 0.3e1 * t129;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->out_dim->v3rho3 + 0] += tv3rho30;

  tv3rho31 = tv3rho30;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->out_dim->v3rho3 + 1] += tv3rho31;

  tv3rho32 = tv3rho31;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->out_dim->v3rho3 + 2] += tv3rho32;

  tv3rho33 = tv3rho32;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->out_dim->v3rho3 + 3] += tv3rho33;

}

#endif


#ifndef XC_DONT_COMPILE_LXC
GPU_DEVICE_FUNCTION static inline void
func_lxc_pol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t1, t2, t3, t4, t5, t6, t7, t8;
  double t10, t11, t12, t13, t14, t21, t24, t26;
  double t30, t31, t36, t38, t39, t43, t45, t49;
  double tzk0;

  double t53, t56, t57, t67, t68, t69, t77, t78;
  double t82, tvrho0, tvrho1;

  double t85, t86, t89, t90, t100, t103, t104, t105;
  double t106, t108, t109, t112, t122, t123, t129, tv2rho20;
  double tv2rho21, tv2rho22;

  double t132, t133, t136, t137, t148, t149, t150, t151;
  double t155, t156, t159, t169, t170, t178, tv3rho30, tv3rho31;
  double tv3rho32, tv3rho33;

  double t183, t184, t197, t201, t205, t210, t213, t218;
  double t235, tv4rho40, tv4rho41, tv4rho42, tv4rho43, tv4rho44;


  t1 = M_CBRT3;
  t2 = 0.1e1 / M_PI;
  t3 = POW_1_3(t2);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = rho[0] + rho[1];
  t8 = POW_1_3(t7);
  t10 = t6 / t8;
  t11 = t4 * t10;
  t12 = t11 / 0.4e1;
  t13 = t12 < 0.7e0;
  t14 = log(t12);
  t21 = t12 < 0.1e2;
  t24 = t1 * t1;
  t26 = t24 / t3;
  t30 = sqrt(0.4e1);
  t31 = sqrt(t11);
  t36 = t3 * t3;
  t38 = t1 / t36;
  t39 = t8 * t8;
  t43 = t24 * t36;
  t45 = t5 / t39;
  t49 = 0.1e1 / t31 / t43 / t45 / 0.4e1;
  tzk0 = my_piecewise5(t13, 0.311e-1 * t14 - 0.48e-1 + 0.225e-2 * t4 * t10 * t14 - 0.425e-2 * t11, t21, -0.6156e-1 + 0.1898e-1 * t14, 0.146e0 * t26 * t5 * t8 + 0.53e1 * t30 / t31 / t11 - 0.49e0 * t38 * t6 * t39 - 0.64e1 * t30 * t49);

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->out_dim->zk + 0] += tzk0;

  t53 = 0.1e1 / t7;
  t56 = 0.1e1 / t8 / t7;
  t57 = t6 * t56;
  t67 = pow(0.4e1, 0.1e1 / 0.6e1);
  t68 = t67 * t49;
  t69 = t4 * t56;
  t77 = 0.1e1 / t31 / t2 / t53 / 0.48e2;
  t78 = t67 * t77;
  t82 = my_piecewise5(t13, -0.10366666666666666667e-1 * t53 - 0.75e-3 * t4 * t57 * t14 + 0.6666666666666666667e-3 * t4 * t57, t21, -0.63266666666666666667e-2 * t53, 0.48666666666666666667e-1 * t26 * t45 + 0.106e2 * t68 * t69 - 0.32666666666666666667e0 * t38 * t10 - 0.21333333333333333333e2 * t78 * t69);
  tvrho0 = t7 * t82 + tzk0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->out_dim->vrho + 0] += tvrho0;

  tvrho1 = tvrho0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->out_dim->vrho + 1] += tvrho1;

  t85 = t7 * t7;
  t86 = 0.1e1 / t85;
  t89 = 0.1e1 / t8 / t85;
  t90 = t6 * t89;
  t100 = t5 / t39 / t7;
  t103 = t67 * t67;
  t104 = t103 * t103;
  t105 = t104 * t67;
  t106 = t105 * t77;
  t108 = 0.1e1 / t39 / t85;
  t109 = t43 * t108;
  t112 = t4 * t89;
  t122 = 0.1e1 / t31 / t1 / t3 / t2 / t57 / 0.48e2;
  t123 = t105 * t122;
  t129 = my_piecewise5(t13, 0.10366666666666666667e-1 * t86 + 0.1e-2 * t4 * t90 * t14 - 0.63888888888888888893e-3 * t4 * t90, t21, 0.63266666666666666667e-2 * t86, -0.32444444444444444445e-1 * t26 * t100 + 0.88333333333333333333e1 * t106 * t109 - 0.14133333333333333333e2 * t68 * t112 + 0.10888888888888888889e0 * t38 * t57 - 0.24888888888888888889e2 * t123 * t109 + 0.28444444444444444444e2 * t78 * t112);
  tv2rho20 = t7 * t129 + 0.2e1 * t82;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->out_dim->v2rho2 + 0] += tv2rho20;

  tv2rho21 = tv2rho20;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->out_dim->v2rho2 + 1] += tv2rho21;

  tv2rho22 = tv2rho21;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->out_dim->v2rho2 + 2] += tv2rho22;

  t132 = t85 * t7;
  t133 = 0.1e1 / t132;
  t136 = 0.1e1 / t8 / t132;
  t137 = t6 * t136;
  t148 = t30 * t122;
  t149 = t85 * t85;
  t150 = 0.1e1 / t149;
  t151 = t2 * t150;
  t155 = 0.1e1 / t39 / t132;
  t156 = t43 * t155;
  t159 = t4 * t136;
  t169 = 0.1e1 / t31 / t24 / t36 / t2 / t100 / 0.192e3;
  t170 = t30 * t169;
  t178 = my_piecewise5(t13, -0.20733333333333333334e-1 * t133 - 0.23333333333333333333e-2 * t4 * t137 * t14 + 0.11574074074074074075e-2 * t4 * t137, t21, -0.12653333333333333333e-1 * t133, 0.54074074074074074075e-1 * t26 * t5 * t108 + 0.12366666666666666667e3 * t148 * t151 - 0.35333333333333333332e2 * t106 * t156 + 0.32977777777777777777e2 * t68 * t159 - 0.14518518518518518519e0 * t38 * t90 - 0.448e3 * t170 * t151 + 0.99555555555555555556e2 * t123 * t156 - 0.66370370370370370369e2 * t78 * t159);
  tv3rho30 = t7 * t178 + 0.3e1 * t129;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->out_dim->v3rho3 + 0] += tv3rho30;

  tv3rho31 = tv3rho30;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->out_dim->v3rho3 + 1] += tv3rho31;

  tv3rho32 = tv3rho31;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->out_dim->v3rho3 + 2] += tv3rho32;

  tv3rho33 = tv3rho32;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->out_dim->v3rho3 + 3] += tv3rho33;

  t183 = 0.1e1 / t8 / t149;
  t184 = t6 * t183;
  t197 = t149 * t7;
  t201 = 0.1e1 / t8 / t197 * t1 * t3;
  t205 = t2 / t197;
  t210 = t43 / t39 / t149;
  t213 = t4 * t183;
  t218 = M_PI * M_PI;
  t235 = my_piecewise5(t13, 0.62200000000000000002e-1 * t150 + 0.77777777777777777777e-2 * t4 * t184 * t14 - 0.30802469135802469139e-2 * t4 * t184, t21, 0.37959999999999999999e-1 * t150, -0.14419753086419753087e0 * t26 * t5 * t155 + 0.74200000000000000002e3 * t67 * t169 * t2 * t201 - 0.98933333333333333333e3 * t148 * t205 + 0.15703703703703703703e3 * t106 * t210 - 0.10992592592592592592e3 * t68 * t213 + 0.33876543209876543211e0 * t38 * t137 - 0.14259259259259259259e1 * t67 / t31 * t218 / t86 * t2 * t201 + 0.3584e4 * t170 * t205 - 0.4424691358024691358e3 * t123 * t210 + 0.2212345679012345679e3 * t78 * t213);
  tv4rho40 = t7 * t235 + 0.4e1 * t178;

  if(out->v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    out->v4rho4[ip*p->out_dim->v4rho4 + 0] += tv4rho40;

  tv4rho41 = tv4rho40;

  if(out->v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    out->v4rho4[ip*p->out_dim->v4rho4 + 1] += tv4rho41;

  tv4rho42 = tv4rho41;

  if(out->v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    out->v4rho4[ip*p->out_dim->v4rho4 + 2] += tv4rho42;

  tv4rho43 = tv4rho42;

  if(out->v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    out->v4rho4[ip*p->out_dim->v4rho4 + 3] += tv4rho43;

  tv4rho44 = tv4rho43;

  if(out->v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    out->v4rho4[ip*p->out_dim->v4rho4 + 4] += tv4rho44;

}

#endif

