/*
  This file was generated automatically with scripts/maple2c.py.
  Do not edit this file directly as it can be overwritten!!

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.

  Maple version     : Maple 2020 (X86 64 LINUX)
  Maple source      : ./maple/lda_exc/lda_c_rpa.mpl
  Type of functional: lda_exc
*/

#define maple2c_order 4
#define MAPLE2C_FLAGS (XC_FLAGS_I_HAVE_EXC | XC_FLAGS_I_HAVE_VXC | XC_FLAGS_I_HAVE_FXC | XC_FLAGS_I_HAVE_KXC | XC_FLAGS_I_HAVE_LXC)


#ifndef XC_DONT_COMPILE_EXC
GPU_DEVICE_FUNCTION static inline void
func_exc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t1, t3, t4, t5, t6, t7, t9, t10;
  double t12, t13, t16, t17, tzk0;


  t1 = M_CBRT3;
  t3 = POW_1_3(0.1e1 / M_PI);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = POW_1_3(rho[0]);
  t9 = t6 / t7;
  t10 = t4 * t9;
  t12 = log(t10 / 0.4e1);
  t13 = 0.311e-1 * t12;
  t16 = 0.225e-2 * t4 * t9 * t12;
  t17 = 0.425e-2 * t10;
  tzk0 = t13 - 0.48e-1 + t16 - t17;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim->zk + 0] += tzk0;

}

#endif


#ifndef XC_DONT_COMPILE_VXC
GPU_DEVICE_FUNCTION static inline void
func_vxc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t1, t3, t4, t5, t6, t7, t9, t10;
  double t12, t13, t16, t17, tzk0;

  double t18, t22, t24, t26, tvrho0;


  t1 = M_CBRT3;
  t3 = POW_1_3(0.1e1 / M_PI);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = POW_1_3(rho[0]);
  t9 = t6 / t7;
  t10 = t4 * t9;
  t12 = log(t10 / 0.4e1);
  t13 = 0.311e-1 * t12;
  t16 = 0.225e-2 * t4 * t9 * t12;
  t17 = 0.425e-2 * t10;
  tzk0 = t13 - 0.48e-1 + t16 - t17;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim->zk + 0] += tzk0;

  t18 = 0.1e1 / rho[0];
  t22 = t6 / t7 / rho[0];
  t24 = t4 * t22 * t12;
  t26 = t4 * t22;
  tvrho0 = t13 - 0.48e-1 + t16 - t17 + rho[0] * (-0.10366666666666666667e-1 * t18 - 0.75e-3 * t24 + 0.6666666666666666667e-3 * t26);

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 0] += tvrho0;

}

#endif


#ifndef XC_DONT_COMPILE_FXC
GPU_DEVICE_FUNCTION static inline void
func_fxc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t1, t3, t4, t5, t6, t7, t9, t10;
  double t12, t13, t16, t17, tzk0;

  double t18, t22, t24, t26, tvrho0;

  double t33, t34, t38, t40, t42, tv2rho20;


  t1 = M_CBRT3;
  t3 = POW_1_3(0.1e1 / M_PI);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = POW_1_3(rho[0]);
  t9 = t6 / t7;
  t10 = t4 * t9;
  t12 = log(t10 / 0.4e1);
  t13 = 0.311e-1 * t12;
  t16 = 0.225e-2 * t4 * t9 * t12;
  t17 = 0.425e-2 * t10;
  tzk0 = t13 - 0.48e-1 + t16 - t17;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim->zk + 0] += tzk0;

  t18 = 0.1e1 / rho[0];
  t22 = t6 / t7 / rho[0];
  t24 = t4 * t22 * t12;
  t26 = t4 * t22;
  tvrho0 = t13 - 0.48e-1 + t16 - t17 + rho[0] * (-0.10366666666666666667e-1 * t18 - 0.75e-3 * t24 + 0.6666666666666666667e-3 * t26);

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 0] += tvrho0;

  t33 = rho[0] * rho[0];
  t34 = 0.1e1 / t33;
  t38 = t6 / t7 / t33;
  t40 = t4 * t38 * t12;
  t42 = t4 * t38;
  tv2rho20 = -0.20733333333333333334e-1 * t18 - 0.15e-2 * t24 + 0.13333333333333333334e-2 * t26 + rho[0] * (0.10366666666666666667e-1 * t34 + 0.1e-2 * t40 - 0.63888888888888888893e-3 * t42);

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim->v2rho2 + 0] += tv2rho20;

}

#endif


#ifndef XC_DONT_COMPILE_KXC
GPU_DEVICE_FUNCTION static inline void
func_kxc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t1, t3, t4, t5, t6, t7, t9, t10;
  double t12, t13, t16, t17, tzk0;

  double t18, t22, t24, t26, tvrho0;

  double t33, t34, t38, t40, t42, tv2rho20;

  double t49, t50, t54, t56, t58, tv3rho30;


  t1 = M_CBRT3;
  t3 = POW_1_3(0.1e1 / M_PI);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = POW_1_3(rho[0]);
  t9 = t6 / t7;
  t10 = t4 * t9;
  t12 = log(t10 / 0.4e1);
  t13 = 0.311e-1 * t12;
  t16 = 0.225e-2 * t4 * t9 * t12;
  t17 = 0.425e-2 * t10;
  tzk0 = t13 - 0.48e-1 + t16 - t17;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim->zk + 0] += tzk0;

  t18 = 0.1e1 / rho[0];
  t22 = t6 / t7 / rho[0];
  t24 = t4 * t22 * t12;
  t26 = t4 * t22;
  tvrho0 = t13 - 0.48e-1 + t16 - t17 + rho[0] * (-0.10366666666666666667e-1 * t18 - 0.75e-3 * t24 + 0.6666666666666666667e-3 * t26);

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 0] += tvrho0;

  t33 = rho[0] * rho[0];
  t34 = 0.1e1 / t33;
  t38 = t6 / t7 / t33;
  t40 = t4 * t38 * t12;
  t42 = t4 * t38;
  tv2rho20 = -0.20733333333333333334e-1 * t18 - 0.15e-2 * t24 + 0.13333333333333333334e-2 * t26 + rho[0] * (0.10366666666666666667e-1 * t34 + 0.1e-2 * t40 - 0.63888888888888888893e-3 * t42);

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim->v2rho2 + 0] += tv2rho20;

  t49 = t33 * rho[0];
  t50 = 0.1e1 / t49;
  t54 = t6 / t7 / t49;
  t56 = t4 * t54 * t12;
  t58 = t4 * t54;
  tv3rho30 = 0.31100000000000000001e-1 * t34 + 0.3e-2 * t40 - 0.19166666666666666668e-2 * t42 + rho[0] * (-0.20733333333333333334e-1 * t50 - 0.23333333333333333333e-2 * t56 + 0.11574074074074074075e-2 * t58);

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->dim->v3rho3 + 0] += tv3rho30;

}

#endif


#ifndef XC_DONT_COMPILE_LXC
GPU_DEVICE_FUNCTION static inline void
func_lxc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t1, t3, t4, t5, t6, t7, t9, t10;
  double t12, t13, t16, t17, tzk0;

  double t18, t22, t24, t26, tvrho0;

  double t33, t34, t38, t40, t42, tv2rho20;

  double t49, t50, t54, t56, t58, tv3rho30;

  double t65, t70, tv4rho40;


  t1 = M_CBRT3;
  t3 = POW_1_3(0.1e1 / M_PI);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = POW_1_3(rho[0]);
  t9 = t6 / t7;
  t10 = t4 * t9;
  t12 = log(t10 / 0.4e1);
  t13 = 0.311e-1 * t12;
  t16 = 0.225e-2 * t4 * t9 * t12;
  t17 = 0.425e-2 * t10;
  tzk0 = t13 - 0.48e-1 + t16 - t17;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim->zk + 0] += tzk0;

  t18 = 0.1e1 / rho[0];
  t22 = t6 / t7 / rho[0];
  t24 = t4 * t22 * t12;
  t26 = t4 * t22;
  tvrho0 = t13 - 0.48e-1 + t16 - t17 + rho[0] * (-0.10366666666666666667e-1 * t18 - 0.75e-3 * t24 + 0.6666666666666666667e-3 * t26);

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 0] += tvrho0;

  t33 = rho[0] * rho[0];
  t34 = 0.1e1 / t33;
  t38 = t6 / t7 / t33;
  t40 = t4 * t38 * t12;
  t42 = t4 * t38;
  tv2rho20 = -0.20733333333333333334e-1 * t18 - 0.15e-2 * t24 + 0.13333333333333333334e-2 * t26 + rho[0] * (0.10366666666666666667e-1 * t34 + 0.1e-2 * t40 - 0.63888888888888888893e-3 * t42);

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim->v2rho2 + 0] += tv2rho20;

  t49 = t33 * rho[0];
  t50 = 0.1e1 / t49;
  t54 = t6 / t7 / t49;
  t56 = t4 * t54 * t12;
  t58 = t4 * t54;
  tv3rho30 = 0.31100000000000000001e-1 * t34 + 0.3e-2 * t40 - 0.19166666666666666668e-2 * t42 + rho[0] * (-0.20733333333333333334e-1 * t50 - 0.23333333333333333333e-2 * t56 + 0.11574074074074074075e-2 * t58);

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->dim->v3rho3 + 0] += tv3rho30;

  t65 = t33 * t33;
  t70 = t6 / t7 / t65;
  tv4rho40 = -0.82933333333333333336e-1 * t50 - 0.93333333333333333333e-2 * t56 + 0.462962962962962963e-2 * t58 + rho[0] * (0.62200000000000000002e-1 / t65 + 0.77777777777777777777e-2 * t4 * t70 * t12 - 0.30802469135802469139e-2 * t4 * t70);

  if(out->v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    out->v4rho4[ip*p->dim->v4rho4 + 0] += tv4rho40;

}

#endif


#ifndef XC_DONT_COMPILE_EXC
GPU_DEVICE_FUNCTION static inline void
func_exc_pol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t1, t3, t4, t5, t6, t7, t8, t10;
  double t11, t13, t14, t17, t18, tzk0;


  t1 = M_CBRT3;
  t3 = POW_1_3(0.1e1 / M_PI);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = rho[0] + rho[1];
  t8 = POW_1_3(t7);
  t10 = t6 / t8;
  t11 = t4 * t10;
  t13 = log(t11 / 0.4e1);
  t14 = 0.311e-1 * t13;
  t17 = 0.225e-2 * t4 * t10 * t13;
  t18 = 0.425e-2 * t11;
  tzk0 = t14 - 0.48e-1 + t17 - t18;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim->zk + 0] += tzk0;

}

#endif


#ifndef XC_DONT_COMPILE_VXC
GPU_DEVICE_FUNCTION static inline void
func_vxc_pol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t1, t3, t4, t5, t6, t7, t8, t10;
  double t11, t13, t14, t17, t18, tzk0;

  double t19, t23, t25, t27, tvrho0, tvrho1;


  t1 = M_CBRT3;
  t3 = POW_1_3(0.1e1 / M_PI);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = rho[0] + rho[1];
  t8 = POW_1_3(t7);
  t10 = t6 / t8;
  t11 = t4 * t10;
  t13 = log(t11 / 0.4e1);
  t14 = 0.311e-1 * t13;
  t17 = 0.225e-2 * t4 * t10 * t13;
  t18 = 0.425e-2 * t11;
  tzk0 = t14 - 0.48e-1 + t17 - t18;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim->zk + 0] += tzk0;

  t19 = 0.1e1 / t7;
  t23 = t6 / t8 / t7;
  t25 = t4 * t23 * t13;
  t27 = t4 * t23;
  tvrho0 = t14 - 0.48e-1 + t17 - t18 + t7 * (-0.10366666666666666667e-1 * t19 - 0.75e-3 * t25 + 0.6666666666666666667e-3 * t27);

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 0] += tvrho0;

  tvrho1 = tvrho0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 1] += tvrho1;

}

#endif


#ifndef XC_DONT_COMPILE_FXC
GPU_DEVICE_FUNCTION static inline void
func_fxc_pol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t1, t3, t4, t5, t6, t7, t8, t10;
  double t11, t13, t14, t17, t18, tzk0;

  double t19, t23, t25, t27, tvrho0, tvrho1;

  double t34, t35, t39, t41, t43, tv2rho20, tv2rho21, tv2rho22;


  t1 = M_CBRT3;
  t3 = POW_1_3(0.1e1 / M_PI);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = rho[0] + rho[1];
  t8 = POW_1_3(t7);
  t10 = t6 / t8;
  t11 = t4 * t10;
  t13 = log(t11 / 0.4e1);
  t14 = 0.311e-1 * t13;
  t17 = 0.225e-2 * t4 * t10 * t13;
  t18 = 0.425e-2 * t11;
  tzk0 = t14 - 0.48e-1 + t17 - t18;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim->zk + 0] += tzk0;

  t19 = 0.1e1 / t7;
  t23 = t6 / t8 / t7;
  t25 = t4 * t23 * t13;
  t27 = t4 * t23;
  tvrho0 = t14 - 0.48e-1 + t17 - t18 + t7 * (-0.10366666666666666667e-1 * t19 - 0.75e-3 * t25 + 0.6666666666666666667e-3 * t27);

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 0] += tvrho0;

  tvrho1 = tvrho0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 1] += tvrho1;

  t34 = t7 * t7;
  t35 = 0.1e1 / t34;
  t39 = t6 / t8 / t34;
  t41 = t4 * t39 * t13;
  t43 = t4 * t39;
  tv2rho20 = -0.20733333333333333334e-1 * t19 - 0.15e-2 * t25 + 0.13333333333333333334e-2 * t27 + t7 * (0.10366666666666666667e-1 * t35 + 0.1e-2 * t41 - 0.63888888888888888893e-3 * t43);

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim->v2rho2 + 0] += tv2rho20;

  tv2rho21 = tv2rho20;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim->v2rho2 + 1] += tv2rho21;

  tv2rho22 = tv2rho21;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim->v2rho2 + 2] += tv2rho22;

}

#endif


#ifndef XC_DONT_COMPILE_KXC
GPU_DEVICE_FUNCTION static inline void
func_kxc_pol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t1, t3, t4, t5, t6, t7, t8, t10;
  double t11, t13, t14, t17, t18, tzk0;

  double t19, t23, t25, t27, tvrho0, tvrho1;

  double t34, t35, t39, t41, t43, tv2rho20, tv2rho21, tv2rho22;

  double t50, t51, t55, t57, t59, tv3rho30, tv3rho31, tv3rho32;
  double tv3rho33;


  t1 = M_CBRT3;
  t3 = POW_1_3(0.1e1 / M_PI);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = rho[0] + rho[1];
  t8 = POW_1_3(t7);
  t10 = t6 / t8;
  t11 = t4 * t10;
  t13 = log(t11 / 0.4e1);
  t14 = 0.311e-1 * t13;
  t17 = 0.225e-2 * t4 * t10 * t13;
  t18 = 0.425e-2 * t11;
  tzk0 = t14 - 0.48e-1 + t17 - t18;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim->zk + 0] += tzk0;

  t19 = 0.1e1 / t7;
  t23 = t6 / t8 / t7;
  t25 = t4 * t23 * t13;
  t27 = t4 * t23;
  tvrho0 = t14 - 0.48e-1 + t17 - t18 + t7 * (-0.10366666666666666667e-1 * t19 - 0.75e-3 * t25 + 0.6666666666666666667e-3 * t27);

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 0] += tvrho0;

  tvrho1 = tvrho0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 1] += tvrho1;

  t34 = t7 * t7;
  t35 = 0.1e1 / t34;
  t39 = t6 / t8 / t34;
  t41 = t4 * t39 * t13;
  t43 = t4 * t39;
  tv2rho20 = -0.20733333333333333334e-1 * t19 - 0.15e-2 * t25 + 0.13333333333333333334e-2 * t27 + t7 * (0.10366666666666666667e-1 * t35 + 0.1e-2 * t41 - 0.63888888888888888893e-3 * t43);

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim->v2rho2 + 0] += tv2rho20;

  tv2rho21 = tv2rho20;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim->v2rho2 + 1] += tv2rho21;

  tv2rho22 = tv2rho21;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim->v2rho2 + 2] += tv2rho22;

  t50 = t34 * t7;
  t51 = 0.1e1 / t50;
  t55 = t6 / t8 / t50;
  t57 = t4 * t55 * t13;
  t59 = t4 * t55;
  tv3rho30 = 0.31100000000000000001e-1 * t35 + 0.3e-2 * t41 - 0.19166666666666666668e-2 * t43 + t7 * (-0.20733333333333333334e-1 * t51 - 0.23333333333333333333e-2 * t57 + 0.11574074074074074075e-2 * t59);

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->dim->v3rho3 + 0] += tv3rho30;

  tv3rho31 = tv3rho30;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->dim->v3rho3 + 1] += tv3rho31;

  tv3rho32 = tv3rho31;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->dim->v3rho3 + 2] += tv3rho32;

  tv3rho33 = tv3rho32;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->dim->v3rho3 + 3] += tv3rho33;

}

#endif


#ifndef XC_DONT_COMPILE_LXC
GPU_DEVICE_FUNCTION static inline void
func_lxc_pol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t1, t3, t4, t5, t6, t7, t8, t10;
  double t11, t13, t14, t17, t18, tzk0;

  double t19, t23, t25, t27, tvrho0, tvrho1;

  double t34, t35, t39, t41, t43, tv2rho20, tv2rho21, tv2rho22;

  double t50, t51, t55, t57, t59, tv3rho30, tv3rho31, tv3rho32;
  double tv3rho33;

  double t66, t71, tv4rho40, tv4rho41, tv4rho42, tv4rho43, tv4rho44;


  t1 = M_CBRT3;
  t3 = POW_1_3(0.1e1 / M_PI);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = rho[0] + rho[1];
  t8 = POW_1_3(t7);
  t10 = t6 / t8;
  t11 = t4 * t10;
  t13 = log(t11 / 0.4e1);
  t14 = 0.311e-1 * t13;
  t17 = 0.225e-2 * t4 * t10 * t13;
  t18 = 0.425e-2 * t11;
  tzk0 = t14 - 0.48e-1 + t17 - t18;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim->zk + 0] += tzk0;

  t19 = 0.1e1 / t7;
  t23 = t6 / t8 / t7;
  t25 = t4 * t23 * t13;
  t27 = t4 * t23;
  tvrho0 = t14 - 0.48e-1 + t17 - t18 + t7 * (-0.10366666666666666667e-1 * t19 - 0.75e-3 * t25 + 0.6666666666666666667e-3 * t27);

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 0] += tvrho0;

  tvrho1 = tvrho0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 1] += tvrho1;

  t34 = t7 * t7;
  t35 = 0.1e1 / t34;
  t39 = t6 / t8 / t34;
  t41 = t4 * t39 * t13;
  t43 = t4 * t39;
  tv2rho20 = -0.20733333333333333334e-1 * t19 - 0.15e-2 * t25 + 0.13333333333333333334e-2 * t27 + t7 * (0.10366666666666666667e-1 * t35 + 0.1e-2 * t41 - 0.63888888888888888893e-3 * t43);

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim->v2rho2 + 0] += tv2rho20;

  tv2rho21 = tv2rho20;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim->v2rho2 + 1] += tv2rho21;

  tv2rho22 = tv2rho21;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim->v2rho2 + 2] += tv2rho22;

  t50 = t34 * t7;
  t51 = 0.1e1 / t50;
  t55 = t6 / t8 / t50;
  t57 = t4 * t55 * t13;
  t59 = t4 * t55;
  tv3rho30 = 0.31100000000000000001e-1 * t35 + 0.3e-2 * t41 - 0.19166666666666666668e-2 * t43 + t7 * (-0.20733333333333333334e-1 * t51 - 0.23333333333333333333e-2 * t57 + 0.11574074074074074075e-2 * t59);

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->dim->v3rho3 + 0] += tv3rho30;

  tv3rho31 = tv3rho30;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->dim->v3rho3 + 1] += tv3rho31;

  tv3rho32 = tv3rho31;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->dim->v3rho3 + 2] += tv3rho32;

  tv3rho33 = tv3rho32;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->dim->v3rho3 + 3] += tv3rho33;

  t66 = t34 * t34;
  t71 = t6 / t8 / t66;
  tv4rho40 = -0.82933333333333333336e-1 * t51 - 0.93333333333333333333e-2 * t57 + 0.462962962962962963e-2 * t59 + t7 * (0.62200000000000000002e-1 / t66 + 0.77777777777777777777e-2 * t4 * t71 * t13 - 0.30802469135802469139e-2 * t4 * t71);

  if(out->v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    out->v4rho4[ip*p->dim->v4rho4 + 0] += tv4rho40;

  tv4rho41 = tv4rho40;

  if(out->v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    out->v4rho4[ip*p->dim->v4rho4 + 1] += tv4rho41;

  tv4rho42 = tv4rho41;

  if(out->v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    out->v4rho4[ip*p->dim->v4rho4 + 2] += tv4rho42;

  tv4rho43 = tv4rho42;

  if(out->v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    out->v4rho4[ip*p->dim->v4rho4 + 3] += tv4rho43;

  tv4rho44 = tv4rho43;

  if(out->v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    out->v4rho4[ip*p->dim->v4rho4 + 4] += tv4rho44;

}

#endif

