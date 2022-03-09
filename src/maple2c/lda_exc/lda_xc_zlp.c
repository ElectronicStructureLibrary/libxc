/*
  This file was generated automatically with scripts/maple2c.py.
  Do not edit this file directly as it can be overwritten!!

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.

  Maple version     : Maple 2020 (X86 64 LINUX)
  Maple source      : ./maple/lda_exc/lda_xc_zlp.mpl
  Type of functional: lda_exc
*/

#define maple2c_order 4
#define MAPLE2C_FLAGS (XC_FLAGS_I_HAVE_EXC | XC_FLAGS_I_HAVE_VXC | XC_FLAGS_I_HAVE_FXC | XC_FLAGS_I_HAVE_KXC | XC_FLAGS_I_HAVE_LXC)


#ifndef XC_DONT_COMPILE_EXC
GPU_DEVICE_FUNCTION static inline void
func_exc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_lda_out_params *out)
{
  double t1, t4, t5, t8, t9, tzk0;


  t1 = POW_1_3(rho[0]);
  t4 = 0.1e1 + 0.10555627099250339363e3 / t1;
  t5 = log(t4);
  t8 = 0.1e1 - 0.947362e-2 * t5 * t1;
  t9 = t8 * t1;
  tzk0 = -0.93222e0 * t9;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim->zk + 0] += tzk0;

}

#endif


#ifndef XC_DONT_COMPILE_VXC
GPU_DEVICE_FUNCTION static inline void
func_vxc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_lda_out_params *out)
{
  double t1, t4, t5, t8, t9, tzk0;

  double t12, t14, t17, t18, t21, tvrho0;


  t1 = POW_1_3(rho[0]);
  t4 = 0.1e1 + 0.10555627099250339363e3 / t1;
  t5 = log(t4);
  t8 = 0.1e1 - 0.947362e-2 * t5 * t1;
  t9 = t8 * t1;
  tzk0 = -0.93222e0 * t9;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim->zk + 0] += tzk0;

  t12 = t1 * rho[0];
  t14 = 0.1e1 / t4;
  t17 = t1 * t1;
  t18 = 0.1e1 / t17;
  t21 = 0.33333333333333333332e0 / rho[0] * t14 - 0.31578733333333333333e-2 * t5 * t18;
  tvrho0 = -0.124296e1 * t9 - 0.93222e0 * t12 * t21;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 0] += tvrho0;

}

#endif


#ifndef XC_DONT_COMPILE_FXC
GPU_DEVICE_FUNCTION static inline void
func_fxc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_lda_out_params *out)
{
  double t1, t4, t5, t8, t9, tzk0;

  double t12, t14, t17, t18, t21, tvrho0;

  double t28, t34, t35, t39, t42, tv2rho20;


  t1 = POW_1_3(rho[0]);
  t4 = 0.1e1 + 0.10555627099250339363e3 / t1;
  t5 = log(t4);
  t8 = 0.1e1 - 0.947362e-2 * t5 * t1;
  t9 = t8 * t1;
  tzk0 = -0.93222e0 * t9;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim->zk + 0] += tzk0;

  t12 = t1 * rho[0];
  t14 = 0.1e1 / t4;
  t17 = t1 * t1;
  t18 = 0.1e1 / t17;
  t21 = 0.33333333333333333332e0 / rho[0] * t14 - 0.31578733333333333333e-2 * t5 * t18;
  tvrho0 = -0.124296e1 * t9 - 0.93222e0 * t12 * t21;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 0] += tvrho0;

  t28 = rho[0] * rho[0];
  t34 = t4 * t4;
  t35 = 0.1e1 / t34;
  t39 = 0.1e1 / t17 / rho[0];
  t42 = -0.22222222222222222221e0 / t28 * t14 + 0.11728474554722599292e2 / t1 / t28 * t35 + 0.21052488888888888889e-2 * t5 * t39;
  tv2rho20 = -0.248592e1 * t21 * t1 - 0.41432e0 * t8 * t18 - 0.93222e0 * t12 * t42;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim->v2rho2 + 0] += tv2rho20;

}

#endif


#ifndef XC_DONT_COMPILE_KXC
GPU_DEVICE_FUNCTION static inline void
func_kxc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_lda_out_params *out)
{
  double t1, t4, t5, t8, t9, tzk0;

  double t12, t14, t17, t18, t21, tvrho0;

  double t28, t34, t35, t39, t42, tv2rho20;

  double t51, t60, t62, t66, t69, tv3rho30;


  t1 = POW_1_3(rho[0]);
  t4 = 0.1e1 + 0.10555627099250339363e3 / t1;
  t5 = log(t4);
  t8 = 0.1e1 - 0.947362e-2 * t5 * t1;
  t9 = t8 * t1;
  tzk0 = -0.93222e0 * t9;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim->zk + 0] += tzk0;

  t12 = t1 * rho[0];
  t14 = 0.1e1 / t4;
  t17 = t1 * t1;
  t18 = 0.1e1 / t17;
  t21 = 0.33333333333333333332e0 / rho[0] * t14 - 0.31578733333333333333e-2 * t5 * t18;
  tvrho0 = -0.124296e1 * t9 - 0.93222e0 * t12 * t21;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 0] += tvrho0;

  t28 = rho[0] * rho[0];
  t34 = t4 * t4;
  t35 = 0.1e1 / t34;
  t39 = 0.1e1 / t17 / rho[0];
  t42 = -0.22222222222222222221e0 / t28 * t14 + 0.11728474554722599292e2 / t1 / t28 * t35 + 0.21052488888888888889e-2 * t5 * t39;
  tv2rho20 = -0.248592e1 * t21 * t1 - 0.41432e0 * t8 * t18 - 0.93222e0 * t12 * t42;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim->v2rho2 + 0] += tv2rho20;

  t51 = t28 * rho[0];
  t60 = 0.1e1 / t17 / t51;
  t62 = 0.1e1 / t34 / t4;
  t66 = 0.1e1 / t17 / t28;
  t69 = 0.37037037037037037035e0 / t51 * t14 - 0.35185423664167797876e2 / t1 / t51 * t35 + 0.82534269228465284243e3 * t60 * t62 - 0.35087481481481481482e-2 * t5 * t66;
  tv3rho30 = -0.372888e1 * t42 * t1 - 0.124296e1 * t21 * t18 + 0.27621333333333333333e0 * t8 * t39 - 0.93222e0 * t12 * t69;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->dim->v3rho3 + 0] += tv3rho30;

}

#endif


#ifndef XC_DONT_COMPILE_LXC
GPU_DEVICE_FUNCTION static inline void
func_lxc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_lda_out_params *out)
{
  double t1, t4, t5, t8, t9, tzk0;

  double t12, t14, t17, t18, t21, tvrho0;

  double t28, t34, t35, t39, t42, tv2rho20;

  double t51, t60, t62, t66, t69, tv3rho30;

  double t80, t94, tv4rho40;


  t1 = POW_1_3(rho[0]);
  t4 = 0.1e1 + 0.10555627099250339363e3 / t1;
  t5 = log(t4);
  t8 = 0.1e1 - 0.947362e-2 * t5 * t1;
  t9 = t8 * t1;
  tzk0 = -0.93222e0 * t9;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim->zk + 0] += tzk0;

  t12 = t1 * rho[0];
  t14 = 0.1e1 / t4;
  t17 = t1 * t1;
  t18 = 0.1e1 / t17;
  t21 = 0.33333333333333333332e0 / rho[0] * t14 - 0.31578733333333333333e-2 * t5 * t18;
  tvrho0 = -0.124296e1 * t9 - 0.93222e0 * t12 * t21;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 0] += tvrho0;

  t28 = rho[0] * rho[0];
  t34 = t4 * t4;
  t35 = 0.1e1 / t34;
  t39 = 0.1e1 / t17 / rho[0];
  t42 = -0.22222222222222222221e0 / t28 * t14 + 0.11728474554722599292e2 / t1 / t28 * t35 + 0.21052488888888888889e-2 * t5 * t39;
  tv2rho20 = -0.248592e1 * t21 * t1 - 0.41432e0 * t8 * t18 - 0.93222e0 * t12 * t42;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim->v2rho2 + 0] += tv2rho20;

  t51 = t28 * rho[0];
  t60 = 0.1e1 / t17 / t51;
  t62 = 0.1e1 / t34 / t4;
  t66 = 0.1e1 / t17 / t28;
  t69 = 0.37037037037037037035e0 / t51 * t14 - 0.35185423664167797876e2 / t1 / t51 * t35 + 0.82534269228465284243e3 * t60 * t62 - 0.35087481481481481482e-2 * t5 * t66;
  tv3rho30 = -0.372888e1 * t42 * t1 - 0.124296e1 * t21 * t18 + 0.27621333333333333333e0 * t8 * t39 - 0.93222e0 * t12 * t69;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->dim->v3rho3 + 0] += tv3rho30;

  t80 = t28 * t28;
  t94 = t34 * t34;
  tv4rho40 = -0.497184e1 * t69 * t1 - 0.248592e1 * t42 * t18 + 0.11048533333333333333e1 * t21 * t39 - 0.46035555555555555555e0 * t8 * t66 - 0.93222e0 * t12 * (-0.98765432098765432088e0 / t80 * t14 + 0.13031638394136221436e3 / t1 / t80 * t35 - 0.55022846152310189495e4 / t17 / t80 * t62 + 0.87120096888481155292e5 / t80 / rho[0] / t94 + 0.93566617283950617285e-2 * t5 * t60);

  if(out->v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    out->v4rho4[ip*p->dim->v4rho4 + 0] += tv4rho40;

}

#endif


#ifndef XC_DONT_COMPILE_EXC
GPU_DEVICE_FUNCTION static inline void
func_exc_pol(const xc_func_type *p, size_t ip, const double *rho, xc_lda_out_params *out)
{
  double t1, t2, t5, t6, t9, t10, tzk0;


  t1 = rho[0] + rho[1];
  t2 = POW_1_3(t1);
  t5 = 0.1e1 + 0.10555627099250339363e3 / t2;
  t6 = log(t5);
  t9 = 0.1e1 - 0.947362e-2 * t6 * t2;
  t10 = t9 * t2;
  tzk0 = -0.93222e0 * t10;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim->zk + 0] += tzk0;

}

#endif


#ifndef XC_DONT_COMPILE_VXC
GPU_DEVICE_FUNCTION static inline void
func_vxc_pol(const xc_func_type *p, size_t ip, const double *rho, xc_lda_out_params *out)
{
  double t1, t2, t5, t6, t9, t10, tzk0;

  double t13, t15, t18, t19, t22, tvrho0, tvrho1;


  t1 = rho[0] + rho[1];
  t2 = POW_1_3(t1);
  t5 = 0.1e1 + 0.10555627099250339363e3 / t2;
  t6 = log(t5);
  t9 = 0.1e1 - 0.947362e-2 * t6 * t2;
  t10 = t9 * t2;
  tzk0 = -0.93222e0 * t10;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim->zk + 0] += tzk0;

  t13 = t2 * t1;
  t15 = 0.1e1 / t5;
  t18 = t2 * t2;
  t19 = 0.1e1 / t18;
  t22 = 0.33333333333333333332e0 / t1 * t15 - 0.31578733333333333333e-2 * t6 * t19;
  tvrho0 = -0.124296e1 * t10 - 0.93222e0 * t13 * t22;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 0] += tvrho0;

  tvrho1 = tvrho0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 1] += tvrho1;

}

#endif


#ifndef XC_DONT_COMPILE_FXC
GPU_DEVICE_FUNCTION static inline void
func_fxc_pol(const xc_func_type *p, size_t ip, const double *rho, xc_lda_out_params *out)
{
  double t1, t2, t5, t6, t9, t10, tzk0;

  double t13, t15, t18, t19, t22, tvrho0, tvrho1;

  double t29, t35, t36, t40, t43, tv2rho20, tv2rho21, tv2rho22;


  t1 = rho[0] + rho[1];
  t2 = POW_1_3(t1);
  t5 = 0.1e1 + 0.10555627099250339363e3 / t2;
  t6 = log(t5);
  t9 = 0.1e1 - 0.947362e-2 * t6 * t2;
  t10 = t9 * t2;
  tzk0 = -0.93222e0 * t10;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim->zk + 0] += tzk0;

  t13 = t2 * t1;
  t15 = 0.1e1 / t5;
  t18 = t2 * t2;
  t19 = 0.1e1 / t18;
  t22 = 0.33333333333333333332e0 / t1 * t15 - 0.31578733333333333333e-2 * t6 * t19;
  tvrho0 = -0.124296e1 * t10 - 0.93222e0 * t13 * t22;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 0] += tvrho0;

  tvrho1 = tvrho0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 1] += tvrho1;

  t29 = t1 * t1;
  t35 = t5 * t5;
  t36 = 0.1e1 / t35;
  t40 = 0.1e1 / t18 / t1;
  t43 = -0.22222222222222222221e0 / t29 * t15 + 0.11728474554722599292e2 / t2 / t29 * t36 + 0.21052488888888888889e-2 * t6 * t40;
  tv2rho20 = -0.248592e1 * t22 * t2 - 0.41432e0 * t9 * t19 - 0.93222e0 * t13 * t43;

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
func_kxc_pol(const xc_func_type *p, size_t ip, const double *rho, xc_lda_out_params *out)
{
  double t1, t2, t5, t6, t9, t10, tzk0;

  double t13, t15, t18, t19, t22, tvrho0, tvrho1;

  double t29, t35, t36, t40, t43, tv2rho20, tv2rho21, tv2rho22;

  double t52, t61, t63, t67, t70, tv3rho30, tv3rho31, tv3rho32;
  double tv3rho33;


  t1 = rho[0] + rho[1];
  t2 = POW_1_3(t1);
  t5 = 0.1e1 + 0.10555627099250339363e3 / t2;
  t6 = log(t5);
  t9 = 0.1e1 - 0.947362e-2 * t6 * t2;
  t10 = t9 * t2;
  tzk0 = -0.93222e0 * t10;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim->zk + 0] += tzk0;

  t13 = t2 * t1;
  t15 = 0.1e1 / t5;
  t18 = t2 * t2;
  t19 = 0.1e1 / t18;
  t22 = 0.33333333333333333332e0 / t1 * t15 - 0.31578733333333333333e-2 * t6 * t19;
  tvrho0 = -0.124296e1 * t10 - 0.93222e0 * t13 * t22;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 0] += tvrho0;

  tvrho1 = tvrho0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 1] += tvrho1;

  t29 = t1 * t1;
  t35 = t5 * t5;
  t36 = 0.1e1 / t35;
  t40 = 0.1e1 / t18 / t1;
  t43 = -0.22222222222222222221e0 / t29 * t15 + 0.11728474554722599292e2 / t2 / t29 * t36 + 0.21052488888888888889e-2 * t6 * t40;
  tv2rho20 = -0.248592e1 * t22 * t2 - 0.41432e0 * t9 * t19 - 0.93222e0 * t13 * t43;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim->v2rho2 + 0] += tv2rho20;

  tv2rho21 = tv2rho20;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim->v2rho2 + 1] += tv2rho21;

  tv2rho22 = tv2rho21;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim->v2rho2 + 2] += tv2rho22;

  t52 = t29 * t1;
  t61 = 0.1e1 / t18 / t52;
  t63 = 0.1e1 / t35 / t5;
  t67 = 0.1e1 / t18 / t29;
  t70 = 0.37037037037037037035e0 / t52 * t15 - 0.35185423664167797876e2 / t2 / t52 * t36 + 0.82534269228465284243e3 * t61 * t63 - 0.35087481481481481482e-2 * t6 * t67;
  tv3rho30 = -0.372888e1 * t43 * t2 - 0.124296e1 * t22 * t19 + 0.27621333333333333333e0 * t9 * t40 - 0.93222e0 * t13 * t70;

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
func_lxc_pol(const xc_func_type *p, size_t ip, const double *rho, xc_lda_out_params *out)
{
  double t1, t2, t5, t6, t9, t10, tzk0;

  double t13, t15, t18, t19, t22, tvrho0, tvrho1;

  double t29, t35, t36, t40, t43, tv2rho20, tv2rho21, tv2rho22;

  double t52, t61, t63, t67, t70, tv3rho30, tv3rho31, tv3rho32;
  double tv3rho33;

  double t81, t95, tv4rho40, tv4rho41, tv4rho42, tv4rho43, tv4rho44;


  t1 = rho[0] + rho[1];
  t2 = POW_1_3(t1);
  t5 = 0.1e1 + 0.10555627099250339363e3 / t2;
  t6 = log(t5);
  t9 = 0.1e1 - 0.947362e-2 * t6 * t2;
  t10 = t9 * t2;
  tzk0 = -0.93222e0 * t10;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim->zk + 0] += tzk0;

  t13 = t2 * t1;
  t15 = 0.1e1 / t5;
  t18 = t2 * t2;
  t19 = 0.1e1 / t18;
  t22 = 0.33333333333333333332e0 / t1 * t15 - 0.31578733333333333333e-2 * t6 * t19;
  tvrho0 = -0.124296e1 * t10 - 0.93222e0 * t13 * t22;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 0] += tvrho0;

  tvrho1 = tvrho0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim->vrho + 1] += tvrho1;

  t29 = t1 * t1;
  t35 = t5 * t5;
  t36 = 0.1e1 / t35;
  t40 = 0.1e1 / t18 / t1;
  t43 = -0.22222222222222222221e0 / t29 * t15 + 0.11728474554722599292e2 / t2 / t29 * t36 + 0.21052488888888888889e-2 * t6 * t40;
  tv2rho20 = -0.248592e1 * t22 * t2 - 0.41432e0 * t9 * t19 - 0.93222e0 * t13 * t43;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim->v2rho2 + 0] += tv2rho20;

  tv2rho21 = tv2rho20;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim->v2rho2 + 1] += tv2rho21;

  tv2rho22 = tv2rho21;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim->v2rho2 + 2] += tv2rho22;

  t52 = t29 * t1;
  t61 = 0.1e1 / t18 / t52;
  t63 = 0.1e1 / t35 / t5;
  t67 = 0.1e1 / t18 / t29;
  t70 = 0.37037037037037037035e0 / t52 * t15 - 0.35185423664167797876e2 / t2 / t52 * t36 + 0.82534269228465284243e3 * t61 * t63 - 0.35087481481481481482e-2 * t6 * t67;
  tv3rho30 = -0.372888e1 * t43 * t2 - 0.124296e1 * t22 * t19 + 0.27621333333333333333e0 * t9 * t40 - 0.93222e0 * t13 * t70;

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

  t81 = t29 * t29;
  t95 = t35 * t35;
  tv4rho40 = -0.497184e1 * t70 * t2 - 0.248592e1 * t43 * t19 + 0.11048533333333333333e1 * t22 * t40 - 0.46035555555555555555e0 * t9 * t67 - 0.93222e0 * t13 * (-0.98765432098765432088e0 / t81 * t15 + 0.13031638394136221436e3 / t2 / t81 * t36 - 0.55022846152310189495e4 / t18 / t81 * t63 + 0.87120096888481155292e5 / t81 / t1 / t95 + 0.93566617283950617285e-2 * t6 * t61);

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

