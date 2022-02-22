/*
  This file was generated automatically with ./scripts/maple2c.py.
  Do not edit this file directly as it can be overwritten!!

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.

  Maple version     : Maple 2016 (X86 64 LINUX)
  Maple source      : ./maple/lda_exc/lda_xc_1d_ehwlrg.mpl
  Type of functional: lda_exc
*/

#define maple2c_order 4
#define MAPLE2C_FLAGS (XC_FLAGS_I_HAVE_EXC | XC_FLAGS_I_HAVE_VXC | XC_FLAGS_I_HAVE_FXC | XC_FLAGS_I_HAVE_KXC | XC_FLAGS_I_HAVE_LXC)


#ifndef XC_DONT_COMPILE_EXC
GPU_DEVICE_FUNCTION static inline void
func_exc_pol(const xc_func_type *p, size_t ip, const double *rho, xc_lda_out_params *out)
{
  double t1, t3, t5, t6, tzk0;

  lda_xc_1d_ehwlrg_params *params;

  assert(p->params != NULL);
  params = (lda_xc_1d_ehwlrg_params * )(p->params);

  t1 = rho[0] + rho[1];
  t3 = t1 * t1;
  t5 = params->a2 * t1 + params->a3 * t3 + params->a1;
  t6 = pow(t1, params->alpha);
  tzk0 = t5 * t6;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim.zk + 0] += tzk0;

}

#endif


#ifndef XC_DONT_COMPILE_VXC
GPU_DEVICE_FUNCTION static inline void
func_vxc_pol(const xc_func_type *p, size_t ip, const double *rho, xc_lda_out_params *out)
{
  double t1, t3, t5, t6, tzk0;

  double t7, t9, tvrho0, tvrho1;

  lda_xc_1d_ehwlrg_params *params;

  assert(p->params != NULL);
  params = (lda_xc_1d_ehwlrg_params * )(p->params);

  t1 = rho[0] + rho[1];
  t3 = t1 * t1;
  t5 = params->a2 * t1 + params->a3 * t3 + params->a1;
  t6 = pow(t1, params->alpha);
  tzk0 = t5 * t6;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim.zk + 0] += tzk0;

  t7 = params->a3 * t1;
  t9 = params->a2 + 0.2e1 * t7;
  tvrho0 = t1 * t9 * t6 + t5 * t6 * params->alpha + tzk0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim.vrho + 0] += tvrho0;

  tvrho1 = tvrho0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim.vrho + 1] += tvrho1;

}

#endif


#ifndef XC_DONT_COMPILE_FXC
GPU_DEVICE_FUNCTION static inline void
func_fxc_pol(const xc_func_type *p, size_t ip, const double *rho, xc_lda_out_params *out)
{
  double t1, t3, t5, t6, tzk0;

  double t7, t9, tvrho0, tvrho1;

  double t14, t16, t17, t23, t24, tv2rho20, tv2rho21, tv2rho22;

  lda_xc_1d_ehwlrg_params *params;

  assert(p->params != NULL);
  params = (lda_xc_1d_ehwlrg_params * )(p->params);

  t1 = rho[0] + rho[1];
  t3 = t1 * t1;
  t5 = params->a2 * t1 + params->a3 * t3 + params->a1;
  t6 = pow(t1, params->alpha);
  tzk0 = t5 * t6;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim.zk + 0] += tzk0;

  t7 = params->a3 * t1;
  t9 = params->a2 + 0.2e1 * t7;
  tvrho0 = t1 * t9 * t6 + t5 * t6 * params->alpha + tzk0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim.vrho + 0] += tvrho0;

  tvrho1 = tvrho0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim.vrho + 1] += tvrho1;

  t14 = t9 * t6;
  t16 = 0.1e1 / t1;
  t17 = params->alpha * t16;
  t23 = params->alpha * params->alpha;
  t24 = t23 * t16;
  tv2rho20 = 0.2e1 * t14 * params->alpha + tzk0 * t17 + tzk0 * t24 + 0.2e1 * t7 * t6 + 0.2e1 * t14;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim.v2rho2 + 0] += tv2rho20;

  tv2rho21 = tv2rho20;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim.v2rho2 + 1] += tv2rho21;

  tv2rho22 = tv2rho21;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim.v2rho2 + 2] += tv2rho22;

}

#endif


#ifndef XC_DONT_COMPILE_KXC
GPU_DEVICE_FUNCTION static inline void
func_kxc_pol(const xc_func_type *p, size_t ip, const double *rho, xc_lda_out_params *out)
{
  double t1, t3, t5, t6, tzk0;

  double t7, t9, tvrho0, tvrho1;

  double t14, t16, t17, t23, t24, tv2rho20, tv2rho21, tv2rho22;

  double t26, t30, t31, t37, t38, tv3rho30, tv3rho31, tv3rho32;
  double tv3rho33;

  lda_xc_1d_ehwlrg_params *params;

  assert(p->params != NULL);
  params = (lda_xc_1d_ehwlrg_params * )(p->params);

  t1 = rho[0] + rho[1];
  t3 = t1 * t1;
  t5 = params->a2 * t1 + params->a3 * t3 + params->a1;
  t6 = pow(t1, params->alpha);
  tzk0 = t5 * t6;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim.zk + 0] += tzk0;

  t7 = params->a3 * t1;
  t9 = params->a2 + 0.2e1 * t7;
  tvrho0 = t1 * t9 * t6 + t5 * t6 * params->alpha + tzk0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim.vrho + 0] += tvrho0;

  tvrho1 = tvrho0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim.vrho + 1] += tvrho1;

  t14 = t9 * t6;
  t16 = 0.1e1 / t1;
  t17 = params->alpha * t16;
  t23 = params->alpha * params->alpha;
  t24 = t23 * t16;
  tv2rho20 = 0.2e1 * t14 * params->alpha + tzk0 * t17 + tzk0 * t24 + 0.2e1 * t7 * t6 + 0.2e1 * t14;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim.v2rho2 + 0] += tv2rho20;

  tv2rho21 = tv2rho20;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim.v2rho2 + 1] += tv2rho21;

  tv2rho22 = tv2rho21;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim.v2rho2 + 2] += tv2rho22;

  t26 = params->a3 * t6;
  t30 = 0.1e1 / t3;
  t31 = params->alpha * t30;
  t37 = t23 * params->alpha;
  t38 = t37 * t30;
  tv3rho30 = 0.3e1 * t14 * t17 + 0.3e1 * t14 * t24 + 0.6e1 * t26 * params->alpha - tzk0 * t31 + tzk0 * t38 + 0.6e1 * t26;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->dim.v3rho3 + 0] += tv3rho30;

  tv3rho31 = tv3rho30;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->dim.v3rho3 + 1] += tv3rho31;

  tv3rho32 = tv3rho31;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->dim.v3rho3 + 2] += tv3rho32;

  tv3rho33 = tv3rho32;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->dim.v3rho3 + 3] += tv3rho33;

}

#endif


#ifndef XC_DONT_COMPILE_LXC
GPU_DEVICE_FUNCTION static inline void
func_lxc_pol(const xc_func_type *p, size_t ip, const double *rho, xc_lda_out_params *out)
{
  double t1, t3, t5, t6, tzk0;

  double t7, t9, tvrho0, tvrho1;

  double t14, t16, t17, t23, t24, tv2rho20, tv2rho21, tv2rho22;

  double t26, t30, t31, t37, t38, tv3rho30, tv3rho31, tv3rho32;
  double tv3rho33;

  double t45, t55, tv4rho40, tv4rho41, tv4rho42, tv4rho43, tv4rho44;

  lda_xc_1d_ehwlrg_params *params;

  assert(p->params != NULL);
  params = (lda_xc_1d_ehwlrg_params * )(p->params);

  t1 = rho[0] + rho[1];
  t3 = t1 * t1;
  t5 = params->a2 * t1 + params->a3 * t3 + params->a1;
  t6 = pow(t1, params->alpha);
  tzk0 = t5 * t6;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim.zk + 0] += tzk0;

  t7 = params->a3 * t1;
  t9 = params->a2 + 0.2e1 * t7;
  tvrho0 = t1 * t9 * t6 + t5 * t6 * params->alpha + tzk0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim.vrho + 0] += tvrho0;

  tvrho1 = tvrho0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim.vrho + 1] += tvrho1;

  t14 = t9 * t6;
  t16 = 0.1e1 / t1;
  t17 = params->alpha * t16;
  t23 = params->alpha * params->alpha;
  t24 = t23 * t16;
  tv2rho20 = 0.2e1 * t14 * params->alpha + tzk0 * t17 + tzk0 * t24 + 0.2e1 * t7 * t6 + 0.2e1 * t14;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim.v2rho2 + 0] += tv2rho20;

  tv2rho21 = tv2rho20;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim.v2rho2 + 1] += tv2rho21;

  tv2rho22 = tv2rho21;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim.v2rho2 + 2] += tv2rho22;

  t26 = params->a3 * t6;
  t30 = 0.1e1 / t3;
  t31 = params->alpha * t30;
  t37 = t23 * params->alpha;
  t38 = t37 * t30;
  tv3rho30 = 0.3e1 * t14 * t17 + 0.3e1 * t14 * t24 + 0.6e1 * t26 * params->alpha - tzk0 * t31 + tzk0 * t38 + 0.6e1 * t26;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->dim.v3rho3 + 0] += tv3rho30;

  tv3rho31 = tv3rho30;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->dim.v3rho3 + 1] += tv3rho31;

  tv3rho32 = tv3rho31;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->dim.v3rho3 + 2] += tv3rho32;

  tv3rho33 = tv3rho32;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->dim.v3rho3 + 3] += tv3rho33;

  t45 = 0.1e1 / t3 / t1;
  t55 = t23 * t23;
  tv4rho40 = -tzk0 * t23 * t45 - 0.2e1 * tzk0 * t37 * t45 + tzk0 * t55 * t45 + 0.2e1 * tzk0 * params->alpha * t45 - 0.4e1 * t14 * t31 + 0.4e1 * t14 * t38 + 0.12e2 * t26 * t17 + 0.12e2 * t26 * t24;

  if(out->v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    out->v4rho4[ip*p->dim.v4rho4 + 0] += tv4rho40;

  tv4rho41 = tv4rho40;

  if(out->v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    out->v4rho4[ip*p->dim.v4rho4 + 1] += tv4rho41;

  tv4rho42 = tv4rho41;

  if(out->v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    out->v4rho4[ip*p->dim.v4rho4 + 2] += tv4rho42;

  tv4rho43 = tv4rho42;

  if(out->v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    out->v4rho4[ip*p->dim.v4rho4 + 3] += tv4rho43;

  tv4rho44 = tv4rho43;

  if(out->v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    out->v4rho4[ip*p->dim.v4rho4 + 4] += tv4rho44;

}

#endif


#ifndef XC_DONT_COMPILE_EXC
GPU_DEVICE_FUNCTION static inline void
func_exc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_lda_out_params *out)
{
  double t1, t4, t5, tzk0;

  lda_xc_1d_ehwlrg_params *params;

  assert(p->params != NULL);
  params = (lda_xc_1d_ehwlrg_params * )(p->params);

  t1 = rho[0] * rho[0];
  t4 = rho[0] * params->a2 + t1 * params->a3 + params->a1;
  t5 = pow(rho[0], params->alpha);
  tzk0 = t4 * t5;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim.zk + 0] += tzk0;

}

#endif


#ifndef XC_DONT_COMPILE_VXC
GPU_DEVICE_FUNCTION static inline void
func_vxc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_lda_out_params *out)
{
  double t1, t4, t5, tzk0;

  double t6, t8, tvrho0;

  lda_xc_1d_ehwlrg_params *params;

  assert(p->params != NULL);
  params = (lda_xc_1d_ehwlrg_params * )(p->params);

  t1 = rho[0] * rho[0];
  t4 = rho[0] * params->a2 + t1 * params->a3 + params->a1;
  t5 = pow(rho[0], params->alpha);
  tzk0 = t4 * t5;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim.zk + 0] += tzk0;

  t6 = rho[0] * params->a3;
  t8 = 0.2e1 * t6 + params->a2;
  tvrho0 = rho[0] * t8 * t5 + t4 * t5 * params->alpha + tzk0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim.vrho + 0] += tvrho0;

}

#endif


#ifndef XC_DONT_COMPILE_FXC
GPU_DEVICE_FUNCTION static inline void
func_fxc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_lda_out_params *out)
{
  double t1, t4, t5, tzk0;

  double t6, t8, tvrho0;

  double t13, t15, t16, t22, t23, tv2rho20;

  lda_xc_1d_ehwlrg_params *params;

  assert(p->params != NULL);
  params = (lda_xc_1d_ehwlrg_params * )(p->params);

  t1 = rho[0] * rho[0];
  t4 = rho[0] * params->a2 + t1 * params->a3 + params->a1;
  t5 = pow(rho[0], params->alpha);
  tzk0 = t4 * t5;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim.zk + 0] += tzk0;

  t6 = rho[0] * params->a3;
  t8 = 0.2e1 * t6 + params->a2;
  tvrho0 = rho[0] * t8 * t5 + t4 * t5 * params->alpha + tzk0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim.vrho + 0] += tvrho0;

  t13 = t8 * t5;
  t15 = 0.1e1 / rho[0];
  t16 = params->alpha * t15;
  t22 = params->alpha * params->alpha;
  t23 = t22 * t15;
  tv2rho20 = 0.2e1 * t13 * params->alpha + tzk0 * t16 + tzk0 * t23 + 0.2e1 * t6 * t5 + 0.2e1 * t13;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim.v2rho2 + 0] += tv2rho20;

}

#endif


#ifndef XC_DONT_COMPILE_KXC
GPU_DEVICE_FUNCTION static inline void
func_kxc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_lda_out_params *out)
{
  double t1, t4, t5, tzk0;

  double t6, t8, tvrho0;

  double t13, t15, t16, t22, t23, tv2rho20;

  double t25, t29, t30, t36, t37, tv3rho30;

  lda_xc_1d_ehwlrg_params *params;

  assert(p->params != NULL);
  params = (lda_xc_1d_ehwlrg_params * )(p->params);

  t1 = rho[0] * rho[0];
  t4 = rho[0] * params->a2 + t1 * params->a3 + params->a1;
  t5 = pow(rho[0], params->alpha);
  tzk0 = t4 * t5;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim.zk + 0] += tzk0;

  t6 = rho[0] * params->a3;
  t8 = 0.2e1 * t6 + params->a2;
  tvrho0 = rho[0] * t8 * t5 + t4 * t5 * params->alpha + tzk0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim.vrho + 0] += tvrho0;

  t13 = t8 * t5;
  t15 = 0.1e1 / rho[0];
  t16 = params->alpha * t15;
  t22 = params->alpha * params->alpha;
  t23 = t22 * t15;
  tv2rho20 = 0.2e1 * t13 * params->alpha + tzk0 * t16 + tzk0 * t23 + 0.2e1 * t6 * t5 + 0.2e1 * t13;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim.v2rho2 + 0] += tv2rho20;

  t25 = params->a3 * t5;
  t29 = 0.1e1 / t1;
  t30 = params->alpha * t29;
  t36 = t22 * params->alpha;
  t37 = t36 * t29;
  tv3rho30 = 0.3e1 * t13 * t16 + 0.3e1 * t13 * t23 + 0.6e1 * t25 * params->alpha - tzk0 * t30 + tzk0 * t37 + 0.6e1 * t25;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->dim.v3rho3 + 0] += tv3rho30;

}

#endif


#ifndef XC_DONT_COMPILE_LXC
GPU_DEVICE_FUNCTION static inline void
func_lxc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_lda_out_params *out)
{
  double t1, t4, t5, tzk0;

  double t6, t8, tvrho0;

  double t13, t15, t16, t22, t23, tv2rho20;

  double t25, t29, t30, t36, t37, tv3rho30;

  double t44, t54, tv4rho40;

  lda_xc_1d_ehwlrg_params *params;

  assert(p->params != NULL);
  params = (lda_xc_1d_ehwlrg_params * )(p->params);

  t1 = rho[0] * rho[0];
  t4 = rho[0] * params->a2 + t1 * params->a3 + params->a1;
  t5 = pow(rho[0], params->alpha);
  tzk0 = t4 * t5;

  if(out->zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    out->zk[ip*p->dim.zk + 0] += tzk0;

  t6 = rho[0] * params->a3;
  t8 = 0.2e1 * t6 + params->a2;
  tvrho0 = rho[0] * t8 * t5 + t4 * t5 * params->alpha + tzk0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->dim.vrho + 0] += tvrho0;

  t13 = t8 * t5;
  t15 = 0.1e1 / rho[0];
  t16 = params->alpha * t15;
  t22 = params->alpha * params->alpha;
  t23 = t22 * t15;
  tv2rho20 = 0.2e1 * t13 * params->alpha + tzk0 * t16 + tzk0 * t23 + 0.2e1 * t6 * t5 + 0.2e1 * t13;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->dim.v2rho2 + 0] += tv2rho20;

  t25 = params->a3 * t5;
  t29 = 0.1e1 / t1;
  t30 = params->alpha * t29;
  t36 = t22 * params->alpha;
  t37 = t36 * t29;
  tv3rho30 = 0.3e1 * t13 * t16 + 0.3e1 * t13 * t23 + 0.6e1 * t25 * params->alpha - tzk0 * t30 + tzk0 * t37 + 0.6e1 * t25;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->dim.v3rho3 + 0] += tv3rho30;

  t44 = 0.1e1 / t1 / rho[0];
  t54 = t22 * t22;
  tv4rho40 = -tzk0 * t22 * t44 - 0.2e1 * tzk0 * t36 * t44 + tzk0 * t54 * t44 + 0.2e1 * tzk0 * params->alpha * t44 - 0.4e1 * t13 * t30 + 0.4e1 * t13 * t37 + 0.12e2 * t25 * t16 + 0.12e2 * t25 * t23;

  if(out->v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    out->v4rho4[ip*p->dim.v4rho4 + 0] += tv4rho40;

}

#endif

