/*
  This file was generated automatically with scripts/maple2c.py.
  Do not edit this file directly as it can be overwritten!!

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.

  Maple version     : Maple 2020 (X86 64 LINUX)
  Maple source      : ./maple/lda_vxc/lda_xc_tih.mpl
  Type of functional: lda_vxc
*/

#define maple2c_order 4
#define MAPLE2C_FLAGS (XC_FLAGS_I_HAVE_VXC | XC_FLAGS_I_HAVE_FXC | XC_FLAGS_I_HAVE_KXC | XC_FLAGS_I_HAVE_LXC)


#ifndef XC_DONT_COMPILE_VXC
GPU_DEVICE_FUNCTION static inline void
func_vxc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t3, t7, t11, t15, t19, t23, t27, t31;
  double tvrho0;


  t3 = tanh(0.10953e1 + 0.334789e-1 * rho[0]);
  t7 = tanh(-0.414661e0 + 0.152399e0 * rho[0]);
  t11 = tanh(-0.354691e0 + 0.390837e-1 * rho[0]);
  t15 = tanh(0.748531e-1 + 0.136598e0 * rho[0]);
  t19 = tanh(-0.141063e1 + 0.496577e-2 * rho[0]);
  t23 = tanh(0.48315e0 + 0.402905e1 * rho[0]);
  t27 = tanh(-0.420166e0 + 0.104352e-1 * rho[0]);
  t31 = tanh(0.147409e1 + 0.442455e0 * rho[0]);
  tvrho0 = 0.625039e0 - 0.130351e1 * t3 - 0.137026e1 * t7 - 0.129598e1 * t11 + 0.104305e1 * t15 - 0.909651e0 * t19 - 0.991782e0 * t23 - 0.915745e0 * t27 - 0.195026e1 * t31;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->out_dim->vrho + 0] += tvrho0;

}

#endif


#ifndef XC_DONT_COMPILE_FXC
GPU_DEVICE_FUNCTION static inline void
func_fxc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t3, t7, t11, t15, t19, t23, t27, t31;
  double tvrho0;

  double t33, t35, t37, t39, t41, t43, t45, t47;
  double tv2rho20;


  t3 = tanh(0.10953e1 + 0.334789e-1 * rho[0]);
  t7 = tanh(-0.414661e0 + 0.152399e0 * rho[0]);
  t11 = tanh(-0.354691e0 + 0.390837e-1 * rho[0]);
  t15 = tanh(0.748531e-1 + 0.136598e0 * rho[0]);
  t19 = tanh(-0.141063e1 + 0.496577e-2 * rho[0]);
  t23 = tanh(0.48315e0 + 0.402905e1 * rho[0]);
  t27 = tanh(-0.420166e0 + 0.104352e-1 * rho[0]);
  t31 = tanh(0.147409e1 + 0.442455e0 * rho[0]);
  tvrho0 = 0.625039e0 - 0.130351e1 * t3 - 0.137026e1 * t7 - 0.129598e1 * t11 + 0.104305e1 * t15 - 0.909651e0 * t19 - 0.991782e0 * t23 - 0.915745e0 * t27 - 0.195026e1 * t31;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->out_dim->vrho + 0] += tvrho0;

  t33 = t3 * t3;
  t35 = t7 * t7;
  t37 = t11 * t11;
  t39 = t15 * t15;
  t41 = t19 * t19;
  t43 = t23 * t23;
  t45 = t27 * t27;
  t47 = t31 * t31;
  tv2rho20 = -0.503355413957527e1 + 0.43640080939e-1 * t33 + 0.20882625374e0 * t35 + 0.50651693526e-1 * t37 - 0.1424785439e0 * t39 + 0.451711764627e-2 * t41 + 0.39959392671e1 * t43 + 0.9555982224e-2 * t45 + 0.8629022883e0 * t47;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->out_dim->v2rho2 + 0] += tv2rho20;

}

#endif


#ifndef XC_DONT_COMPILE_KXC
GPU_DEVICE_FUNCTION static inline void
func_kxc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t3, t7, t11, t15, t19, t23, t27, t31;
  double tvrho0;

  double t33, t35, t37, t39, t41, t43, t45, t47;
  double tv2rho20;

  double t50, t54, t58, t62, t66, t70, t74, t78;
  double tv3rho30;


  t3 = tanh(0.10953e1 + 0.334789e-1 * rho[0]);
  t7 = tanh(-0.414661e0 + 0.152399e0 * rho[0]);
  t11 = tanh(-0.354691e0 + 0.390837e-1 * rho[0]);
  t15 = tanh(0.748531e-1 + 0.136598e0 * rho[0]);
  t19 = tanh(-0.141063e1 + 0.496577e-2 * rho[0]);
  t23 = tanh(0.48315e0 + 0.402905e1 * rho[0]);
  t27 = tanh(-0.420166e0 + 0.104352e-1 * rho[0]);
  t31 = tanh(0.147409e1 + 0.442455e0 * rho[0]);
  tvrho0 = 0.625039e0 - 0.130351e1 * t3 - 0.137026e1 * t7 - 0.129598e1 * t11 + 0.104305e1 * t15 - 0.909651e0 * t19 - 0.991782e0 * t23 - 0.915745e0 * t27 - 0.195026e1 * t31;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->out_dim->vrho + 0] += tvrho0;

  t33 = t3 * t3;
  t35 = t7 * t7;
  t37 = t11 * t11;
  t39 = t15 * t15;
  t41 = t19 * t19;
  t43 = t23 * t23;
  t45 = t27 * t27;
  t47 = t31 * t31;
  tv2rho20 = -0.503355413957527e1 + 0.43640080939e-1 * t33 + 0.20882625374e0 * t35 + 0.50651693526e-1 * t37 - 0.1424785439e0 * t39 + 0.451711764627e-2 * t41 + 0.39959392671e1 * t43 + 0.9555982224e-2 * t45 + 0.8629022883e0 * t47;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->out_dim->v2rho2 + 0] += tv2rho20;

  t50 = 0.334789e-1 - 0.334789e-1 * t33;
  t54 = 0.152399e0 - 0.152399e0 * t35;
  t58 = 0.390837e-1 - 0.390837e-1 * t37;
  t62 = 0.136598e0 - 0.136598e0 * t39;
  t66 = 0.496577e-2 - 0.496577e-2 * t41;
  t70 = 0.402905e1 - 0.402905e1 * t43;
  t74 = 0.104352e-1 - 0.104352e-1 * t45;
  t78 = 0.442455e0 - 0.442455e0 * t47;
  tv3rho30 = 0.87280161878e-1 * t3 * t50 + 0.41765250748e0 * t7 * t54 + 0.101303387052e0 * t11 * t58 - 0.2849570878e0 * t15 * t62 + 0.903423529254e-2 * t19 * t66 + 0.79918785342e1 * t23 * t70 + 0.19111964448e-1 * t27 * t74 + 0.17258045766e1 * t31 * t78;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->out_dim->v3rho3 + 0] += tv3rho30;

}

#endif


#ifndef XC_DONT_COMPILE_LXC
GPU_DEVICE_FUNCTION static inline void
func_lxc_unpol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t3, t7, t11, t15, t19, t23, t27, t31;
  double tvrho0;

  double t33, t35, t37, t39, t41, t43, t45, t47;
  double tv2rho20;

  double t50, t54, t58, t62, t66, t70, t74, t78;
  double tv3rho30;

  double t81, t85, t89, t93, t97, t101, t105, t109;
  double tv4rho40;


  t3 = tanh(0.10953e1 + 0.334789e-1 * rho[0]);
  t7 = tanh(-0.414661e0 + 0.152399e0 * rho[0]);
  t11 = tanh(-0.354691e0 + 0.390837e-1 * rho[0]);
  t15 = tanh(0.748531e-1 + 0.136598e0 * rho[0]);
  t19 = tanh(-0.141063e1 + 0.496577e-2 * rho[0]);
  t23 = tanh(0.48315e0 + 0.402905e1 * rho[0]);
  t27 = tanh(-0.420166e0 + 0.104352e-1 * rho[0]);
  t31 = tanh(0.147409e1 + 0.442455e0 * rho[0]);
  tvrho0 = 0.625039e0 - 0.130351e1 * t3 - 0.137026e1 * t7 - 0.129598e1 * t11 + 0.104305e1 * t15 - 0.909651e0 * t19 - 0.991782e0 * t23 - 0.915745e0 * t27 - 0.195026e1 * t31;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->out_dim->vrho + 0] += tvrho0;

  t33 = t3 * t3;
  t35 = t7 * t7;
  t37 = t11 * t11;
  t39 = t15 * t15;
  t41 = t19 * t19;
  t43 = t23 * t23;
  t45 = t27 * t27;
  t47 = t31 * t31;
  tv2rho20 = -0.503355413957527e1 + 0.43640080939e-1 * t33 + 0.20882625374e0 * t35 + 0.50651693526e-1 * t37 - 0.1424785439e0 * t39 + 0.451711764627e-2 * t41 + 0.39959392671e1 * t43 + 0.9555982224e-2 * t45 + 0.8629022883e0 * t47;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->out_dim->v2rho2 + 0] += tv2rho20;

  t50 = 0.334789e-1 - 0.334789e-1 * t33;
  t54 = 0.152399e0 - 0.152399e0 * t35;
  t58 = 0.390837e-1 - 0.390837e-1 * t37;
  t62 = 0.136598e0 - 0.136598e0 * t39;
  t66 = 0.496577e-2 - 0.496577e-2 * t41;
  t70 = 0.402905e1 - 0.402905e1 * t43;
  t74 = 0.104352e-1 - 0.104352e-1 * t45;
  t78 = 0.442455e0 - 0.442455e0 * t47;
  tv3rho30 = 0.87280161878e-1 * t3 * t50 + 0.41765250748e0 * t7 * t54 + 0.101303387052e0 * t11 * t58 - 0.2849570878e0 * t15 * t62 + 0.903423529254e-2 * t19 * t66 + 0.79918785342e1 * t23 * t70 + 0.19111964448e-1 * t27 * t74 + 0.17258045766e1 * t31 * t78;

  if(out->v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    out->v3rho3[ip*p->out_dim->v3rho3 + 0] += tv3rho30;

  t81 = t50 * t50;
  t85 = t54 * t54;
  t89 = t58 * t58;
  t93 = t62 * t62;
  t97 = t66 * t66;
  t101 = t70 * t70;
  t105 = t74 * t74;
  t109 = t78 * t78;
  tv4rho40 = 0.87280161878e-1 * t81 - 0.58440876229947484e-2 * t33 * t50 + 0.41765250748e0 * t85 - 0.12729964897488904e0 * t35 * t54 + 0.101303387052e0 * t89 - 0.79186223770485048e-2 * t37 * t58 - 0.2849570878e0 * t93 + 0.778491365586088e-1 * t39 * t62 + 0.903423529254e-2 * t97 - 0.897238691772727116e-4 * t41 * t66 + 0.79918785342e1 * t101 - 0.6439935641643702e2 * t43 * t70 + 0.19111964448e-1 * t105 - 0.3988743428155392e-3 * t45 * t74 + 0.17258045766e1 * t109 - 0.1527181727879106e1 * t47 * t78;

  if(out->v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    out->v4rho4[ip*p->out_dim->v4rho4 + 0] += tv4rho40;

}

#endif


#ifndef XC_DONT_COMPILE_VXC
GPU_DEVICE_FUNCTION static inline void
func_vxc_pol(const xc_func_type *p, size_t ip, const double *rho, xc_output_variables *out)
{
  double t4, t9, t14, t19, t24, t29, t34, t39;
  double tvrho0, tvrho1;


  t4 = tanh(0.10953e1 + 0.334789e-1 * rho[0] + 0.334789e-1 * rho[1]);
  t9 = tanh(-0.414661e0 + 0.152399e0 * rho[0] + 0.152399e0 * rho[1]);
  t14 = tanh(-0.354691e0 + 0.390837e-1 * rho[0] + 0.390837e-1 * rho[1]);
  t19 = tanh(0.748531e-1 + 0.136598e0 * rho[0] + 0.136598e0 * rho[1]);
  t24 = tanh(-0.141063e1 + 0.496577e-2 * rho[0] + 0.496577e-2 * rho[1]);
  t29 = tanh(0.48315e0 + 0.402905e1 * rho[0] + 0.402905e1 * rho[1]);
  t34 = tanh(-0.420166e0 + 0.104352e-1 * rho[0] + 0.104352e-1 * rho[1]);
  t39 = tanh(0.147409e1 + 0.442455e0 * rho[0] + 0.442455e0 * rho[1]);
  tvrho0 = 0.625039e0 - 0.130351e1 * t4 - 0.137026e1 * t9 - 0.129598e1 * t14 + 0.104305e1 * t19 - 0.909651e0 * t24 - 0.991782e0 * t29 - 0.915745e0 * t34 - 0.195026e1 * t39;

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
  double t4, t9, t14, t19, t24, t29, t34, t39;
  double tvrho0, tvrho1;

  double t41, t43, t45, t47, t49, t51, t53, t55;
  double tv2rho20, tv2rho21, tv2rho22;


  t4 = tanh(0.10953e1 + 0.334789e-1 * rho[0] + 0.334789e-1 * rho[1]);
  t9 = tanh(-0.414661e0 + 0.152399e0 * rho[0] + 0.152399e0 * rho[1]);
  t14 = tanh(-0.354691e0 + 0.390837e-1 * rho[0] + 0.390837e-1 * rho[1]);
  t19 = tanh(0.748531e-1 + 0.136598e0 * rho[0] + 0.136598e0 * rho[1]);
  t24 = tanh(-0.141063e1 + 0.496577e-2 * rho[0] + 0.496577e-2 * rho[1]);
  t29 = tanh(0.48315e0 + 0.402905e1 * rho[0] + 0.402905e1 * rho[1]);
  t34 = tanh(-0.420166e0 + 0.104352e-1 * rho[0] + 0.104352e-1 * rho[1]);
  t39 = tanh(0.147409e1 + 0.442455e0 * rho[0] + 0.442455e0 * rho[1]);
  tvrho0 = 0.625039e0 - 0.130351e1 * t4 - 0.137026e1 * t9 - 0.129598e1 * t14 + 0.104305e1 * t19 - 0.909651e0 * t24 - 0.991782e0 * t29 - 0.915745e0 * t34 - 0.195026e1 * t39;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->out_dim->vrho + 0] += tvrho0;

  tvrho1 = tvrho0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->out_dim->vrho + 1] += tvrho1;

  t41 = t4 * t4;
  t43 = t9 * t9;
  t45 = t14 * t14;
  t47 = t19 * t19;
  t49 = t24 * t24;
  t51 = t29 * t29;
  t53 = t34 * t34;
  t55 = t39 * t39;
  tv2rho20 = -0.503355413957527e1 + 0.43640080939e-1 * t41 + 0.20882625374e0 * t43 + 0.50651693526e-1 * t45 - 0.1424785439e0 * t47 + 0.451711764627e-2 * t49 + 0.39959392671e1 * t51 + 0.9555982224e-2 * t53 + 0.8629022883e0 * t55;

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
  double t4, t9, t14, t19, t24, t29, t34, t39;
  double tvrho0, tvrho1;

  double t41, t43, t45, t47, t49, t51, t53, t55;
  double tv2rho20, tv2rho21, tv2rho22;

  double t58, t62, t66, t70, t74, t78, t82, t86;
  double tv3rho30, tv3rho31, tv3rho32, tv3rho33;


  t4 = tanh(0.10953e1 + 0.334789e-1 * rho[0] + 0.334789e-1 * rho[1]);
  t9 = tanh(-0.414661e0 + 0.152399e0 * rho[0] + 0.152399e0 * rho[1]);
  t14 = tanh(-0.354691e0 + 0.390837e-1 * rho[0] + 0.390837e-1 * rho[1]);
  t19 = tanh(0.748531e-1 + 0.136598e0 * rho[0] + 0.136598e0 * rho[1]);
  t24 = tanh(-0.141063e1 + 0.496577e-2 * rho[0] + 0.496577e-2 * rho[1]);
  t29 = tanh(0.48315e0 + 0.402905e1 * rho[0] + 0.402905e1 * rho[1]);
  t34 = tanh(-0.420166e0 + 0.104352e-1 * rho[0] + 0.104352e-1 * rho[1]);
  t39 = tanh(0.147409e1 + 0.442455e0 * rho[0] + 0.442455e0 * rho[1]);
  tvrho0 = 0.625039e0 - 0.130351e1 * t4 - 0.137026e1 * t9 - 0.129598e1 * t14 + 0.104305e1 * t19 - 0.909651e0 * t24 - 0.991782e0 * t29 - 0.915745e0 * t34 - 0.195026e1 * t39;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->out_dim->vrho + 0] += tvrho0;

  tvrho1 = tvrho0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->out_dim->vrho + 1] += tvrho1;

  t41 = t4 * t4;
  t43 = t9 * t9;
  t45 = t14 * t14;
  t47 = t19 * t19;
  t49 = t24 * t24;
  t51 = t29 * t29;
  t53 = t34 * t34;
  t55 = t39 * t39;
  tv2rho20 = -0.503355413957527e1 + 0.43640080939e-1 * t41 + 0.20882625374e0 * t43 + 0.50651693526e-1 * t45 - 0.1424785439e0 * t47 + 0.451711764627e-2 * t49 + 0.39959392671e1 * t51 + 0.9555982224e-2 * t53 + 0.8629022883e0 * t55;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->out_dim->v2rho2 + 0] += tv2rho20;

  tv2rho21 = tv2rho20;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->out_dim->v2rho2 + 1] += tv2rho21;

  tv2rho22 = tv2rho21;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->out_dim->v2rho2 + 2] += tv2rho22;

  t58 = 0.334789e-1 - 0.334789e-1 * t41;
  t62 = 0.152399e0 - 0.152399e0 * t43;
  t66 = 0.390837e-1 - 0.390837e-1 * t45;
  t70 = 0.136598e0 - 0.136598e0 * t47;
  t74 = 0.496577e-2 - 0.496577e-2 * t49;
  t78 = 0.402905e1 - 0.402905e1 * t51;
  t82 = 0.104352e-1 - 0.104352e-1 * t53;
  t86 = 0.442455e0 - 0.442455e0 * t55;
  tv3rho30 = 0.87280161878e-1 * t4 * t58 + 0.41765250748e0 * t9 * t62 + 0.101303387052e0 * t14 * t66 - 0.2849570878e0 * t19 * t70 + 0.903423529254e-2 * t24 * t74 + 0.79918785342e1 * t29 * t78 + 0.19111964448e-1 * t34 * t82 + 0.17258045766e1 * t39 * t86;

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
  double t4, t9, t14, t19, t24, t29, t34, t39;
  double tvrho0, tvrho1;

  double t41, t43, t45, t47, t49, t51, t53, t55;
  double tv2rho20, tv2rho21, tv2rho22;

  double t58, t62, t66, t70, t74, t78, t82, t86;
  double tv3rho30, tv3rho31, tv3rho32, tv3rho33;

  double t89, t93, t97, t101, t105, t109, t113, t117;
  double tv4rho40, tv4rho41, tv4rho42, tv4rho43, tv4rho44;


  t4 = tanh(0.10953e1 + 0.334789e-1 * rho[0] + 0.334789e-1 * rho[1]);
  t9 = tanh(-0.414661e0 + 0.152399e0 * rho[0] + 0.152399e0 * rho[1]);
  t14 = tanh(-0.354691e0 + 0.390837e-1 * rho[0] + 0.390837e-1 * rho[1]);
  t19 = tanh(0.748531e-1 + 0.136598e0 * rho[0] + 0.136598e0 * rho[1]);
  t24 = tanh(-0.141063e1 + 0.496577e-2 * rho[0] + 0.496577e-2 * rho[1]);
  t29 = tanh(0.48315e0 + 0.402905e1 * rho[0] + 0.402905e1 * rho[1]);
  t34 = tanh(-0.420166e0 + 0.104352e-1 * rho[0] + 0.104352e-1 * rho[1]);
  t39 = tanh(0.147409e1 + 0.442455e0 * rho[0] + 0.442455e0 * rho[1]);
  tvrho0 = 0.625039e0 - 0.130351e1 * t4 - 0.137026e1 * t9 - 0.129598e1 * t14 + 0.104305e1 * t19 - 0.909651e0 * t24 - 0.991782e0 * t29 - 0.915745e0 * t34 - 0.195026e1 * t39;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->out_dim->vrho + 0] += tvrho0;

  tvrho1 = tvrho0;

  if(out->vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    out->vrho[ip*p->out_dim->vrho + 1] += tvrho1;

  t41 = t4 * t4;
  t43 = t9 * t9;
  t45 = t14 * t14;
  t47 = t19 * t19;
  t49 = t24 * t24;
  t51 = t29 * t29;
  t53 = t34 * t34;
  t55 = t39 * t39;
  tv2rho20 = -0.503355413957527e1 + 0.43640080939e-1 * t41 + 0.20882625374e0 * t43 + 0.50651693526e-1 * t45 - 0.1424785439e0 * t47 + 0.451711764627e-2 * t49 + 0.39959392671e1 * t51 + 0.9555982224e-2 * t53 + 0.8629022883e0 * t55;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->out_dim->v2rho2 + 0] += tv2rho20;

  tv2rho21 = tv2rho20;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->out_dim->v2rho2 + 1] += tv2rho21;

  tv2rho22 = tv2rho21;

  if(out->v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    out->v2rho2[ip*p->out_dim->v2rho2 + 2] += tv2rho22;

  t58 = 0.334789e-1 - 0.334789e-1 * t41;
  t62 = 0.152399e0 - 0.152399e0 * t43;
  t66 = 0.390837e-1 - 0.390837e-1 * t45;
  t70 = 0.136598e0 - 0.136598e0 * t47;
  t74 = 0.496577e-2 - 0.496577e-2 * t49;
  t78 = 0.402905e1 - 0.402905e1 * t51;
  t82 = 0.104352e-1 - 0.104352e-1 * t53;
  t86 = 0.442455e0 - 0.442455e0 * t55;
  tv3rho30 = 0.87280161878e-1 * t4 * t58 + 0.41765250748e0 * t9 * t62 + 0.101303387052e0 * t14 * t66 - 0.2849570878e0 * t19 * t70 + 0.903423529254e-2 * t24 * t74 + 0.79918785342e1 * t29 * t78 + 0.19111964448e-1 * t34 * t82 + 0.17258045766e1 * t39 * t86;

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

  t89 = t58 * t58;
  t93 = t62 * t62;
  t97 = t66 * t66;
  t101 = t70 * t70;
  t105 = t74 * t74;
  t109 = t78 * t78;
  t113 = t82 * t82;
  t117 = t86 * t86;
  tv4rho40 = 0.87280161878e-1 * t89 - 0.58440876229947484e-2 * t41 * t58 + 0.41765250748e0 * t93 - 0.12729964897488904e0 * t43 * t62 + 0.101303387052e0 * t97 - 0.79186223770485048e-2 * t45 * t66 - 0.2849570878e0 * t101 + 0.778491365586088e-1 * t47 * t70 + 0.903423529254e-2 * t105 - 0.897238691772727116e-4 * t49 * t74 + 0.79918785342e1 * t109 - 0.6439935641643702e2 * t51 * t78 + 0.19111964448e-1 * t113 - 0.3988743428155392e-3 * t53 * t82 + 0.17258045766e1 * t117 - 0.1527181727879106e1 * t55 * t86;

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

