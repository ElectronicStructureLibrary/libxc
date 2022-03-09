/*
  This file was generated automatically with scripts/maple2c.py.
  Do not edit this file directly as it can be overwritten!!

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.

  Maple version     : Maple 2020 (X86 64 LINUX)
  Maple source      : ./maple/lda_exc/lda_x_2d.mpl
  Type of functional: lda_exc
*/

#define maple2c_order 4
#define MAPLE2C_FLAGS (XC_FLAGS_I_HAVE_EXC | XC_FLAGS_I_HAVE_VXC | XC_FLAGS_I_HAVE_FXC | XC_FLAGS_I_HAVE_KXC | XC_FLAGS_I_HAVE_LXC)


#ifndef XC_DONT_COMPILE_EXC
GPU_DEVICE_FUNCTION static inline void
func_exc_unpol(const xc_func_type * const p, size_t ip, const double * const rho, xc_lda_out_params *out)
{
  double t1, t2, t4, t6, t8, t9, t11, tzk0;


  t1 = M_SQRT2;
  t2 = sqrt(M_PI);
  t4 = t1 / t2;
  t6 = sqrt(p->zeta_threshold);
  t8 = my_piecewise3(0.1e1 <= p->zeta_threshold, t6 * p->zeta_threshold, 1);
  t9 = sqrt(rho[0]);
  t11 = t4 * t8 * t9;
  tzk0 = -0.4e1 / 0.3e1 * t11;

    out->zk[ip*p->dim.zk + 0] += tzk0;

}

#endif


#ifndef XC_DONT_COMPILE_VXC
GPU_DEVICE_FUNCTION static inline void
func_vxc_unpol(const xc_func_type * const p, size_t ip, const double * const rho, xc_lda_out_params *out)
{
  double t1, t2, t4, t6, t8, t9, t11, tzk0;

  double tvrho0;


  t1 = M_SQRT2;
  t2 = sqrt(M_PI);
  t4 = t1 / t2;
  t6 = sqrt(p->zeta_threshold);
  t8 = my_piecewise3(0.1e1 <= p->zeta_threshold, t6 * p->zeta_threshold, 1);
  t9 = sqrt(rho[0]);
  t11 = t4 * t8 * t9;
  tzk0 = -0.4e1 / 0.3e1 * t11;

    out->zk[ip*p->dim.zk + 0] += tzk0;

  tvrho0 = -0.2e1 * t11;

    out->vrho[ip*p->dim.vrho + 0] += tvrho0;

}

#endif


#ifndef XC_DONT_COMPILE_FXC
GPU_DEVICE_FUNCTION static inline void
func_fxc_unpol(const xc_func_type * const p, size_t ip, const double * const rho, xc_lda_out_params *out)
{
  double t1, t2, t4, t6, t8, t9, t11, tzk0;

  double tvrho0;

  double tv2rho20;


  t1 = M_SQRT2;
  t2 = sqrt(M_PI);
  t4 = t1 / t2;
  t6 = sqrt(p->zeta_threshold);
  t8 = my_piecewise3(0.1e1 <= p->zeta_threshold, t6 * p->zeta_threshold, 1);
  t9 = sqrt(rho[0]);
  t11 = t4 * t8 * t9;
  tzk0 = -0.4e1 / 0.3e1 * t11;

    out->zk[ip*p->dim.zk + 0] += tzk0;

  tvrho0 = -0.2e1 * t11;

    out->vrho[ip*p->dim.vrho + 0] += tvrho0;

  tv2rho20 = -t4 * t8 / t9;

    out->v2rho2[ip*p->dim.v2rho2 + 0] += tv2rho20;

}

#endif


#ifndef XC_DONT_COMPILE_KXC
GPU_DEVICE_FUNCTION static inline void
func_kxc_unpol(const xc_func_type * const p, size_t ip, const double * const rho, xc_lda_out_params *out)
{
  double t1, t2, t4, t6, t8, t9, t11, tzk0;

  double tvrho0;

  double tv2rho20;

  double tv3rho30;


  t1 = M_SQRT2;
  t2 = sqrt(M_PI);
  t4 = t1 / t2;
  t6 = sqrt(p->zeta_threshold);
  t8 = my_piecewise3(0.1e1 <= p->zeta_threshold, t6 * p->zeta_threshold, 1);
  t9 = sqrt(rho[0]);
  t11 = t4 * t8 * t9;
  tzk0 = -0.4e1 / 0.3e1 * t11;

    out->zk[ip*p->dim.zk + 0] += tzk0;

  tvrho0 = -0.2e1 * t11;

    out->vrho[ip*p->dim.vrho + 0] += tvrho0;

  tv2rho20 = -t4 * t8 / t9;

    out->v2rho2[ip*p->dim.v2rho2 + 0] += tv2rho20;

  tv3rho30 = t4 * t8 / t9 / rho[0] / 0.2e1;

    out->v3rho3[ip*p->dim.v3rho3 + 0] += tv3rho30;

}

#endif


#ifndef XC_DONT_COMPILE_LXC
GPU_DEVICE_FUNCTION static inline void
func_lxc_unpol(const xc_func_type * const p, size_t ip, const double * const rho, xc_lda_out_params *out)
{
  double t1, t2, t4, t6, t8, t9, t11, tzk0;

  double tvrho0;

  double tv2rho20;

  double tv3rho30;

  double t21, tv4rho40;


  t1 = M_SQRT2;
  t2 = sqrt(M_PI);
  t4 = t1 / t2;
  t6 = sqrt(p->zeta_threshold);
  t8 = my_piecewise3(0.1e1 <= p->zeta_threshold, t6 * p->zeta_threshold, 1);
  t9 = sqrt(rho[0]);
  t11 = t4 * t8 * t9;
  tzk0 = -0.4e1 / 0.3e1 * t11;

    out->zk[ip*p->dim.zk + 0] += tzk0;

  tvrho0 = -0.2e1 * t11;

    out->vrho[ip*p->dim.vrho + 0] += tvrho0;

  tv2rho20 = -t4 * t8 / t9;

    out->v2rho2[ip*p->dim.v2rho2 + 0] += tv2rho20;

  tv3rho30 = t4 * t8 / t9 / rho[0] / 0.2e1;

    out->v3rho3[ip*p->dim.v3rho3 + 0] += tv3rho30;

  t21 = rho[0] * rho[0];
  tv4rho40 = -0.3e1 / 0.4e1 * t4 * t8 / t9 / t21;

    out->v4rho4[ip*p->dim.v4rho4 + 0] += tv4rho40;

}

#endif


#ifndef XC_DONT_COMPILE_EXC
GPU_DEVICE_FUNCTION static inline void
func_exc_pol(const xc_func_type * const p, size_t ip, const double * const rho, xc_lda_out_params *out)
{
  double t1, t2, t3, t4, t5, t6, t7, t8;
  double t9, t10, t11, t12, t13, t14, t15, t16;
  double t17, t18, t19, t20, t22, t23, t25, tzk0;


  t1 = M_SQRT2;
  t2 = sqrt(M_PI);
  t3 = 0.1e1 / t2;
  t4 = t1 * t3;
  t5 = rho[0] - rho[1];
  t6 = rho[0] + rho[1];
  t7 = 0.1e1 / t6;
  t8 = t5 * t7;
  t9 = 0.1e1 + t8;
  t10 = t9 <= p->zeta_threshold;
  t11 = sqrt(p->zeta_threshold);
  t12 = t11 * p->zeta_threshold;
  t13 = sqrt(t9);
  t14 = t13 * t9;
  t15 = my_piecewise3(t10, t12, t14);
  t16 = 0.1e1 - t8;
  t17 = t16 <= p->zeta_threshold;
  t18 = sqrt(t16);
  t19 = t18 * t16;
  t20 = my_piecewise3(t17, t12, t19);
  t22 = t15 / 0.2e1 + t20 / 0.2e1;
  t23 = sqrt(t6);
  t25 = t4 * t22 * t23;
  tzk0 = -0.4e1 / 0.3e1 * t25;

    out->zk[ip*p->dim.zk + 0] += tzk0;

}

#endif


#ifndef XC_DONT_COMPILE_VXC
GPU_DEVICE_FUNCTION static inline void
func_vxc_pol(const xc_func_type * const p, size_t ip, const double * const rho, xc_lda_out_params *out)
{
  double t1, t2, t3, t4, t5, t6, t7, t8;
  double t9, t10, t11, t12, t13, t14, t15, t16;
  double t17, t18, t19, t20, t22, t23, t25, tzk0;

  double t27, t28, t29, t30, t31, t32, t33, t36;
  double t37, t40, t42, tvrho0, t46, t49, t50, t53;
  double t56, tvrho1;


  t1 = M_SQRT2;
  t2 = sqrt(M_PI);
  t3 = 0.1e1 / t2;
  t4 = t1 * t3;
  t5 = rho[0] - rho[1];
  t6 = rho[0] + rho[1];
  t7 = 0.1e1 / t6;
  t8 = t5 * t7;
  t9 = 0.1e1 + t8;
  t10 = t9 <= p->zeta_threshold;
  t11 = sqrt(p->zeta_threshold);
  t12 = t11 * p->zeta_threshold;
  t13 = sqrt(t9);
  t14 = t13 * t9;
  t15 = my_piecewise3(t10, t12, t14);
  t16 = 0.1e1 - t8;
  t17 = t16 <= p->zeta_threshold;
  t18 = sqrt(t16);
  t19 = t18 * t16;
  t20 = my_piecewise3(t17, t12, t19);
  t22 = t15 / 0.2e1 + t20 / 0.2e1;
  t23 = sqrt(t6);
  t25 = t4 * t22 * t23;
  tzk0 = -0.4e1 / 0.3e1 * t25;

    out->zk[ip*p->dim.zk + 0] += tzk0;

  t27 = 0.2e1 * t25;
  t28 = t23 * t6;
  t29 = t28 * t1;
  t30 = t6 * t6;
  t31 = 0.1e1 / t30;
  t32 = t5 * t31;
  t33 = t7 - t32;
  t36 = my_piecewise3(t10, 0, 0.3e1 / 0.2e1 * t13 * t33);
  t37 = -t33;
  t40 = my_piecewise3(t17, 0, 0.3e1 / 0.2e1 * t18 * t37);
  t42 = t36 / 0.2e1 + t40 / 0.2e1;
  tvrho0 = -t27 - 0.4e1 / 0.3e1 * t29 * t3 * t42;

    out->vrho[ip*p->dim.vrho + 0] += tvrho0;

  t46 = -t7 - t32;
  t49 = my_piecewise3(t10, 0, 0.3e1 / 0.2e1 * t13 * t46);
  t50 = -t46;
  t53 = my_piecewise3(t17, 0, 0.3e1 / 0.2e1 * t18 * t50);
  t56 = t3 * (t49 / 0.2e1 + t53 / 0.2e1);
  tvrho1 = -t27 - 0.4e1 / 0.3e1 * t29 * t56;

    out->vrho[ip*p->dim.vrho + 1] += tvrho1;

}

#endif


#ifndef XC_DONT_COMPILE_FXC
GPU_DEVICE_FUNCTION static inline void
func_fxc_pol(const xc_func_type * const p, size_t ip, const double * const rho, xc_lda_out_params *out)
{
  double t1, t2, t3, t4, t5, t6, t7, t8;
  double t9, t10, t11, t12, t13, t14, t15, t16;
  double t17, t18, t19, t20, t22, t23, t25, tzk0;

  double t27, t28, t29, t30, t31, t32, t33, t36;
  double t37, t40, t42, tvrho0, t46, t49, t50, t53;
  double t56, tvrho1;

  double t60, t62, t64, t65, t66, t70, t71, t73;
  double t77, t78, t79, t82, t86, t88, tv2rho20, t93;
  double t94, t96, t99, t103, t104, t107, t111, t114;
  double tv2rho21, t118, t122, t126, t127, t130, t134, t137;
  double tv2rho22;


  t1 = M_SQRT2;
  t2 = sqrt(M_PI);
  t3 = 0.1e1 / t2;
  t4 = t1 * t3;
  t5 = rho[0] - rho[1];
  t6 = rho[0] + rho[1];
  t7 = 0.1e1 / t6;
  t8 = t5 * t7;
  t9 = 0.1e1 + t8;
  t10 = t9 <= p->zeta_threshold;
  t11 = sqrt(p->zeta_threshold);
  t12 = t11 * p->zeta_threshold;
  t13 = sqrt(t9);
  t14 = t13 * t9;
  t15 = my_piecewise3(t10, t12, t14);
  t16 = 0.1e1 - t8;
  t17 = t16 <= p->zeta_threshold;
  t18 = sqrt(t16);
  t19 = t18 * t16;
  t20 = my_piecewise3(t17, t12, t19);
  t22 = t15 / 0.2e1 + t20 / 0.2e1;
  t23 = sqrt(t6);
  t25 = t4 * t22 * t23;
  tzk0 = -0.4e1 / 0.3e1 * t25;

    out->zk[ip*p->dim.zk + 0] += tzk0;

  t27 = 0.2e1 * t25;
  t28 = t23 * t6;
  t29 = t28 * t1;
  t30 = t6 * t6;
  t31 = 0.1e1 / t30;
  t32 = t5 * t31;
  t33 = t7 - t32;
  t36 = my_piecewise3(t10, 0, 0.3e1 / 0.2e1 * t13 * t33);
  t37 = -t33;
  t40 = my_piecewise3(t17, 0, 0.3e1 / 0.2e1 * t18 * t37);
  t42 = t36 / 0.2e1 + t40 / 0.2e1;
  tvrho0 = -t27 - 0.4e1 / 0.3e1 * t29 * t3 * t42;

    out->vrho[ip*p->dim.vrho + 0] += tvrho0;

  t46 = -t7 - t32;
  t49 = my_piecewise3(t10, 0, 0.3e1 / 0.2e1 * t13 * t46);
  t50 = -t46;
  t53 = my_piecewise3(t17, 0, 0.3e1 / 0.2e1 * t18 * t50);
  t56 = t3 * (t49 / 0.2e1 + t53 / 0.2e1);
  tvrho1 = -t27 - 0.4e1 / 0.3e1 * t29 * t56;

    out->vrho[ip*p->dim.vrho + 1] += tvrho1;

  t60 = t4 * t42 * t23;
  t62 = 0.1e1 / t23;
  t64 = t4 * t22 * t62;
  t65 = 0.1e1 / t13;
  t66 = t33 * t33;
  t70 = 0.1e1 / t30 / t6;
  t71 = t5 * t70;
  t73 = -0.2e1 * t31 + 0.2e1 * t71;
  t77 = my_piecewise3(t10, 0, 0.3e1 / 0.4e1 * t65 * t66 + 0.3e1 / 0.2e1 * t13 * t73);
  t78 = 0.1e1 / t18;
  t79 = t37 * t37;
  t82 = -t73;
  t86 = my_piecewise3(t17, 0, 0.3e1 / 0.4e1 * t78 * t79 + 0.3e1 / 0.2e1 * t18 * t82);
  t88 = t77 / 0.2e1 + t86 / 0.2e1;
  tv2rho20 = -0.4e1 * t60 - t64 - 0.4e1 / 0.3e1 * t29 * t3 * t88;

    out->v2rho2[ip*p->dim.v2rho2 + 0] += tv2rho20;

  t93 = t23 * t1;
  t94 = t93 * t56;
  t96 = t65 * t46;
  t99 = t13 * t5;
  t103 = my_piecewise3(t10, 0, 0.3e1 / 0.4e1 * t96 * t33 + 0.3e1 * t99 * t70);
  t104 = t78 * t50;
  t107 = t18 * t5;
  t111 = my_piecewise3(t17, 0, 0.3e1 / 0.4e1 * t104 * t37 - 0.3e1 * t107 * t70);
  t114 = t3 * (t103 / 0.2e1 + t111 / 0.2e1);
  tv2rho21 = -0.2e1 * t60 - t64 - 0.2e1 * t94 - 0.4e1 / 0.3e1 * t29 * t114;

    out->v2rho2[ip*p->dim.v2rho2 + 1] += tv2rho21;

  t118 = t46 * t46;
  t122 = 0.2e1 * t31 + 0.2e1 * t71;
  t126 = my_piecewise3(t10, 0, 0.3e1 / 0.4e1 * t65 * t118 + 0.3e1 / 0.2e1 * t13 * t122);
  t127 = t50 * t50;
  t130 = -t122;
  t134 = my_piecewise3(t17, 0, 0.3e1 / 0.4e1 * t78 * t127 + 0.3e1 / 0.2e1 * t18 * t130);
  t137 = t3 * (t126 / 0.2e1 + t134 / 0.2e1);
  tv2rho22 = -0.4e1 * t94 - t64 - 0.4e1 / 0.3e1 * t29 * t137;

    out->v2rho2[ip*p->dim.v2rho2 + 2] += tv2rho22;

}

#endif


#ifndef XC_DONT_COMPILE_KXC
GPU_DEVICE_FUNCTION static inline void
func_kxc_pol(const xc_func_type * const p, size_t ip, const double * const rho, xc_lda_out_params *out)
{
  double t1, t2, t3, t4, t5, t6, t7, t8;
  double t9, t10, t11, t12, t13, t14, t15, t16;
  double t17, t18, t19, t20, t22, t23, t25, tzk0;

  double t27, t28, t29, t30, t31, t32, t33, t36;
  double t37, t40, t42, tvrho0, t46, t49, t50, t53;
  double t56, tvrho1;

  double t60, t62, t64, t65, t66, t70, t71, t73;
  double t77, t78, t79, t82, t86, t88, tv2rho20, t93;
  double t94, t96, t99, t103, t104, t107, t111, t114;
  double tv2rho21, t118, t122, t126, t127, t130, t134, t137;
  double tv2rho22;

  double t141, t144, t146, t149, t150, t151, t154, t157;
  double t158, t159, t161, t165, t166, t167, t170, t173;
  double t177, t179, tv3rho30, t185, t186, t188, t189, t192;
  double t203, t204, t207, t218, t221, tv3rho31, t225, t227;
  double t232, t237, t241, t242, t247, t250, t254, t257;
  double tv3rho32, t262, t268, t272, t273, t278, t282, t285;
  double tv3rho33;


  t1 = M_SQRT2;
  t2 = sqrt(M_PI);
  t3 = 0.1e1 / t2;
  t4 = t1 * t3;
  t5 = rho[0] - rho[1];
  t6 = rho[0] + rho[1];
  t7 = 0.1e1 / t6;
  t8 = t5 * t7;
  t9 = 0.1e1 + t8;
  t10 = t9 <= p->zeta_threshold;
  t11 = sqrt(p->zeta_threshold);
  t12 = t11 * p->zeta_threshold;
  t13 = sqrt(t9);
  t14 = t13 * t9;
  t15 = my_piecewise3(t10, t12, t14);
  t16 = 0.1e1 - t8;
  t17 = t16 <= p->zeta_threshold;
  t18 = sqrt(t16);
  t19 = t18 * t16;
  t20 = my_piecewise3(t17, t12, t19);
  t22 = t15 / 0.2e1 + t20 / 0.2e1;
  t23 = sqrt(t6);
  t25 = t4 * t22 * t23;
  tzk0 = -0.4e1 / 0.3e1 * t25;

    out->zk[ip*p->dim.zk + 0] += tzk0;

  t27 = 0.2e1 * t25;
  t28 = t23 * t6;
  t29 = t28 * t1;
  t30 = t6 * t6;
  t31 = 0.1e1 / t30;
  t32 = t5 * t31;
  t33 = t7 - t32;
  t36 = my_piecewise3(t10, 0, 0.3e1 / 0.2e1 * t13 * t33);
  t37 = -t33;
  t40 = my_piecewise3(t17, 0, 0.3e1 / 0.2e1 * t18 * t37);
  t42 = t36 / 0.2e1 + t40 / 0.2e1;
  tvrho0 = -t27 - 0.4e1 / 0.3e1 * t29 * t3 * t42;

    out->vrho[ip*p->dim.vrho + 0] += tvrho0;

  t46 = -t7 - t32;
  t49 = my_piecewise3(t10, 0, 0.3e1 / 0.2e1 * t13 * t46);
  t50 = -t46;
  t53 = my_piecewise3(t17, 0, 0.3e1 / 0.2e1 * t18 * t50);
  t56 = t3 * (t49 / 0.2e1 + t53 / 0.2e1);
  tvrho1 = -t27 - 0.4e1 / 0.3e1 * t29 * t56;

    out->vrho[ip*p->dim.vrho + 1] += tvrho1;

  t60 = t4 * t42 * t23;
  t62 = 0.1e1 / t23;
  t64 = t4 * t22 * t62;
  t65 = 0.1e1 / t13;
  t66 = t33 * t33;
  t70 = 0.1e1 / t30 / t6;
  t71 = t5 * t70;
  t73 = -0.2e1 * t31 + 0.2e1 * t71;
  t77 = my_piecewise3(t10, 0, 0.3e1 / 0.4e1 * t65 * t66 + 0.3e1 / 0.2e1 * t13 * t73);
  t78 = 0.1e1 / t18;
  t79 = t37 * t37;
  t82 = -t73;
  t86 = my_piecewise3(t17, 0, 0.3e1 / 0.4e1 * t78 * t79 + 0.3e1 / 0.2e1 * t18 * t82);
  t88 = t77 / 0.2e1 + t86 / 0.2e1;
  tv2rho20 = -0.4e1 * t60 - t64 - 0.4e1 / 0.3e1 * t29 * t3 * t88;

    out->v2rho2[ip*p->dim.v2rho2 + 0] += tv2rho20;

  t93 = t23 * t1;
  t94 = t93 * t56;
  t96 = t65 * t46;
  t99 = t13 * t5;
  t103 = my_piecewise3(t10, 0, 0.3e1 / 0.4e1 * t96 * t33 + 0.3e1 * t99 * t70);
  t104 = t78 * t50;
  t107 = t18 * t5;
  t111 = my_piecewise3(t17, 0, 0.3e1 / 0.4e1 * t104 * t37 - 0.3e1 * t107 * t70);
  t114 = t3 * (t103 / 0.2e1 + t111 / 0.2e1);
  tv2rho21 = -0.2e1 * t60 - t64 - 0.2e1 * t94 - 0.4e1 / 0.3e1 * t29 * t114;

    out->v2rho2[ip*p->dim.v2rho2 + 1] += tv2rho21;

  t118 = t46 * t46;
  t122 = 0.2e1 * t31 + 0.2e1 * t71;
  t126 = my_piecewise3(t10, 0, 0.3e1 / 0.4e1 * t65 * t118 + 0.3e1 / 0.2e1 * t13 * t122);
  t127 = t50 * t50;
  t130 = -t122;
  t134 = my_piecewise3(t17, 0, 0.3e1 / 0.4e1 * t78 * t127 + 0.3e1 / 0.2e1 * t18 * t130);
  t137 = t3 * (t126 / 0.2e1 + t134 / 0.2e1);
  tv2rho22 = -0.4e1 * t94 - t64 - 0.4e1 / 0.3e1 * t29 * t137;

    out->v2rho2[ip*p->dim.v2rho2 + 2] += tv2rho22;

  t141 = t4 * t88 * t23;
  t144 = t4 * t42 * t62;
  t146 = 0.1e1 / t28;
  t149 = t4 * t22 * t146 / 0.2e1;
  t150 = 0.1e1 / t14;
  t151 = t66 * t33;
  t154 = t65 * t33;
  t157 = t30 * t30;
  t158 = 0.1e1 / t157;
  t159 = t5 * t158;
  t161 = 0.6e1 * t70 - 0.6e1 * t159;
  t165 = my_piecewise3(t10, 0, -0.3e1 / 0.8e1 * t150 * t151 + 0.9e1 / 0.4e1 * t154 * t73 + 0.3e1 / 0.2e1 * t13 * t161);
  t166 = 0.1e1 / t19;
  t167 = t79 * t37;
  t170 = t78 * t37;
  t173 = -t161;
  t177 = my_piecewise3(t17, 0, -0.3e1 / 0.8e1 * t166 * t167 + 0.9e1 / 0.4e1 * t170 * t82 + 0.3e1 / 0.2e1 * t18 * t173);
  t179 = t165 / 0.2e1 + t177 / 0.2e1;
  tv3rho30 = -0.6e1 * t141 - 0.3e1 * t144 + t149 - 0.4e1 / 0.3e1 * t29 * t3 * t179;

    out->v3rho3[ip*p->dim.v3rho3 + 0] += tv3rho30;

  t185 = t62 * t1;
  t186 = t185 * t56;
  t188 = 0.4e1 * t93 * t114;
  t189 = t150 * t46;
  t192 = t65 * t5;
  t203 = my_piecewise3(t10, 0, -0.3e1 / 0.8e1 * t189 * t66 + 0.3e1 * t192 * t70 * t33 + 0.3e1 / 0.4e1 * t96 * t73 + 0.3e1 * t13 * t70 - 0.9e1 * t99 * t158);
  t204 = t166 * t50;
  t207 = t78 * t5;
  t218 = my_piecewise3(t17, 0, -0.3e1 / 0.8e1 * t204 * t79 - 0.3e1 * t207 * t70 * t37 + 0.3e1 / 0.4e1 * t104 * t82 - 0.3e1 * t18 * t70 + 0.9e1 * t107 * t158);
  t221 = t3 * (t203 / 0.2e1 + t218 / 0.2e1);
  tv3rho31 = -0.2e1 * t141 - 0.2e1 * t144 + t149 - t186 - t188 - 0.4e1 / 0.3e1 * t29 * t221;

    out->v3rho3[ip*p->dim.v3rho3 + 1] += tv3rho31;

  t225 = t93 * t137;
  t227 = t150 * t118;
  t232 = t65 * t122;
  t237 = -0.2e1 * t70 - 0.6e1 * t159;
  t241 = my_piecewise3(t10, 0, -0.3e1 / 0.8e1 * t227 * t33 + 0.3e1 * t96 * t71 + 0.3e1 / 0.4e1 * t232 * t33 + 0.3e1 / 0.2e1 * t13 * t237);
  t242 = t166 * t127;
  t247 = t78 * t130;
  t250 = -t237;
  t254 = my_piecewise3(t17, 0, -0.3e1 / 0.8e1 * t242 * t37 - 0.3e1 * t104 * t71 + 0.3e1 / 0.4e1 * t247 * t37 + 0.3e1 / 0.2e1 * t18 * t250);
  t257 = t3 * (t241 / 0.2e1 + t254 / 0.2e1);
  tv3rho32 = -0.2e1 * t186 - t188 - t144 + t149 - 0.2e1 * t225 - 0.4e1 / 0.3e1 * t29 * t257;

    out->v3rho3[ip*p->dim.v3rho3 + 2] += tv3rho32;

  t262 = t118 * t46;
  t268 = -0.6e1 * t70 - 0.6e1 * t159;
  t272 = my_piecewise3(t10, 0, -0.3e1 / 0.8e1 * t150 * t262 + 0.9e1 / 0.4e1 * t96 * t122 + 0.3e1 / 0.2e1 * t13 * t268);
  t273 = t127 * t50;
  t278 = -t268;
  t282 = my_piecewise3(t17, 0, -0.3e1 / 0.8e1 * t166 * t273 + 0.9e1 / 0.4e1 * t104 * t130 + 0.3e1 / 0.2e1 * t18 * t278);
  t285 = t3 * (t272 / 0.2e1 + t282 / 0.2e1);
  tv3rho33 = -0.3e1 * t186 - 0.6e1 * t225 + t149 - 0.4e1 / 0.3e1 * t29 * t285;

    out->v3rho3[ip*p->dim.v3rho3 + 3] += tv3rho33;

}

#endif


#ifndef XC_DONT_COMPILE_LXC
GPU_DEVICE_FUNCTION static inline void
func_lxc_pol(const xc_func_type * const p, size_t ip, const double * const rho, xc_lda_out_params *out)
{
  double t1, t2, t3, t4, t5, t6, t7, t8;
  double t9, t10, t11, t12, t13, t14, t15, t16;
  double t17, t18, t19, t20, t22, t23, t25, tzk0;

  double t27, t28, t29, t30, t31, t32, t33, t36;
  double t37, t40, t42, tvrho0, t46, t49, t50, t53;
  double t56, tvrho1;

  double t60, t62, t64, t65, t66, t70, t71, t73;
  double t77, t78, t79, t82, t86, t88, tv2rho20, t93;
  double t94, t96, t99, t103, t104, t107, t111, t114;
  double tv2rho21, t118, t122, t126, t127, t130, t134, t137;
  double tv2rho22;

  double t141, t144, t146, t149, t150, t151, t154, t157;
  double t158, t159, t161, t165, t166, t167, t170, t173;
  double t177, t179, tv3rho30, t185, t186, t188, t189, t192;
  double t203, t204, t207, t218, t221, tv3rho31, t225, t227;
  double t232, t237, t241, t242, t247, t250, t254, t257;
  double tv3rho32, t262, t268, t272, t273, t278, t282, t285;
  double tv3rho33;

  double t289, t292, t295, t301, t302, t304, t305, t311;
  double t317, t318, t320, t324, t325, t327, t328, t334;
  double t343, tv4rho40, t353, t355, t356, t357, t383, t385;
  double t410, t412, tv4rho41, t420, t421, t432, t435, t451;
  double t477, tv4rho42, t487, t507, t511, t532, tv4rho43, t541;
  double t546, t552, t556, t557, t562, t571, tv4rho44;


  t1 = M_SQRT2;
  t2 = sqrt(M_PI);
  t3 = 0.1e1 / t2;
  t4 = t1 * t3;
  t5 = rho[0] - rho[1];
  t6 = rho[0] + rho[1];
  t7 = 0.1e1 / t6;
  t8 = t5 * t7;
  t9 = 0.1e1 + t8;
  t10 = t9 <= p->zeta_threshold;
  t11 = sqrt(p->zeta_threshold);
  t12 = t11 * p->zeta_threshold;
  t13 = sqrt(t9);
  t14 = t13 * t9;
  t15 = my_piecewise3(t10, t12, t14);
  t16 = 0.1e1 - t8;
  t17 = t16 <= p->zeta_threshold;
  t18 = sqrt(t16);
  t19 = t18 * t16;
  t20 = my_piecewise3(t17, t12, t19);
  t22 = t15 / 0.2e1 + t20 / 0.2e1;
  t23 = sqrt(t6);
  t25 = t4 * t22 * t23;
  tzk0 = -0.4e1 / 0.3e1 * t25;

    out->zk[ip*p->dim.zk + 0] += tzk0;

  t27 = 0.2e1 * t25;
  t28 = t23 * t6;
  t29 = t28 * t1;
  t30 = t6 * t6;
  t31 = 0.1e1 / t30;
  t32 = t5 * t31;
  t33 = t7 - t32;
  t36 = my_piecewise3(t10, 0, 0.3e1 / 0.2e1 * t13 * t33);
  t37 = -t33;
  t40 = my_piecewise3(t17, 0, 0.3e1 / 0.2e1 * t18 * t37);
  t42 = t36 / 0.2e1 + t40 / 0.2e1;
  tvrho0 = -t27 - 0.4e1 / 0.3e1 * t29 * t3 * t42;

    out->vrho[ip*p->dim.vrho + 0] += tvrho0;

  t46 = -t7 - t32;
  t49 = my_piecewise3(t10, 0, 0.3e1 / 0.2e1 * t13 * t46);
  t50 = -t46;
  t53 = my_piecewise3(t17, 0, 0.3e1 / 0.2e1 * t18 * t50);
  t56 = t3 * (t49 / 0.2e1 + t53 / 0.2e1);
  tvrho1 = -t27 - 0.4e1 / 0.3e1 * t29 * t56;

    out->vrho[ip*p->dim.vrho + 1] += tvrho1;

  t60 = t4 * t42 * t23;
  t62 = 0.1e1 / t23;
  t64 = t4 * t22 * t62;
  t65 = 0.1e1 / t13;
  t66 = t33 * t33;
  t70 = 0.1e1 / t30 / t6;
  t71 = t5 * t70;
  t73 = -0.2e1 * t31 + 0.2e1 * t71;
  t77 = my_piecewise3(t10, 0, 0.3e1 / 0.4e1 * t65 * t66 + 0.3e1 / 0.2e1 * t13 * t73);
  t78 = 0.1e1 / t18;
  t79 = t37 * t37;
  t82 = -t73;
  t86 = my_piecewise3(t17, 0, 0.3e1 / 0.4e1 * t78 * t79 + 0.3e1 / 0.2e1 * t18 * t82);
  t88 = t77 / 0.2e1 + t86 / 0.2e1;
  tv2rho20 = -0.4e1 * t60 - t64 - 0.4e1 / 0.3e1 * t29 * t3 * t88;

    out->v2rho2[ip*p->dim.v2rho2 + 0] += tv2rho20;

  t93 = t23 * t1;
  t94 = t93 * t56;
  t96 = t65 * t46;
  t99 = t13 * t5;
  t103 = my_piecewise3(t10, 0, 0.3e1 / 0.4e1 * t96 * t33 + 0.3e1 * t99 * t70);
  t104 = t78 * t50;
  t107 = t18 * t5;
  t111 = my_piecewise3(t17, 0, 0.3e1 / 0.4e1 * t104 * t37 - 0.3e1 * t107 * t70);
  t114 = t3 * (t103 / 0.2e1 + t111 / 0.2e1);
  tv2rho21 = -0.2e1 * t60 - t64 - 0.2e1 * t94 - 0.4e1 / 0.3e1 * t29 * t114;

    out->v2rho2[ip*p->dim.v2rho2 + 1] += tv2rho21;

  t118 = t46 * t46;
  t122 = 0.2e1 * t31 + 0.2e1 * t71;
  t126 = my_piecewise3(t10, 0, 0.3e1 / 0.4e1 * t65 * t118 + 0.3e1 / 0.2e1 * t13 * t122);
  t127 = t50 * t50;
  t130 = -t122;
  t134 = my_piecewise3(t17, 0, 0.3e1 / 0.4e1 * t78 * t127 + 0.3e1 / 0.2e1 * t18 * t130);
  t137 = t3 * (t126 / 0.2e1 + t134 / 0.2e1);
  tv2rho22 = -0.4e1 * t94 - t64 - 0.4e1 / 0.3e1 * t29 * t137;

    out->v2rho2[ip*p->dim.v2rho2 + 2] += tv2rho22;

  t141 = t4 * t88 * t23;
  t144 = t4 * t42 * t62;
  t146 = 0.1e1 / t28;
  t149 = t4 * t22 * t146 / 0.2e1;
  t150 = 0.1e1 / t14;
  t151 = t66 * t33;
  t154 = t65 * t33;
  t157 = t30 * t30;
  t158 = 0.1e1 / t157;
  t159 = t5 * t158;
  t161 = 0.6e1 * t70 - 0.6e1 * t159;
  t165 = my_piecewise3(t10, 0, -0.3e1 / 0.8e1 * t150 * t151 + 0.9e1 / 0.4e1 * t154 * t73 + 0.3e1 / 0.2e1 * t13 * t161);
  t166 = 0.1e1 / t19;
  t167 = t79 * t37;
  t170 = t78 * t37;
  t173 = -t161;
  t177 = my_piecewise3(t17, 0, -0.3e1 / 0.8e1 * t166 * t167 + 0.9e1 / 0.4e1 * t170 * t82 + 0.3e1 / 0.2e1 * t18 * t173);
  t179 = t165 / 0.2e1 + t177 / 0.2e1;
  tv3rho30 = -0.6e1 * t141 - 0.3e1 * t144 + t149 - 0.4e1 / 0.3e1 * t29 * t3 * t179;

    out->v3rho3[ip*p->dim.v3rho3 + 0] += tv3rho30;

  t185 = t62 * t1;
  t186 = t185 * t56;
  t188 = 0.4e1 * t93 * t114;
  t189 = t150 * t46;
  t192 = t65 * t5;
  t203 = my_piecewise3(t10, 0, -0.3e1 / 0.8e1 * t189 * t66 + 0.3e1 * t192 * t70 * t33 + 0.3e1 / 0.4e1 * t96 * t73 + 0.3e1 * t13 * t70 - 0.9e1 * t99 * t158);
  t204 = t166 * t50;
  t207 = t78 * t5;
  t218 = my_piecewise3(t17, 0, -0.3e1 / 0.8e1 * t204 * t79 - 0.3e1 * t207 * t70 * t37 + 0.3e1 / 0.4e1 * t104 * t82 - 0.3e1 * t18 * t70 + 0.9e1 * t107 * t158);
  t221 = t3 * (t203 / 0.2e1 + t218 / 0.2e1);
  tv3rho31 = -0.2e1 * t141 - 0.2e1 * t144 + t149 - t186 - t188 - 0.4e1 / 0.3e1 * t29 * t221;

    out->v3rho3[ip*p->dim.v3rho3 + 1] += tv3rho31;

  t225 = t93 * t137;
  t227 = t150 * t118;
  t232 = t65 * t122;
  t237 = -0.2e1 * t70 - 0.6e1 * t159;
  t241 = my_piecewise3(t10, 0, -0.3e1 / 0.8e1 * t227 * t33 + 0.3e1 * t96 * t71 + 0.3e1 / 0.4e1 * t232 * t33 + 0.3e1 / 0.2e1 * t13 * t237);
  t242 = t166 * t127;
  t247 = t78 * t130;
  t250 = -t237;
  t254 = my_piecewise3(t17, 0, -0.3e1 / 0.8e1 * t242 * t37 - 0.3e1 * t104 * t71 + 0.3e1 / 0.4e1 * t247 * t37 + 0.3e1 / 0.2e1 * t18 * t250);
  t257 = t3 * (t241 / 0.2e1 + t254 / 0.2e1);
  tv3rho32 = -0.2e1 * t186 - t188 - t144 + t149 - 0.2e1 * t225 - 0.4e1 / 0.3e1 * t29 * t257;

    out->v3rho3[ip*p->dim.v3rho3 + 2] += tv3rho32;

  t262 = t118 * t46;
  t268 = -0.6e1 * t70 - 0.6e1 * t159;
  t272 = my_piecewise3(t10, 0, -0.3e1 / 0.8e1 * t150 * t262 + 0.9e1 / 0.4e1 * t96 * t122 + 0.3e1 / 0.2e1 * t13 * t268);
  t273 = t127 * t50;
  t278 = -t268;
  t282 = my_piecewise3(t17, 0, -0.3e1 / 0.8e1 * t166 * t273 + 0.9e1 / 0.4e1 * t104 * t130 + 0.3e1 / 0.2e1 * t18 * t278);
  t285 = t3 * (t272 / 0.2e1 + t282 / 0.2e1);
  tv3rho33 = -0.3e1 * t186 - 0.6e1 * t225 + t149 - 0.4e1 / 0.3e1 * t29 * t285;

    out->v3rho3[ip*p->dim.v3rho3 + 3] += tv3rho33;

  t289 = t4 * t179 * t23;
  t292 = t4 * t88 * t62;
  t295 = t4 * t42 * t146;
  t301 = 0.3e1 / 0.4e1 * t4 * t22 / t23 / t30;
  t302 = t9 * t9;
  t304 = 0.1e1 / t13 / t302;
  t305 = t66 * t66;
  t311 = t73 * t73;
  t317 = 0.1e1 / t157 / t6;
  t318 = t5 * t317;
  t320 = -0.24e2 * t158 + 0.24e2 * t318;
  t324 = my_piecewise3(t10, 0, 0.9e1 / 0.16e2 * t304 * t305 - 0.9e1 / 0.4e1 * t150 * t66 * t73 + 0.9e1 / 0.4e1 * t65 * t311 + 0.3e1 * t154 * t161 + 0.3e1 / 0.2e1 * t13 * t320);
  t325 = t16 * t16;
  t327 = 0.1e1 / t18 / t325;
  t328 = t79 * t79;
  t334 = t82 * t82;
  t343 = my_piecewise3(t17, 0, 0.9e1 / 0.16e2 * t327 * t328 - 0.9e1 / 0.4e1 * t166 * t79 * t82 + 0.9e1 / 0.4e1 * t78 * t334 + 0.3e1 * t170 * t173 - 0.3e1 / 0.2e1 * t18 * t320);
  tv4rho40 = -0.8e1 * t289 - 0.6e1 * t292 + 0.2e1 * t295 - t301 - 0.4e1 / 0.3e1 * t29 * t3 * (t324 / 0.2e1 + t343 / 0.2e1);

    out->v4rho4[ip*p->dim.v4rho4 + 0] += tv4rho40;

  t353 = t146 * t1 * t56;
  t355 = t185 * t114;
  t356 = 0.3e1 * t355;
  t357 = t93 * t221;
  t383 = 0.36e2 * t99 * t317;
  t385 = my_piecewise3(t10, 0, 0.9e1 / 0.16e2 * t304 * t46 * t151 - 0.9e1 / 0.4e1 * t150 * t5 * t70 * t66 - 0.9e1 / 0.8e1 * t189 * t33 * t73 + 0.9e1 / 0.2e1 * t65 * t70 * t33 - 0.27e2 / 0.2e1 * t192 * t158 * t33 + 0.9e1 / 0.2e1 * t192 * t70 * t73 + 0.3e1 / 0.4e1 * t96 * t161 - 0.18e2 * t13 * t158 + t383);
  t410 = 0.36e2 * t107 * t317;
  t412 = my_piecewise3(t17, 0, 0.9e1 / 0.16e2 * t327 * t50 * t167 + 0.9e1 / 0.4e1 * t166 * t5 * t70 * t79 - 0.9e1 / 0.8e1 * t204 * t37 * t82 - 0.9e1 / 0.2e1 * t78 * t70 * t37 + 0.27e2 / 0.2e1 * t207 * t158 * t37 - 0.9e1 / 0.2e1 * t207 * t70 * t82 + 0.3e1 / 0.4e1 * t104 * t173 + 0.18e2 * t18 * t158 - t410);
  tv4rho41 = -0.2e1 * t289 - 0.3e1 * t292 + 0.3e1 / 0.2e1 * t295 - t301 + t353 / 0.2e1 - t356 - 0.6e1 * t357 - 0.4e1 / 0.3e1 * t29 * t3 * (t385 / 0.2e1 + t412 / 0.2e1);

    out->v4rho4[ip*p->dim.v4rho4 + 1] += tv4rho41;

  t420 = t185 * t137;
  t421 = t93 * t257;
  t432 = t5 * t5;
  t435 = 0.1e1 / t157 / t30;
  t451 = my_piecewise3(t10, 0, 0.9e1 / 0.16e2 * t304 * t118 * t66 - 0.3e1 * t189 * t33 * t5 * t70 - 0.3e1 / 0.8e1 * t227 * t73 + 0.6e1 * t65 * t432 * t435 + 0.3e1 * t96 * t70 - 0.9e1 * t96 * t159 - 0.3e1 / 0.8e1 * t150 * t122 * t66 + 0.3e1 / 0.2e1 * t65 * t237 * t33 + 0.3e1 / 0.4e1 * t232 * t73 + t383);
  t477 = my_piecewise3(t17, 0, 0.9e1 / 0.16e2 * t327 * t127 * t79 + 0.3e1 * t204 * t37 * t5 * t70 - 0.3e1 / 0.8e1 * t242 * t82 + 0.6e1 * t78 * t432 * t435 - 0.3e1 * t104 * t70 + 0.9e1 * t104 * t159 - 0.3e1 / 0.8e1 * t166 * t130 * t79 + 0.3e1 / 0.2e1 * t78 * t250 * t37 + 0.3e1 / 0.4e1 * t247 * t82 - t410);
  tv4rho42 = t353 - 0.4e1 * t355 - 0.4e1 * t357 - t292 + t295 - t301 - t420 - 0.4e1 * t421 - 0.4e1 / 0.3e1 * t29 * t3 * (t451 / 0.2e1 + t477 / 0.2e1);

    out->v4rho4[ip*p->dim.v4rho4 + 2] += tv4rho42;

  t487 = t93 * t285;
  t507 = 0.12e2 * t158 + 0.24e2 * t318;
  t511 = my_piecewise3(t10, 0, 0.9e1 / 0.16e2 * t304 * t262 * t33 - 0.9e1 / 0.4e1 * t227 * t71 - 0.9e1 / 0.8e1 * t189 * t122 * t33 + 0.9e1 / 0.2e1 * t192 * t70 * t122 + 0.9e1 / 0.4e1 * t96 * t237 + 0.3e1 / 0.4e1 * t65 * t268 * t33 + 0.3e1 / 0.2e1 * t13 * t507);
  t532 = my_piecewise3(t17, 0, 0.9e1 / 0.16e2 * t327 * t273 * t37 + 0.9e1 / 0.4e1 * t242 * t71 - 0.9e1 / 0.8e1 * t204 * t130 * t37 - 0.9e1 / 0.2e1 * t207 * t70 * t130 + 0.9e1 / 0.4e1 * t104 * t250 + 0.3e1 / 0.4e1 * t78 * t278 * t37 - 0.3e1 / 0.2e1 * t18 * t507);
  tv4rho43 = 0.3e1 / 0.2e1 * t353 - t356 - 0.3e1 * t420 - 0.6e1 * t421 + t295 / 0.2e1 - t301 - 0.2e1 * t487 - 0.4e1 / 0.3e1 * t29 * t3 * (t511 / 0.2e1 + t532 / 0.2e1);

    out->v4rho4[ip*p->dim.v4rho4 + 3] += tv4rho43;

  t541 = t118 * t118;
  t546 = t122 * t122;
  t552 = 0.24e2 * t158 + 0.24e2 * t318;
  t556 = my_piecewise3(t10, 0, 0.9e1 / 0.16e2 * t304 * t541 - 0.9e1 / 0.4e1 * t227 * t122 + 0.9e1 / 0.4e1 * t65 * t546 + 0.3e1 * t96 * t268 + 0.3e1 / 0.2e1 * t13 * t552);
  t557 = t127 * t127;
  t562 = t130 * t130;
  t571 = my_piecewise3(t17, 0, 0.9e1 / 0.16e2 * t327 * t557 - 0.9e1 / 0.4e1 * t242 * t130 + 0.9e1 / 0.4e1 * t78 * t562 + 0.3e1 * t104 * t278 - 0.3e1 / 0.2e1 * t18 * t552);
  tv4rho44 = 0.2e1 * t353 - 0.6e1 * t420 - 0.8e1 * t487 - t301 - 0.4e1 / 0.3e1 * t29 * t3 * (t556 / 0.2e1 + t571 / 0.2e1);

    out->v4rho4[ip*p->dim.v4rho4 + 4] += tv4rho44;

}

#endif

