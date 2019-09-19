/* 
  This file was generated automatically with ./maple2c_new.pl.
  Do not edit this file directly as it can be overwritten!!

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.

  Maple version     : Maple 2018 (X86 64 LINUX)
  Maple source      : ../maple/lda_exc/lda_c_chachiyo_mod.mpl
  Type of functional: lda_exc
*/

#define maple2c_order 3

static inline void
func_unpol(const xc_func_type *p, int order, const double *rho, double *zk, double *vrho, double *v2rho2, double *v3rho3)
{
  double t1, t2, t3, t5, t7, t8, t9, t13;
  double t14, t16, t17, t18, t22, t23, t24, t33;
  double t34, t49, t52, t53, t54, t63;

  lda_c_chachiyo_mod_params *params;

  assert(p->params != NULL);
  params = (lda_c_chachiyo_mod_params * )(p->params);

  t1 = M_CBRT3;
  t2 = t1 * t1;
  t3 = params->bp * t2;
  t5 = POW_1_3(0.1e1 / M_PI);
  t7 = M_CBRT4;
  t8 = 0.1e1 / t5 * t7;
  t9 = POW_1_3(rho[0]);
  t13 = params->bp * t1;
  t14 = t5 * t5;
  t16 = t7 * t7;
  t17 = 0.1e1 / t14 * t16;
  t18 = t9 * t9;
  t22 = 0.1e1 + t3 * t8 * t9 / 0.3e1 + t13 * t17 * t18 / 0.3e1;
  t23 = log(t22);
  if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    *zk = params->ap * t23;

#ifndef XC_DONT_COMPILE_VXC

  if(order < 1) return;


  t24 = rho[0] * params->ap;
  t33 = t3 * t8 / t18 / 0.9e1 + 0.2e1 / 0.9e1 * t13 * t17 / t9;
  t34 = 0.1e1 / t22;
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[0] = t24 * t33 * t34 + (params->ap * t23);

#ifndef XC_DONT_COMPILE_FXC

  if(order < 2) return;


  t49 = -0.2e1 / 0.27e2 * t3 * t8 / t18 / rho[0] - 0.2e1 / 0.27e2 * t13 * t17 / t9 / rho[0];
  t52 = t33 * t33;
  t53 = t22 * t22;
  t54 = 0.1e1 / t53;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[0] = t24 * t49 * t34 - t24 * t52 * t54 + 0.2e1 * params->ap * t33 * t34;

#ifndef XC_DONT_COMPILE_KXC

  if(order < 3) return;


  t63 = rho[0] * rho[0];
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[0] = 0.3e1 * params->ap * t49 * t34 - 0.3e1 * params->ap * t52 * t54 + t24 * (0.10e2 / 0.81e2 * t3 * t8 / t18 / t63 + 0.8e1 / 0.81e2 * t13 * t17 / t9 / t63) * t34 - 0.3e1 * t24 * t49 * t54 * t33 + 0.2e1 * t24 * t52 * t33 / t53 / t22;

#ifndef XC_DONT_COMPILE_LXC

  if(order < 4) return;


#endif

#endif

#endif

#endif


}


static inline void
func_ferr(const xc_func_type *p, int order, const double *rho, double *zk, double *vrho, double *v2rho2, double *v3rho3)
{
  double t1, t2, t3, t5, t7, t8, t9, t13;
  double t14, t16, t17, t18, t22, t23, t24, t33;
  double t34, t49, t52, t53, t54, t63;

  lda_c_chachiyo_mod_params *params;

  assert(p->params != NULL);
  params = (lda_c_chachiyo_mod_params * )(p->params);

  t1 = M_CBRT3;
  t2 = t1 * t1;
  t3 = params->bf * t2;
  t5 = POW_1_3(0.1e1 / M_PI);
  t7 = M_CBRT4;
  t8 = 0.1e1 / t5 * t7;
  t9 = POW_1_3(rho[0]);
  t13 = params->bf * t1;
  t14 = t5 * t5;
  t16 = t7 * t7;
  t17 = 0.1e1 / t14 * t16;
  t18 = t9 * t9;
  t22 = 0.1e1 + t3 * t8 * t9 / 0.3e1 + t13 * t17 * t18 / 0.3e1;
  t23 = log(t22);
  if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    *zk = params->af * t23;

#ifndef XC_DONT_COMPILE_VXC

  if(order < 1) return;


  t24 = rho[0] * params->af;
  t33 = t3 * t8 / t18 / 0.9e1 + 0.2e1 / 0.9e1 * t13 * t17 / t9;
  t34 = 0.1e1 / t22;
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[0] = t24 * t33 * t34 + (params->af * t23);

#ifndef XC_DONT_COMPILE_FXC

  if(order < 2) return;


  t49 = -0.2e1 / 0.27e2 * t3 * t8 / t18 / rho[0] - 0.2e1 / 0.27e2 * t13 * t17 / t9 / rho[0];
  t52 = t33 * t33;
  t53 = t22 * t22;
  t54 = 0.1e1 / t53;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[0] = t24 * t49 * t34 - t24 * t52 * t54 + 0.2e1 * params->af * t33 * t34;

#ifndef XC_DONT_COMPILE_KXC

  if(order < 3) return;


  t63 = rho[0] * rho[0];
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[0] = 0.3e1 * params->af * t49 * t34 - 0.3e1 * params->af * t52 * t54 + t24 * (0.10e2 / 0.81e2 * t3 * t8 / t18 / t63 + 0.8e1 / 0.81e2 * t13 * t17 / t9 / t63) * t34 - 0.3e1 * t24 * t49 * t54 * t33 + 0.2e1 * t24 * t52 * t33 / t53 / t22;

#ifndef XC_DONT_COMPILE_LXC

  if(order < 4) return;


#endif

#endif

#endif

#endif


}


static inline void
func_pol(const xc_func_type *p, int order, const double *rho, double *zk, double *vrho, double *v2rho2, double *v3rho3)
{
  double t1, t2, t3, t5, t7, t8, t9, t10;
  double t11, t14, t15, t17, t18, t19, t20, t23;
  double t24, t25, t26, t29, t32, t33, t35, t36;
  double t37, t38, t39, t40, t41, t42, t43, t44;
  double t46, t47, t50, t51, t53, t57, t60, t62;
  double t63, t68, t70, t72, t73, t74, t75, t76;
  double t77, t78, t79, t81, t82, t85, t86, t87;
  double t90, t92, t95, t96, t97, t100, t101, t105;
  double t109, t113, t114, t115, t117, t118, t119, t124;
  double t126, t128, t129, t131, t132, t133, t134, t135;
  double t136, t137, t138, t139, t141, t142, t146, t147;
  double t149, t153, t154, t157, t160, t161, t162, t166;
  double t168, t169, t171, t174, t177, t180, t183, t184;
  double t189, t190, t191, t192, t193, t197, t200, t203;
  double t206, t207, t208, t211, t212, t213, t219, t224;
  double t229, t232, t238, t256, t257, t258, t260, t261;
  double t263, t271, t273, t280, t281, t282, t284, t287;
  double t289, t305, t306, t310, t313, t315, t354, t359;
  double t367, t382, t399, t419;

  lda_c_chachiyo_mod_params *params;

  assert(p->params != NULL);
  params = (lda_c_chachiyo_mod_params * )(p->params);

  t1 = M_CBRT3;
  t2 = t1 * t1;
  t3 = params->bp * t2;
  t5 = POW_1_3(0.1e1 / M_PI);
  t7 = M_CBRT4;
  t8 = 0.1e1 / t5 * t7;
  t9 = rho[0] + rho[1];
  t10 = POW_1_3(t9);
  t11 = t8 * t10;
  t14 = params->bp * t1;
  t15 = t5 * t5;
  t17 = t7 * t7;
  t18 = 0.1e1 / t15 * t17;
  t19 = t10 * t10;
  t20 = t18 * t19;
  t23 = 0.1e1 + t3 * t11 / 0.3e1 + t14 * t20 / 0.3e1;
  t24 = log(t23);
  t25 = params->ap * t24;
  t26 = params->bf * t2;
  t29 = params->bf * t1;
  t32 = 0.1e1 + t26 * t11 / 0.3e1 + t29 * t20 / 0.3e1;
  t33 = log(t32);
  t35 = params->af * t33 - t25;
  t36 = rho[0] - rho[1];
  t37 = 0.1e1 / t9;
  t38 = t36 * t37;
  t39 = 0.1e1 + t38;
  t40 = POW_1_3(t39);
  t41 = t40 * t40;
  t42 = 0.1e1 - t38;
  t43 = POW_1_3(t42);
  t44 = t43 * t43;
  t46 = t41 / 0.2e1 + t44 / 0.2e1;
  t47 = t46 * t46;
  t50 = -0.2e1 * t47 * t46 + 0.2e1;
  t51 = t35 * t50;
  if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    *zk = t25 + t51;

#ifndef XC_DONT_COMPILE_VXC

  if(order < 1) return;


  t53 = t8 / t19;
  t57 = t18 / t10;
  t60 = t3 * t53 / 0.9e1 + 0.2e1 / 0.9e1 * t14 * t57;
  t62 = 0.1e1 / t23;
  t63 = params->ap * t60 * t62;
  t68 = t26 * t53 / 0.9e1 + 0.2e1 / 0.9e1 * t29 * t57;
  t70 = 0.1e1 / t32;
  t72 = params->af * t68 * t70 - t63;
  t73 = t72 * t50;
  t74 = t35 * t47;
  t75 = 0.1e1 / t40;
  t76 = t9 * t9;
  t77 = 0.1e1 / t76;
  t78 = t36 * t77;
  t79 = t37 - t78;
  t81 = 0.1e1 / t43;
  t82 = -t79;
  t85 = t75 * t79 / 0.3e1 + t81 * t82 / 0.3e1;
  t86 = t74 * t85;
  t87 = 0.6e1 * t86;
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[0] = t25 + t51 + t9 * (t63 + t73 - t87);

  t90 = -t37 - t78;
  t92 = -t90;
  t95 = t75 * t90 / 0.3e1 + t81 * t92 / 0.3e1;
  t96 = t74 * t95;
  t97 = 0.6e1 * t96;
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[1] = t25 + t51 + t9 * (t63 + t73 - t97);

#ifndef XC_DONT_COMPILE_FXC

  if(order < 2) return;


  t100 = 0.2e1 * t63;
  t101 = 0.2e1 * t73;
  t105 = t8 / t19 / t9;
  t109 = t18 / t10 / t9;
  t113 = params->ap * (-0.2e1 / 0.27e2 * t3 * t105 - 0.2e1 / 0.27e2 * t14 * t109);
  t114 = t113 * t62;
  t115 = t60 * t60;
  t117 = t23 * t23;
  t118 = 0.1e1 / t117;
  t119 = params->ap * t115 * t118;
  t124 = params->af * (-0.2e1 / 0.27e2 * t26 * t105 - 0.2e1 / 0.27e2 * t29 * t109);
  t126 = t68 * t68;
  t128 = t32 * t32;
  t129 = 0.1e1 / t128;
  t131 = -params->af * t126 * t129 + t124 * t70 - t114 + t119;
  t132 = t131 * t50;
  t133 = t72 * t47;
  t134 = t133 * t85;
  t135 = 0.12e2 * t134;
  t136 = t35 * t46;
  t137 = t85 * t85;
  t138 = t136 * t137;
  t139 = 0.12e2 * t138;
  t141 = 0.1e1 / t40 / t39;
  t142 = t79 * t79;
  t146 = 0.1e1 / t76 / t9;
  t147 = t36 * t146;
  t149 = -0.2e1 * t77 + 0.2e1 * t147;
  t153 = 0.1e1 / t43 / t42;
  t154 = t82 * t82;
  t157 = -t149;
  t160 = -t141 * t142 / 0.9e1 + t75 * t149 / 0.3e1 - t153 * t154 / 0.9e1 + t81 * t157 / 0.3e1;
  t161 = t74 * t160;
  t162 = 0.6e1 * t161;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[0] = t100 + t101 - 0.12e2 * t86 + t9 * (t114 - t119 + t132 - t135 - t139 - t162);

  t166 = t133 * t95;
  t168 = t95 * t85;
  t169 = t136 * t168;
  t171 = t141 * t90;
  t174 = t75 * t36;
  t177 = t153 * t92;
  t180 = t81 * t36;
  t183 = -t171 * t79 / 0.9e1 + 0.2e1 / 0.3e1 * t174 * t146 - t177 * t82 / 0.9e1 - 0.2e1 / 0.3e1 * t180 * t146;
  t184 = t74 * t183;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[1] = t100 + t101 - t87 - t97 + t9 * (t114 - t119 + t132 - 0.6e1 * t134 - 0.6e1 * t166 - 0.12e2 * t169 - 0.6e1 * t184);

  t189 = 0.12e2 * t166;
  t190 = t95 * t95;
  t191 = t136 * t190;
  t192 = 0.12e2 * t191;
  t193 = t90 * t90;
  t197 = 0.2e1 * t77 + 0.2e1 * t147;
  t200 = t92 * t92;
  t203 = -t197;
  t206 = -t141 * t193 / 0.9e1 + t75 * t197 / 0.3e1 - t153 * t200 / 0.9e1 + t81 * t203 / 0.3e1;
  t207 = t74 * t206;
  t208 = 0.6e1 * t207;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[2] = t100 + t101 - 0.12e2 * t96 + t9 * (t114 - t119 + t132 - t189 - t192 - t208);

#ifndef XC_DONT_COMPILE_KXC

  if(order < 3) return;


  t211 = 0.3e1 * t114;
  t212 = 0.3e1 * t119;
  t213 = 0.3e1 * t132;
  t219 = t8 / t19 / t76;
  t224 = t18 / t10 / t76;
  t229 = params->ap * (0.10e2 / 0.81e2 * t3 * t219 + 0.8e1 / 0.81e2 * t14 * t224) * t62;
  t232 = 0.3e1 * t113 * t118 * t60;
  t238 = 0.2e1 * params->ap * t115 * t60 / t117 / t23;
  t256 = (params->af * (0.10e2 / 0.81e2 * t26 * t219 + 0.8e1 / 0.81e2 * t29 * t224) * t70 - 0.3e1 * t124 * t129 * t68 + 0.2e1 * params->af * t126 * t68 / t128 / t32 - t229 + t232 - t238) * t50;
  t257 = t131 * t47;
  t258 = t257 * t85;
  t260 = t72 * t46;
  t261 = t260 * t137;
  t263 = t133 * t160;
  t271 = t39 * t39;
  t273 = 0.1e1 / t40 / t271;
  t280 = t76 * t76;
  t281 = 0.1e1 / t280;
  t282 = t36 * t281;
  t284 = 0.6e1 * t146 - 0.6e1 * t282;
  t287 = t42 * t42;
  t289 = 0.1e1 / t43 / t287;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[0] = t211 - t212 + t213 - 0.36e2 * t134 - 0.36e2 * t138 - 0.18e2 * t161 + t9 * (t229 - t232 + t238 + t256 - 0.18e2 * t258 - 0.36e2 * t261 - 0.18e2 * t263 - 0.12e2 * t35 * t137 * t85 - 0.36e2 * t136 * t85 * t160 - 0.6e1 * t74 * (0.4e1 / 0.27e2 * t273 * t142 * t79 - t141 * t79 * t149 / 0.3e1 + t75 * t284 / 0.3e1 + 0.4e1 / 0.27e2 * t289 * t154 * t82 - t153 * t82 * t157 / 0.3e1 - t81 * t284 / 0.3e1));

  t305 = 0.24e2 * t169;
  t306 = 0.12e2 * t184;
  t310 = t257 * t95;
  t313 = 0.24e2 * t260 * t168;
  t315 = 0.12e2 * t133 * t183;
  t354 = t229 - t232 + t238 + t256 - 0.12e2 * t258 - 0.12e2 * t261 - 0.6e1 * t263 - 0.6e1 * t310 - t313 - t315 - 0.12e2 * t35 * t137 * t95 - 0.24e2 * t136 * t183 * t85 - 0.12e2 * t136 * t95 * t160 - 0.6e1 * t74 * (0.4e1 / 0.27e2 * t273 * t90 * t142 - 0.4e1 / 0.9e1 * t141 * t36 * t146 * t79 - t171 * t149 / 0.9e1 + 0.2e1 / 0.3e1 * t75 * t146 - 0.2e1 * t174 * t281 + 0.4e1 / 0.27e2 * t289 * t92 * t154 + 0.4e1 / 0.9e1 * t153 * t36 * t146 * t82 - t177 * t157 / 0.9e1 - 0.2e1 / 0.3e1 * t81 * t146 + 0.2e1 * t180 * t281);
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[1] = t9 * t354 - 0.24e2 * t134 - t139 - t162 - t189 + t211 - t212 + t213 - t305 - t306;

  t359 = t260 * t190;
  t367 = t133 * t206;
  t382 = -0.2e1 * t146 - 0.6e1 * t282;
  t399 = t229 - t232 + t238 + t256 - 0.6e1 * t258 - 0.12e2 * t310 - t313 - t315 - 0.12e2 * t359 - 0.12e2 * t35 * t85 * t190 - 0.24e2 * t136 * t95 * t183 - 0.6e1 * t367 - 0.12e2 * t136 * t206 * t85 - 0.6e1 * t74 * (0.4e1 / 0.27e2 * t273 * t193 * t79 - 0.4e1 / 0.9e1 * t171 * t147 - t141 * t197 * t79 / 0.9e1 + t75 * t382 / 0.3e1 + 0.4e1 / 0.27e2 * t289 * t200 * t82 + 0.4e1 / 0.9e1 * t177 * t147 - t153 * t203 * t82 / 0.9e1 - t81 * t382 / 0.3e1);
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[2] = t9 * t399 - t135 - 0.24e2 * t166 - t192 - t208 + t211 - t212 + t213 - t305 - t306;

  t419 = -0.6e1 * t146 - 0.6e1 * t282;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[3] = t211 - t212 + t213 - 0.36e2 * t166 - 0.36e2 * t191 - 0.18e2 * t207 + t9 * (t229 - t232 + t238 + t256 - 0.18e2 * t310 - 0.36e2 * t359 - 0.18e2 * t367 - 0.12e2 * t35 * t190 * t95 - 0.36e2 * t136 * t95 * t206 - 0.6e1 * t74 * (0.4e1 / 0.27e2 * t273 * t193 * t90 - t171 * t197 / 0.3e1 + t75 * t419 / 0.3e1 + 0.4e1 / 0.27e2 * t289 * t200 * t92 - t177 * t203 / 0.3e1 - t81 * t419 / 0.3e1));

#ifndef XC_DONT_COMPILE_LXC

  if(order < 4) return;


#endif

#endif

#endif

#endif


}
