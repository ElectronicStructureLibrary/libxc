/* 
  This file was generated automatically with ./scripts/maple2c_new.pl.
  Do not edit this file directly as it can be overwritten!!

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.

  Maple version     : Maple 2016 (X86 64 LINUX)
  Maple source      : ./maple/lda_c_vwn_1.mpl
  Type of functional: lda_exc
*/

#define maple2c_order 3

static inline void
func_unpol(const xc_func_type *p, int order, const double *rho, double *zk, double *vrho, double *v2rho2, double *v3rho3)
{
  double t1, t2, t3, t4, t5, t6, t7, t8;
  double t9, t10, t12, t14, t15, t19, t20, t21;
  double t24, t25, t27, t28, t30, t31, t33, t34;
  double t38, t39, t40, t41, t42, t44, t45, t46;
  double t50, t55, t57, t58, t59, t60, t61, t63;
  double t64, t66, t68, t69, t72, t74, t75, t78;
  double t80, t81, t82, t83, t90, t92, t93, t95;
  double t97, t102, t103, t104, t108, t111, t112, t113;
  double t114, t115, t117, t124, t130, t131, t134, t135;
  double t136, t138, t139, t143, t144, t146, t149, t152;
  double t156, t158, t161, t162, t163, t164, t166, t169;
  double t170, t175, t176, t178, t183, t187, t188, t189;
  double t192, t193, t194, t195, t196, t198, t212, t214;
  double t220, t225, t227, t231, t251, t261, t263, t264;
  double t265, t268, t272, t275, t276, t278, t281, t292;
  double t294, t335, t342, t372, t382, t407;


  t1 = M_CBRT3;
  t2 = 0.1e1 / M_PI;
  t3 = POW_1_3(t2);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = POW_1_3(rho[0]);
  t8 = 0.1e1 / t7;
  t9 = t6 * t8;
  t10 = t4 * t9;
  t12 = sqrt(t10);
  t14 = t10 / 0.4e1 + 0.18637200000000000000e1 * t12 + 0.129352e2;
  t15 = 0.1e1 / t14;
  t19 = log(t4 * t9 * t15 / 0.4e1);
  t20 = 0.310907e-1 * t19;
  t21 = t12 + 0.372744e1;
  t24 = atan(0.61519908197590802322e1 / t21);
  t25 = 0.38783294878113014393e-1 * t24;
  t27 = t12 / 0.2e1 + 0.10498e0;
  t28 = t27 * t27;
  t30 = log(t28 * t15);
  t31 = 0.96902277115443742139e-3 * t30;
  if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    *zk = t20 + t25 + t31;

  if(order < 1) return;


  t33 = 0.1e1 / t7 / rho[0];
  t34 = t6 * t33;
  t38 = t4 * t6;
  t39 = t14 * t14;
  t40 = 0.1e1 / t39;
  t41 = t8 * t40;
  t42 = t4 * t34;
  t44 = 0.1e1 / t12;
  t45 = t44 * t1;
  t46 = t3 * t6;
  t50 = -t42 / 0.12e2 - 0.31062000000000000000e0 * t45 * t46 * t33;
  t55 = t1 * t1;
  t57 = 0.1e1 / t3;
  t58 = (-t4 * t34 * t15 / 0.12e2 - t38 * t41 * t50 / 0.4e1) * t55 * t57;
  t59 = t5 * t7;
  t60 = t59 * t14;
  t61 = t58 * t60;
  t63 = t21 * t21;
  t64 = 0.1e1 / t63;
  t66 = t64 * t44 * t1;
  t68 = 0.37846991046400000000e2 * t64 + 0.1e1;
  t69 = 0.1e1 / t68;
  t72 = t66 * t46 * t33 * t69;
  t74 = t27 * t15;
  t75 = t74 * t44;
  t78 = t28 * t40;
  t80 = -t75 * t42 / 0.6e1 - t78 * t50;
  t81 = 0.1e1 / t28;
  t82 = t80 * t81;
  t83 = t82 * t14;
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[0] = t20 + t25 + t31 + rho[0] * (0.10363566666666666667e-1 * t61 + 0.39765745675026770179e-1 * t72 + 0.96902277115443742139e-3 * t83);

  if(order < 2) return;


  t90 = rho[0] * rho[0];
  t92 = 0.1e1 / t7 / t90;
  t93 = t6 * t92;
  t95 = t4 * t93 * t15;
  t97 = t33 * t40;
  t102 = 0.1e1 / t39 / t14;
  t103 = t8 * t102;
  t104 = t50 * t50;
  t108 = t4 * t93;
  t111 = 0.1e1 / t12 / t10;
  t112 = t111 * t55;
  t113 = t3 * t3;
  t114 = t113 * t5;
  t115 = t7 * t7;
  t117 = 0.1e1 / t115 / t90;
  t124 = t108 / 0.9e1 - 0.20708000000000000000e0 * t112 * t114 * t117 + 0.41416000000000000000e0 * t45 * t46 * t92;
  t130 = (t95 / 0.9e1 + t38 * t97 * t50 / 0.6e1 + t38 * t103 * t104 / 0.2e1 - t38 * t41 * t124 / 0.4e1) * t55 * t57;
  t131 = t130 * t60;
  t134 = t5 / t115;
  t135 = t134 * t14;
  t136 = t58 * t135;
  t138 = t59 * t50;
  t139 = t58 * t138;
  t143 = 0.1e1 / t63 / t21 * t1;
  t144 = t143 * t3;
  t146 = t144 * t93 * t69;
  t149 = t64 * t111 * t55;
  t152 = t149 * t114 * t117 * t69;
  t156 = t66 * t46 * t92 * t69;
  t158 = t63 * t63;
  t161 = 0.1e1 / t158 / t21 * t1;
  t162 = t161 * t3;
  t163 = t68 * t68;
  t164 = 0.1e1 / t163;
  t166 = t162 * t93 * t164;
  t169 = t27 * t40;
  t170 = t169 * t45;
  t175 = t74 * t111;
  t176 = t55 * t113;
  t178 = t176 * t5 * t117;
  t183 = t28 * t102;
  t187 = t95 / 0.72e2 + t170 * t46 * t33 * t50 / 0.3e1 - t175 * t178 / 0.9e1 + 0.2e1 / 0.9e1 * t75 * t108 + 0.2e1 * t183 * t104 - t78 * t124;
  t188 = t187 * t81;
  t189 = t188 * t14;
  t192 = 0.1e1 / t28 / t27;
  t193 = t80 * t192;
  t194 = t14 * t44;
  t195 = t193 * t194;
  t196 = t195 * t42;
  t198 = t82 * t50;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[0] = 0.20727133333333333334e-1 * t61 + 0.79531491350053540358e-1 * t72 + 0.19380455423088748428e-2 * t83 + rho[0] * (0.10363566666666666667e-1 * t131 + 0.34545222222222222223e-2 * t136 + 0.10363566666666666667e-1 * t139 + 0.13255248558342256726e-1 * t146 + 0.26510497116684513453e-1 * t152 - 0.53020994233369026905e-1 * t156 - 0.50167127350538589834e0 * t166 + 0.96902277115443742139e-3 * t189 + 0.16150379519240623690e-3 * t196 + 0.96902277115443742139e-3 * t198);

  if(order < 3) return;


  t212 = t90 * rho[0];
  t214 = 0.1e1 / t115 / t212;
  t220 = 0.1e1 / t7 / t212;
  t225 = t6 * t220;
  t227 = t4 * t225 * t15;
  t231 = t38 * t92 * t40 * t50;
  t251 = 0.1e1 / t113;
  t261 = 0.1e1 / t12 / t176 / t134 / 0.4e1;
  t263 = t90 * t90;
  t264 = 0.1e1 / t263;
  t265 = t261 * t2 * t264;
  t268 = t5 * t214;
  t272 = t4 * t225;
  t275 = t39 * t39;
  t276 = 0.1e1 / t275;
  t278 = t104 * t50;
  t281 = t50 * t124;
  t292 = -0.7e1 / 0.27e2 * t272 - 0.12424800000000000000e1 * t265 + 0.82832000000000000000e0 * t112 * t114 * t214 - 0.96637333333333333333e0 * t45 * t46 * t220;
  t294 = -0.11e2 / 0.216e3 * t227 - t231 / 0.24e2 - t27 * t102 * t45 * t46 * t33 * t104 + t169 * t112 * t114 * t117 * t50 / 0.3e1 - 0.2e1 / 0.3e1 * t170 * t46 * t92 * t50 + t170 * t46 * t33 * t124 / 0.2e1 + t1 * t251 * t6 * t220 * t2 * t15 / 0.432e3 - 0.2e1 / 0.3e1 * t74 * t265 + 0.4e1 / 0.9e1 * t175 * t176 * t268 - 0.14e2 / 0.27e2 * t75 * t272 - 0.6e1 * t28 * t276 * t278 + 0.6e1 * t183 * t281 - t78 * t292;
  t335 = t225 * t164;
  t342 = -0.10604198846673805381e0 * t149 * t114 * t214 * t69 + 0.12371565321119439611e0 * t66 * t46 * t220 * t69 + 0.96902277115443742139e-3 * t294 * t81 * t14 + 0.19380455423088748428e-2 * t188 * t50 + 0.96902277115443742139e-3 * t82 * t124 + 0.10363566666666666667e-1 * (-0.7e1 / 0.27e2 * t227 - t231 / 0.3e1 - t38 * t33 * t102 * t104 / 0.2e1 + t38 * t97 * t124 / 0.4e1 - 0.3e1 / 0.2e1 * t38 * t8 * t276 * t278 + 0.3e1 / 0.2e1 * t38 * t103 * t281 - t38 * t41 * t292 / 0.4e1) * t55 * t57 * t60 + 0.20727133333333333334e-1 * t130 * t138 + 0.69090444444444444446e-2 * t58 * t134 * t50 + 0.10363566666666666667e-1 * t58 * t59 * t124 - 0.26614487661862784322e-1 * t161 * t251 * t335 + 0.18394613361864149606e1 * t162 * t335 + 0.69090444444444444446e-2 * t130 * t135;
  t372 = t158 * t158;
  t382 = t28 * t28;
  t407 = -0.23030148148148148149e-2 * t58 * t5 / t115 / rho[0] * t14 - 0.48602578047254941329e-1 * t144 * t225 * t69 + 0.26510497116684513452e-1 / t158 * t55 * t113 * t268 * t69 * t44 - 0.23411326096918008589e1 / t158 / t63 * t55 * t113 * t268 * t164 * t44 + 0.22092080930570427878e-2 * t143 * t251 * t225 * t2 * t69 + 0.50631328524251801694e2 / t372 * t55 * t113 * t268 / t163 / t68 * t44 + 0.40375948798101559225e-4 * t80 / t382 * t14 * t108 + 0.15906298270010708072e0 * t64 * t261 * t2 * t264 * t69 - 0.21533839358987498253e-3 * t195 * t108 + 0.32300759038481247380e-3 * t187 * t192 * t194 * t42 + 0.32300759038481247380e-3 * t193 * t50 * t44 * t42 + 0.10766919679493749127e-3 * t193 * t14 * t111 * t178;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[0] = 0.31090700000000000001e-1 * t131 + 0.10363566666666666667e-1 * t136 + 0.31090700000000000001e-1 * t139 + 0.39765745675026770179e-1 * t146 + 0.79531491350053540358e-1 * t152 - 0.15906298270010708072e0 * t156 - 0.15050138205161576950e1 * t166 + 0.29070683134633122642e-2 * t189 + 0.48451138557721871070e-3 * t196 + 0.29070683134633122642e-2 * t198 + rho[0] * (t342 + t407);

  if(order < 4) return;



}


static inline void
func_ferr(const xc_func_type *p, int order, const double *rho, double *zk, double *vrho, double *v2rho2, double *v3rho3)
{
  double t1, t2, t3, t4, t5, t6, t7, t8;
  double t9, t10, t12, t14, t15, t19, t20, t21;
  double t24, t25, t27, t28, t30, t31, t33, t34;
  double t38, t39, t40, t41, t42, t44, t45, t46;
  double t50, t55, t57, t58, t59, t60, t61, t63;
  double t64, t66, t68, t69, t72, t74, t75, t78;
  double t80, t81, t82, t83, t90, t92, t93, t95;
  double t97, t102, t103, t104, t108, t111, t112, t113;
  double t114, t115, t117, t124, t130, t131, t134, t135;
  double t136, t138, t139, t143, t144, t146, t149, t152;
  double t156, t158, t161, t162, t163, t164, t166, t169;
  double t170, t175, t176, t178, t183, t187, t188, t189;
  double t192, t193, t194, t195, t196, t198, t212, t214;
  double t220, t225, t227, t231, t251, t261, t263, t264;
  double t265, t268, t272, t275, t276, t278, t281, t292;
  double t294, t335, t342, t372, t382, t407;


  t1 = M_CBRT3;
  t2 = 0.1e1 / M_PI;
  t3 = POW_1_3(t2);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = POW_1_3(rho[0]);
  t8 = 0.1e1 / t7;
  t9 = t6 * t8;
  t10 = t4 * t9;
  t12 = sqrt(t10);
  t14 = t10 / 0.4e1 + 0.35302100000000000000e1 * t12 + 0.180578e2;
  t15 = 0.1e1 / t14;
  t19 = log(t4 * t9 * t15 / 0.4e1);
  t20 = 0.1554535e-1 * t19;
  t21 = t12 + 0.706042e1;
  t24 = atan(0.47309269095601128300e1 / t21);
  t25 = 0.52491393169780936218e-1 * t24;
  t27 = t12 / 0.2e1 + 0.32500e0;
  t28 = t27 * t27;
  t30 = log(t28 * t15);
  t31 = 0.22478670955426118383e-2 * t30;
  if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    *zk = t20 + t25 + t31;

  if(order < 1) return;


  t33 = 0.1e1 / t7 / rho[0];
  t34 = t6 * t33;
  t38 = t4 * t6;
  t39 = t14 * t14;
  t40 = 0.1e1 / t39;
  t41 = t8 * t40;
  t42 = t4 * t34;
  t44 = 0.1e1 / t12;
  t45 = t44 * t1;
  t46 = t3 * t6;
  t50 = -t42 / 0.12e2 - 0.58836833333333333333e0 * t45 * t46 * t33;
  t55 = t1 * t1;
  t57 = 0.1e1 / t3;
  t58 = (-t4 * t34 * t15 / 0.12e2 - t38 * t41 * t50 / 0.4e1) * t55 * t57;
  t59 = t5 * t7;
  t60 = t59 * t14;
  t61 = t58 * t60;
  t63 = t21 * t21;
  t64 = 0.1e1 / t63;
  t66 = t64 * t44 * t1;
  t68 = 0.22381669423600000000e2 * t64 + 0.1e1;
  t69 = 0.1e1 / t68;
  t72 = t66 * t46 * t33 * t69;
  t74 = t27 * t15;
  t75 = t74 * t44;
  t78 = t28 * t40;
  t80 = -t75 * t42 / 0.6e1 - t78 * t50;
  t81 = 0.1e1 / t28;
  t82 = t80 * t81;
  t83 = t82 * t14;
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[0] = t20 + t25 + t31 + rho[0] * (0.51817833333333333333e-2 * t61 + 0.41388824077869423261e-1 * t72 + 0.22478670955426118383e-2 * t83);

  if(order < 2) return;


  t90 = rho[0] * rho[0];
  t92 = 0.1e1 / t7 / t90;
  t93 = t6 * t92;
  t95 = t4 * t93 * t15;
  t97 = t33 * t40;
  t102 = 0.1e1 / t39 / t14;
  t103 = t8 * t102;
  t104 = t50 * t50;
  t108 = t4 * t93;
  t111 = 0.1e1 / t12 / t10;
  t112 = t111 * t55;
  t113 = t3 * t3;
  t114 = t113 * t5;
  t115 = t7 * t7;
  t117 = 0.1e1 / t115 / t90;
  t124 = t108 / 0.9e1 - 0.39224555555555555555e0 * t112 * t114 * t117 + 0.78449111111111111110e0 * t45 * t46 * t92;
  t130 = (t95 / 0.9e1 + t38 * t97 * t50 / 0.6e1 + t38 * t103 * t104 / 0.2e1 - t38 * t41 * t124 / 0.4e1) * t55 * t57;
  t131 = t130 * t60;
  t134 = t5 / t115;
  t135 = t134 * t14;
  t136 = t58 * t135;
  t138 = t59 * t50;
  t139 = t58 * t138;
  t143 = 0.1e1 / t63 / t21 * t1;
  t144 = t143 * t3;
  t146 = t144 * t93 * t69;
  t149 = t64 * t111 * t55;
  t152 = t149 * t114 * t117 * t69;
  t156 = t66 * t46 * t92 * t69;
  t158 = t63 * t63;
  t161 = 0.1e1 / t158 / t21 * t1;
  t162 = t161 * t3;
  t163 = t68 * t68;
  t164 = 0.1e1 / t163;
  t166 = t162 * t93 * t164;
  t169 = t27 * t40;
  t170 = t169 * t45;
  t175 = t74 * t111;
  t176 = t55 * t113;
  t178 = t176 * t5 * t117;
  t183 = t28 * t102;
  t187 = t95 / 0.72e2 + t170 * t46 * t33 * t50 / 0.3e1 - t175 * t178 / 0.9e1 + 0.2e1 / 0.9e1 * t75 * t108 + 0.2e1 * t183 * t104 - t78 * t124;
  t188 = t187 * t81;
  t189 = t188 * t14;
  t192 = 0.1e1 / t28 / t27;
  t193 = t80 * t192;
  t194 = t14 * t44;
  t195 = t193 * t194;
  t196 = t195 * t42;
  t198 = t82 * t50;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[0] = 0.10363566666666666667e-1 * t61 + 0.82777648155738846522e-1 * t72 + 0.44957341910852236766e-2 * t83 + rho[0] * (0.51817833333333333333e-2 * t131 + 0.17272611111111111111e-2 * t136 + 0.51817833333333333333e-2 * t139 + 0.13796274692623141087e-1 * t146 + 0.27592549385246282174e-1 * t152 - 0.55185098770492564348e-1 * t156 - 0.30878365944746984535e0 * t166 + 0.22478670955426118383e-2 * t189 + 0.37464451592376863972e-3 * t196 + 0.22478670955426118383e-2 * t198);

  if(order < 3) return;


  t212 = t90 * rho[0];
  t214 = 0.1e1 / t115 / t212;
  t220 = 0.1e1 / t7 / t212;
  t225 = t6 * t220;
  t227 = t4 * t225 * t15;
  t231 = t38 * t92 * t40 * t50;
  t251 = 0.1e1 / t113;
  t261 = 0.1e1 / t12 / t176 / t134 / 0.4e1;
  t263 = t90 * t90;
  t264 = 0.1e1 / t263;
  t265 = t261 * t2 * t264;
  t268 = t5 * t214;
  t272 = t4 * t225;
  t275 = t39 * t39;
  t276 = 0.1e1 / t275;
  t278 = t104 * t50;
  t281 = t50 * t124;
  t292 = -0.7e1 / 0.27e2 * t272 - 0.23534733333333333333e1 * t265 + 0.15689822222222222222e1 * t112 * t114 * t214 - 0.18304792592592592592e1 * t45 * t46 * t220;
  t294 = -0.11e2 / 0.216e3 * t227 - t231 / 0.24e2 - t27 * t102 * t45 * t46 * t33 * t104 + t169 * t112 * t114 * t117 * t50 / 0.3e1 - 0.2e1 / 0.3e1 * t170 * t46 * t92 * t50 + t170 * t46 * t33 * t124 / 0.2e1 + t1 * t251 * t6 * t220 * t2 * t15 / 0.432e3 - 0.2e1 / 0.3e1 * t74 * t265 + 0.4e1 / 0.9e1 * t175 * t176 * t268 - 0.14e2 / 0.27e2 * t75 * t272 - 0.6e1 * t28 * t276 * t278 + 0.6e1 * t183 * t281 - t78 * t292;
  t335 = t225 * t164;
  t342 = -0.11037019754098512870e0 * t149 * t114 * t214 * t69 + 0.12876523046448265015e0 * t66 * t46 * t220 * t69 + 0.22478670955426118383e-2 * t294 * t81 * t14 + 0.44957341910852236766e-2 * t188 * t50 + 0.22478670955426118383e-2 * t82 * t124 + 0.51817833333333333333e-2 * (-0.7e1 / 0.27e2 * t227 - t231 / 0.3e1 - t38 * t33 * t102 * t104 / 0.2e1 + t38 * t97 * t124 / 0.4e1 - 0.3e1 / 0.2e1 * t38 * t8 * t276 * t278 + 0.3e1 / 0.2e1 * t38 * t103 * t281 - t38 * t41 * t292 / 0.4e1) * t55 * t57 * t60 + 0.10363566666666666667e-1 * t130 * t138 + 0.34545222222222222222e-2 * t58 * t134 * t50 + 0.51817833333333333333e-2 * t58 * t59 * t124 - 0.16381481915689750931e-1 * t161 * t251 * t335 + 0.11322067513073894330e1 * t162 * t335 + 0.34545222222222222222e-2 * t130 * t135;
  t372 = t158 * t158;
  t382 = t28 * t28;
  t407 = -0.11515074074074074074e-2 * t58 * t5 / t115 / rho[0] * t14 - 0.50586340539618183985e-1 * t144 * t225 * t69 + 0.27592549385246282174e-1 / t158 * t55 * t113 * t268 * t69 * t44 - 0.14409904107548592783e1 / t158 / t63 * t55 * t113 * t268 * t164 * t44 + 0.22993791154371901812e-2 * t143 * t251 * t225 * t2 * t69 + 0.18429583437767336288e2 / t372 * t55 * t113 * t268 / t163 / t68 * t44 + 0.93661128980942159930e-4 * t80 / t382 * t14 * t108 + 0.16555529631147769304e0 * t64 * t261 * t2 * t264 * t69 - 0.49952602123169151963e-3 * t195 * t108 + 0.74928903184753727944e-3 * t187 * t192 * t194 * t42 + 0.74928903184753727944e-3 * t193 * t50 * t44 * t42 + 0.24976301061584575981e-3 * t193 * t14 * t111 * t178;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[0] = 0.15545350000000000000e-1 * t131 + 0.51817833333333333334e-2 * t136 + 0.15545350000000000000e-1 * t139 + 0.41388824077869423261e-1 * t146 + 0.82777648155738846522e-1 * t152 - 0.16555529631147769305e0 * t156 - 0.92635097834240953604e0 * t166 + 0.67436012866278355149e-2 * t189 + 0.11239335477713059192e-2 * t196 + 0.67436012866278355149e-2 * t198 + rho[0] * (t342 + t407);

  if(order < 4) return;



}


static inline void
func_pol(const xc_func_type *p, int order, const double *rho, double *zk, double *vrho, double *v2rho2, double *v3rho3)
{
  double t1, t2, t3, t4, t5, t6, t7, t8;
  double t9, t10, t11, t12, t13, t15, t16, t20;
  double t22, t25, t27, t28, t29, t31, t33, t34;
  double t35, t36, t37, t38, t40, t41, t43, t44;
  double t47, t49, t50, t52, t53, t57, t59, t62;
  double t64, t65, t67, t69, t71, t73, t74, t78;
  double t79, t80, t81, t82, t83, t84, t85, t86;
  double t88, t90, t95, t97, t98, t99, t100, t103;
  double t104, t106, t108, t109, t114, t115, t118, t120;
  double t121, t122, t125, t126, t127, t128, t129, t130;
  double t132, t135, t137, t141, t142, t143, t145, t151;
  double t152, t155, t156, t158, t160, t161, t166, t167;
  double t170, t172, t173, t174, t177, t179, t181, t184;
  double t186, t189, t191, t193, t196, t198, t201, t202;
  double t204, t206, t211, t212, t213, t217, t218, t220;
  double t221, t222, t223, t224, t226, t228, t231, t233;
  double t239, t243, t244, t247, t252, t253, t258, t267;
  double t270, t271, t272, t273, t278, t279, t284, t285;
  double t287, t292, t296, t297, t301, t302, t303, t304;
  double t309, t310, t312, t313, t314, t315, t316, t319;
  double t320, t321, t323, t326, t327, t328, t331, t334;
  double t336, t338, t340, t345, t346, t347, t353, t359;
  double t362, t365, t370, t371, t376, t385, t388, t389;
  double t390, t391, t396, t397, t402, t407, t411, t412;
  double t416, t417, t418, t419, t424, t426, t428, t429;
  double t431, t435, t436, t439, t442, t445, t448, t450;
  double t452, t454, t459, t460, t464, t467, t470, t473;
  double t475, t476, t478, t481, t484, t491, t496, t497;
  double t505, t506, t519, t525, t535, t542, t546, t567;
  double t568, t575, t577, t578, t579, t582, t585, t588;
  double t589, t591, t594, t597, t600, t603, t605, t607;
  double t617, t671, t677, t681, t683, t685, t688, t691;
  double t698, t700, t704, t714, t718, t722, t751, t752;
  double t754, t757, t763, t765, t804, t828, t838, t877;
  double t882, t885, t887, t890, t898, t899, t902, t905;
  double t932, t937, t940, t943, t949, t960, t974, t979;
  double t982, t996, t1007;


  t1 = M_CBRT3;
  t2 = 0.1e1 / M_PI;
  t3 = POW_1_3(t2);
  t4 = t1 * t3;
  t5 = M_CBRT4;
  t6 = t5 * t5;
  t7 = rho[0] + rho[1];
  t8 = POW_1_3(t7);
  t9 = 0.1e1 / t8;
  t10 = t6 * t9;
  t11 = t4 * t10;
  t12 = t11 / 0.4e1;
  t13 = sqrt(t11);
  t15 = t12 + 0.18637200000000000000e1 * t13 + 0.129352e2;
  t16 = 0.1e1 / t15;
  t20 = log(t4 * t10 * t16 / 0.4e1);
  t22 = t13 + 0.372744e1;
  t25 = atan(0.61519908197590802322e1 / t22);
  t27 = t13 / 0.2e1;
  t28 = t27 + 0.10498e0;
  t29 = t28 * t28;
  t31 = log(t29 * t16);
  t33 = 0.310907e-1 * t20 + 0.38783294878113014393e-1 * t25 + 0.96902277115443742139e-3 * t31;
  t34 = rho[0] - rho[1];
  t35 = 0.1e1 / t7;
  t36 = t34 * t35;
  t37 = 0.1e1 + t36;
  t38 = POW_1_3(t37);
  t40 = 0.1e1 - t36;
  t41 = POW_1_3(t40);
  t43 = t38 * t37 + t41 * t40 - 0.2e1;
  t44 = M_CBRT2;
  t47 = 0.1e1 / (0.2e1 * t44 - 0.2e1);
  t49 = -t43 * t47 + 0.1e1;
  t50 = t33 * t49;
  t52 = t12 + 0.35302100000000000000e1 * t13 + 0.180578e2;
  t53 = 0.1e1 / t52;
  t57 = log(t4 * t10 * t53 / 0.4e1);
  t59 = t13 + 0.706042e1;
  t62 = atan(0.47309269095601128300e1 / t59);
  t64 = t27 + 0.32500e0;
  t65 = t64 * t64;
  t67 = log(t65 * t53);
  t69 = 0.1554535e-1 * t57 + 0.52491393169780936218e-1 * t62 + 0.22478670955426118383e-2 * t67;
  t71 = t69 * t43 * t47;
  if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    *zk = t50 + t71;

  if(order < 1) return;


  t73 = 0.1e1 / t8 / t7;
  t74 = t6 * t73;
  t78 = t4 * t6;
  t79 = t15 * t15;
  t80 = 0.1e1 / t79;
  t81 = t9 * t80;
  t82 = t4 * t74;
  t83 = t82 / 0.12e2;
  t84 = 0.1e1 / t13;
  t85 = t84 * t1;
  t86 = t3 * t6;
  t88 = t85 * t86 * t73;
  t90 = -t83 - 0.31062000000000000000e0 * t88;
  t95 = t1 * t1;
  t97 = 0.1e1 / t3;
  t98 = (-t4 * t74 * t16 / 0.12e2 - t78 * t81 * t90 / 0.4e1) * t95 * t97;
  t99 = t5 * t8;
  t100 = t99 * t15;
  t103 = t22 * t22;
  t104 = 0.1e1 / t103;
  t106 = t104 * t84 * t1;
  t108 = 0.37846991046400000000e2 * t104 + 0.1e1;
  t109 = 0.1e1 / t108;
  t114 = t28 * t16;
  t115 = t114 * t84;
  t118 = t29 * t80;
  t120 = -t115 * t82 / 0.6e1 - t118 * t90;
  t121 = 0.1e1 / t29;
  t122 = t120 * t121;
  t125 = 0.10363566666666666667e-1 * t98 * t100 + 0.39765745675026770179e-1 * t106 * t86 * t73 * t109 + 0.96902277115443742139e-3 * t122 * t15;
  t126 = t125 * t49;
  t127 = t7 * t7;
  t128 = 0.1e1 / t127;
  t129 = t34 * t128;
  t130 = t35 - t129;
  t132 = -t130;
  t135 = 0.4e1 / 0.3e1 * t38 * t130 + 0.4e1 / 0.3e1 * t41 * t132;
  t137 = t33 * t135 * t47;
  t141 = t52 * t52;
  t142 = 0.1e1 / t141;
  t143 = t9 * t142;
  t145 = -t83 - 0.58836833333333333333e0 * t88;
  t151 = (-t4 * t74 * t53 / 0.12e2 - t78 * t143 * t145 / 0.4e1) * t95 * t97;
  t152 = t99 * t52;
  t155 = t59 * t59;
  t156 = 0.1e1 / t155;
  t158 = t156 * t84 * t1;
  t160 = 0.22381669423600000000e2 * t156 + 0.1e1;
  t161 = 0.1e1 / t160;
  t166 = t64 * t53;
  t167 = t166 * t84;
  t170 = t65 * t142;
  t172 = -t167 * t82 / 0.6e1 - t170 * t145;
  t173 = 0.1e1 / t65;
  t174 = t172 * t173;
  t177 = 0.51817833333333333333e-2 * t151 * t152 + 0.41388824077869423261e-1 * t158 * t86 * t73 * t161 + 0.22478670955426118383e-2 * t174 * t52;
  t179 = t177 * t43 * t47;
  t181 = t69 * t135 * t47;
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[0] = t50 + t71 + t7 * (t126 - t137 + t179 + t181);

  t184 = -t35 - t129;
  t186 = -t184;
  t189 = 0.4e1 / 0.3e1 * t38 * t184 + 0.4e1 / 0.3e1 * t41 * t186;
  t191 = t33 * t189 * t47;
  t193 = t69 * t189 * t47;
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[1] = t50 + t71 + t7 * (t126 - t191 + t179 + t193);

  if(order < 2) return;


  t196 = 0.2e1 * t126;
  t198 = 0.2e1 * t179;
  t201 = 0.1e1 / t8 / t127;
  t202 = t6 * t201;
  t204 = t4 * t202 * t16;
  t206 = t73 * t80;
  t211 = 0.1e1 / t79 / t15;
  t212 = t9 * t211;
  t213 = t90 * t90;
  t217 = t4 * t202;
  t218 = t217 / 0.9e1;
  t220 = 0.1e1 / t13 / t11;
  t221 = t220 * t95;
  t222 = t3 * t3;
  t223 = t222 * t5;
  t224 = t8 * t8;
  t226 = 0.1e1 / t224 / t127;
  t228 = t221 * t223 * t226;
  t231 = t85 * t86 * t201;
  t233 = t218 - 0.20708000000000000000e0 * t228 + 0.41416000000000000000e0 * t231;
  t239 = (t204 / 0.9e1 + t78 * t206 * t90 / 0.6e1 + t78 * t212 * t213 / 0.2e1 - t78 * t81 * t233 / 0.4e1) * t95 * t97;
  t243 = t5 / t224;
  t244 = t243 * t15;
  t247 = t99 * t90;
  t252 = 0.1e1 / t103 / t22 * t1;
  t253 = t252 * t3;
  t258 = t104 * t220 * t95;
  t267 = t103 * t103;
  t270 = 0.1e1 / t267 / t22 * t1;
  t271 = t270 * t3;
  t272 = t108 * t108;
  t273 = 0.1e1 / t272;
  t278 = t28 * t80;
  t279 = t278 * t85;
  t284 = t114 * t220;
  t285 = t95 * t222;
  t287 = t285 * t5 * t226;
  t292 = t29 * t211;
  t296 = t204 / 0.72e2 + t279 * t86 * t73 * t90 / 0.3e1 - t284 * t287 / 0.9e1 + 0.2e1 / 0.9e1 * t115 * t217 + 0.2e1 * t292 * t213 - t118 * t233;
  t297 = t296 * t121;
  t301 = 0.1e1 / t29 / t28;
  t302 = t120 * t301;
  t303 = t15 * t84;
  t304 = t302 * t303;
  t309 = 0.10363566666666666667e-1 * t239 * t100 + 0.34545222222222222223e-2 * t98 * t244 + 0.10363566666666666667e-1 * t98 * t247 + 0.13255248558342256726e-1 * t253 * t202 * t109 + 0.26510497116684513453e-1 * t258 * t223 * t226 * t109 - 0.53020994233369026905e-1 * t106 * t86 * t201 * t109 - 0.50167127350538589834e0 * t271 * t202 * t273 + 0.96902277115443742139e-3 * t297 * t15 + 0.16150379519240623690e-3 * t304 * t82 + 0.96902277115443742139e-3 * t122 * t90;
  t310 = t309 * t49;
  t312 = t125 * t135 * t47;
  t313 = 0.2e1 * t312;
  t314 = t38 * t38;
  t315 = 0.1e1 / t314;
  t316 = t130 * t130;
  t319 = t127 * t7;
  t320 = 0.1e1 / t319;
  t321 = t34 * t320;
  t323 = -0.2e1 * t128 + 0.2e1 * t321;
  t326 = t41 * t41;
  t327 = 0.1e1 / t326;
  t328 = t132 * t132;
  t331 = -t323;
  t334 = 0.4e1 / 0.9e1 * t315 * t316 + 0.4e1 / 0.3e1 * t38 * t323 + 0.4e1 / 0.9e1 * t327 * t328 + 0.4e1 / 0.3e1 * t41 * t331;
  t336 = t33 * t334 * t47;
  t338 = t4 * t202 * t53;
  t340 = t73 * t142;
  t345 = 0.1e1 / t141 / t52;
  t346 = t9 * t345;
  t347 = t145 * t145;
  t353 = t218 - 0.39224555555555555555e0 * t228 + 0.78449111111111111110e0 * t231;
  t359 = (t338 / 0.9e1 + t78 * t340 * t145 / 0.6e1 + t78 * t346 * t347 / 0.2e1 - t78 * t143 * t353 / 0.4e1) * t95 * t97;
  t362 = t243 * t52;
  t365 = t99 * t145;
  t370 = 0.1e1 / t155 / t59 * t1;
  t371 = t370 * t3;
  t376 = t156 * t220 * t95;
  t385 = t155 * t155;
  t388 = 0.1e1 / t385 / t59 * t1;
  t389 = t388 * t3;
  t390 = t160 * t160;
  t391 = 0.1e1 / t390;
  t396 = t64 * t142;
  t397 = t396 * t85;
  t402 = t166 * t220;
  t407 = t65 * t345;
  t411 = t338 / 0.72e2 + t397 * t86 * t73 * t145 / 0.3e1 - t402 * t287 / 0.9e1 + 0.2e1 / 0.9e1 * t167 * t217 + 0.2e1 * t407 * t347 - t170 * t353;
  t412 = t411 * t173;
  t416 = 0.1e1 / t65 / t64;
  t417 = t172 * t416;
  t418 = t52 * t84;
  t419 = t417 * t418;
  t424 = 0.51817833333333333333e-2 * t359 * t152 + 0.17272611111111111111e-2 * t151 * t362 + 0.51817833333333333333e-2 * t151 * t365 + 0.13796274692623141087e-1 * t371 * t202 * t161 + 0.27592549385246282174e-1 * t376 * t223 * t226 * t161 - 0.55185098770492564348e-1 * t158 * t86 * t201 * t161 - 0.30878365944746984535e0 * t389 * t202 * t391 + 0.22478670955426118383e-2 * t412 * t52 + 0.37464451592376863972e-3 * t419 * t82 + 0.22478670955426118383e-2 * t174 * t145;
  t426 = t424 * t43 * t47;
  t428 = t177 * t135 * t47;
  t429 = 0.2e1 * t428;
  t431 = t69 * t334 * t47;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[0] = t196 - 0.2e1 * t137 + t198 + 0.2e1 * t181 + t7 * (t310 - t313 - t336 + t426 + t429 + t431);

  t435 = t125 * t189 * t47;
  t436 = t315 * t184;
  t439 = t38 * t34;
  t442 = t327 * t186;
  t445 = t41 * t34;
  t448 = 0.4e1 / 0.9e1 * t436 * t130 + 0.8e1 / 0.3e1 * t439 * t320 + 0.4e1 / 0.9e1 * t442 * t132 - 0.8e1 / 0.3e1 * t445 * t320;
  t450 = t33 * t448 * t47;
  t452 = t177 * t189 * t47;
  t454 = t69 * t448 * t47;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[1] = t196 - t137 + t198 + t181 - t191 + t193 + t7 * (t310 - t312 - t435 - t450 + t426 + t428 + t452 + t454);

  t459 = 0.2e1 * t435;
  t460 = t184 * t184;
  t464 = 0.2e1 * t128 + 0.2e1 * t321;
  t467 = t186 * t186;
  t470 = -t464;
  t473 = 0.4e1 / 0.9e1 * t315 * t460 + 0.4e1 / 0.3e1 * t38 * t464 + 0.4e1 / 0.9e1 * t327 * t467 + 0.4e1 / 0.3e1 * t41 * t470;
  t475 = t33 * t473 * t47;
  t476 = 0.2e1 * t452;
  t478 = t69 * t473 * t47;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[2] = t196 - 0.2e1 * t191 + t198 + 0.2e1 * t193 + t7 * (t310 - t459 - t475 + t426 + t476 + t478);

  if(order < 3) return;


  t481 = 0.3e1 * t310;
  t484 = 0.3e1 * t426;
  t491 = t5 / t224 / t7;
  t496 = 0.1e1 / t8 / t319;
  t497 = t6 * t496;
  t505 = 0.1e1 / t224 / t319;
  t506 = t5 * t505;
  t519 = 0.1e1 / t222;
  t525 = t267 * t267;
  t535 = t29 * t29;
  t542 = t4 * t497 * t16;
  t546 = t78 * t201 * t80 * t90;
  t567 = t1 * t519 * t6;
  t568 = t496 * t2;
  t575 = 0.1e1 / t13 / t285 / t243 / 0.4e1;
  t577 = t127 * t127;
  t578 = 0.1e1 / t577;
  t579 = t575 * t2 * t578;
  t582 = t285 * t506;
  t585 = t4 * t497;
  t588 = t79 * t79;
  t589 = 0.1e1 / t588;
  t591 = t213 * t90;
  t594 = t90 * t233;
  t597 = 0.7e1 / 0.27e2 * t585;
  t600 = t221 * t223 * t505;
  t603 = t85 * t86 * t496;
  t605 = -t597 - 0.12424800000000000000e1 * t579 + 0.82832000000000000000e0 * t600 - 0.96637333333333333333e0 * t603;
  t607 = -0.11e2 / 0.216e3 * t542 - t546 / 0.24e2 - t28 * t211 * t85 * t86 * t73 * t213 + t278 * t221 * t223 * t226 * t90 / 0.3e1 - 0.2e1 / 0.3e1 * t279 * t86 * t201 * t90 + t279 * t86 * t73 * t233 / 0.2e1 + t567 * t568 * t16 / 0.432e3 - 0.2e1 / 0.3e1 * t114 * t579 + 0.4e1 / 0.9e1 * t284 * t582 - 0.14e2 / 0.27e2 * t115 * t585 - 0.6e1 * t29 * t589 * t591 + 0.6e1 * t292 * t594 - t118 * t605;
  t617 = 0.69090444444444444446e-2 * t239 * t244 - 0.23030148148148148149e-2 * t98 * t491 * t15 - 0.48602578047254941329e-1 * t253 * t497 * t109 + 0.26510497116684513452e-1 / t267 * t95 * t222 * t506 * t109 * t84 - 0.23411326096918008589e1 / t267 / t103 * t95 * t222 * t506 * t273 * t84 + 0.22092080930570427878e-2 * t252 * t519 * t497 * t2 * t109 + 0.50631328524251801694e2 / t525 * t95 * t222 * t506 / t272 / t108 * t84 + 0.40375948798101559225e-4 * t120 / t535 * t15 * t217 + 0.96902277115443742139e-3 * t607 * t121 * t15 + 0.19380455423088748428e-2 * t297 * t90 + 0.96902277115443742139e-3 * t122 * t233 - 0.21533839358987498253e-3 * t304 * t217;
  t671 = t497 * t273;
  t677 = t2 * t578;
  t681 = -0.10604198846673805381e0 * t258 * t223 * t505 * t109 + 0.12371565321119439611e0 * t106 * t86 * t496 * t109 + 0.32300759038481247380e-3 * t296 * t301 * t303 * t82 + 0.32300759038481247380e-3 * t302 * t90 * t84 * t82 + 0.10766919679493749127e-3 * t302 * t15 * t220 * t287 + 0.10363566666666666667e-1 * (-0.7e1 / 0.27e2 * t542 - t546 / 0.3e1 - t78 * t73 * t211 * t213 / 0.2e1 + t78 * t206 * t233 / 0.4e1 - 0.3e1 / 0.2e1 * t78 * t9 * t589 * t591 + 0.3e1 / 0.2e1 * t78 * t212 * t594 - t78 * t81 * t605 / 0.4e1) * t95 * t97 * t100 + 0.20727133333333333334e-1 * t239 * t247 + 0.69090444444444444446e-2 * t98 * t243 * t90 + 0.10363566666666666667e-1 * t98 * t99 * t233 - 0.26614487661862784322e-1 * t270 * t519 * t671 + 0.18394613361864149606e1 * t271 * t671 + 0.15906298270010708072e0 * t104 * t575 * t677 * t109;
  t683 = (t617 + t681) * t49;
  t685 = t309 * t135 * t47;
  t688 = t125 * t334 * t47;
  t691 = 0.1e1 / t314 / t37;
  t698 = t34 * t578;
  t700 = 0.6e1 * t320 - 0.6e1 * t698;
  t704 = 0.1e1 / t326 / t40;
  t714 = -0.8e1 / 0.27e2 * t691 * t316 * t130 + 0.4e1 / 0.3e1 * t315 * t130 * t323 + 0.4e1 / 0.3e1 * t38 * t700 - 0.8e1 / 0.27e2 * t704 * t328 * t132 + 0.4e1 / 0.3e1 * t327 * t132 * t331 - 0.4e1 / 0.3e1 * t41 * t700;
  t718 = t4 * t497 * t53;
  t722 = t78 * t201 * t142 * t145;
  t751 = t141 * t141;
  t752 = 0.1e1 / t751;
  t754 = t347 * t145;
  t757 = t145 * t353;
  t763 = -t597 - 0.23534733333333333333e1 * t579 + 0.15689822222222222222e1 * t600 - 0.18304792592592592592e1 * t603;
  t765 = -0.11e2 / 0.216e3 * t718 - t722 / 0.24e2 - t64 * t345 * t85 * t86 * t73 * t347 + t396 * t221 * t223 * t226 * t145 / 0.3e1 - 0.2e1 / 0.3e1 * t397 * t86 * t201 * t145 + t397 * t86 * t73 * t353 / 0.2e1 + t567 * t568 * t53 / 0.432e3 - 0.2e1 / 0.3e1 * t166 * t579 + 0.4e1 / 0.9e1 * t402 * t582 - 0.14e2 / 0.27e2 * t167 * t585 - 0.6e1 * t65 * t752 * t754 + 0.6e1 * t407 * t757 - t170 * t763;
  t804 = 0.22478670955426118383e-2 * t765 * t173 * t52 + 0.44957341910852236766e-2 * t412 * t145 + 0.22478670955426118383e-2 * t174 * t353 - 0.49952602123169151963e-3 * t419 * t217 - 0.11037019754098512870e0 * t376 * t223 * t505 * t161 + 0.12876523046448265015e0 * t158 * t86 * t496 * t161 + 0.74928903184753727944e-3 * t411 * t416 * t418 * t82 + 0.74928903184753727944e-3 * t417 * t145 * t84 * t82 + 0.24976301061584575981e-3 * t417 * t52 * t220 * t287 + 0.16555529631147769304e0 * t156 * t575 * t677 * t161 + 0.34545222222222222222e-2 * t359 * t362 - 0.11515074074074074074e-2 * t151 * t491 * t52;
  t828 = t385 * t385;
  t838 = t65 * t65;
  t877 = t497 * t391;
  t882 = -0.50586340539618183985e-1 * t371 * t497 * t161 + 0.27592549385246282174e-1 / t385 * t95 * t222 * t506 * t161 * t84 - 0.14409904107548592783e1 / t385 / t155 * t95 * t222 * t506 * t391 * t84 + 0.22993791154371901812e-2 * t370 * t519 * t497 * t2 * t161 + 0.18429583437767336288e2 / t828 * t95 * t222 * t506 / t390 / t160 * t84 + 0.93661128980942159930e-4 * t172 / t838 * t52 * t217 + 0.51817833333333333333e-2 * (-0.7e1 / 0.27e2 * t718 - t722 / 0.3e1 - t78 * t73 * t345 * t347 / 0.2e1 + t78 * t340 * t353 / 0.4e1 - 0.3e1 / 0.2e1 * t78 * t9 * t752 * t754 + 0.3e1 / 0.2e1 * t78 * t346 * t757 - t78 * t143 * t763 / 0.4e1) * t95 * t97 * t152 + 0.10363566666666666667e-1 * t359 * t365 + 0.34545222222222222222e-2 * t151 * t243 * t145 + 0.51817833333333333333e-2 * t151 * t99 * t353 - 0.16381481915689750931e-1 * t388 * t519 * t877 + 0.11322067513073894330e1 * t389 * t877;
  t885 = (t804 + t882) * t43 * t47;
  t887 = t424 * t135 * t47;
  t890 = t177 * t334 * t47;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[0] = t481 - 0.6e1 * t312 - 0.3e1 * t336 + t484 + 0.6e1 * t428 + 0.3e1 * t431 + t7 * (-t33 * t714 * t47 + t69 * t714 * t47 + t683 - 0.3e1 * t685 - 0.3e1 * t688 + t885 + 0.3e1 * t887 + 0.3e1 * t890);

  t898 = 0.2e1 * t450;
  t899 = 0.2e1 * t454;
  t902 = t309 * t189 * t47;
  t905 = 0.2e1 * t125 * t448 * t47;
  t932 = -0.8e1 / 0.27e2 * t691 * t184 * t316 + 0.16e2 / 0.9e1 * t315 * t34 * t320 * t130 + 0.4e1 / 0.9e1 * t436 * t323 + 0.8e1 / 0.3e1 * t38 * t320 - 0.8e1 * t439 * t578 - 0.8e1 / 0.27e2 * t704 * t186 * t328 - 0.16e2 / 0.9e1 * t327 * t34 * t320 * t132 + 0.4e1 / 0.9e1 * t442 * t331 - 0.8e1 / 0.3e1 * t41 * t320 + 0.8e1 * t445 * t578;
  t937 = t424 * t189 * t47;
  t940 = 0.2e1 * t177 * t448 * t47;
  t943 = -t33 * t932 * t47 + t69 * t932 * t47 + t683 - 0.2e1 * t685 - t688 + t885 + 0.2e1 * t887 + t890 - t902 - t905 + t937 + t940;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[1] = t7 * t943 - 0.4e1 * t312 - t336 + 0.4e1 * t428 + t431 - t459 + t476 + t481 + t484 - t898 + t899;

  t949 = t125 * t473 * t47;
  t960 = -0.2e1 * t320 - 0.6e1 * t698;
  t974 = -0.8e1 / 0.27e2 * t691 * t460 * t130 + 0.16e2 / 0.9e1 * t436 * t321 + 0.4e1 / 0.9e1 * t315 * t464 * t130 + 0.4e1 / 0.3e1 * t38 * t960 - 0.8e1 / 0.27e2 * t704 * t467 * t132 - 0.16e2 / 0.9e1 * t442 * t321 + 0.4e1 / 0.9e1 * t327 * t470 * t132 - 0.4e1 / 0.3e1 * t41 * t960;
  t979 = t177 * t473 * t47;
  t982 = -t33 * t974 * t47 + t69 * t974 * t47 + t683 - t685 + t885 + t887 - 0.2e1 * t902 - t905 + 0.2e1 * t937 + t940 - t949 + t979;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[2] = t7 * t982 - t313 + t429 - 0.4e1 * t435 + 0.4e1 * t452 - t475 + t478 + t481 + t484 - t898 + t899;

  t996 = -0.6e1 * t320 - 0.6e1 * t698;
  t1007 = -0.8e1 / 0.27e2 * t691 * t460 * t184 + 0.4e1 / 0.3e1 * t436 * t464 + 0.4e1 / 0.3e1 * t38 * t996 - 0.8e1 / 0.27e2 * t704 * t467 * t186 + 0.4e1 / 0.3e1 * t442 * t470 - 0.4e1 / 0.3e1 * t41 * t996;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[3] = t481 - 0.6e1 * t435 - 0.3e1 * t475 + t484 + 0.6e1 * t452 + 0.3e1 * t478 + t7 * (-t33 * t1007 * t47 + t69 * t1007 * t47 + t683 + t885 - 0.3e1 * t902 + 0.3e1 * t937 - 0.3e1 * t949 + 0.3e1 * t979);

  if(order < 4) return;



}

