/* 
  This file was generated automatically with ../scripts/maple2c.pl.
  Do not edit this file directly as it can be overwritten!!

  Maple source      : ../maple/lda_c_2d_prm.mpl
  Type of functional: work_lda
*/

static void
func0(const XC(func_type) *p, XC(lda_work_t) *r)
{
  double t1, t3, t4, t5, t7, t8, t9, t10;
  double t11, t14, t15, t18, t19, t21, t24, t25;
  double t26, t29, t32, t33, t34, t37, t39, t40;
  double t41, t43, t44, t51, t58, t70, t73, t78;
  double t79, t80, t86, t87, t93, t94, t103, t106;
  double t112, t113, t125, t132, t135, t142, t145, t151;
  double t152, t155, t156, t158, t159, t172, t200;

  lda_c_2d_prm_params *params;

  assert(p->params != NULL);
  params = (lda_c_2d_prm_params * )(p->params);

  assert(params->N > 1.0);

  t1 = 0.1e1 / r->rs;
  t3 = 0.22157981704254580414e1 * t1 + 0.88622692545275801365e0;
  t4 = 0.1e1 / t3;
  t5 = t1 * t4;
  t7 = 0.22157981704254580414e1 * t5 - 0.10e1;
  t8 = t1 * t7;
  t9 = 0.20e1 + params->c;
  t10 = sqrt(t9);
  t11 = 0.1e1 / t10;
  t14 = 0.1e1 / t9;
  t15 = t7 * t14;
  t18 = t3 * t3;
  t19 = 0.1e1 / t18;
  t21 = pow(t9, -0.15e1);
  t24 = 0.10e1 + params->c;
  t25 = sqrt(t24);
  t26 = 0.1e1 / t25;
  t29 = 0.1e1 / t24;
  r->f = 0.19997916265148655845e0 * t8 * t11 + 0.22565232098914243869e0 * t5 * t15 + 0.99989581325743279224e-1 * t1 * t19 * t21 + 0.39995832530297311691e0 * t8 * t26 + 0.22565232098914243869e0 * t5 * t29;

  if(r->order < 1) return;

  t32 = r->rs * r->rs;
  t33 = 0.1e1 / t32;
  t34 = t33 * t7;
  t37 = t33 * t4;
  t39 = t32 * r->rs;
  t40 = 0.1e1 / t39;
  t41 = t40 * t19;
  t43 = -0.22157981704254580414e1 * t37 + 0.49097615320608071993e1 * t41;
  t44 = t1 * t43;
  t51 = t43 * t14;
  t58 = 0.1e1 / t18 / t3;
  r->dfdrs = -0.19997916265148655845e0 * t34 * t11 + 0.19997916265148655845e0 * t44 * t11 - 0.22565232098914243869e0 * t37 * t15 + 0.50000000000000000004e0 * t41 * t15 + 0.22565232098914243869e0 * t5 * t51 - 0.99989581325743279224e-1 * t33 * t19 * t21 + 0.44311346272637900685e0 * t40 * t58 * t21 - 0.39995832530297311691e0 * t34 * t26 + 0.39995832530297311691e0 * t44 * t26 - 0.22565232098914243869e0 * t37 * t29 + 0.50000000000000000004e0 * t41 * t29;

  if(r->order < 2) return;

  t70 = t40 * t4;
  t73 = t40 * t7;
  t78 = t32 * t32;
  t79 = 0.1e1 / t78;
  t80 = t79 * t19;
  t86 = 0.1e1 / t78 / r->rs;
  t87 = t86 * t58;
  t93 = 0.44315963408509160828e1 * t70 - 0.19639046128243228797e2 * t80 + 0.21758081239931260892e2 * t87;
  t94 = t93 * t14;
  t103 = t33 * t43;
  t106 = t1 * t93;
  t112 = t18 * t18;
  t113 = 0.1e1 / t112;
  r->d2fdrs2 = 0.45130464197828487738e0 * t70 * t15 + 0.39995832530297311690e0 * t73 * t11 - 0.45130464197828487738e0 * t37 * t51 - 0.20000000000000000001e1 * t80 * t15 + 0.10000000000000000001e1 * t41 * t51 + 0.22157981704254580416e1 * t87 * t15 + 0.22565232098914243869e0 * t5 * t94 + 0.19997916265148655845e0 * t41 * t21 + 0.79991665060594623382e0 * t73 * t26 + 0.45130464197828487738e0 * t70 * t29 - 0.39995832530297311690e0 * t103 * t11 + 0.19997916265148655845e0 * t106 * t11 - 0.17724538509055160274e1 * t79 * t58 * t21 + 0.29455500000000000002e1 * t86 * t113 * t21 - 0.79991665060594623382e0 * t103 * t26 + 0.39995832530297311691e0 * t106 * t26 - 0.20000000000000000001e1 * t80 * t29 + 0.22157981704254580416e1 * t87 * t29;

  if(r->order < 3) return;

  t125 = t79 * t7;
  t132 = t79 * t4;
  t135 = t40 * t43;
  t142 = t86 * t19;
  t145 = t33 * t93;
  t151 = 0.1e1 / t78 / t32;
  t152 = t151 * t58;
  t155 = 0.1e1 / t78 / t39;
  t156 = t155 * t113;
  t158 = -0.13294789022552748248e2 * t132 + 0.88375707577094529586e2 * t142 - 0.19582273115938134803e3 * t152 + 0.14463454981022450832e3 * t156;
  t159 = t1 * t158;
  t172 = -0.11998749759089193507e1 * t125 * t11 - 0.59993748795445967535e0 * t80 * t21 - 0.23997499518178387015e1 * t125 * t26 - 0.13539139259348546321e1 * t132 * t29 + 0.11998749759089193507e1 * t135 * t11 + 0.79760423290748221233e1 * t87 * t21 + 0.23997499518178387014e1 * t135 * t26 + 0.90000000000000000005e1 * t142 * t29 - 0.59993748795445967535e0 * t145 * t11 + 0.19997916265148655845e0 * t159 * t11 - 0.26509950000000000002e2 * t151 * t113 * t21 + 0.26106977203586831737e2 * t155 / t112 / t3 * t21 - 0.11998749759089193507e1 * t145 * t26;
  t200 = 0.39995832530297311691e0 * t159 * t26 - 0.19942183533829122374e2 * t152 * t29 + 0.14729284596182421599e2 * t156 * t29 - 0.13539139259348546321e1 * t132 * t15 + 0.13539139259348546321e1 * t70 * t51 + 0.90000000000000000005e1 * t142 * t15 - 0.67695696296742731607e0 * t37 * t94 - 0.60000000000000000005e1 * t80 * t51 - 0.19942183533829122374e2 * t152 * t15 + 0.15000000000000000001e1 * t41 * t94 + 0.66473945112763741248e1 * t87 * t51 + 0.14729284596182421599e2 * t156 * t15 + 0.22565232098914243869e0 * t5 * t158 * t14;
  r->d3fdrs3 = t172 + t200;

  if(r->order < 4) return;


}

static void
func1(const XC(func_type) *p, XC(lda_work_t) *r)
{
  double t1, t3, t4, t5, t7, t8, t9, t10;
  double t11, t14, t15, t18, t19, t21, t24, t25;
  double t26, t29, t32, t33, t34, t37, t39, t40;
  double t41, t43, t44, t51, t58, t70, t73, t78;
  double t79, t80, t86, t87, t93, t94, t103, t106;
  double t112, t113, t125, t130, t138, t139, t147, t148;
  double t155, t159, t166, t169, t178, t181, t200;

  lda_c_2d_prm_params *params;

  assert(p->params != NULL);
  params = (lda_c_2d_prm_params * )(p->params);

  assert(params->N > 1.0);

  t1 = 0.1e1 / r->rs;
  t3 = 0.22157981704254580414e1 * t1 + 0.88622692545275801365e0;
  t4 = 0.1e1 / t3;
  t5 = t1 * t4;
  t7 = 0.22157981704254580414e1 * t5 - 0.10e1;
  t8 = t1 * t7;
  t9 = 0.20e1 + params->c;
  t10 = sqrt(t9);
  t11 = 0.1e1 / t10;
  t14 = 0.1e1 / t9;
  t15 = t7 * t14;
  t18 = t3 * t3;
  t19 = 0.1e1 / t18;
  t21 = pow(t9, -0.15e1);
  t24 = 0.10e1 + params->c;
  t25 = sqrt(t24);
  t26 = 0.1e1 / t25;
  t29 = 0.1e1 / t24;
  r->f = 0.19997916265148655845e0 * t8 * t11 + 0.22565232098914243869e0 * t5 * t15 + 0.99989581325743279224e-1 * t1 * t19 * t21 + 0.39995832530297311691e0 * t8 * t26 + 0.22565232098914243869e0 * t5 * t29;

  if(r->order < 1) return;

  t32 = r->rs * r->rs;
  t33 = 0.1e1 / t32;
  t34 = t33 * t7;
  t37 = t33 * t4;
  t39 = t32 * r->rs;
  t40 = 0.1e1 / t39;
  t41 = t40 * t19;
  t43 = -0.22157981704254580414e1 * t37 + 0.49097615320608071993e1 * t41;
  t44 = t1 * t43;
  t51 = t43 * t14;
  t58 = 0.1e1 / t18 / t3;
  r->dfdrs = -0.19997916265148655845e0 * t34 * t11 + 0.19997916265148655845e0 * t44 * t11 - 0.22565232098914243869e0 * t37 * t15 + 0.50000000000000000004e0 * t41 * t15 + 0.22565232098914243869e0 * t5 * t51 - 0.99989581325743279224e-1 * t33 * t19 * t21 + 0.44311346272637900685e0 * t40 * t58 * t21 - 0.39995832530297311691e0 * t34 * t26 + 0.39995832530297311691e0 * t44 * t26 - 0.22565232098914243869e0 * t37 * t29 + 0.50000000000000000004e0 * t41 * t29;
  r->dfdz = 0;

  if(r->order < 2) return;

  t70 = t40 * t4;
  t73 = t40 * t7;
  t78 = t32 * t32;
  t79 = 0.1e1 / t78;
  t80 = t79 * t19;
  t86 = 0.1e1 / t78 / r->rs;
  t87 = t86 * t58;
  t93 = 0.44315963408509160828e1 * t70 - 0.19639046128243228797e2 * t80 + 0.21758081239931260892e2 * t87;
  t94 = t93 * t14;
  t103 = t33 * t43;
  t106 = t1 * t93;
  t112 = t18 * t18;
  t113 = 0.1e1 / t112;
  r->d2fdrs2 = 0.45130464197828487738e0 * t70 * t15 + 0.39995832530297311690e0 * t73 * t11 - 0.45130464197828487738e0 * t37 * t51 - 0.20000000000000000001e1 * t80 * t15 + 0.10000000000000000001e1 * t41 * t51 + 0.22157981704254580416e1 * t87 * t15 + 0.22565232098914243869e0 * t5 * t94 + 0.19997916265148655845e0 * t41 * t21 + 0.79991665060594623382e0 * t73 * t26 + 0.45130464197828487738e0 * t70 * t29 - 0.39995832530297311690e0 * t103 * t11 + 0.19997916265148655845e0 * t106 * t11 - 0.17724538509055160274e1 * t79 * t58 * t21 + 0.29455500000000000002e1 * t86 * t113 * t21 - 0.79991665060594623382e0 * t103 * t26 + 0.39995832530297311691e0 * t106 * t26 - 0.20000000000000000001e1 * t80 * t29 + 0.22157981704254580416e1 * t87 * t29;
  r->d2fdrsz = 0;
  r->d2fdz2 = 0;

  if(r->order < 3) return;

  t125 = t79 * t4;
  t130 = t86 * t19;
  t138 = 0.1e1 / t78 / t32;
  t139 = t138 * t58;
  t147 = 0.1e1 / t78 / t39;
  t148 = t147 * t113;
  t155 = -0.13294789022552748248e2 * t125 + 0.88375707577094529586e2 * t130 - 0.19582273115938134803e3 * t139 + 0.14463454981022450832e3 * t148;
  t159 = t79 * t7;
  t166 = -0.13539139259348546321e1 * t125 * t15 + 0.13539139259348546321e1 * t70 * t51 + 0.90000000000000000005e1 * t130 * t15 - 0.67695696296742731607e0 * t37 * t94 - 0.60000000000000000005e1 * t80 * t51 - 0.19942183533829122374e2 * t139 * t15 + 0.15000000000000000001e1 * t41 * t94 + 0.66473945112763741248e1 * t87 * t51 + 0.14729284596182421599e2 * t148 * t15 + 0.22565232098914243869e0 * t5 * t155 * t14 - 0.11998749759089193507e1 * t159 * t11 - 0.59993748795445967535e0 * t80 * t21 - 0.23997499518178387015e1 * t159 * t26;
  t169 = t40 * t43;
  t178 = t33 * t93;
  t181 = t1 * t155;
  t200 = -0.13539139259348546321e1 * t125 * t29 + 0.11998749759089193507e1 * t169 * t11 + 0.79760423290748221233e1 * t87 * t21 + 0.23997499518178387014e1 * t169 * t26 + 0.90000000000000000005e1 * t130 * t29 - 0.59993748795445967535e0 * t178 * t11 + 0.19997916265148655845e0 * t181 * t11 - 0.26509950000000000002e2 * t138 * t113 * t21 + 0.26106977203586831737e2 * t147 / t112 / t3 * t21 - 0.11998749759089193507e1 * t178 * t26 + 0.39995832530297311691e0 * t181 * t26 - 0.19942183533829122374e2 * t139 * t29 + 0.14729284596182421599e2 * t148 * t29;
  r->d3fdrs3 = t166 + t200;
  r->d3fdrs2z = 0;
  r->d3fdrsz2 = 0;
  r->d3fdz3 = 0;

  if(r->order < 4) return;


}

void 
XC(lda_c_2d_prm_func)(const XC(func_type) *p, XC(lda_work_t) *r)
{
  if(p->nspin == XC_UNPOLARIZED)
    func0(p, r);
  else
    func1(p, r);
}

#define maple2c_order 3
#define maple2c_func  XC(lda_c_2d_prm_func)