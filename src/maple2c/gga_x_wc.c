
void XC(gga_x_wc_enhance)
  (const XC(func_type) *p, int order, 
   FLOAT x, FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{

  double t1, t2, t3, t5, t6, t8, t9, t11;
  double t12, t13, t15, t16, t17, t22, t23, t25;
  double t28, t31, t32, t33, t36, t40, t43, t44;
  double t53, t60, t61, t62, t64, t65, t69, t71;
  double t85, t98;

  if(order < 0) return;

  t1 = X2S * X2S;
  t2 = MU_GE * t1;
  t3 = x * x;
  t5 = mu - MU_GE;
  t6 = t5 * t1;
  t8 = exp(-t1 * t3);
  t9 = t3 * t8;
  t11 = t1 * t1;
  t12 = c * t11;
  t13 = t3 * t3;
  t15 = 0.10e1 + t12 * t13;
  t16 = log(t15);
  t17 = t2 * t3 + t6 * t9 + kappa + t16;

  *f = 0.10e1 + kappa * (0.10e1 - kappa / t17);

  if(order < 0+1) return;

  t22 = kappa * kappa;
  t23 = t17 * t17;
  t25 = t22 / t23;
  t28 = x * t8;
  t31 = t5 * t11;
  t32 = t3 * x;
  t33 = t32 * t8;
  t36 = 0.1e1 / t15;
  t40 = 0.4e1 * t12 * t32 * t36 + 0.2e1 * t2 * x + 0.2e1 * t6 * t28 - 0.2e1 * t31 * t33;

  *dfdx = t25 * t40;

  if(order < 1+1) return;

  t43 = t22 / t23 / t17;
  t44 = t40 * t40;
  t53 = t5 * t11 * t1;
  t60 = c * c;
  t61 = t11 * t11;
  t62 = t60 * t61;
  t64 = t15 * t15;
  t65 = 0.1e1 / t64;
  t69 = -0.16e2 * t62 * t13 * t3 * t65 + 0.12e2 * t12 * t3 * t36 + 0.4e1 * t53 * t13 * t8 - 0.10e2 * t31 * t9 + 0.2e1 * t6 * t8 + 0.2e1 * t2;

  *d2fdx2 = t25 * t69 - 0.2e1 * t43 * t44;

  if(order < 2+1) return;

  t71 = t23 * t23;
  t85 = t13 * x;
  t98 = t13 * t13;

  *d3fdx3 = 0.6e1 * t22 / t71 * t44 * t40 - 0.6e1 * t43 * t40 * t69 + t25 * (-0.24e2 * t31 * t28 + 0.36e2 * t53 * t33 - 0.8e1 * t5 * t61 * t85 * t8 + 0.24e2 * t12 * x * t36 - 0.144e3 * t62 * t85 * t65 + 0.128e3 * t60 * c * t61 * t11 * t98 * x / t64 / t15);

  if(order < 3+1) return;

}


#define maple2c_order 3
#define maple2c_func  XC(gga_x_wc_enhance)

