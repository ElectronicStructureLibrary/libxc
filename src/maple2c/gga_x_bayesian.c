
void XC(gga_x_bayesian_enhance)
  (const XC(func_type) *p, int order, 
   FLOAT x, FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  double t1, t2, t3, t5, t6, t7, t10, t11;
  double t13, t16, t17, t19, t20, t27, t28, t30;
  double t33, t38, t39, t40, t41, t42, t45, t53;
  double t54, t56, t61, t69, t71;


  if(order < 0) return;

  t1 = X2S * X2S;
  t2 = x * x;
  t3 = t1 * t2;
  t5 = 0.10e1 + X2S * x;
  t6 = t5 * t5;
  t7 = 0.1e1 / t6;
  t10 = 0.1926e0 + 0.18962e1 * t3 * t7;
  t11 = t7 * t10;

  *f = 0.10008e1 + t3 * t11;

  if(order < 0+1) return;

  t13 = t1 * x;
  t16 = t1 * X2S;
  t17 = t16 * t2;
  t19 = 0.1e1 / t6 / t5;
  t20 = t19 * t10;
  t27 = 0.37924e1 * t13 * t7 - 0.37924e1 * t17 * t19;
  t28 = t7 * t27;

  *dfdx = 0.2e1 * t13 * t11 - 0.2e1 * t17 * t20 + t3 * t28;

  if(order < 1+1) return;

  t30 = t1 * t7;
  t33 = t16 * x;
  t38 = t1 * t1;
  t39 = t38 * t2;
  t40 = t6 * t6;
  t41 = 0.1e1 / t40;
  t42 = t41 * t10;
  t45 = t19 * t27;
  t53 = 0.37924e1 * t30 - 0.151696e2 * t33 * t19 + 0.113772e2 * t39 * t41;
  t54 = t7 * t53;

  *d2fdx2 = 0.2e1 * t30 * t10 + 0.4e1 * t13 * t28 - 0.4e1 * t17 * t45 - 0.8e1 * t33 * t20 + t3 * t54 + 0.6e1 * t39 * t42;

  if(order < 2+1) return;

  t56 = t16 * t19;
  t61 = t38 * x;
  t69 = t38 * X2S * t2;
  t71 = 0.1e1 / t40 / t5;

  *d3fdx3 = -0.12e2 * t56 * t10 + 0.6e1 * t30 * t27 + 0.36e2 * t61 * t42 - 0.24e2 * t33 * t45 + 0.6e1 * t13 * t54 - 0.24e2 * t69 * t71 * t10 + 0.18e2 * t39 * t41 * t27 - 0.6e1 * t17 * t19 * t53 + t3 * t7 * (-0.227544e2 * t56 + 0.682632e2 * t61 * t41 - 0.455088e2 * t69 * t71);

  if(order < 3+1) return;

}


#define maple2c_order 3
#define maple2c_func  XC(gga_x_bayesian_enhance)

