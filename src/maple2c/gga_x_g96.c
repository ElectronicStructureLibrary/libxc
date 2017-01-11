
void XC(gga_x_g96_enhance)
  (const XC(func_type) *p, int order, 
   FLOAT x, FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  double t1, t2, t5, t7, t9;


  if(order < 0) return;

  t1 = 0.1e1 / X_FACTOR_C;
  t2 = pow(x, 0.15000000000000000000e1);

  *f = 0.10e1 + 0.72992700729927007299e-2 * t1 * t2;

  if(order < 0+1) return;

  t5 = pow(x, 0.5000000000000000000e0);

  *dfdx = 0.10948905109489051095e-1 * t1 * t5;

  if(order < 1+1) return;

  t7 = pow(x, -0.5000000000000000000e0);

  *d2fdx2 = 0.54744525547445255475e-2 * t1 * t7;

  if(order < 2+1) return;

  t9 = pow(x, -0.15000000000000000000e1);

  *d3fdx3 = -0.27372262773722627738e-2 * t1 * t9;

  if(order < 3+1) return;

}


#define maple2c_order 3
#define maple2c_func  XC(gga_x_g96_enhance)

