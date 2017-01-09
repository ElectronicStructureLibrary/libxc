
void XC(gga_x_herman_enhance)
  (const XC(func_type) *p, int order, 
   FLOAT x, FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
  double t1, t2;


  if(order < 0) return;

  t1 = 0.1e1 / X_FACTOR_C;
  t2 = x * x;

  *f = 0.10e1 + 0.3e-2 * t1 * t2;

  if(order < 0+1) return;


  *dfdx = 0.6e-2 * t1 * x;

  if(order < 1+1) return;


  *d2fdx2 = 0.6e-2 * t1;

  if(order < 2+1) return;


  *d3fdx3 = 0;

  if(order < 3+1) return;

}


#define maple2c_order 3
#define maple2c_func  XC(gga_x_herman_enhance)

