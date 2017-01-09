
void XC(gga_x_wc_enhance)
  (const XC(func_type) *p, int order, 
   FLOAT x, FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{

  double t1, t2, t3, t5, t6, t8, t9, t11;
  double t12, t15, t16, t17, t20, t21, t24, t27;
  double t28, t29, t33, t36, t39, t40, t49, t56;
  double t59, t60, t63, t66, t79, t90;

  if(order < 0) return;

  t1 = X2S * X2S;
  t2 = MU_GE * t1;
  t3 = x * x;
  t5 = 0.2195149727645171e0 - MU_GE;
  t6 = t5 * t1;
  t8 = exp(-t1 * t3);
  t9 = t3 * t8;
  t11 = t1 * t1;
  t12 = t3 * t3;
  t15 = 0.10e1 + 0.793746933516181879e-2 * t11 * t12;
  t16 = log(t15);
  t17 = 0.8040e0 + t2 * t3 + t6 * t9 + t16;

  *f = 0.180400e1 - 0.64641600e0 / t17;

  if(order < 0+1) return;

  t20 = t17 * t17;
  t21 = 0.1e1 / t20;
  t24 = x * t8;
  t27 = t5 * t11;
  t28 = t3 * x;
  t29 = t28 * t8;
  t33 = 0.1e1 / t15;
  t36 = 0.2e1 * t2 * x + 0.2e1 * t6 * t24 - 0.2e1 * t27 * t29 + 0.3174987734064727516e-1 * t11 * t28 * t33;

  *dfdx = 0.64641600e0 * t21 * t36;

  if(order < 1+1) return;

  t39 = 0.1e1 / t20 / t17;
  t40 = t36 * t36;
  t49 = t5 * t11 * t1;
  t56 = t11 * t11;
  t59 = t15 * t15;
  t60 = 0.1e1 / t59;
  t63 = 0.2e1 * t2 + 0.2e1 * t6 * t8 - 0.10e2 * t27 * t9 + 0.4e1 * t49 * t12 * t8 + 0.9524963202194182548e-1 * t11 * t3 * t33 - 0.10080547111461472895e-2 * t56 * t12 * t3 * t60;

  *d2fdx2 = -0.129283200e1 * t39 * t40 + 0.64641600e0 * t21 * t63;

  if(order < 2+1) return;

  t66 = t20 * t20;
  t79 = t12 * x;
  t90 = t12 * t12;

  *d3fdx3 = 0.387849600e1 / t66 * t40 * t36 - 0.387849600e1 * t39 * t36 * t63 + 0.64641600e0 * t21 * (-0.24e2 * t27 * t24 + 0.36e2 * t49 * t29 - 0.8e1 * t5 * t56 * t79 * t8 + 0.19049926404388365096e0 * t11 * x * t33 - 0.90724924003153256054e-2 * t56 * t79 * t60 + 0.64011226863103592059e-4 * t56 * t11 * t90 * x / t59 / t15);

  if(order < 3+1) return;

}


#define maple2c_order 3
#define maple2c_func  XC(gga_x_wc_enhance)

