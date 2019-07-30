/*
 Copyright (C) 2019 Daniel Mejia-Rodriguez

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

static inline void
ked_unpol(int order, const double *rho, const double *sigma, const double *lapl, double *taus, double *vrho, double *vsigma, double *vlapl)
{
  double kf, kf2, p, q;
  double ge2, dge2dp, dge2dq;
  double ge4, dge4dp, dge4dq;
  double ftw, oneftw, factor, dfactordp, dfactordq, dftwdp;
  double mge4, dmge4dp, dmge4dq;
  double z, dzdp, dzdq;
  double exp1, exp2, fz, fzb, dfzb, exp1b;
  double fs, dfsdp, dfsdq;
  double factor2;
  //
  const double a = 1.784720e0;
  const double b = 0.258304e0;
  const double thresh1 = 0.00066e0;
  const double thresh2 = 1.73289e0;
  //
  const double f53 = 5.0/3.0;
  const double ckf = CBRT(3.0*M_PI*M_PI);
  const double ckf2 = POW_2(ckf);
  const double ctf = 0.3*ckf2;
  //
  kf = ckf * CBRT(*rho);
  kf2 = kf*kf;
  //
  p = *sigma / (4.0*kf2*POW_2( *rho ) );
  q = *lapl / (4.0*kf2* *rho);
  //
  ge2 = 5.0/27.0*p + 20.0/9.0*q;
  ge4 = 8.0/81.0*q*q - p*q/9.0 + 8.0/243.0*p*p;
  ftw = 5.0/3.0*p;
  oneftw = 1.0 + ftw;
  factor = ge4/oneftw;
  factor2 = factor*factor;
  mge4 = (1.0 + ge2 + ge4)/sqrt(1.0 + factor2);
  z = mge4 - ftw;
  //
  if ( z > thresh1 && z < thresh2 )
  {
    exp1 = exp(-a/z);
    exp1b = exp(-a*b/z);
    exp2 = exp(-a/(a-z));
    fzb = exp1b*pow((1.0+exp2)/(exp1+exp2),b);
    dfzb = a*b*exp2*fzb/(exp1+exp2)*(1.0/pow(z,2) + (1.0-exp1)/((1.0+exp2)*pow(a-z,2)));
  } 
  else if ( z >= thresh2 )
  {
    fzb = 1.0;
    dfzb = 0.0;
  } 
  else 
  {
    fzb = 0.0;
    dfzb = 0.0;
  };
  //
  fs = z*fzb + ftw;
  *taus = ctf*fs*POW_5_3(*rho);
  //
  if(order < 1) return;
  //
  dge2dp = 5.0/27.0;
  dge2dq = 20.0/9.0;
  dge4dp = -q/9.0 + 16.0/243.0*p;
  dge4dq = 16.0/81.0*q - p/9.0;
  dftwdp = 5.0/3.0;
  dfactordp = (dge4dp - dftwdp*factor)/oneftw;
  dfactordq = dge4dq/oneftw;
  //
  dmge4dp = (dge2dp + dge4dp)/sqrt(1.0 + factor2) - 
    mge4*factor*dfactordp/(1.0 + factor2);
  dmge4dq = (dge2dq + dge4dq)/sqrt(1.0 + factor2) -
    mge4*factor*dfactordq/(1.0 + factor2);
  //
  dzdp = dmge4dp - dftwdp;
  dzdq = dmge4dq;
  //
  dfsdp = dzdp*fzb + z*dfzb*dzdp + dftwdp;
  dfsdq = dzdq*fzb + z*dfzb*dzdq;
  *vrho = 5.0/3.0*ctf*POW_2_3( *rho ) * (fs - 1.6*p*dfsdp - q*dfsdq);
  *vsigma = 0.075*dfsdp / *rho;
  *vlapl =  0.075*dfsdq;
};
