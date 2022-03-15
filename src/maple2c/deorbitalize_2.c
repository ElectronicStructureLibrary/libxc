out->VAR(v2rho2, ip, 0) = ked1->VAR(v2rho2, ip, 0)*mgga->VAR(vtau, ip, 0) + ked1->VAR(vrho, ip,
  0)*ked1->VAR(vrho, ip, 0)*mgga->VAR(v2tau2, ip, 0) + 2*ked1->VAR(vrho, ip, 0)*mgga->VAR(v2rhotau,
  ip, 0) + mgga->VAR(v2rho2, ip, 0);
out->VAR(v2rhosigma, ip, 0) = ked1->VAR(v2rhosigma, ip, 0)*mgga->VAR(vtau, ip, 0) + ked1->VAR(vrho,
  ip, 0)*mgga->VAR(v2sigmatau, ip, 0) + ked1->VAR(vsigma, ip, 0)*(ked1->VAR(vrho, ip,
  0)*mgga->VAR(v2tau2, ip, 0) + mgga->VAR(v2rhotau, ip, 0)) + mgga->VAR(v2rhosigma, ip, 0);
out->VAR(v2rholapl, ip, 0) = ked1->VAR(v2rholapl, ip, 0)*mgga->VAR(vtau, ip, 0) + ked1->VAR(vrho,
  ip, 0)*mgga->VAR(v2lapltau, ip, 0) + ked1->VAR(vlapl, ip, 0)*(ked1->VAR(vrho, ip,
  0)*mgga->VAR(v2tau2, ip, 0) + mgga->VAR(v2rhotau, ip, 0)) + mgga->VAR(v2rholapl, ip, 0);
out->VAR(v2sigma2, ip, 0) = ked1->VAR(v2sigma2, ip, 0)*mgga->VAR(vtau, ip, 0) + ked1->VAR(vsigma,
  ip, 0)*ked1->VAR(vsigma, ip, 0)*mgga->VAR(v2tau2, ip, 0) + 2*ked1->VAR(vsigma, ip,
  0)*mgga->VAR(v2sigmatau, ip, 0) + mgga->VAR(v2sigma2, ip, 0);
out->VAR(v2sigmalapl, ip, 0) = ked1->VAR(v2sigmalapl, ip, 0)*mgga->VAR(vtau, ip, 0) +
  ked1->VAR(vsigma, ip, 0)*mgga->VAR(v2lapltau, ip, 0) + ked1->VAR(vlapl, ip, 0)*(ked1->VAR(vsigma,
  ip, 0)*mgga->VAR(v2tau2, ip, 0) + mgga->VAR(v2sigmatau, ip, 0)) + mgga->VAR(v2sigmalapl, ip, 0);
out->VAR(v2lapl2, ip, 0) = ked1->VAR(v2lapl2, ip, 0)*mgga->VAR(vtau, ip, 0) + ked1->VAR(vlapl, ip,
  0)*ked1->VAR(vlapl, ip, 0)*mgga->VAR(v2tau2, ip, 0) + 2*ked1->VAR(vlapl, ip,
  0)*mgga->VAR(v2lapltau, ip, 0) + mgga->VAR(v2lapl2, ip, 0);

if(p->nspin == XC_POLARIZED){
  out->VAR(v2rho2, ip, 1) = ked1->VAR(vrho, ip, 0)*(ked2->VAR(vrho, ip, 0)*mgga->VAR(v2tau2, ip, 1)
    + mgga->VAR(v2rhotau, ip, 2)) + ked2->VAR(vrho, ip, 0)*mgga->VAR(v2rhotau, ip, 1) +
    mgga->VAR(v2rho2, ip, 1);
  out->VAR(v2rho2, ip, 2) = ked2->VAR(v2rho2, ip, 0)*mgga->VAR(vtau, ip, 1) + ked2->VAR(vrho, ip,
    0)*ked2->VAR(vrho, ip, 0)*mgga->VAR(v2tau2, ip, 2) + 2*ked2->VAR(vrho, ip,
    0)*mgga->VAR(v2rhotau, ip, 3) + mgga->VAR(v2rho2, ip, 2);
  out->VAR(v2rhosigma, ip, 1) = ked1->VAR(vrho, ip, 0)*mgga->VAR(v2sigmatau, ip, 2) +
    mgga->VAR(v2rhosigma, ip, 1);
  out->VAR(v2rhosigma, ip, 2) = ked1->VAR(vrho, ip, 0)*mgga->VAR(v2sigmatau, ip, 4) +
    ked2->VAR(vsigma, ip, 0)*(ked1->VAR(vrho, ip, 0)*mgga->VAR(v2tau2, ip, 1) + mgga->VAR(v2rhotau,
    ip, 1)) + mgga->VAR(v2rhosigma, ip, 2);
  out->VAR(v2rhosigma, ip, 3) = ked2->VAR(vrho, ip, 0)*mgga->VAR(v2sigmatau, ip, 1) +
    ked1->VAR(vsigma, ip, 0)*(ked2->VAR(vrho, ip, 0)*mgga->VAR(v2tau2, ip, 1) + mgga->VAR(v2rhotau,
    ip, 2)) + mgga->VAR(v2rhosigma, ip, 3);
  out->VAR(v2rhosigma, ip, 4) = ked2->VAR(vrho, ip, 0)*mgga->VAR(v2sigmatau, ip, 3) +
    mgga->VAR(v2rhosigma, ip, 4);
  out->VAR(v2rhosigma, ip, 5) = ked2->VAR(v2rhosigma, ip, 0)*mgga->VAR(vtau, ip, 1) +
    ked2->VAR(vrho, ip, 0)*mgga->VAR(v2sigmatau, ip, 5) + ked2->VAR(vsigma, ip, 0)*(ked2->VAR(vrho,
    ip, 0)*mgga->VAR(v2tau2, ip, 2) + mgga->VAR(v2rhotau, ip, 3)) + mgga->VAR(v2rhosigma, ip, 5);
  out->VAR(v2rholapl, ip, 1) = ked1->VAR(vrho, ip, 0)*mgga->VAR(v2lapltau, ip, 2) + ked2->VAR(vlapl,
    ip, 0)*(ked1->VAR(vrho, ip, 0)*mgga->VAR(v2tau2, ip, 1) + mgga->VAR(v2rhotau, ip, 1)) +
    mgga->VAR(v2rholapl, ip, 1);
  out->VAR(v2rholapl, ip, 2) = ked2->VAR(vrho, ip, 0)*mgga->VAR(v2lapltau, ip, 1) + ked1->VAR(vlapl,
    ip, 0)*(ked2->VAR(vrho, ip, 0)*mgga->VAR(v2tau2, ip, 1) + mgga->VAR(v2rhotau, ip, 2)) +
    mgga->VAR(v2rholapl, ip, 2);
  out->VAR(v2rholapl, ip, 3) = ked2->VAR(v2rholapl, ip, 0)*mgga->VAR(vtau, ip, 1) + ked2->VAR(vrho,
    ip, 0)*mgga->VAR(v2lapltau, ip, 3) + ked2->VAR(vlapl, ip, 0)*(ked2->VAR(vrho, ip,
    0)*mgga->VAR(v2tau2, ip, 2) + mgga->VAR(v2rhotau, ip, 3)) + mgga->VAR(v2rholapl, ip, 3);
  out->VAR(v2sigma2, ip, 1) = ked1->VAR(vsigma, ip, 0)*mgga->VAR(v2sigmatau, ip, 2) +
    mgga->VAR(v2sigma2, ip, 1);
  out->VAR(v2sigma2, ip, 2) = ked1->VAR(vsigma, ip, 0)*(ked2->VAR(vsigma, ip, 0)*mgga->VAR(v2tau2,
    ip, 1) + mgga->VAR(v2sigmatau, ip, 4)) + ked2->VAR(vsigma, ip, 0)*mgga->VAR(v2sigmatau, ip, 1) +
    mgga->VAR(v2sigma2, ip, 2);
  out->VAR(v2sigma2, ip, 3) = mgga->VAR(v2sigma2, ip, 3);
  out->VAR(v2sigma2, ip, 4) = ked2->VAR(vsigma, ip, 0)*mgga->VAR(v2sigmatau, ip, 3) +
    mgga->VAR(v2sigma2, ip, 4);
  out->VAR(v2sigma2, ip, 5) = ked2->VAR(v2sigma2, ip, 0)*mgga->VAR(vtau, ip, 1) + ked2->VAR(vsigma,
    ip, 0)*ked2->VAR(vsigma, ip, 0)*mgga->VAR(v2tau2, ip, 2) + 2*ked2->VAR(vsigma, ip,
    0)*mgga->VAR(v2sigmatau, ip, 5) + mgga->VAR(v2sigma2, ip, 5);
  out->VAR(v2sigmalapl, ip, 1) = ked1->VAR(vsigma, ip, 0)*mgga->VAR(v2lapltau, ip, 2) +
    ked2->VAR(vlapl, ip, 0)*(ked1->VAR(vsigma, ip, 0)*mgga->VAR(v2tau2, ip, 1) +
    mgga->VAR(v2sigmatau, ip, 1)) + mgga->VAR(v2sigmalapl, ip, 1);
  out->VAR(v2sigmalapl, ip, 2) = ked1->VAR(vlapl, ip, 0)*mgga->VAR(v2sigmatau, ip, 2) +
    mgga->VAR(v2sigmalapl, ip, 2);
  out->VAR(v2sigmalapl, ip, 3) = ked2->VAR(vlapl, ip, 0)*mgga->VAR(v2sigmatau, ip, 3) +
    mgga->VAR(v2sigmalapl, ip, 3);
  out->VAR(v2sigmalapl, ip, 4) = ked2->VAR(vsigma, ip, 0)*mgga->VAR(v2lapltau, ip, 1) +
    ked1->VAR(vlapl, ip, 0)*(ked2->VAR(vsigma, ip, 0)*mgga->VAR(v2tau2, ip, 1) +
    mgga->VAR(v2sigmatau, ip, 4)) + mgga->VAR(v2sigmalapl, ip, 4);
  out->VAR(v2sigmalapl, ip, 5) = ked2->VAR(v2sigmalapl, ip, 0)*mgga->VAR(vtau, ip, 1) +
    ked2->VAR(vsigma, ip, 0)*mgga->VAR(v2lapltau, ip, 3) + ked2->VAR(vlapl, ip,
    0)*(ked2->VAR(vsigma, ip, 0)*mgga->VAR(v2tau2, ip, 2) + mgga->VAR(v2sigmatau, ip, 5)) +
    mgga->VAR(v2sigmalapl, ip, 5);
  out->VAR(v2lapl2, ip, 1) = ked1->VAR(vlapl, ip, 0)*(ked2->VAR(vlapl, ip, 0)*mgga->VAR(v2tau2, ip,
    1) + mgga->VAR(v2lapltau, ip, 2)) + ked2->VAR(vlapl, ip, 0)*mgga->VAR(v2lapltau, ip, 1) +
    mgga->VAR(v2lapl2, ip, 1);
  out->VAR(v2lapl2, ip, 2) = ked2->VAR(v2lapl2, ip, 0)*mgga->VAR(vtau, ip, 1) + ked2->VAR(vlapl, ip,
    0)*ked2->VAR(vlapl, ip, 0)*mgga->VAR(v2tau2, ip, 2) + 2*ked2->VAR(vlapl, ip,
    0)*mgga->VAR(v2lapltau, ip, 3) + mgga->VAR(v2lapl2, ip, 2);
}

