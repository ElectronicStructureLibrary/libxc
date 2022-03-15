out->VAR(vrho, ip, 0) = ked1->VAR(vrho, ip, 0)*mgga->VAR(vtau, ip, 0) + mgga->VAR(vrho, ip, 0);
out->VAR(vsigma, ip, 0) = ked1->VAR(vsigma, ip, 0)*mgga->VAR(vtau, ip, 0) + mgga->VAR(vsigma, ip,
  0);
out->VAR(vlapl, ip, 0) = ked1->VAR(vlapl, ip, 0)*mgga->VAR(vtau, ip, 0) + mgga->VAR(vlapl, ip, 0);

if(p->nspin == XC_POLARIZED){
  out->VAR(vrho, ip, 1) = ked2->VAR(vrho, ip, 0)*mgga->VAR(vtau, ip, 1) + mgga->VAR(vrho, ip, 1);
  out->VAR(vsigma, ip, 1) = mgga->VAR(vsigma, ip, 1);
  out->VAR(vsigma, ip, 2) = ked2->VAR(vsigma, ip, 0)*mgga->VAR(vtau, ip, 1) + mgga->VAR(vsigma, ip,
    2);
  out->VAR(vlapl, ip, 1) = ked2->VAR(vlapl, ip, 0)*mgga->VAR(vtau, ip, 1) + mgga->VAR(vlapl, ip, 1);
}

