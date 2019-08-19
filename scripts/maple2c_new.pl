#!/usr/bin/env perl

# Copyright (C) 2017 M.A.L. Marques
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

use Data::Dumper;
use FindBin;                 # locate this script
use lib "$FindBin::Bin/.";   # use current directory

use maple2c_work;

die "Usage: $0 srcdir functional max_order (optional args) (simplify)\n"
    if ($#ARGV < 2);

$config{"srcdir"}     = $ARGV[0];
$config{"functional"} = $ARGV[1];
$config{"max_order"}  = $ARGV[2];
$config{"simplify"}   = ($#ARGV >= 3 && $ARGV[3] eq "yes") ? 1 : 0;
$config{"prefix"}     = "";

$config{"simplify_begin"} = ($config{'simplify'} == 1) ? "simplify(" : "";
$config{"simplify_end"}   = ($config{'simplify'} == 1) ? ", symbolic)" : "";

$config{"replace"} = [];

# find out where the maple file resides
if(-f $config{"functional"}.".mpl"){
  $config{"mathfile"} = $config{"functional"}.".mpl";
}elsif(-f $config{"srcdir"}."/maple/".$config{"functional"}.".mpl"){
  $config{"mathfile"} = $config{"srcdir"}."/maple/".$config{"functional"}.".mpl";
}else{
  $temp = $config{"functional"};
  $temp =~ s/hyb_//;
  $temp =~ s/_.*$//;

  if(-f $config{"srcdir"}."/maple/".$temp."_exc/".$config{"functional"}.".mpl"){
    $config{"mathfile"} =  $config{"srcdir"}."/maple/".$temp."_exc/".$config{"functional"}.".mpl";
  }elsif(-f $config{"srcdir"}."/maple/".$temp."_vxc/".$config{"functional"}.".mpl"){
    $config{"mathfile"} =  $config{"srcdir"}."/maple/".$temp."_vxc/".$config{"functional"}.".mpl";
  }
}

open my $in, '<', $config{"mathfile"} or die "File $mathfile does not exist\n";

# Find out the type of functional
while($_ = <$in>){
  if(/^\(\* type:\s(\S*)\s/){
    $config{"functype"} = $1;
  };
  if(/^\(\* replace:\s*"([^"]*)"\s*->\s*"([^"]*)"/){
    push @{$config{"replace"}}, "$1";
    push @{$config{"replace"}}, "$2";
  };
  if(/^\(\* prefix:/){
    while( ($_ = <$in>) && ! /^\*\)/ ){
      $config{"prefix"} .= $_;
    }
  }
}
close($in);

my %commands = (
  "lda_exc"   => \&work_lda_exc,
  "lda_vxc"   => \&work_lda_vxc,
  "gga_exc"   => \&work_gga_exc,
  "mgga_exc"  => \&work_mgga_exc,
    );

if ($commands{$config{"functype"}}) {
  $commands{$config{"functype"}}->();
} else {
  die "No such type: $string\n";
} 

sub work_lda_exc {
  # these are the variables that the functional depends on
  my @variables = ("rho_0_", "rho_1_");

  # the definition of the derivatives that libxc transmits to the calling program
  my @derivatives = (
    [[[0,0], "_s_zk"]],
    [[[1,0], "vrho_0_"],   [[0,1], "vrho_1_"]],
    [[[2,0], "v2rho2_0_"], [[1,1], "v2rho2_1_"], [[0,2], "v2rho2_2_"]],
    [[[3,0], "v3rho3_0_"], [[2,1], "v3rho3_1_"], [[1,2], "v3rho3_2_"], [[0,3], "v3rho3_3_"]],
    );

  # get arguments of the functions
  my ($input_args, $output_args) = maple2c_construct_arguments(\@variables, \@derivatives);

  # honor max_order
  splice @derivatives, $config{"max_order"}+1, $#derivatives, ;

  my ($der_def_unpol, @out_c_unpol) = 
      maple2c_create_derivatives(\@variables, \@derivatives, "mf", "unpol");
  my $out_c_unpol = join(", ", @out_c_unpol);
  my ($der_def_pol, @out_c_pol) = 
      maple2c_create_derivatives(\@variables, \@derivatives, "mf", "pol");
  my $out_c_pol = join(", ", @out_c_pol);

  # we join all the pieces
  my $maple_code = "
# zk is energy per unit particle
mzk  := (r0, r1) -> $config{'simplify_begin'} f(r_ws(dens(r0, r1)), zeta(r0, r1)) $config{'simplify_end'}:

(* mf is energy per unit volume *)
mf   := (r0, r1) -> eval(dens(r0, r1)*mzk(r0, r1)):

\$include <util.mpl>
";
  my $maple_zk = " _s_zk = mzk(".join(", ", @variables).")";

  # we build 3 variants of the functional, for unpolarized, ferromagnetic, and polarized densities
  @variants = (
    "unpol", "
dens := (r0, r1) -> r0:
zeta := (r0, r1) -> 0:

$der_def_unpol

$maple_code
C([$maple_zk, $out_c_unpol], optimize, deducetypes=false):
",

    "ferr", "
dens := (r0, r1) -> r0:
zeta := (r0, r1) -> 1:

$der_def_unpol

$maple_code
C([$maple_zk, $out_c_unpol], optimize, deducetypes=false):
\n",

    "pol", "
dens := (r0, r1) -> r0 + r1:
zeta := (r0, r1) -> (r0 - r1)/(r0 + r1):

$der_def_pol

$maple_code
C([$maple_zk, $out_c_pol], optimize, deducetypes=false):
\n",
      );

  maple2c_run(\@variables, \@derivatives, \@variants, 0, $input_args, $output_args);
}


sub work_lda_vxc {
  # these are the variables that the functional depends on
  my @variables = ("rho_0_", "rho_1_");

  # the definition of the derivatives that libxc transmits to the calling program
  my @derivatives1 = (
    [[[0,0], "vrho_0_"]],
    [[[1,0], "v2rho2_0_"], [[0,1], "v2rho2_1_"]],
    [[[2,0], "v3rho3_0_"], [[1,1], "v3rho3_1_"], [[0,2], "v3rho3_2_"]],
      );
  my @derivatives2 = (
    [[[0,0], "vrho_1_"]],
    [[[0,1], "v2rho2_2_"]],
    [[[0,2], "v3rho3_3_"]],    
      );
  
  my @derivatives = ();
  for(my $i=0; $i<=$#derivatives1; $i++){
    @{$derivatives[$i]} = @{$derivatives1[$i]};
    push(@{$derivatives[$i]}, @{$derivatives2[$i]});
  }

  my ($input_args, $output_args) = maple2c_construct_arguments(\@variables, \@derivatives);

  # honor max_order
  splice @derivatives1, $config{"max_order"}+1, $#derivatives1, ;
  splice @derivatives2, $config{"max_order"}+1, $#derivatives2, ;
  splice @derivatives,  $config{"max_order"}+1, $#derivatives,  ;

  # we obtain the missing pieces for maple
  # unpolarized calculation
  my ($der_def_unpol, @out_c_unpol) = 
      maple2c_create_derivatives(\@variables, \@derivatives1, "mf0", "unpol");
  my $out_c_unpol = join(", ", @out_c_unpol);

  # polarized calculation
  my ($der_def_pol, @out_c_pol1) = 
      maple2c_create_derivatives(\@variables, \@derivatives1, "mf0", "pol");
  my ($der_def_pol2, @out_c_pol2) = 
      maple2c_create_derivatives(\@variables, \@derivatives2, "mf1", "pol");

  $der_def_pol .= $der_def_pol2;
  
  push(@out_c_pol1, @out_c_pol2);
  my $out_c_pol = join(", ", sort(@out_c_pol1));

  # we join all the pieces
  my $maple_code1 = "
(* mf is the up potential *)
mzk   := (r0, r1) -> $config{'simplify_begin'} f(r_ws(dens(r0, r1)), zeta(r0, r1)) $config{'simplify_end'}:
mf0   := (r0, r1) -> eval(mzk(r0, r1)):
mf1   := (r0, r1) -> eval(mzk(r1, r0)):

\$include <util.mpl>
";

  my $maple_vrho0 = "vrho_0_ = mf0(".join(", ", @variables).")";
  my $maple_vrho1 = "vrho_1_ = mf1(".join(", ", @variables).")"; 

  # we build 3 variants of the functional, for unpolarized, ferromagnetic, and polarized densities
  @variants = (
    "unpol", "
dens := (r0, r1) -> r0:
zeta := (r0, r1) -> 0:

$der_def_unpol

$maple_code1
C([$maple_vrho0, $out_c_unpol], optimize, deducetypes=false):
",

    "ferr", "
dens := (r0, r1) -> r0:
zeta := (r0, r1) -> 1:

$der_def_unpol

$maple_code1
C([$maple_vrho0, $out_c_unpol], optimize, deducetypes=false):
",

    "pol", "
dens := (r0, r1) -> r0 + r1:
zeta := (r0, r1) -> (r0 - r1)/(r0 + r1):

$der_def_pol

$maple_code1
C([$maple_vrho0, $maple_vrho1, $out_c_pol], optimize, deducetypes=false):
"
      );

  maple2c_run(\@variables, \@derivatives, \@variants, 1, $input_args, $output_args);
}

sub work_gga_exc {
  # these are the variables that the functional depends on
  my @variables = ("rho_0_", "rho_1_", "sigma_0_", "sigma_1_", "sigma_2_");

  # the definition of the derivatives that libxc transmits to the calling program
  my @derivatives = (
    [
     [[0,0,0,0,0], "_s_zk"]
    ], [
      [[1,0,0,0,0], "vrho_0_"],   [[0,1,0,0,0], "vrho_1_"],
      [[0,0,1,0,0], "vsigma_0_"], [[0,0,0,1,0], "vsigma_1_"], [[0,0,0,0,1], "vsigma_2_"]
    ], [
      [[2,0,0,0,0], "v2rho2_0_"], [[1,1,0,0,0], "v2rho2_1_"], [[0,2,0,0,0], "v2rho2_2_"],
      [[1,0,1,0,0], "v2rhosigma_0_"], [[1,0,0,1,0], "v2rhosigma_1_"],  [[1,0,0,0,1], "v2rhosigma_2_"],
      [[0,1,1,0,0], "v2rhosigma_3_"], [[0,1,0,1,0], "v2rhosigma_4_"],  [[0,1,0,0,1], "v2rhosigma_5_"],
      [[0,0,2,0,0], "v2sigma2_0_"], [[0,0,1,1,0], "v2sigma2_1_"], [[0,0,1,0,1], "v2sigma2_2_"],
      [[0,0,0,2,0], "v2sigma2_3_"], [[0,0,0,1,1], "v2sigma2_4_"], [[0,0,0,0,2], "v2sigma2_5_"], 
    ], [
      [[3,0,0,0,0], "v3rho3_0_"], [[2,1,0,0,0], "v3rho3_1_"], 
      [[1,2,0,0,0], "v3rho3_2_"], [[0,3,0,0,0], "v3rho3_3_"],
      [[2,0,1,0,0], "v3rho2sigma_0_"], [[2,0,0,1,0], "v3rho2sigma_1_"],  [[2,0,0,0,1], "v3rho2sigma_2_"],
      [[1,1,1,0,0], "v3rho2sigma_3_"], [[1,1,0,1,0], "v3rho2sigma_4_"],  [[1,1,0,0,1], "v3rho2sigma_5_"],
      [[0,2,1,0,0], "v3rho2sigma_6_"], [[0,2,0,1,0], "v3rho2sigma_7_"],  [[0,2,0,0,1], "v3rho2sigma_8_"],
      [[1,0,2,0,0], "v3rhosigma2_0_"], [[1,0,1,1,0], "v3rhosigma2_1_"],  [[1,0,1,0,1], "v3rhosigma2_2_"],
      [[1,0,0,2,0], "v3rhosigma2_3_"], [[1,0,0,1,1], "v3rhosigma2_4_"],  [[1,0,0,0,2], "v3rhosigma2_5_"],
      [[0,1,2,0,0], "v3rhosigma2_6_"], [[0,1,1,1,0], "v3rhosigma2_7_"],  [[0,1,1,0,1], "v3rhosigma2_8_"],
      [[0,1,0,2,0], "v3rhosigma2_9_"], [[0,1,0,1,1], "v3rhosigma2_10_"], [[0,1,0,0,2], "v3rhosigma2_11_"],
      [[0,0,3,0,0], "v3sigma3_0_"], [[0,0,2,1,0], "v3sigma3_1_"], [[0,0,2,0,1], "v3sigma3_2_"],
      [[0,0,1,2,0], "v3sigma3_3_"], [[0,0,1,1,1], "v3sigma3_4_"], [[0,0,1,0,2], "v3sigma3_5_"],
      [[0,0,0,3,0], "v3sigma3_6_"], [[0,0,0,2,1], "v3sigma3_7_"], [[0,0,0,1,2], "v3sigma3_8_"],
      [[0,0,0,0,3], "v3sigma3_9_"],
    ],
      );

  # get arguments of the functions
  my ($input_args, $output_args) = maple2c_construct_arguments(\@variables, \@derivatives);

  # honor max_order
  splice @derivatives, $config{"max_order"}+1, $#derivatives, ;
  
  my ($der_def, @out_c) = 
      maple2c_create_derivatives(\@variables, \@derivatives, "mf");
  my $out_c = join(", ", @out_c);
  $out_c = ", $out_c" if ($out_c ne "");

  # we join all the pieces
  my $maple_code = "
# zk is energy per unit particle
mzk  := (r0, r1, s0, s1, s2) -> \\
  $config{'simplify_begin'} \\
    f(r_ws(dens(r0, r1)), zeta(r0, r1), xt(r0, r1, s0, s1, s2), xs0(r0, r1, s0, s2), xs1(r0, r1, s0, s2)) \\
  $config{'simplify_end'}:

(* mf is energy per unit volume *)
mf   := (r0, r1, s0, s1, s2) -> eval(dens(r0, r1)*mzk(r0, r1, s0, s1, s2)):

\$include <util.mpl>
";
  my $maple_zk = " _s_zk = mzk(".join(", ", @variables).")";

  # we build 3 variants of the functional, for unpolarized, ferromagnetic, and polarized densities
  @variants = (
    "unpol", "
dens := (r0, r1) -> r0:
zeta := (r0, r1) -> 0:
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):

$der_def

$maple_code
C([$maple_zk$out_c], optimize, deducetypes=false):
",

    "ferr", "
dens := (r0, r1) -> r0:
zeta := (r0, r1) -> 1:
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
xs1  := (r0, r1, sigma0, sigma2) -> 0:
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):

$der_def

$maple_code
C([$maple_zk$out_c], optimize, deducetypes=false):
\n",

    "pol", "
dens := (r0, r1) -> r0 + r1:
zeta := (r0, r1) -> (r0 - r1)/(r0 + r1):
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma2)/r1^(1 + 1/DIMENSIONS):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0 + 2*sigma1 + sigma2)/(r0 + r1)^(1 + 1/DIMENSIONS):

$der_def

$maple_code
C([$maple_zk$out_c], optimize, deducetypes=false):
\n",
      );

  maple2c_run(\@variables, \@derivatives, \@variants, 0, $input_args, $output_args);
}

sub work_mgga_exc {
  # these are the variables that the functional depends on
  my @variables = ("rho_0_", "rho_1_", "sigma_0_", "sigma_1_", "sigma_2_", "lapl_0_", "lapl_1_", "tau_0_", "tau_1_");

  # the definition of the derivatives that libxc transmits to the calling program
  my @derivatives = (
    [
     [[0,0,0,0,0,0,0,0,0], "_s_zk"]
    ], [
      [[1,0,0,0,0,0,0,0,0], "vrho_0_"],   
      [[0,1,0,0,0,0,0,0,0], "vrho_1_"],
      [[0,0,1,0,0,0,0,0,0], "vsigma_0_"],
      [[0,0,0,1,0,0,0,0,0], "vsigma_1_"],
      [[0,0,0,0,1,0,0,0,0], "vsigma_2_"],
      [[0,0,0,0,0,1,0,0,0], "vlapl_0_"],
      [[0,0,0,0,0,0,1,0,0], "vlapl_1_"],
      [[0,0,0,0,0,0,0,1,0], "vtau_0_"],
      [[0,0,0,0,0,0,0,0,1], "vtau_1_"],
    ], [
      [[2,0,0,0,0,0,0,0,0], "v2rho2_0_"],
      [[1,1,0,0,0,0,0,0,0], "v2rho2_1_"],
      [[0,2,0,0,0,0,0,0,0], "v2rho2_2_"],
      [[1,0,1,0,0,0,0,0,0], "v2rhosigma_0_"],
      [[1,0,0,1,0,0,0,0,0], "v2rhosigma_1_"],
      [[1,0,0,0,1,0,0,0,0], "v2rhosigma_2_"],
      [[0,1,1,0,0,0,0,0,0], "v2rhosigma_3_"],
      [[0,1,0,1,0,0,0,0,0], "v2rhosigma_4_"],
      [[0,1,0,0,1,0,0,0,0], "v2rhosigma_5_"],
      [[1,0,0,0,0,1,0,0,0], "v2rholapl_0_"],
      [[1,0,0,0,0,0,1,0,0], "v2rholapl_1_"],
      [[0,1,0,0,0,1,0,0,0], "v2rholapl_2_"],
      [[0,1,0,0,0,0,1,0,0], "v2rholapl_3_"],
      [[1,0,0,0,0,0,0,1,0], "v2rhotau_0_"],
      [[1,0,0,0,0,0,0,0,1], "v2rhotau_1_"],
      [[0,1,0,0,0,0,0,1,0], "v2rhotau_2_"],
      [[0,1,0,0,0,0,0,0,1], "v2rhotau_3_"],
      [[0,0,2,0,0,0,0,0,0], "v2sigma2_0_"],
      [[0,0,1,1,0,0,0,0,0], "v2sigma2_1_"],
      [[0,0,1,0,1,0,0,0,0], "v2sigma2_2_"],
      [[0,0,0,2,0,0,0,0,0], "v2sigma2_3_"],
      [[0,0,0,1,1,0,0,0,0], "v2sigma2_4_"],
      [[0,0,0,0,2,0,0,0,0], "v2sigma2_5_"],
      [[0,0,1,0,0,1,0,0,0], "v2sigmalapl_0_"],
      [[0,0,1,0,0,0,1,0,0], "v2sigmalapl_1_"],
      [[0,0,0,1,0,1,0,0,0], "v2sigmalapl_2_"],
      [[0,0,0,1,0,0,1,0,0], "v2sigmalapl_3_"],
      [[0,0,0,0,1,1,0,0,0], "v2sigmalapl_4_"],
      [[0,0,0,0,1,0,1,0,0], "v2sigmalapl_5_"],
      [[0,0,1,0,0,0,0,1,0], "v2sigmatau_0_"],
      [[0,0,1,0,0,0,0,0,1], "v2sigmatau_1_"],
      [[0,0,0,1,0,0,0,1,0], "v2sigmatau_2_"],
      [[0,0,0,1,0,0,0,0,1], "v2sigmatau_3_"],
      [[0,0,0,0,1,0,0,1,0], "v2sigmatau_4_"],
      [[0,0,0,0,1,0,0,0,1], "v2sigmatau_5_"],
      [[0,0,0,0,0,2,0,0,0], "v2lapl2_0_"],
      [[0,0,0,0,0,1,1,0,0], "v2lapl2_1_"],
      [[0,0,0,0,0,0,2,0,0], "v2lapl2_2_"],
      [[0,0,0,0,0,1,0,1,0], "v2lapltau_0_"],
      [[0,0,0,0,0,1,0,0,1], "v2lapltau_1_"],
      [[0,0,0,0,0,0,1,1,0], "v2lapltau_2_"],
      [[0,0,0,0,0,0,1,0,1], "v2lapltau_3_"],
      [[0,0,0,0,0,0,0,2,0], "v2tau2_0_"],
      [[0,0,0,0,0,0,0,1,1], "v2tau2_1_"],
      [[0,0,0,0,0,0,0,0,2], "v2tau2_2_"],
    ], [
      [[3,0,0,0,0,0,0,0,0], "v3rho3_0_"],
      [[2,1,0,0,0,0,0,0,0], "v3rho3_1_"], 
      [[1,2,0,0,0,0,0,0,0], "v3rho3_2_"],
      [[0,3,0,0,0,0,0,0,0], "v3rho3_3_"],
      [[2,0,1,0,0,0,0,0,0], "v3rho2sigma_0_"],
      [[2,0,0,1,0,0,0,0,0], "v3rho2sigma_1_"], 
      [[2,0,0,0,1,0,0,0,0], "v3rho2sigma_2_"],
      [[1,1,1,0,0,0,0,0,0], "v3rho2sigma_3_"],
      [[1,1,0,1,0,0,0,0,0], "v3rho2sigma_4_"], 
      [[1,1,0,0,1,0,0,0,0], "v3rho2sigma_5_"],
      [[0,2,1,0,0,0,0,0,0], "v3rho2sigma_6_"],
      [[0,2,0,1,0,0,0,0,0], "v3rho2sigma_7_"],
      [[0,2,0,0,1,0,0,0,0], "v3rho2sigma_8_"],
      [[2,0,0,0,0,1,0,0,0], "v3rho2lapl_0_"],
      [[2,0,0,0,0,0,1,0,0], "v3rho2lapl_1_"],
      [[1,1,0,0,0,1,0,0,0], "v3rho2lapl_2_"],
      [[1,1,0,0,0,0,1,0,0], "v3rho2lapl_3_"],
      [[0,2,0,0,0,1,0,0,0], "v3rho2lapl_4_"],
      [[0,2,0,0,0,0,1,0,0], "v3rho2lapl_5_"],
      [[2,0,0,0,0,0,0,1,0], "v3rho2tau_0_"],
      [[2,0,0,0,0,0,0,0,1], "v3rho2tau_1_"],
      [[1,1,0,0,0,0,0,1,0], "v3rho2tau_2_"],
      [[1,1,0,0,0,0,0,0,1], "v3rho2tau_3_"],
      [[0,2,0,0,0,0,0,1,0], "v3rho2tau_4_"],
      [[0,2,0,0,0,0,0,0,1], "v3rho2tau_5_"],
      [[1,0,2,0,0,0,0,0,0], "v3rhosigma2_0_"],
      [[1,0,1,1,0,0,0,0,0], "v3rhosigma2_1_"],
      [[1,0,1,0,1,0,0,0,0], "v3rhosigma2_2_"],
      [[1,0,0,2,0,0,0,0,0], "v3rhosigma2_3_"],
      [[1,0,0,1,1,0,0,0,0], "v3rhosigma2_4_"],
      [[1,0,0,0,2,0,0,0,0], "v3rhosigma2_5_"],
      [[0,1,2,0,0,0,0,0,0], "v3rhosigma2_6_"],
      [[0,1,1,1,0,0,0,0,0], "v3rhosigma2_7_"],
      [[0,1,1,0,1,0,0,0,0], "v3rhosigma2_8_"],
      [[0,1,0,2,0,0,0,0,0], "v3rhosigma2_9_"],
      [[0,1,0,1,1,0,0,0,0], "v3rhosigma2_10_"],
      [[0,1,0,0,2,0,0,0,0], "v3rhosigma2_11_"],
      [[1,0,1,0,0,1,0,0,0], "v3rhosigmalapl_0_"],
      [[1,0,1,0,0,0,1,0,0], "v3rhosigmalapl_1_"],
      [[1,0,0,1,0,1,0,0,0], "v3rhosigmalapl_2_"],
      [[1,0,0,1,0,0,1,0,0], "v3rhosigmalapl_3_"],
      [[1,0,0,0,1,1,0,0,0], "v3rhosigmalapl_4_"],
      [[1,0,0,0,1,0,1,0,0], "v3rhosigmalapl_5_"],
      [[0,1,1,0,0,1,0,0,0], "v3rhosigmalapl_6_"],
      [[0,1,1,0,0,0,1,0,0], "v3rhosigmalapl_7_"],
      [[0,1,0,1,0,1,0,0,0], "v3rhosigmalapl_8_"],
      [[0,1,0,1,0,0,1,0,0], "v3rhosigmalapl_9_"],
      [[0,1,0,0,1,1,0,0,0], "v3rhosigmalapl_10_"],
      [[0,1,0,0,1,0,1,0,0], "v3rhosigmalapl_11_"],
      [[1,0,1,0,0,0,0,1,0], "v3rhosigmatau_0_"],
      [[1,0,1,0,0,0,0,0,1], "v3rhosigmatau_1_"],
      [[1,0,0,1,0,0,0,1,0], "v3rhosigmatau_2_"],
      [[1,0,0,1,0,0,0,0,1], "v3rhosigmatau_3_"],
      [[1,0,0,0,1,0,0,1,0], "v3rhosigmatau_4_"],
      [[1,0,0,0,1,0,0,0,1], "v3rhosigmatau_5_"],
      [[0,1,1,0,0,0,0,1,0], "v3rhosigmatau_6_"],
      [[0,1,1,0,0,0,0,0,1], "v3rhosigmatau_7_"],
      [[0,1,0,1,0,0,0,1,0], "v3rhosigmatau_8_"],
      [[0,1,0,1,0,0,0,0,1], "v3rhosigmatau_9_"],
      [[0,1,0,0,1,0,0,1,0], "v3rhosigmatau_10_"],
      [[0,1,0,0,1,0,0,0,1], "v3rhosigmatau_11_"],
      [[1,0,0,0,0,2,0,0,0], "v3rholapl2_0_"],
      [[1,0,0,0,0,1,1,0,0], "v3rholapl2_1_"],
      [[1,0,0,0,0,0,2,0,0], "v3rholapl2_2_"],
      [[0,1,0,0,0,2,0,0,0], "v3rholapl2_3_"],
      [[0,1,0,0,0,1,1,0,0], "v3rholapl2_4_"],
      [[0,1,0,0,0,0,2,0,0], "v3rholapl2_5_"],
      [[1,0,0,0,0,1,0,1,0], "v3rholapltau_0_"],
      [[1,0,0,0,0,1,0,0,1], "v3rholapltau_1_"],
      [[1,0,0,0,0,0,1,1,0], "v3rholapltau_2_"],
      [[1,0,0,0,0,0,1,0,1], "v3rholapltau_3_"],
      [[0,1,0,0,0,1,0,1,0], "v3rholapltau_4_"],
      [[0,1,0,0,0,1,0,0,1], "v3rholapltau_5_"],
      [[0,1,0,0,0,0,1,1,0], "v3rholapltau_6_"],
      [[0,1,0,0,0,0,1,0,1], "v3rholapltau_7_"],
      [[1,0,0,0,0,0,0,2,0], "v3rhotau2_0_"],
      [[1,0,0,0,0,0,0,1,1], "v3rhotau2_1_"],
      [[1,0,0,0,0,0,0,0,2], "v3rhotau2_2_"],
      [[0,1,0,0,0,0,0,2,0], "v3rhotau2_3_"],
      [[0,1,0,0,0,0,0,1,1], "v3rhotau2_4_"],
      [[0,1,0,0,0,0,0,0,2], "v3rhotau2_5_"],
      [[0,0,3,0,0,0,0,0,0], "v3sigma3_0_"],
      [[0,0,2,1,0,0,0,0,0], "v3sigma3_1_"],
      [[0,0,2,0,1,0,0,0,0], "v3sigma3_2_"],
      [[0,0,1,2,0,0,0,0,0], "v3sigma3_3_"],
      [[0,0,1,1,1,0,0,0,0], "v3sigma3_4_"],
      [[0,0,1,0,2,0,0,0,0], "v3sigma3_5_"],
      [[0,0,0,3,0,0,0,0,0], "v3sigma3_6_"],
      [[0,0,0,2,1,0,0,0,0], "v3sigma3_7_"],
      [[0,0,0,1,2,0,0,0,0], "v3sigma3_8_"],
      [[0,0,0,0,3,0,0,0,0], "v3sigma3_9_"],
      [[0,0,2,0,0,1,0,0,0], "v3sigma2lapl_0_"],
      [[0,0,2,0,0,0,1,0,0], "v3sigma2lapl_1_"],
      [[0,0,1,1,0,1,0,0,0], "v3sigma2lapl_2_"],
      [[0,0,1,1,0,0,1,0,0], "v3sigma2lapl_3_"],
      [[0,0,1,0,1,1,0,0,0], "v3sigma2lapl_4_"],
      [[0,0,1,0,1,0,1,0,0], "v3sigma2lapl_5_"],
      [[0,0,0,2,0,1,0,0,0], "v3sigma2lapl_6_"],
      [[0,0,0,2,0,0,1,0,0], "v3sigma2lapl_7_"],
      [[0,0,0,1,1,1,0,0,0], "v3sigma2lapl_8_"],
      [[0,0,0,1,1,0,1,0,0], "v3sigma2lapl_9_"],
      [[0,0,0,0,2,1,0,0,0], "v3sigma2lapl_10_"],
      [[0,0,0,0,2,0,1,0,0], "v3sigma2lapl_11_"],
      [[0,0,2,0,0,0,0,1,0], "v3sigma2tau_0_"],
      [[0,0,2,0,0,0,0,0,1], "v3sigma2tau_1_"],
      [[0,0,1,1,0,0,0,1,0], "v3sigma2tau_2_"],
      [[0,0,1,1,0,0,0,0,1], "v3sigma2tau_3_"],
      [[0,0,1,0,1,0,0,1,0], "v3sigma2tau_4_"],
      [[0,0,1,0,1,0,0,0,1], "v3sigma2tau_5_"],
      [[0,0,0,2,0,0,0,1,0], "v3sigma2tau_6_"],
      [[0,0,0,2,0,0,0,0,1], "v3sigma2tau_7_"],
      [[0,0,0,1,1,0,0,1,0], "v3sigma2tau_8_"],
      [[0,0,0,1,1,0,0,0,1], "v3sigma2tau_9_"],
      [[0,0,0,0,2,0,0,1,0], "v3sigma2tau_10_"],
      [[0,0,0,0,2,0,0,0,1], "v3sigma2tau_11_"],
      [[0,0,1,0,0,2,0,0,0], "v3sigmalapl2_0_"],
      [[0,0,1,0,0,1,1,0,0], "v3sigmalapl2_1_"],
      [[0,0,1,0,0,0,2,0,0], "v3sigmalapl2_2_"],
      [[0,0,0,1,0,2,0,0,0], "v3sigmalapl2_3_"],
      [[0,0,0,1,0,1,1,0,0], "v3sigmalapl2_4_"],
      [[0,0,0,1,0,0,2,0,0], "v3sigmalapl2_5_"],
      [[0,0,0,0,1,2,0,0,0], "v3sigmalapl2_6_"],
      [[0,0,0,0,1,1,1,0,0], "v3sigmalapl2_7_"],
      [[0,0,0,0,1,0,2,0,0], "v3sigmalapl2_8_"],
      [[0,0,1,0,0,1,0,1,0], "v3sigmalapltau_0_"],
      [[0,0,1,0,0,1,0,0,1], "v3sigmalapltau_1_"],
      [[0,0,1,0,0,0,1,1,0], "v3sigmalapltau_2_"],
      [[0,0,1,0,0,0,1,0,1], "v3sigmalapltau_3_"],
      [[0,0,0,1,0,1,0,1,0], "v3sigmalapltau_4_"],
      [[0,0,0,1,0,1,0,0,1], "v3sigmalapltau_5_"],
      [[0,0,0,1,0,0,1,1,0], "v3sigmalapltau_6_"],
      [[0,0,0,1,0,0,1,0,1], "v3sigmalapltau_7_"],
      [[0,0,0,0,1,1,0,1,0], "v3sigmalapltau_8_"],
      [[0,0,0,0,1,1,0,0,1], "v3sigmalapltau_9_"],
      [[0,0,0,0,1,0,1,1,0], "v3sigmalapltau_10_"],
      [[0,0,0,0,1,0,1,0,1], "v3sigmalapltau_11_"],
      [[0,0,1,0,0,0,0,2,0], "v3sigmatau2_0_"],
      [[0,0,1,0,0,0,0,1,1], "v3sigmatau2_1_"],
      [[0,0,1,0,0,0,0,0,2], "v3sigmatau2_2_"],
      [[0,0,0,1,0,0,0,2,0], "v3sigmatau2_3_"],
      [[0,0,0,1,0,0,0,1,1], "v3sigmatau2_4_"],
      [[0,0,0,1,0,0,0,0,2], "v3sigmatau2_5_"],
      [[0,0,0,0,1,0,0,2,0], "v3sigmatau2_6_"],
      [[0,0,0,0,1,0,0,1,1], "v3sigmatau2_7_"],
      [[0,0,0,0,1,0,0,0,2], "v3sigmatau2_8_"],
      [[0,0,0,0,0,3,0,0,0], "v3lapl3_0_"],
      [[0,0,0,0,0,2,1,0,0], "v3lapl3_1_"],
      [[0,0,0,0,0,1,2,0,0], "v3lapl3_2_"],
      [[0,0,0,0,0,0,3,0,0], "v3lapl3_3_"],
      [[0,0,0,0,0,2,0,1,0], "v3lapl2tau_0_"],
      [[0,0,0,0,0,2,0,0,1], "v3lapl2tau_1_"],
      [[0,0,0,0,0,1,1,1,0], "v3lapl2tau_2_"],
      [[0,0,0,0,0,1,1,0,1], "v3lapl2tau_3_"],
      [[0,0,0,0,0,0,2,1,0], "v3lapl2tau_4_"],
      [[0,0,0,0,0,0,2,0,1], "v3lapl2tau_5_"],
      [[0,0,0,0,0,1,0,2,0], "v3lapltau2_0_"],
      [[0,0,0,0,0,1,0,1,1], "v3lapltau2_1_"],
      [[0,0,0,0,0,1,0,0,2], "v3lapltau2_2_"],
      [[0,0,0,0,0,0,1,2,0], "v3lapltau2_3_"],
      [[0,0,0,0,0,0,1,1,1], "v3lapltau2_4_"],
      [[0,0,0,0,0,0,1,0,2], "v3lapltau2_5_"],
      [[0,0,0,0,0,0,0,3,0], "v3tau3_0_"],
      [[0,0,0,0,0,0,0,2,1], "v3tau3_1_"],
      [[0,0,0,0,0,0,0,1,2], "v3tau3_2_"],
      [[0,0,0,0,0,0,0,0,3], "v3tau3_3_"],
    ],
      );

  # get arguments of the functions
  my ($input_args, $output_args) = maple2c_construct_arguments(\@variables, \@derivatives);

  # honor max_order
  splice @derivatives, $config{"max_order"}+1, $#derivatives, ;
  
  my ($der_def, @out_c) = 
      maple2c_create_derivatives(\@variables, \@derivatives, "mf");
  my $out_c = join(", ", @out_c);
  $out_c = ", $out_c" if ($out_c ne "");

  # we join all the pieces
  my $maple_code = "
# zk is energy per unit particle
mzk  := (r0, r1, s0, s1, s2, l0, l1, tau0, tau1) -> \\
  $config{'simplify_begin'} \\
    f(r_ws(dens(r0, r1)), zeta(r0, r1), xt(r0, r1, s0, s1, s2), xs0(r0, r1, s0, s2), xs1(r0, r1, s0, s2), u0(r0, r1, l0, l1), u1(r0, r1, l0, l1), t0(r0, r1, tau0, tau1), t1(r0, r1, tau0, tau1)) \\
  $config{'simplify_end'}:

(* mf is energy per unit volume *)
mf   := (r0, r1, s0, s1, s2, l0, l1, tau0, tau1) -> eval(dens(r0, r1)*mzk(r0, r1, s0, s1, s2, l0, l1, tau0, tau1)):

\$include <util.mpl>
";
  my $maple_zk = " _s_zk = mzk(".join(", ", @variables).")";

  # we build 3 variants of the functional, for unpolarized, ferromagnetic, and polarized densities
  @variants = (
    "unpol", "
dens := (r0, r1) -> r0:
zeta := (r0, r1) -> 0:
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0/4)/((r0/2)^(1 + 1/DIMENSIONS)):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
u0   := (r0, r1, l0, l1) -> (l0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
u1   := (r0, r1, l0, l1) -> (l0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
t0   := (r0, r1, tau0, tau1) -> (tau0/2)/((r0/2)^(1 + 2/DIMENSIONS)):
t1   := (r0, r1, tau0, tau1) -> (tau0/2)/((r0/2)^(1 + 2/DIMENSIONS)):

$der_def

$maple_code
C([$maple_zk$out_c], optimize, deducetypes=false):
",

    "ferr", "
dens := (r0, r1) -> r0:
zeta := (r0, r1) -> 1:
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
xs1  := (r0, r1, sigma0, sigma2) -> 0:
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
u0   := (r0, r1, l0, l1) -> l0/(r0^(1 + 2/DIMENSIONS)):
u1   := (r0, r1, l0, l1) -> 0:
t0   := (r0, r1, tau0, tau1) -> tau0/(r0^(1 + 2/DIMENSIONS)):
t1   := (r0, r1, tau0, tau1) -> 0:

$der_def

$maple_code
C([$maple_zk$out_c], optimize, deducetypes=false):
\n",

    "pol", "
dens := (r0, r1) -> r0 + r1:
zeta := (r0, r1) -> (r0 - r1)/(r0 + r1):
xs0  := (r0, r1, sigma0, sigma2) -> sqrt(sigma0)/r0^(1 + 1/DIMENSIONS):
xs1  := (r0, r1, sigma0, sigma2) -> sqrt(sigma2)/r1^(1 + 1/DIMENSIONS):
xt   := (r0, r1, sigma0, sigma1, sigma2) -> sqrt(sigma0 + 2*sigma1 + sigma2)/(r0 + r1)^(1 + 1/DIMENSIONS):
u0   := (r0, r1, l0, l1) -> l0/(r0^(1 + 2/DIMENSIONS)):
u1   := (r0, r1, l0, l1) -> l1/(r1^(1 + 2/DIMENSIONS)):
t0   := (r0, r1, tau0, tau1) -> tau0/(r0^(1 + 2/DIMENSIONS)):
t1   := (r0, r1, tau0, tau1) -> tau1/(r1^(1 + 2/DIMENSIONS)):

$der_def

$maple_code
C([$maple_zk$out_c], optimize, deducetypes=false):
\n",
      );

  maple2c_run(\@variables, \@derivatives, \@variants, 0, $input_args, $output_args);
}
