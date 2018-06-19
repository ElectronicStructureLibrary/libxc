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
$config{"mathfile"}   = $config{"srcdir"}."/maple/".$config{"functional"}.".mpl";
$config{"prefix"}     = "";

$config{"simplify_begin"} = ($config{'simplify'} == 1) ? "simplify(" : "";
$config{"simplify_end"}   = ($config{'simplify'} == 1) ? ", symbolic)" : "";
    
# Find out the type of functional
open my $in, '<', $config{"mathfile"} or die "File $mathfile does not exist\n";
while($_ = <$in>){
  if(/^\(\* type:\s(\S*)\s/){
    $config{"functype"} = $1;
  };
  if(/^\(\* prefix:/){
    while( ($_ = <$in>) && ! /^\*\)/ ){
      $config{"prefix"} .= $_;
    }
  }
}
close($in);

my %commands = (
  "lda_exc"     => \&work_lda_exc,
  "lda_vxc"     => \&work_lda_vxc,
#  "work_gga_x"   => \&work_gga_x,
#  "work_gga_c"   => \&work_gga_c,
#  "work_mgga_x"  => \&work_mgga_x,
#  "work_mgga_c"  => \&work_mgga_c,
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
  
  # we build 3 variants of the functional, for unpolarized, ferromagnetic, and polarized densities
  @variants = (
    "unpol", "
dens := (r0, r1) -> r0:
zeta := (r0, r1) -> 0:
\n",
    "ferr", "
dens := (r0, r1) -> r0:
zeta := (r0, r1) -> 1:
\n",
    "pol", "
dens := (r0, r1) -> r0 + r1:
zeta := (r0, r1) -> (r0 - r1)/(r0 + r1):
\n",
      );

  # we obtain the missing pieces for maple
  my ($out1, $out2) = maple2c_create_derivatives(\@variables, \@derivatives);

  my $maple_code = "
# zk is energy per unit particle
mzk  := (r0, r1) -> $config{'simplify_begin'} f(r_ws(dens(r0, r1)), zeta(r0, r1)) $config{'simplify_end'}:

(* mf is energy per unit volume *)
mf   := (r0, r1) -> $config{'simplify_begin'} dens(r0, r1)*mzk(r0, r1) $config{'simplify_end'}:

\$include <util.mpl>

$out1

C([_s_zk = mzk(".join(", ", @{$variables})."), $out2], optimize, deducetypes=false):
";
  
  maple2c_run(\@variables, \@derivatives, \@variants, $maple_code, 0);
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

  # we obtain the missing pieces for maple
  # unpolarized calculation
  my ($der_def_unpol, @out_c_unpol) = 
      maple2c_create_derivatives(\@variables, \@derivatives1, "mf0", "unpol");

  # polarized calculation
  my ($der_def_pol, @out_c_unpol1) = 
      maple2c_create_derivatives(\@variables, \@derivatives1, "mf0", "pol");
  my ($der_def_pol2, @out_c_unpol2) = 
      maple2c_create_derivatives(\@variables, \@derivatives2, "mf1", "pol");

  $der_def_pol .= $der_def_pol2;
  
  push(@out_c_unpol1, @out_c_unpol2);
  my $out_c_unpol = join(", ", sort(@out_c_unpol1));

  # we join all the pieces
  my $maple_code1 = "
(* mf is the up potential *)
mzk   := (r0, r1) -> $config{'simplify_begin'} f(r_ws(dens(r0, r1)), zeta(r0, r1)) $config{'simplify_end'}:
mf0   := (r0, r1) -> mzk(r0, r1):
mf1   := (r0, r1) -> mzk(r1, r0):

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
C([$maple_vrho0, $maple_vrho1, $out_c_unpol], optimize, deducetypes=false):
"
      );

  maple2c_run(\@variables, \@derivatives, \@variants, 1);
}
