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

  my $maple_code = "
# zk is energy per unit particle
mzk  := (r0, r1) -> f(r_ws(dens(r0, r1)), zeta(r0, r1)):

(* mf is energy per unit volume *)
mf   := (r0, r1) -> simplify(dens(r0, r1)*mzk(r0, r1), symbolic):
";
  
  math2c_run(\@variables, \@derivatives, \@variants, $maple_code);
}
