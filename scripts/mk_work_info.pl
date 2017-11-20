#!/usr/bin/env perl

# Copyright (C) 2017 M.A.L. Marques
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# this provides a first approximation to the info structure
# we need to define the derivatives. This structure then has
# to be edited by hand.

use Data::Dumper;

sub mk_info {
  ($vars, $f, $order) = @_;

  $info = [];
  return $info if($order < 0);

  @{$info}[0] = [["r_a_f", "$f"]];
  return $info if($order < 1);

  @{$info}[1] = [];
  for my $i (0 .. $#{$vars}){
    my $v1 = ${$vars}[$i];
    push(@{$info}[1], ["r_a_dfd${v1}", "diff($f, r_a_$v1)"]);
  }
  return $info if($order < 2);

  @{$info}[2] = [];
  for my $i (0 .. $#{$vars}){
    my $v1 = ${$vars}[$i];
    push(@{$info}[2], ["r_a_d2fd${v1}2", "diff($f, r_a_$v1\$2)"]);
    for my $j ($i+1 .. $#{$vars}){
      my $v2 = ${$vars}[$j];
      push(@{$info}[2], ["r_a_d2fd${v1}${v2}", "diff($f, r_a_$v1, r_a_$v2)"]);
    }
  }
  return $info if($order < 3);
  
  @{$info}[3] = [];
  for my $i (0 .. $#{$vars}){
    my $v1 = ${$vars}[$i];
    push(@{$info}[3], ["r_a_d3fd${v1}3", "diff($f, r_a_$v1\$3)"]);
    for my $j ($i+1 .. $#{$vars}){
      my $v2 = ${$vars}[$j];
      push(@{$info}[3], ["r_a_d3fd${v1}2${v2}", "diff($f, r_a_$v1\$2, r_a_$v2)"]);
      push(@{$info}[3], ["r_a_d3fd${v1}${v2}2", "diff($f, r_a_$v1, r_a_$v2\$2)"]);
      for my $k ($j+1 .. $#{$vars}){
        my $v3 = ${$vars}[$k];
        push(@{$info}[3], ["r_a_d3fd${v1}${v2}${v3}", "diff($f, r_a_$v1, r_a_$v2, r_a_$v3)"]);
      }
    }
  }
  return $info if($order < 4);
}

# for work_lda_0
#my $f    = "f(r_a_rs, 0.0)";
#my $vars = ["rs"];

# for work_lda_1
#my $f    = "f(r_a_rs, r_a_z)";
#my $vars = ["rs", "z"];

# for work_gga_x
#my $f    = "f(r_a_x)";
#my $vars = ["x"];

# for work_gga_c
#my $f    = "f(r_a_rs, r_a_z, r_a_xt, r_a_xs_0_, r_a_xs_1_)";
#my $vars = ["rs", "z", "xt", "xs_0_", "xs_1_"];

# for work_mgga_c
my $f    = "f(r_a_rs, r_a_z, r_a_xt, r_a_xs_0_, r_a_xs_1_, r_a_ts_0, r_a_ts_1, r_a_us_0_, r_a_us_1_)";
my $vars = ["rs", "z", "xt", "xs_0_", "xs_1_", "ts_0_", "ts_1_", "us_0_", "us_1_"];

my $info =  mk_info($vars, $f, 3);

$Data::Dumper::Indent = 0;
print Dumper($info);
