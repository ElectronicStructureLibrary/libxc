#!/usr/bin/env perl

# Copyright (C) 2006-2007 M.A.L. Marques
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

if(@ARGV < 2) {
    print STDERR "Usage: get_funcs.pl srcdir builddir\n";
    exit(1);
}

$srcdir = shift;
$builddir = shift;

my @funcs = ("lda", "gga", "hyb_gga", "mgga", "hyb_mgga");
my %all_ids;

open(DOCS, ">$builddir/libxc_docs.txt") or die("Could not open '$builddir/libxc_docs.txt.'\n");

$s0 = ""; $s3 = ""; $s4 = ""; $s5 = "";
$publiclist = ""; $xclist = ""; $fxclist = ""; $xcf90list = ""; $xcfclist = "";

foreach $func (@funcs){
  undef %deflist_f;
  undef %deflist_c;
  undef %num;

  read_file($srcdir, $func);

  $s1 = ""; $s2 = "";
  foreach $key (sort { $a <=> $b } keys %deflist_f) {
    $s0 .= sprintf "%s %-30s %3s  /*%-70s*/\n", "#define ",
      $deflist_f{$key}, $key, $deflist_c{$key};

    $t = $deflist_f{$key};
    $t =~ s/XC_(.*)/\L$1/;

    $s4 .= ",\n" if($s4);
    $s4 .= sprintf "{\"%s\", %d}", $t, $key;

    $s3 .= sprintf "  %s %-30s = %3s  ! %s\n", "integer, parameter ::",
      $deflist_f{$key}, $key, $deflist_c{$key};

    $s5 .= sprintf "  %s %-30s = %3s  ! %s\n", "integer(c_int), parameter, public ::",
      $deflist_f{$key}, $key, $deflist_c{$key};

    $s1 .= "extern xc_func_info_type xc_func_info_$t;\n";
    $s2 .= "  &xc_func_info_$t,\n";
  }

  open(OUT, ">$builddir/funcs_$func.c") or die("Could not open '$builddir/funcs_$func.c'.\n");
  print OUT <<EOF
#include "util.h"

$s1

const xc_func_info_type *xc_${func}_known_funct[] = {
$s2  NULL
};
EOF
    ;
  close OUT;
}

close DOCS;

open(OUT, ">$builddir/funcs_key.c") or die("Could not open '$builddir/funcs_key.c'.\n");
print OUT <<EOF
#include "util.h"

xc_functional_key_t xc_functional_keys[] = {
$s4,
{"", -1}
};
EOF
;

open(OUT, ">$builddir/xc_funcs.h") or die("Could not open '$builddir/xc_funcs.h'.\n");
print OUT $s0;
print $so;
close OUT;

open(OUT, ">$builddir/libxc_funcs.f90") or die("Could not open '$builddir/libxc_funcs.f90'.\n");
print OUT <<EOF
!! Copyright (C) 2003-2015 Miguel Marques
!! All rights reserved.
!!
!! This file is dual-licensed under a GPL and a BSD license
!!
!! MPL License:
!!
!! This Source Code Form is subject to the terms of the Mozilla Public
!! License, v. 2.0. If a copy of the MPL was not distributed with this
!! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!!
!! BSD License:
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!! 1. Redistributions of source code must retain the above copyright
!! notice, this list of conditions and the following disclaimer.
!!
!! 2. Redistributions in binary form must reproduce the above
!! copyright notice, this list of conditions and the following
!! disclaimer in the documentation and/or other materials provided
!! with the distribution.
!!
!! 3. Neither the name of the copyright holder nor the names of its
!! contributors may be used to endorse or promote products derived
!! from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
!! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
!! HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
!! STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
!! OF THE POSSIBILITY OF SUCH DAMAGE.

module libxc_funcs_m
  implicit none

  public

$s3
end module libxc_funcs_m
EOF
  ;
close OUT;


open(OUT, ">$builddir/libxc_inc.f03") or die("Could not open '$builddir/libxc_incs.f03'.\n");
print OUT <<EOF
$s5
EOF
  ;
close OUT;

sub read_file() {
  my ($dir, $type) = @_;
  $type =~ s/(.*)/\L$1/;

  my $TYPE = $type;
  $TYPE =~ s/(.*)/\U$1/;

  # we remove the hyb from the filenames
  $save_type = $type;
  $type =~ s/^hyb_//;

  $xc_info_exe = "$builddir/xc-info";

  opendir(DIR, "$dir/") || die "cannot opendir '$dir': $!";
  while($_ = readdir(DIR)){
    next if(!/^${type}_.*\.c$/ && !/^hyb_${type}_.*\.c$/ );

    $file=$_;
    open(IN, "<$dir/$_") or die("Could not open '$dir/$_'.\n");
    while($_=<IN>){
      if(/#define\s+(XC_${TYPE}_\S+)\s+(\S+)\s+\/\*(.*)\*\//){
	$deflist_f{$2} = $1;
	$deflist_c{$2} = $3;
	$num{$1} = $2;

        # check if ID is already in use
        if ( $all_ids{$2} ){
          printf stderr "Error: ID $2 repeated in\n  $1\n  $all_ids{$2}\n";
          exit 1;
        }else{
          $all_ids{$2} = $1;
        }
      }

      if(/^(const |)xc_func_info_type xc_func_info_${save_type}/){
	  $infostr = "";
	  while($_=<IN>){
	      if(/([^}])*};/) { 
		  $infostr .= $1;
		  last;
	      }
	      # remove C comments
	      $_ =~ s|/\*[^\*]*\*/||;
	      # remove braces
	      $_ =~ s/[{}]//g;
	      chomp($_);
	      $infostr .= $_;
	  }
	  @infos = split('"', $infostr);
	  @infos0 = split(',', $infos[0]);
	  $infos0[0] =~ s/^\s*//;
	  $infos0[1] =~ s/^\s*//;
	  @infos2 = split(',', $infos[2]);

	  for($ii = 0; $ii <= $#infos2; $ii++) {
	      # remove leading spaces
	      $infos2[$ii] =~ s/^\s*//;
	  }

	  print DOCS "Number         : $num{$infos0[0]}\n";
	  print DOCS "File           : $file\n";
	  print DOCS "Codename       : $infos0[0]\n";
	  print DOCS "Kind           : $infos0[1]\n";
	  print DOCS "Description 1  : $infos[1]\n";
	  $deflist_c{$num{$infos0[0]}} =~ s/^\s*//;
	  $deflist_c{$num{$infos0[0]}} =~ s/\s*$//;
	  if($deflist_c{$num{$infos0[0]}} ne $infos[1]) {
	      print DOCS "Description 2  : $deflist_c{$num{$infos0[0]}}\n";
	  }
	  #infos2[0] will be blank
	  print DOCS "Family         : $infos2[1]\n";

	  if(-e "$xc_info_exe" && -x "$xc_info_exe") {
	      $xc_info = `$xc_info_exe $num{$infos0[0]}`;
	      @refs = split('\n', $xc_info);
	      if($refs[4] =~ /Reference\(s\)/) {
		  print DOCS "References     : ";
		  print DOCS $refs[5] . "\n";
		  $ref_start = 6;
	      } else {
		  print DOCS $refs[4] . "\n";
		  print DOCS "References     : ";
		  print DOCS $refs[7] . "\n";
		  $ref_start = 8;
	      }
	      for($ii = $ref_start; $ii <= $#refs; $ii++) {
		  print DOCS "                 " . $refs[$ii] . "\n";
	      }
	  } else {
	      # print only the names of the variables in references.c
	      print DOCS "References     :";
	      for($ii = 2; $ii <= 6; $ii++) {
		  if($infos2[$ii] ne "NULL") {
		      $infos2[$ii] =~ s/&xc_ref_//;
		      print DOCS " $infos2[$ii]";
		  }
	      }
	      print DOCS "\n";
	  }

	  if(($infos2[7] =~ /XC_FLAGS_(.)D/) != 1) {
            print STDERR $infos2[7], "\n";
	      print STDERR "$infos0[0]: Must set exactly one dimensionality flag.\n";
	      exit(1);
	  }
          print DOCS "Dimensionality : $1\n";
	   
	  print DOCS "Quantities     : ";
	  @quantities = ($infos2[7] =~ /XC_FLAGS_HAVE_(.XC)/g);
	  print DOCS join(" ", @quantities) . "\n";

	  $infos2[7] =~ s/XC_FLAGS_.D//;
	  $infos2[7] =~ s/XC_FLAGS_HAVE_.XC//g;
	  $infos2[7] =~ s/\|//g;
	  $infos2[7] =~ s/^\s*//;
	  $infos2[7] =~ s/^s*$//;
	  
	  print DOCS "Other flags    : $infos2[7]\n";

          $shortname = lc(substr($infos0[0], 3));
          $set_params = `grep "xc_${shortname}_set_params(xc_func_type" $srcdir/$file`;
	  chomp $set_params;
	  if($set_params ne "") {
	      if($set_params !~ /void/) {
		  $set_params = "void $set_params";
	      }
	      print DOCS $set_params . "\n";
	  }
	  
	  print DOCS "min dens       : $infos2[8]\n";
	  print DOCS "min grad       : $infos2[9]\n";
	  print DOCS "min tau        : $infos2[10]\n";
	  print DOCS "min zeta       : $infos2[11]\n";
	  print DOCS "init           : $infos2[12]\n";
	  # apparently, the end is always NULL
	  print DOCS "end            : $infos2[13]\n";
	  print DOCS "work lda       : $infos2[14]\n";
	  print DOCS "work gga       : $infos2[15]\n";
	  print DOCS "work mgga      : $infos2[16]\n";
	  print DOCS "----------------------------\n";

#	  print "$file $infos0[0] $infos2[12]\n";

	  if($num{$infos0[0]} eq "") {
	      print STDERR "ERROR: missing number\n";
              print STDERR $infos0[0], "\n";
	      exit(1);
	  }

	  if($deflist_f{$num{$infos0[0]}} ne $infos0[0]) {
	      print STDERR $deflist_f{$num{$infos0[0]}} . " " . $infos0[0] . "\n";
	      print STDERR "Mismatch of names.\n";
	      exit(1);
	  }
      }
    }
    close(IN);
  }
  closedir DIR;
}
