#!/usr/bin/perl

my @funcs = ("lda", "gga");

$s0 = ""; $s3 = "";
foreach $func (@funcs){
  undef %deflist_f;
  undef %deflist_c;

  read_file($func);

  $s1 = ""; $s2 = "";
  foreach $key (sort { $a <=> $b } keys %deflist_f) {
    $s0 .= sprintf "%s %-20s %3s  /*%-60s*/\n", "#define ",
      $deflist_f{$key}, $key, $deflist_c{$key};

    $s3 .= sprintf "  %s %-20s = %3s  ! %s\n", "integer, parameter ::",
      $deflist_f{$key}, $key, $deflist_c{$key};

    $t = $deflist_f{$key};
    $t =~ s/XC_(.*)/\L$1/;

    $s1 .= "extern xc_func_info_type func_info_$t;\n";
    $s2 .= "  &func_info_$t,\n";
  }

  open(OUT, ">funcs_$func.c");
  print OUT <<EOF
#include "util.h"

$s1

const xc_func_info_type *${func}_known_funct[] = {
$s2  NULL
};
EOF
    ;
  close OUT;
}

open(OUT, ">xc_funcs.h");
print OUT $s0;
close OUT;

open(OUT, ">libxc_funcs.f90");
print OUT <<EOF
module libxc_funcs
  implicit none

  public

$s3
end module libxc_funcs
EOF
  ;
close OUT;

sub read_file() {
  my $type = $_[0];
  $type =~ s/(.*)/\L$1/;

  my $TYPE = $type;
  $TYPE =~ s/(.*)/\U$1/;

  opendir(DIR, "./") || die "can’t opendir ./: $!";
  while($_ = readdir(DIR)){
    next if(!/^${type}_.*\.c$/);

    open(IN, "<$_");
    while($_=<IN>){
      if(/#define\s+(XC_${TYPE}_\S+)\s+(\S+)\s+\/\*(.*)\*\//){
	$deflist_f{$2} = $1;
	$deflist_c{$2} = $3;
      }
    }
    close(IN);
  }
  closedir DIR;
}


