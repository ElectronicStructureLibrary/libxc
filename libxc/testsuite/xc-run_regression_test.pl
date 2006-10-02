#!/usr/bin/perl
#
# $Id: oct-run_regression_test.pl 2423 2006-09-24 21:25:52Z acastro $

use Getopt::Std;

getopts("hf:");
$opt_h && usage();
$opt_f || usage();

# Handle options
$opt_f =~ s/(.*)/\L$1\E/;

my $data_dir = ".";
my $tmp_file =  "/tmp/xc.tmp.$$";
my $exec_cmd = "./xc-get_data";

# start by reading xc.h to get a list of the defined constants
my %constants;
read_xc_h(\%constants);
$constants{"$opt_f"} || die "Functional '$opt_f' not found";

# check if we have a data file
my $data_file = "$data_dir/$opt_f.data";
(-f $data_file && -r $data_file) || die "Could not read data file '$data_file'";
open DATA, "<$data_file";
my %data;
while(data_read(*DATA, \%data) != 0){
  my $pol;
  $pol = ($data{"rhoa"}    == $data{"rhob"}    &&
	  $data{"sigmaaa"} == $data{"sigmabb"} &&
	  $data{"sigmaab"} == $data{"sigmabb"}) ? 1 : 2;

  for(;$pol<=2; $pol++){
    $cmd1  = "$exec_cmd ".$constants{"$opt_f"};
    $cmd2  = " ".$data{"rhoa"}." ".$data{"rhob"};
    $cmd2 .= " ".$data{"sigmaaa"}." ".$data{"sigmaab"}." ".$data{"sigmabb"};
    `$cmd1 $pol $cmd2 >$tmp_file`;

    open DATA2, "<$tmp_file";
    my %data2;
    data_read(*DATA2, \%data2) || die "Could not read data file '$tmp_file'";
    close DATA2;

    if($pol == 1){
      my $tmp = $data{"vsigmaaa"};
      $data{"vsigmaaa"}  = ($data{"vsigmaaa"} + $data{"vsigmaab"} + $data{"vsigmabb"})/4.0;

      my @cmp = ("zk", "vrhoa", "vsigmaaa");
      cmp_data(\%data, \%data2, \@cmp);

      $data{"vsigmaaa"} = $tmp;
    }else{
      my @cmp = ("zk", "vrhoa", "vrhob", "vsigmaaa", "vsigmaab", "vsigmabb");
      cmp_data(\%data, \%data2, \@cmp);
    }
  }
}
close DATA;
unlink $tmp_file;

###########################################
sub usage {
  print <<EndOfUsage;

 Copyright (C) 2006 by M.A.L. Marques

Usage: $0 [options] -f functional

    -h        this usage
    -f        functional to test

Report bugs to <marques\@tddft.org>.
EndOfUsage
  exit 0;
}


###########################################
sub read_xc_h {
  my $c = shift;

  open FILE, "<../src/xc.h";
  while($_ = <FILE>){
    if(/^#define +(\S*) +(\S*)/){
      my $name = $1;
      my $value = $2;

      $name =~ s/^XC_(.*)/\L$1\E/;
      $$c{$name} = $value;
    }
  }
  close FILE;
}

###########################################
sub data_read {
  my ($FILE, $data) = @_;

  while( ($line = <$FILE>) && !($line =~ /rhoa/) ){}
  $line || return 0;

  $line =~ / rhoa= (\S*) rhob= (\S*) sigmaaa= (\S*) sigmaab= (\S*) sigmabb= (\S*)/;
  $$data{"rhoa"} = $1;
  $$data{"rhob"} = $2;
  $$data{"sigmaaa"} = $3;
  $$data{"sigmaab"} = $4;
  $$data{"sigmabb"} = $5;

  my $n = 0;
  while($n++ < 24){
    $line = <$FILE> || return 0;
    $line =~ /\s*(\S*)\s*=\s*(\S*)/;
    $$data{$1} = $2;
  }
  return 1;
}

sub cmp_data {
  my ($d1, $d2, $what) = @_;
  my $tol = 1e-10;

  foreach $var (@$what){
    $ok = (abs($$d1{$var}) < 1e-15 && abs($$d2{$var}) < 1e-15);
    if(!$ok){
      $ok = (abs($$d1{$var} - $$d2{$var}) <= $tol*abs($$d1{$var}));
    }
    if(!$ok){
      print "$var mismatch: ", $$d1{$var}, " != ", $$d2{$var}, "\n";
    }
  }
}
