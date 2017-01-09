#!/usr/bin/env perl

my $srcdir     = $ARGV[0];
my $functional = $ARGV[1];
my $max_order  = $ARGV[2];

# Find out the type of functional
my $mathfile = "$srcdir/maple/$functional.mpl";
open my $in, '<', $mathfile or die "File $mathfile does not exist\n";

my $functype = "";
my $prefix   = "";
while($_ = <$in>){
  if(/^\(\* type:\s(\S*)\s/){
    $functype = $1;
  };
  if(/^\(\* prefix:/){
    while( ($_ = <$in>) && ! /^\*\)/ ){
      $prefix .= $_;
    }
  }
}
close($in);

my %all_vars;
my @new_c_code;

my @result_f = ("*f", "*dfdx", "*d2fdx2", "*d3fdx3", "*d4fdx4", "*d5fdx5");

# generate the derivatives to generate
my $cmd = "[";
for(my $order=0; $order <= $max_order; $order++){
  $cmd .= "result_f$order = ";
  $cmd .= ($order == 0) ? "f(x)" : "diff(f(x), x\$$order)";
  $cmd .= ", " if($order != $max_order);
}
$cmd .= "]";

# save maple file
open($mfile, ">/tmp/$$.mpl");
printf $mfile "%s", "
Digits := 20:

\$include <$functional.mpl>

with(CodeGeneration):
  C($cmd, optimize, defaulttype=numeric, precision=double, resultname=result_f);
";
close($mfile);

# run maple
$c_code = `maple -q -u /tmp/$$.mpl`;
#unlink "/tmp/$$.mpl";

# find all variables defined
my $new_c_code = "  if(order < 0) return;\n\n";
my $variables  = "";
my $n_var = 0;

for (split /^/, $c_code) {
  s/params/params->/g;

  if(/result_f/){
    s/result_f(\d+)/\n  $result_f[$1]/;
    $new_c_code .= $_."\n  if(order < $1+1) return;\n\n";
  }elsif(/(t\d+) =/){
    if($n_var % 8 == 0){
      $variables .= ";\n" if($n_var != 0);
      $variables .= "  double ";
    }else{
      $variables .= ", ";
    }
    $n_var++;

    $variables .= "$1";
    $new_c_code .= "  ".$_;
  }
}
$variables .= ";\n" if($n_var != 0);


# now we print the c file
  print "
void XC(${functional}_enhance)
  (const XC(func_type) *p, int order, 
   FLOAT x, FLOAT *f, FLOAT *dfdx, FLOAT *d2fdx2, FLOAT *d3fdx3)
{
";

print $variables, "\n", $prefix, "\n", $new_c_code;

  print "}\n

#define maple2c_order $max_order
#define maple2c_func  XC(${functional}_enhance)
\n";
