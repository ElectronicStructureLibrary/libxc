#!/usr/bin/env perl

# Copyright (C) 2006-2007 M.A.L. Marques
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

if(@ARGV < 1) {
    print STDERR "Usage: get_references.pl BIBFILE\n";
    exit(1);
}

$bibfile   = shift;

use Cwd;

# initialize stuff
@reference = ();
%bibtex    = {};
%doi       = {};
%bibitem   = {};

# first we read the libxc.bib file and get the references
open(BIB, "<", $bibfile);
while($_=<BIB>){
  if($_ =~ /@.*\{(.*),/){
    $ref = $1;
    push(@reference, $ref);
    $bibtex{$ref} = $_;

    while(($_=<BIB>) && !($_=~/^\s*\}/)){
      $bibtex{$ref} .= $_;

      if($_ =~ /doi\s*=\s*[{\"]([^}\"]*)[}\"]*/){
        $doi{$ref} = $1;
        $doi{$ref} =~ s/http:\/\/doi.org\///;
      }
    }
    $bibtex{$ref} .= "}";
    $bibtex{$ref} =~ s/\\/\\\\/g; # protect backslashes
    $bibtex{$ref} =~ s/\"/\\\"/g; # protect quotation marks
    $bibtex{$ref} =~ s/\n/\\n/g; # protect new lines
  }
}
close(BIB);

# now we have to make the reference. For that we will run latex
my $cwd = getcwd();

$dir = "/tmp";

open(TEX, ">$dir/$$.tex");
print TEX "\\documentclass[aps,prl]{revtex4-1}
\\usepackage[utf8]{inputenc}
\\begin{document}
\\nocite{*}
\\bibliography{$cwd/../libxc}
\\end{document}
";
close(TEX);

# run latex and bibtex
system "cd $dir && latex $$.tex && bibtex $$.aux && latex $$.tex && latex $$.tex && /bin/mv -f $$.dvi libxc.dvi";

# now we parse the bbl file
open(BBL, "<$dir/$$.bbl");
$item = "";
while($_=<BBL>){
  if($item ne ""){
    if($_ =~ /^\s\s/){
      $_ =~ s/\{NoStop\}%//;
      $item .= $_;
    }else{
      # we now clean and parse the bibitem
      $item =~ s/\n//gm;
      $item =~ s/%//gm;

      # remove accents
      $item =~ s/\\o({})?/o/g;
      $item =~ s/\\l({})?/l/g;
      $item =~ s/\\["'`~ov]{(.?)}/$1/g;
      $item =~ s/\\["'`~]//g;

      $item =~ s/\{\\natexlab\{.\}(.*?)\}/$1/g;

      $item =~ s/^\\bibitem\s*\[.*\]\{(.*?)\}//;
      $label = $1;

      $item =~ s/\\BibitemOpen//;
      $item =~ s/\\BibitemShut.*//;

      $item =~ s/\\bibf*namefont\s*\{(.*?)\}/$1/g;
      $item =~ s/\\bibinfo\s*\{.*?\}\s*\{(.*?)\}/$1/g;
      $item =~ s/\\bibinfo\s*\{editor\s+\{(.*?)\}(.*?)\}/$1$2/g; # result of nested bibinfos for book with editor
      $item =~ s/\\bibfield\s*\{.*?\}\s*\{(.*?)\}/$1/g;

      $item =~ s/\\textbf\s*\{(.*?)\}/$1/g;
      $item =~ s/\\emph\s*\{(.*?)\}/$1/g;
      $item =~ s/\\enquote\s*\{(.*?)\}/$1/g;

      $item =~ s/\\href.*?\{.*?\}\s*\{(.*?)\}/$1/g;
      $item =~ s/,\\ \\Eprint.*?\{.*?\}\s*\{(http:\/\/.*?)\}\s*//g; # wipe URL that is not arxiv
      $item =~ s/\\Eprint.*?\{.*?\}\s*\{(.*?)\}\s*/$1/g; # arxiv
      $item =~ s/\\v\{(.)(.*?)\}/\\v{$1}$2/g; # special rule for haceks \v{ } in names

      $item =~ s/\\ / /g;
      $item =~ s/~/ /g;
      $item =~ s/^\s*//;
      $item =~ s/\s+/ /g;

      # double up remaining backslashes so they print and don't try to escape the next character
      $item =~ s/\\/\\\\/g;

      # check if things seem ok
      if($item =~ /\\/) {
        print STDERR "WARNING: backslashes remain.\n";
        print STDERR $item . "\n";
      }
      if($item =~ /\{/ || $item =~ /\}/) {
        print STDERR "WARNING: braces remain.\n";
        print STDERR $item . "\n";
      }

      # remove remaining braces and backslashes
      $item =~ s/{(.*?)}/$1/g;
      $item =~ s/\\//g;

      $bibitem{$label} = $item;
      $item = "";
    }
  }

  if($_ =~ /^\\bibitem/){
    $item = $_;
  }

  while($item =~ /\%/){
    $item =~ s/\%//g;
    $item .= <BBL>;
  }
};
close(BBL);

# delete garbage
system "cd $dir && /bin/rm -f $$.tex $$.aux $$.bbl $$.blg $$.log $$Notes.bib $$.pol $$.pol_ref $$.unpol $$.unpol_ref";

# now we make a nice output
open(OUT_C, ">references.c");
open(OUT_H, ">references.h");
$str = "/*
   This file is autogenerated. Please do not change it directly. Run instead
   ./get_references.pl
*/

/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <xc.h>

";

print OUT_C $str;
print OUT_H $str;

foreach $r (@reference){
  $r2 = $r;
  $r2 =~ s/[^a-zA-Z0-9]/_/g;

  print OUT_C "func_reference_type xc_ref_$r2 = {
  \"$bibitem{$r}\",
  \"$doi{$r}\",
  \"$bibtex{$r}\"
};
\n\n";

  print OUT_H "extern func_reference_type xc_ref_$r2;\n";
}

close(OUT_C);
close(OUT_H);
