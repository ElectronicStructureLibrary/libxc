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
	$doi{$ref} =~ s/http:\/\/dx.doi.org\///;
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

%journal_abbreviations =
  (
   "Mathematical Proceedings of the Cambridge Philosophical Society" => "Math. Proc. Cambridge Philos. Soc.",
   "Zeitschrift für Physik" => "Z. Phys.",
   "Journal of Physics C: Solid State Physics" => "J. Phys. C: Solid State Phys.",
   "Canadian Journal of Physics" => "Can. J. Phys.",
   "The Journal of Chemical Physics" => "J. Chem. Phys.",
   "Rendiconti dell'Accademia Nazionale dei Lincei" => "Rend. Accad. Naz. Lincei ",
   "Journal of Chemical Theory and Computation" => "J. Chem. Theory Comput.",
   "Journal of Physics: Condensed Matter" => "J. Phys.: Condens. Matter",
   "Chemical Physics Letters" => "Chem. Phys. Lett.",
   "International Journal of Quantum Chemistry" => "Int. J. Quantum Chem.",
   "The Journal of Physical Chemistry A" => "J. Phys. Chem. A",
   "The Journal of Physical Chemistry B" => "J. Phys. Chem. B",
   "The Journal of Physical Chemistry Letters" => "J. Phys. Chem. Lett.",
   "The Journal of Physical Chemistry" => "J. Phys. Chem",
   "Molecular Physics" => "Mol. Phys.",
   "Physica Scripta" => "Phys. Scr.",
   "Journal of Computational Chemistry" => "J. Comput. Chem.",
   "Proceedings of the National Academy of Sciences of the U. S. A." => "Proc. Natl. Acad. Sci. U. S. A.",
   "Theoretical Chemistry Accounts" => "Theor. Chem. Acc.",
   "Theoretica chimica acta" => "Theor. Chim. Acta",
   "Journal of Computational Methods in Science and Engineering" => "J. Comput. Methods Sci. Eng.",
   "Journal of the Physical Society of Japan" => "J. Phys. Soc. Jpn.",
   "Physics Letters A" => "Phys. Lett. A",
   "Journal of Molecular Structure: THEOCHEM" => "J. Mol. Struct.: THEOCHEM",
   "Computational and Theoretical Chemistry" => "Comput. Theor. Chem.",
   "Journal of Photochemistry and Photobiology A: Chemistry" => "J. Photochem. Photobiol., A",
  );

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
      $item =~ s/\{\\natexlab\{.\}(.*?)\}/$1/g;

      $item =~ s/^\\bibitem\s*\[.*\]\{(.*?)\}//;
      $label = $1;

      $item =~ s/\\BibitemOpen//;
      $item =~ s/\\BibitemShut.*//;

      $item =~ s/\\bibf*namefont\s*\{(.*?)\}/$1/g;
      $item =~ s/\\bibinfo\s*\{.*?\}\s*\{(.*?)\}/$1/g;
      $item =~ s/\\bibinfo\s*\{editor \{(.*?)\}(.*?)\}/$1$2/g; # result of nested bibinfos for book with editor
      $item =~ s/\\bibfield\s*\{.*?\}\s*\{(.*?)\}/$1/g;
      $item =~ s/\\href.*?\{.*?\}\s*\{(.*?)\}/$1/g;
      $item =~ s/,\\ \\Eprint.*?\{.*?\}\s*\{(http:\/\/.*?)\}\s*//g; # wipe URL that is not arxiv
      $item =~ s/\\Eprint.*?\{.*?\}\s*\{(.*?)\}\s*/$1/g; # arxiv
      $item =~ s/\\v\{(.)(.*?)\}/\\v{$1}$2/g; # special rule for haceks \v{ } in names

      $item =~ s/\\textbf\s*\{(.*?)\}/$1/g;
      $item =~ s/\\emph\s*\{(.*?)\}/$1/g;

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

      # and now we abbreviate the journal names
      foreach $key ( keys %journal_abbreviations ){
	$item =~ s/$key/$journal_abbreviations{$key}/;
	
      }

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
