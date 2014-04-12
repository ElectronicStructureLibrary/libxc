#!/usr/bin/env perl

# Copyright (C) 2006-2007 M.A.L. Marques
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#  
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

use Cwd;

# initialize stuff
@reference = ();
%bibtex    = {};
%doi       = {};
%bibitem   = {};

# first we read the libxc.bib file and get the references
open(BIB, "<../libxc.bib");
while($_=<BIB>){
  if($_ =~ /@.*\{(.*),/){
    $ref = $1;
    push(@reference, $ref);
    $bibtex{$ref} = $_;

    while(($_=<BIB>) && !($_=~/^\s*\}/)){
      $bibtex{$ref} .= $_;

      if($_ =~ /doi\s*=\s*[{\"]([^}\"]*)[}\"]*/){
	$doi{$ref} = $1;
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

open(TEX, ">/tmp/$$.tex");
print TEX "\\documentclass[prl]{revtex4-1}
\\usepackage[utf8]{inputenc}
\\begin{document}
\\nocite{*}
\\bibliography{$cwd/../libxc}
\\end{document}
";
close(TEX);

# run latex and bibtex
`(cd /tmp && latex $$.tex && bibtex $$.aux && latex $$.tex && latex $$.tex && /bin/mv -f $$.dvi libxc.dvi)`;

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
   "Proceedings of the National Academy of Sciences of the United States of America" => "Proc. Natl. Acad. Sci. U. S. A.",
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
open(BBL, "</tmp/$$.bbl");
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
      $item =~ s/\\bibfield\s*\{.*?\}\s*\{(.*?)\}/$1/g;
      $item =~ s/\\href.*?\{.*?\}\s*\{(.*?)\}/$1/g;

      $item =~ s/\\textbf\s*\{(.*?)\}/$1/g;
      $item =~ s/\\emph\s*\{(.*?)\}/$1/g;

      $item =~ s/\\ / /g;
      $item =~ s/~/ /g;
      $item =~ s/^\s*//;
      $item =~ s/\s+/ /g;

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
`(cd /tmp && /bin/rm -f $$.tex $$.aux $$.bbl $$.blg $$.log)`;

# now we make a nice output
open(OUT_C, ">references.c");
open(OUT_H, ">references.h");
$str = "/*
   This file is autogenerated. Please do not change it directly. Run instead
   ./get_references.pl
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
