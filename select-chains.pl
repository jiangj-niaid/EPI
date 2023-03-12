#!/usr/bin/perl
# file: select-chain.pl contact A B C
# select-chains.pl for editing 'contact.inp'
# "contact.inp" is a CNS template input file 
# default two chain ids are "A" and "B"
# default cutoff distance "C" is 3.8
# replace with two new chain ids: $A and $B; and cutoff=$C
# i.e "select-chains.pl contact A B C"
#
# Jiansheng Jiang, 10-06-21

$f = shift;
$A = shift; if (!$A) {$A = "A";}
$B = shift; if (!$B) {$B = "B";}
$C = shift; if (!$C) {$C = 3.8}
$in = "../Scripts/" . $f . ".inp";
$ot = $f . "-" . $A . $B . ".inp";
$oldA = "segid AAA";
$oldB = "segid BBB";
$oldC = "cutoff=CCC";
$newA = "segid " . $A;
$newB = "segid " . $B;
$newC = "cutoff=" . $C;

open (O,">$ot") || die("cannot open the file $ot");
open (I,"<$in") || die("cannot open the file $in");
while (<I>) {
  chop($_);
  $_ =~ s/$oldA/$newA/g;
  $_ =~ s/$oldB/$newB/g;
  $_ =~ s/$oldC/$newC/g;
  print O "$_\n";
}
close(I);
close(O);


1;


