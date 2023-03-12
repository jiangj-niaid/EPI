#!/usr/bin/perl
# Usage: "readnwrite-fcs.pl fcs-define.txt contact.txt"
# 
# "contact.txt" is the list of contact distance between repitope and paratope
# "fcs-define.txt" is the definition of FCS for epitope residues (may include classification)
# Extract the CDR loops and place TAG to CDR loop residues 
# The purpose of this script is to identify/match the paratope residue with CDR loop tag
#
# Jiansheng Jiang, 10-06-21
# Jiansheng Jiang, 03-10-22, v2
#
#
#
# presettings
#
# Three letter to one code for aa
my %ONE = (
	"ALA" => "A", "PHE" => "F", "ILE" => "I", "LEU" => "L", 
	"MET" => "M", "PRO" => "P", "VAL" => "V", "TRP" => "W", 
	"ASP" => "D", "GLU" => "E", "LYS" => "K", "ARG" => "R", "HIS" => "H",  
	"CYS" => "C", "GLY" => "G", "ASN" => "N", 
	"GLN" => "Q", "SER" => "S", "THR" => "T", "TYR" => "Y", 
	);


$file2 = $ARGV[0]; # fcs.txt 
$file = $ARGV[1];  # contact.txt


my(@W);
my(@NM,@PNM,@PD,@PAR); my(@FCS) = 0;
my($j,$mj,$nj,$len) = 0; my($na,$se,$st,$tg); my(@NAM2,@SEQ,@CLA); my(@SQ);
for ($k=328;$k<=540;$k++) {
  $FCS[$k] = 0;
}

open (II,"<$file2") || die("ERROR: Could not open the file - $file2 ?");
while ($l = <II>)
{
  chop($l);
#  $l =~ s/^\*//;
  next if ($l =~ /^#/);
  $j++;
  ($na,$se,$fc,$cl) = split("\t",$l);
#  print "dbx: $na\t$se\n";
  $NAM2[$j] = $na;
  $SEQ[$j] = $se;
  $FCS[$se] = $fc;
  $CLA[$j] = $cl;
#  (@SQ) = split(",",$se);
#  $len = scalar(@SQ) - 1;
#  for ($k=0;$k<=$len;$k++) {
#    ($st,$tg) = split(/\./,$SQ[$k]);
#    $k1 = $k + 1;
#    $CDR[$j][$k1] = $tg;
#  }
}
$mj = $j;
close(II);

# print FCS array
my($k) = 0; my($s) = "";
for ($j=1;$j<=$mj;$j++) {
  $k = $k + 1;
  $se = $SEQ[$j]; 
  $s .= $FCS[$se] . ",";
  if ($k >= 40) {
#    print "$s\n";
    $k = 0; $s = "";
  }
}
#print "$s\n";


#PDB_NAME	CHAINS	RBD	Rid	FCS	AB	Rid	DIST	CDR 	(9 data iterms)
print "#PDB_NAME	CHAINS	RBD	Rid	FCS	AB	Rid	DIST	CDR\n";

open (I,"<$file") || die("ERROR: Could not open the file - $file ?");
$n = 0;
while ($l = <I>)
{
  chop($l);
  $n++;
#  $l =~ s/^\*//;
  next if ($l =~ /^#/);
  @W = split("\t",$l);
  $e = $W[3];
  $fc = $FCS[$e];
  for ($j=0;$j<=3;$j++) {
    print "$W[$j]\t";
  }
  print "$fc\t";
#  for ($j=4;$j<=6;$j++) {
  for ($j=4;$j<=7;$j++) {
    print "$W[$j]\t";
  }
#  print "$W[7]\n";  
  print "$W[8]\n";  
}
$nn = $n;

close(I);


1;
