#!/usr/bin/perl
# Usage: "pickup-outliers.pl pick-list.txt"
# "pick-list.txt" is a list of names to be pickup.
# 
# Jiansheng Jiang, 05-19-22
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
#
# RBD sequence 332-528
my $rbdseq  = "ITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVS";
   $rbdseq .= "PTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCV";
   $rbdseq .= "IAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGV";
   $rbdseq .= "EGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPK";

my (@RBD) = "";
my ($q0) = 332; 
for ($i=0;$i<=196;$i++) {
  $q = $i + $q0;
  $s = substr($rbdseq,$i,1);
  $sq = $s . sprintf("%3s",$q);
  $RBD[$q] = $sq;
#  print "$q\t$s\t$RBD[$q]\n";
}
#

$file = $ARGV[0];

my (@EPI) = 0;

open (I,"<$file") || die("ERROR: Could not open the file - $file ?");

my(@W);
my($i,$ni,$j,$nj) = 0;
my(@NM,@PNM,@PD,@CA,@CB,@EP,@NEP,@PAR,@DIS,@PT,@CDR) = 0;
my(@NBK) = 0;
my(@PAT) = "";
my(@EPI) = "";
my(@RB) = "";

while ($l = <I>)
{
  chomp($l);
  $j++;
#  $l =~ s/^\*//;
#  print "$l\n";
  @W = split("\t",$l);
#  print "@W\n";
  $PNM[$j] = $W[0];
  ($PD[$j],$na) = split("_",$W[0]);
  print "$PNM[$j]\t$na\n";
  system('/usr/bin/perl ../../Scripts/epitope-statis+fcs+cdr-abs-v5.pl 217abs-contact-dis-H+fcs+cdr-0405.txt 5 $na');
}
$nj = $j;

close(I);

exit;



