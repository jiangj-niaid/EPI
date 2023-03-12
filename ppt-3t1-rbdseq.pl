#!/usr/bin/perl
# file: three-to-one-code+property.pl
# Jiansheng Jiang @NIH, 10/05/2021
#
# Map property of aa
# Hydrophobic: H=(A,F,I,L,M,P,V,W)
# Charge:      C=(D,E,K,R,H)
# Polar:       P=(C,G,N,Q,S,T,Y)
my %ppt = (
	"A" => "H", "F" => "H", "I" => "H", "L" => "H", 
	"M" => "H", "P" => "H", "V" => "H", "W" => "H", 
	"D" => "C", "E" => "C", "K" => "C", "R" => "C", "H" => "C",  
	"C" => "P", "G" => "P", "N" => "P", 
	"Q" => "P", "S" => "P", "T" => "P", "Y" => "P", 
	);
#
# Three letter to one code for aa
my %one = (
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

my (@EPI) = 0;



1;

