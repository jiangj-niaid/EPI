#!/usr/bin/perl
# Usage: "read-cdr-tags-v2.pl nbs|abs align-defined.txt contact.txt"
# 
# "contact.txt" is the list of contact distance between repitope and paratope
# "align.txt" is the list of sequences which aligned with CDR loops by Cluster Omega
# Extract the CDR loops and place TAG to CDR loog residues 
# The purpose of this script is to identify/match the paratope residue with CDR loop tag
#
# If cdr-loop is not defined in "align-defined" file,
# the cdr-tags is asigned just according to the alignment range:
# nbs	CDR1=27-37	CDR2=53-62	CDR3=100-123
# abs   CDR1=28-38	CDR2=54-64	CDR3=102-127
# 
#
# Jiansheng Jiang, 10-06-21
# Jiansheng Jiang, 03-30-22, v2
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

$option = $ARGV[0];
$align = $ARGV[1]; # align.txt 
$file = $ARGV[2];  # contact.txt


my(@W);
my(@NM,@PNM,@PD,@PAR);
my($j,$mj,$nj,$len) = 0; my($na,$se,$st,$tg); my(@NAM2,@SEQ); my(@SQ);


open (II,"<$align") || die("ERROR: Could not open the file - $align ?");
while ($l = <II>)
{
  chomp($l);
#  $l =~ s/^\*//;
  next if ($l =~ /^#/);
  $j++;
  ($na,$se) = split("\t",$l);
#  print "dbx: $na\t$se\n";
  $NAM2[$j] = $na;
  $SEQ[$j] = $se;
}
$mj = $j;
close(II);



# tagging CDR loop id
my ($ks1,$ke1,$ks2,$ke2,$ks3,$ke3) = 0;
# CDR	11111111111	2222222222	333333333333333333333333
# 07-06-22 for 267abs
# nbs	CDR1=28-39	CDR2=55-65	CDR3=103-128
# 07-06-22 for 58nbs
# nbs	CDR1=27-38	CDR2=54-63	CDR3=101-127
# k should be less -1
if ($option =~ "nbs") {
#  $ks1 = 26; $ke1 = 37; $ks2 = 53; $ke2 = 62; $ks3 = 100; $ke3 = 126;
  $ks1 = 27; $ke1 = 38; $ks2 = 54; $ke2 = 63; $ks3 = 101; $ke3 = 126;
}
# CDR	11111111111	2222222222	333333333333333333333333
# abs   CDR1=28-38	CDR2=54-64	CDR3=102-127
# k should be less -1
if ($option =~ "abs") {
#  $ks1 = 27; $ke1 = 38; $ks2 = 54; $ke2 = 64; $ks3 = 101; $ke3 = 125;
  $ks1 = 27; $ke1 = 38; $ks2 = 54; $ke2 = 64; $ks3 = 101; $ke3 = 127;
}
#

#PDB_NAME	CHAINS	RBD	Rid	AB	Rid	DIST	CDR	(8 data terms)
print "#PDB_NAME	CHAINS	RBD	Rid	AB	Rid	DIST	CDR\n";

my(@SST,@TAG); my($pdb,$nam,$resid); my($len,$last,$st,$tg);
open (I,"<$file") || die("ERROR: Could not open the file - $file ?");
my $n = 0;
while ($l = <I>)
{
  chomp($l);
  $n++;
#  $l =~ s/^\*//;
  next if ($l =~ /^#/);
  @W = split("\t",$l);
  ($pdb,$nam) = split(/\_/,$W[0]);
  $p = $ONE{$W[4]} . $W[5];
  $resid = $W[5]; $len = length($W[5]); $last = substr($W[5],$len-1,1);
  if ($last =~ /A|B|C|D/) {
    $r = substr($W[5],0,$len-1); $resid = int($r);
  }
  $jj = 0;
  for ($j=1;$j<=$mj;$j++) {
    $na = $NAM2[$j];
    if ($na =~ $nam) {
      $jj = $j;
    }
  }
#
#
  (@SQ) = split(",",$SEQ[$jj]);
  $len = scalar(@SQ) - 1;
  $SST[$n] = $p; $TAG[$n] = 0;
  for ($k=0;$k<=$len;$k++) {
    ($st,$tg) = split(/\./,$SQ[$k]);
# if match, use the predefined cdr-tag
    if ($p =~ $st) {
      $TAG[$n] = $tg;
      $SST[$n] = $st;
#      if ($nam =~ /H11-H4m/) {
#        print "dbx: k=$k\tp=$p\tst=$st\ttg=$tg\tn=$n\tTAG=$TAG[$n]\n";
#      }
# otherwise check the range
    } elsif ($resid >= $ks1 && $resid <= $ke1) {
      $TAG[$n] = 1;
      $SST[$n] = $p;
    } elsif ($resid >= $ks2 && $resid <= $ke2) {
      $TAG[$n] = 2;
      $SST[$n] = $p;
    } elsif ($resid >= $ks3 && $resid <= $ke3) {
      $TAG[$n] = 3;
      $SST[$n] = $p;
    } else {
      $TAG[$n] = 0;
      $SST[$n] = $p;
    }
  }
  print "$W[0]\t$W[1]\t$W[2]\t$W[3]\t$W[4]\t$W[5]\t$W[6]\t$TAG[$n]\n";
}
$nn = $n;

close(I);



exit;
