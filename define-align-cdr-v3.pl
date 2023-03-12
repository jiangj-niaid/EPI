#!/usr/bin/perl
# Usage: "define-align-cdr.pl $option align-wclusterO.txt"
#
# $option: "nbs", "absH" or "absL"
# "align-wclusterO.txt" is a sequence alignment by Cluster Omega (merged as single line)
# Extract the CDR loops and place TAG to CDR loog residues
#
# Jiansheng Jiang, 10-06-21
# Jiansheng Jiang, 03-30-22, v2
# Jiansheng Jiang, 02-06-23, v3
# 
#

my ($HorL) = "";
my $option = $ARGV[0];
if ($ARGV[0] eq "abs" || $ARGV[0] eq "absH") { 
  $HorL = "H"; $option = "abs";
}
if ($ARGV[0] eq "absL") {
  $HorL = "L";  $option = "abs";
}
$file = $ARGV[1];
($ofile = $file) =~ s/(\.txt)$//; $ofile = $ofile . "-cdr.txt";

my(@NAM,@SEQ);
my($j,$nj)=0; my($na,$se);
open (I,"<$file") || die("ERROR: Could not open the file - $file ?");
while ($l = <I>)
{
  chomp($l);
#  $l =~ s/^\*//;
  next if ($l =~ /^#/);
  $j++;
  ($na,$se) = split("\t",$l);
  $NAM[$j] = $na; $SEQ[$j]=$se;
#
}
$nj = $j;
#print "nj=$nj\n";

# tagging CDR loop id
my ($ks1,$ke1,$ks2,$ke2,$ks3,$ke3) = 0;
# CDR	11111111111	2222222222	333333333333333333333333
# 07-06-22 for 267abs
# nbs	CDR1=28-39	CDR2=55-65	CDR3=103-128
# 07-06-22 for 58nbs
# nbs	CDR1=27-38	CDR2=54-63	CDR3=101-127
# k should be less -1
if ($option =~ "nbs") {
#  $ks1 = 26; $ke1 = 36; $ks2 = 52; $ke2 = 61; $ks3 = 99; $ke3 = 122;
#  $ks1 = 26; $ke1 = 37; $ks2 = 53; $ke2 = 62; $ks3 = 100; $ke3 = 126;
  $ks1 = 26; $ke1 = 37; $ks2 = 54; $ke2 = 64; $ks3 = 101; $ke3 = 127;
}
# CDR	11111111111	2222222222	333333333333333333333333
# abs   CDR1=28-38	CDR2=54-64	CDR3=102-127
# k should be less -1
if ($option =~ "abs") {
  if ($HorL =~ "H") {
#    $ks1 = 27; $ke1 = 38; $ks2 = 54; $ke2 = 64; $ks3 = 101; $ke3 = 125;
    $ks1 = 27; $ke1 = 37; $ks2 = 54; $ke2 = 64; $ks3 = 101; $ke3 = 125;
  }
  if ($HorL =~ "L") {
#    $ks1 = 27; $ke1 = 38; $ks2 = 54; $ke2 = 65; $ks3 = 103; $ke3 = 116;
    $ks1 = 26; $ke1 = 39; $ks2 = 54; $ke2 = 65; $ks3 = 103; $ke3 = 116;
  }
}
#
my(@CDR); my($cdr1,$cdr2,$cdr3); my($str,$len,$s,$i,$tag); 
my(@NI); my(@SQ,@CD,@TG);
#
for ($j=0;$j<=$nj;$j++) {
  $NI[$j] = ""; $TG[$j] = ""; $SQ[$j] = ""; $CD[$j] = "";
}
#
for ($j=1;$j<=$nj;$j++) {
  $str = $SEQ[$j]; $len = length($str) - 1;
  $i = 0;
  ($cdr1,$cdr2,$cdr3)="";
  for ($k=0;$k<=$len;$k++) {
    $s = substr($str,$k,1);
# skip off "-"
    next if ($s eq "-");
#
    $i++; $tag = "0";
# k should be less -1
    if ($k >= $ks1 and $k <=$ke1) {
      $tag = "1";
      $cdr1 = $cdr1 . $s;
    }
    if ($k >= $ks2 and $k <=$ke2) {
      $tag = "2";
      $cdr2 = $cdr2 . $s;
    }
    if ($k >= $ks3 and $k <=$ke3) {
      $tag = "3";      
      $cdr3 = $cdr3 . $s;
    }
    $SQ[$j] = $SQ[$j] . $s;
    $CD[$j] = $CD[$j] . $tag;
    $stag = $s . $i . "." . $tag;
    $TG[$j] = $TG[$j] . "," . $stag;
  }
  $CDR[$j] = $cdr1 . "," . $cdr2 . "," . $cdr3;
  $NI[$j] = length($cdr1) . "," . length($cdr2) . "," . length($cdr3);
  $TG[$j] =~ s/^\,//;
}

# print out
open (O,">$ofile") || die("ERROR: Could not open the file - $ofile ?");
# CDR loops sequences: CDR1, CDR2, CDR3
print O "#NAME      \tCDR1,CDR2,CDR3\tLENGTH\n";
my ($sum1,$sum2,$sum3) = 0; my($l1,$l2,$l3);
for ($j=1;$j<=$nj;$j++) {
  ($l1,$l2,$l3) = split(",",$NI[$j]);
  $sum1 += $l1; $sum2 += $l2; $sum3 += $l3;
  my ($nn) = sprintf("%-12s",$NAM[$j]);
  print O "$nn\t$CDR[$j]\t$NI[$j]\n";
  print "$nn\t$TG[$j]\n";
}
$l1 = sprintf("%3.1f",$sum1/$nj);
$l2 = sprintf("%3.1f",$sum2/$nj);
$l3 = sprintf("%3.1f",$sum3/$nj);
print O "#AVERAGE-Length\tCDR1,CDR2,CDR3\t$l1,$l2,$l3\n";




1; 


