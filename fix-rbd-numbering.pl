#!/usr/bin/perl
# Usage: "fix-rbd-numbering.pl contact-dis.txt"
# RBD residue numbering + 3
# Beta-26,Beta-32,Beta-49,Beta-50
# 
# Jiansheng Jiang, 05-19-22
#

$file = $ARGV[0];
#$selection = $ARGV[1];

my (@PKLS) = ("Beta-26","Beta-32","Beta-49","Beta-50");
my ($npk) = 3;


open (I,"<$file") || die("ERROR: Could not open the file - $file ?");

my(@W);
my($k,$ni,$j,$nj) = 0;
my(@NM,@PNM,@PD,@CA,@CB,@EP,@NEP,@PAR,@DIS,@PT,@CDR) = 0;

while ($l = <I>)
{
  chomp($l);
  $j++;
#  $l =~ s/^\*//;
  @W = split("\t",$l);
#  print "@W\n";
  $NEP[$j] = $W[3];
  ($PD[$j],$na) = split("_",$W[0]);
  for ($k=0;$k<=$npk;$k++) {
    if ($na =~ $PKLS[$k]) {
      $NEP[$j] = $W[3] + 3;
    }
  }
  print "$W[0]\t$W[1]\t$W[2]\t$NEP[$j]\t$W[4]\t$W[5]\t$W[6]\t$W[7]\t$W[8]\n";
}
$nj = $j;

close(I);

exit;


