#!/usr/bin/perl
# Usage: "parse-contact.pl contact.list"
# 
# "contact.list" is the output distance file produced by CNS "contact.inp"
# Jiansheng Jiang, 10-04-21
#

$file = $ARGV[0];
$pdb = $ARGV[1];
$A = $ARGV[2];
$B = $ARGV[3];

open (I,"<$file") || die("ERROR: Could not open the file - $file ?");

my(@W);
my($i,$j,$nj) = 0;
my(@AA,@AI,@BA,@BI,@DD) = 0;

while ($l = <I>)
{
  chomp($l);
  $i++;
  $l =~ s/^\*//;
  if ($l) {
#  print "$l\n";
  if ($i > 14) {
    $l =~ tr/\[\]//d;
    @W = split(" ",$l);
#    print "@W\n";
    $j = $j + 1;
    $AA[$j] = $W[0];
    $AI[$j] = $W[1];
    $BA[$j] = $W[3];
    $BI[$j] = $W[4];
    $DD[$j] = $W[6];
    if ($j > 1) {
      $j1 = $j - 1;
      $aaj = $AA[$j] . $AI[$j];
      $aaj1 = $AA[$j1] . $AI[$j1];
      if ($aaj =~ $aaj1 ) {
        if ($DD[$j] < $DD[$j1]) {
          $AA[$j1] = $AA[$j];
          $AI[$j1] = $AI[$j];
          $BA[$j1] = $BA[$j];
          $BI[$j1] = $BI[$j];
          $DD[$j1] = $DD[$j];
        }
        $j = $j - 1;
      }
    }
  }
  }
  $nj = $j;
}
close(I);

# print out
if ($nj => 1) {
  for ($j=1;$j<=$nj;$j++) {
    print "$pdb\t$A$B\t$AA[$j]\t$AI[$j]\t$BA[$j]\t$BI[$j]\t";
    printf "%4.2f\n",$DD[$j];
  }
}

1; 


