#!/usr/bin/perl
# Usage: "reordering.pl infile list outfile"
# 
# Jiansheng Jiang, 03-30-22, v2
#


$ifile = $ARGV[0];
$list = $ARGV[1];
$ofile = $ARGV[2];


my(@ONM); my($l); my($k,$kn) = 0;

open (I,"<$list") || die("ERROR: Could not open the file - $list ?");
while ($l = <I>)
{
  chomp($l);
  $k++;
#  $l =~ s/^\*//;
  next if ($l =~ /^#/);
  $ONM[$k] = $l;  
}
$nk = $k;
close (I);



my(@W); my ($na);
my($i,$j,$nj,$hh,$ll) = 0;
my($pick,$j1,$j2) = 0;

open (O,">$ofile") || die("ERROR: Could not open the file - $onfile ?");

open (I,"<$ifile") || die("ERROR: Could not open the file - $ifile ?");
while ($l = <I>)
{
  chomp($l);
  $j++;
#  $l =~ s/^\*//;
  next if ($l =~ /^#/);
  @W = split("\t",$l);
#  print "@W\n";
  $na = $W[0];
  $pick = 0;
  for ($k=1;$k<=$nk;$k++) {
    if ($na eq $ONM[$k]) {
      print O "$l\n";
      $j2 = $j2 + 1;
      $pick = 1;
    }
  }
  if (!$pick) {
    print "$l\n";
    $j1 = $j1 + 1;
  }
}
$nj = $j;
close(I);
close(O);

print "nj=$nj\tj1=$j1\tj2=$j2\tnk=$nk\n";

exit; 


