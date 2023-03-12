#!/usr/bin/perl
# Usage: "pick-named-sequences.pl sequence namelist"
#
# read in all sequences that presents in the file of  "*-sequence-aligned-merge.txt"
# pick up the named sequences that presents in "-namelist.txt"
#
# Jiansheng Jiang, 07-19-22
#
#


$seqid = $ARGV[0]; # the file contains both aligned sequences and names
$namefile = $ARGV[1]; # the names of antibody/nanobody to be picked up 

my (@NAM,@NM,@PD,@CA,@CB,@SQ,@LN);
my ($lnam,$pdb,$name,$seq,$len,$A,$B,$cha,$chb,$b1,$b2,$bb);
my ($d1) = "--------------------------"; # 26x -
my ($d2) = "----------------"; # 16x -
my ($d3) = "--------------------------------------"; # 38x -


open (II,"<$namefile") || die("ERROR: Could not open the namelist file - $namefile ?");
chomp(@NAM = <II>);
close (II);
$lnam = scalar(@NAM) - 1;
#print "@NAM\n"; 

#my($j) = 0;
open (I,"<$seqid") || die("ERROR: Could not open the file - $seqid ?");
while ($l = <I>)
{
  chomp($l);
  if ($l =~ /^#/) {
    print "$l\n";
    next;
  }
  ($name,$seq,$len) = split("\t",$l);
  $name =~ s/\t//g; $name =~ s/\s//g;
  ($c1,$c2,$c3) = split("\,",$seq);
#  $sq = $d1 . $c1 . $d2 . $c2 . $d3 . $c3;
  $sq = $c1 . $d2 . $c2 . $d3 . $c3;
  for ($i=0;$i<=$lnam;$i++) {
    next if ($NAM[$i] ne $name);
    if ($NAM[$i] eq $name) {
      $NM[$i] = $name;
      $SQ[$i] = $sq;
      $LN[$i] = $len;
#      $nn = sprintf("%-12s",$name);
#      print "$nn\t$sq\n";
      next;
    }
  }

}
close(I);
#$nj = $j;

# sorting to the list order
for ($i=0;$i<=$lnam;$i++) {
  $nn = sprintf("%-12s",$NM[$i]);
#  print "$nn\t$SQ[$i]\t$LN[$i]\n";
  print "$nn\t$SQ[$i]\n";
}





exit;
