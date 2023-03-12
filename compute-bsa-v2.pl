#!/usr/bin/perl
#use warnings;
#use strict;
# file: compute-bsa.pl nbs|abs pdbid.txt HorL
# 
# Usage: "compute-bsa.pl nbs|abs pdbid.txt HorL"
# 
# Jiansheng Jiang, 12-18-22, v1
#


my $mode = $ARGV[0];
my $file = $ARGV[1];
my $HorL = $ARGV[2];
if (!$HorL) {$HorL = "H";}

my ($l,$hh,$ll) = ""; my (@W,@NM,@PD,@CA,@CH,@CL,@PNM);
my ($j, $nj) = 0;
open (I,"<$file") || die("ERROR: Could not open the file - $file ?");

while ($l = <I>)
{
  chomp($l);
  next if ($l =~ /^#/);
  $j++;
#  $l =~ s/^\*//;
  @W = split("\t",$l);
#  print "@W\n";
  $NM[$j] = $W[0];
  $PD[$j] = $W[1];
  $CA[$j] = $W[2];
  ($hh,$ll) = split(",",$W[3]);
  $CH[$j] = substr($hh,0,1);
  $CL[$j] = substr($ll,0,1);
  $PNM[$j] = $PD[$j] . "_" . $NM[$j];
#  print "$j\t$PNM[$j]\t$CA[$j]\t$CH[$j]\t$CL[$j]\n";
}
$nj = $j;

close(I);


my ($p,$n,$a); 
my ($ab,$pdb,$pdbab,$pbsa,$bsa,$bs); my ($tmp)="./tmp/";
print "#Name\tPDB\tChains\tBSA\n";
for ($j=1;$j<=$nj;$j++) {
  $p = $PD[$j];
  $n = $NM[$j]; 
  $a = $CA[$j]; 
  $hh = $CH[$j]; 
  $ll = $CL[$j];
  $pdb = "./pdb/" . $p . ".pdb";
  if ($HorL eq "H") {
    $ab = "-" . $a . "," . $hh;
    $pdbab = "./pdb/" . $p . "_" . $a . $hh . ".pdb";
  } else {
    $ab = "-" . $a . "," . $ll;
    $pdbab = "./pdb/" . $p . "_" . $a . $ll . ".pdb";
  }
  system ("/usr/local/bin/pdb_selchain $ab $pdb | egrep -v 'HETATM|ANISOU' > $pdbab");
  $pbsa = $tmp . $p . "_" . $n . "_pisa.txt";
  system ("csh ../../Scripts/pisa-bsa.csh $pdbab > $pbsa");
  sleep 1;
#
  $bsa = 0;
  open(my $si,"<",$pbsa) or die $!;
  chomp(my @lines = <$si>);
  close $si;
  print @lines;
  $l = grep("Buried", @lines);
  $bsa = substr($l,50,7) + 0.5;
  $bs = sprintf("%6d",$bsa);
  $ab =~ s/^-//;
  print "$n\t$p\t$ab\t$bs\n";
}


exit;