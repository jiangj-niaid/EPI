#!/usr/bin/perl
# Usage: "evaluate-dis.pl pdb-list.txt HorL"
# 
# "pdb-list.txt" is the list of pdbid and names of nbs or abs
# Format is: "NAME	PDB	Chain-A	Chain-B"
# Chain-A should be RBD chain
# For Abs, Chain-B might have two chains: both chains (HL), Heavy (H) or Light (L)
# default is H
#
# Jiansheng Jiang, 10-06-21
# Jiansheng Jiang, 03-30-22, v2
# Jiansheng Jiang, 06-21-23, v3
#
# this run "run-cns-dis.csh" with the option of distance cutoff $cut
$runcns = "run-cns-dis-v3.csh";


$file = $ARGV[0];
#($ofile = $file) =~ s/(pdbid\.txt)$//; $ofile = $ofile . "contact-dis.txt";

$HorL = $ARGV[1];
if (!$HorL) {$HorL = "H";}

$cut = $ARGV[2];
if (!$cut) {$cut = 5.0;}

my(@W);
my($i,$j,$nj,$hh,$ll) = 0;
my(@NM,@PD,@PNM,@CA,@CH,@CL) = 0;


open (I,"<$file") || die("ERROR: Could not open the file - $file ?");

while ($l = <I>)
{
  chomp($l);
  $j++;
#  $l =~ s/^\*//;
  next if ($l =~ /^#/);
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

# print out
#for ($j=1;$j<=$nj;$j++) {
#  print "$PDB[$j]\t$CA[$j]\t$CB[$j]\n";
#}

#For each $pdb; (assume a "pdb/" folder in the current working directory)
# start or each $pdb to run CNS
#
for ($j=1;$j<=$nj;$j++) {
  $p = $PD[$j]; 
  $n = $NM[$j]; 
  $a = $CA[$j]; 
  $hh = $CH[$j];
  $ll = $CL[$j];
  if ($HorL =~ "L") {
    print "evaluate distance for $p\t$a\t$ll\n";
    system ("csh ../Scripts/$runcns $p $n $a $ll $cut");
  }  
  if ($HorL =~ "H") {
    print "evaluate distance for $p\t$a\t$hh\n";
    system ("csh ../Scripts/$runcns $p $n $a $hh $cut");
  }
  if ($HorL =~ "HL") {
    print "evaluate distance for $p\t$a\t$hh\n";
    system ("csh ../Scripts/$runcns $p $n $a $hh $cut");
    print "evaluate distance for $p\t$a\t$ll\n";
    system ("csh ../Scripts/$runcns $p $n $a $ll $cut");
  }
}


exit; 


