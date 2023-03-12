#!/usr/bin/perl
# Usage: "pick-extract-namedpdb.pl pdbid.list namelist cluster=number"
#
# read all pdb files from the "pdb/" that presents in "*-pdbid.txt"
# pick up the named pdb file that presents in "-namelist.txt" 
#
# (1) write all pdb files to directory "pdb1/"
# (2) extract RBD domain ($A chain), and Variable H-chain ($B chain)
#     default RBD from 328-529; default Abs from 1-128
#     write out pdb file with "-ex.pdb".
# (3) write out associate Pymol script file, "pymol-load.pml"
#
# Jiansheng Jiang, 07-19-22
#
#
# ES definitions
my @ESSEQ = ("--------","GEV","NTRFAS","YAWNRKR","LYNSASF","STFKCY","GVSPTK","RGDEVRQ",
    "APGQTGK","DYNYKLPDD","NSNNLDS","KVGGNY","NYLYR","LFRKSN","KPFERD","ISTEI","YQAGSTP",
	"NGVE","GFN","CYFP","LQSYG","FQPTNG","VGYQPYR","ELLHA");
#my @ESRNGE = ("-------","339:341","343:349","351:357","368:374","375:380","381:386","403:409",
#	"411:417","420:428","437:443","444:449","450:454","455:460","462:467","468:472","473:479",
#	"481:484","485:487","488:491","492:496","497:502","503:509","516:520");
my @ESRNGE = ("-------","339:341","343:349","351:357","368:374","375:380","381:386","403:409",
	"415:417","420:422","437:443","444:449","450:454","455:460","462:467","468:472","473:479",
	"481:484","485:487","488:491","492:496","497:502","503:509","516:520");


my $pdbid = $ARGV[0]; # the file contains both pdbid and names
my $namefile = $ARGV[1]; # the names of antibody/nanobody to be picked up 
# specify which cluster to be pick
my ($fo,$cluster,$type);
if (!$ARGV[2]) {
  $cluster = 0;
  print "cluster=?\n";
  exit;
} else {
  ($fo,$cluster) = split("=",$ARGV[2]);
  if ($cluster =~ /^[1-9]/) {
    $type = "bynumber"; # (list of clusters)
  } else {
    $type = "byname"; # (list of individual names)
  }
}
print "dbx: cluster=$cluster\ttype=$type\n";

my (@NAM,@NM,@PD,@CA,@CB,@EST,@PDB,@ES,@EP,@NN) = "";
my ($cl,$sz,$e,$s,$list,$n,$na,$p);
my ($na,$A,$B,$cha,$chb,$b1,$b2,$bb);

if ($type =~ /bynumber/ && $cluster >= 1) {

  open (II,"<$namefile") || die("ERROR: Could not open the namelist file - $namefile ?");
  while ($l = <II>)
  {
    chomp($l);
    next if ($l =~ /^#/);
    ($cl,$sz,$e,$s,$list) = split("\t",$l);
    if ($cl eq $cluster) {
      $list =~ s/\"+//g;
      @NAM = split(",",$list);
#      print "@NAM\n";
#      print "$cl\t$sz\t$ep\t$es\t$list\n";
      $s =~ s/\"+//g; $s =~ s/\(//g; $s =~ s/\)//g;
      @EST = split(",",$s);
#      print "@EST\n";
    }
  }
}

if ($type =~ /byname/) {
  my ($j) = -1;
  open (II,"<$namefile") || die("ERROR: Could not open the namelist file - $namefile ?");
  while ($l = <II>)
  {
    chomp($l);
    next if ($l =~ /^#/);
    ($n,$na,$p,$e,$s) = split("\t",$l);
     print "$n\t$na\t$p\t$e\t$s\n";
    $j = $j + 1;
    $NAM[$j] = $na; $NN[$j] = $n;
    $PDB[$j] = $p;
    $ES[$j] = $e;
    $EP[$j] = $s;    
  }
  $nj = $j;
  $s = $ES[0]; 
  $s =~ s/\"+//g; $s =~ s/\(//g; $s =~ s/\)//g;
  @EST = split(",",$s);
#  print "@EST\n";
}  
close (II);

my($lnam) = scalar(@NAM) - 1;

$namefile =~ s/.txt//g;
my $pymol = "./pdb1/" . $namefile . "-cl-" . $cluster . ".pml";

open (O,">$pymol") || die("ERROR: Could not open the file - $pymol ?");
print O "/cmd.set('bg_rgb',0,'',0)\n";

$j = 0;
open (I,"<$pdbid") || die("ERROR: Could not open the file - $pdbid ?");
while ($l = <I>)
{
  chomp($l);
#  $l =~ s/^\*//;
#  print "\n" if ($l !~ /\S/);
  next if ($l =~ /^#/);
  ($name,$pid,$A,$B) = split("\t",$l);

  for ($i=0;$i<=$lnam;$i++) {
    next if ($NAM[$i] ne $name);
    if ($NAM[$i] eq $name) {
      $j = $j + 1;
      $pdb = $pid . ".pdb";
      system("cp ./pdb/$pdb ./pdb1/$pdb");
      $NM[$j] = $name;
      $PD[$j] = $pid; 
      $CA[$j] = $A;
      ($b1,$b2) = split("\,",$B); $bb = substr($b1,0,1);
      $CB[$j] = $bb;
      $pex = $pid . "-" . $A . $bb; $po = $pex . ".pdb"; 
      $pnam = $pid . "-" . $name; 
      if ($j == 1) {
        $pnam1 = $pnam; $a1 = $A;
        $rbd = "rbd"; $rbdpdb = $pid . "-rbd" . ".pdb";
        print O "\n";
        print O "load $rbdpdb,$rbd\n";
      }
      print O "\n";
      print O "load $po,$pnam\n";
      print O "align ($pnam and chain $A),($pnam1 and chain $a1)\n";
      print O "color gray, ($pnam and chain $A)\n";
      next;
    }
  }

}
close(I);
$nj = $j;

print O "\ncenter\n";
print O "color gray, rbd\n";
print O "show surface, rbd\n";
# paiting ES footprint
foreach $k (@EST) {
  $e = "ES" . $k;
  $r = $ESRNGE[$k];
  print O "select $e, (rbd and resid $r)\n";
  print O "color magenta, $e\n";
}
close O;

# default begin and end of residue numbers of RBD, as "A" chain
my($sa) = 328; my($ea) = 520;
# default begin and end of residue numbers of Antibody, as "A" chain
my($sb) = 1; my($eb) = 118;

for ($j=1;$j<=$nj;$j++) {

  $pid = $PD[$j];
  $pdb = $pid . ".pdb"; $pi = "./pdb1/" . $pdb; 
  $cha = $CA[$j];
  $chb = $CB[$j];
  $po = "./pdb1/" . $pid . "-" . $cha . $chb . ".pdb"; 

# open output pdb file
  open(O,">$po") || die("cannot open file $po");

# extract A chain
  open(I,"<$pi") || die("cannot open file $pi");
  while ($line = <I>) {
    chomp($line);
    if ($line =~ /^ATOM/) {
      $chn = substr($line,21,1);
      $res = substr($line,23,4);
      if ($chn eq $cha) {
        if ($res >= $sa && $res <= $ea) {
          print O "$line\n";
        }
      }
    }
  }
  close (I);
  print O "TER\n";
  if ($j == 1) {
    $rbdpdb = "./pdb1/" .  $pid . "-rbd" . ".pdb";
    system("cp $po $rbdpdb");
  }

# extract B chain
  open(I,"<$pi") || die("cannot open file $pi");
  while ($line = <I>) {
    chomp($line);
    if ($line =~ /^ATOM/) {
      $chn = substr($line,21,1);
      $res = substr($line,23,4);
      if ($chn eq $chb) {
        if ($res >= $sb && $res <= $eb) {
          print O "$line\n";
        }
      }
    }
  }
  close (I);
  print O "TER\nEND\n";
  close (O);
}



exit;
