#!/usr/bin/perl
# Usage: "evaluate-asa-es.pl key getarea.txt|cns-asa.txt"
# 
# key = "GA" or "CNS"
# "getarea.txt" is the ASA value for individual residue computed by 
# https://curie.utmb.edu/getarea.html
# Fraczkiewicz, R. and Braun, W. (1998) 
# J. Comp. Chem., 19, 319-333.
# or
# "cns-asa.txt" is computed by CNS 1.3 - a list of ASA for individual residue
#
# To estimated ASA for each each ES range (2-9 residues):
# Sum ( lamda * ASAj ), where j from 1 to n, n is the number of residues in the range
# where lamda is an efficiency factor, i.e. 95% (take 0.5% away of overlap area in neighboring residues)
# 
# 
# Jiansheng Jiang, 10-26-22, 12-25-22
#
# FCS definitions
my @FCSCLS = (0,3,3,3,4,4,4,1,0,0,2,2,2,1,3,2,1,2,1,2,2,1,1,0);
my @FCSSEQ = ("--------","GEV","NTRFAS","YAWNRKR","LYNSASF","STFKCY","GVSPTK","RGDEVRQ",
    "APGQTGK","DYNYKLPDD","NSNNLDS","KVGGNY","NYLYR","LFRKSN","KPFERD","ISTEI","YQAGSTP",
	"NGVE","GFN","CYFP","LQSYG","FQPTNG","VGYQPYR","ELLHA");
my @FCSRNGE = ("-------","339-341","343-349","351-357","368-374","375-380","381-386","403-409",
	"411-417","420-428","437-443","444-449","450-454","455-460","462-467","468-472","473-479",
	"481-484","485-487","488-491","492-496","497-502","503-509","516-520");
my @FCSL = (0,3,6,7,7,6,6,7,7,9,7,6,5,6,6,5,7,4,3,4,5,6,7,5);
my @FCSN = (0,0,0,0,0,0,0,1,1,1,
	0,2,2,2,2,2,2,2,0,3,3,3,3,3,3,3,0,0,0,0,0,0,0,0,0,0,4,4,4,4,4,4,4,5,5,5,5,5,5,6,6,
	6,6,6,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,7,7,7,7,7,7,0,8,8,8,8,8,8,8,0,0,9,9,9,
	9,9,9,9,9,9,0,0,0,0,0,0,0,0,10,10,10,10,10,10,10,11,11,11,11,11,11,12,12,12,12,12,13,13,
	13,13,13,13,0,14,14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,16,16,0,17,17,17,17,
	18,18,18,19,19,19,19,20,20,20,20,20,21,21,21,21,21,21,
	22,22,22,22,22,22,22,0,0,0,0,0,0,23,23,23,23,23,0,0,0,0,0,0,0,0,0);
my @FCSCLR = ("--------","olive","aquamarine","forestgreen","limon","marine","hotpink","brightorange",
    "lightblue","cyan","salmon","slate","green","orange","wheat","lightmagenta","pink",
	"palegreen","palecyan","splitpea","red","teal","yellow","gray50");
#
# shift 332 i.e. from 332 to 529
for ($i=0;$i<=197;$i++) {
  $i1 = $i + 332;
  $FCSN[$i1] = $FCSN[$i];
}
#
my @ASA = 0;

# key
my $key = $ARGV[0];
# input file
my $file = $ARGV[1];

my(@W);
my($i,$j,$nj) = 0;
my(@RN,@ID,@TT) = 0;

open (I,"<$file") || die("ERROR: Could not open the file - $file ?");
while ($l = <I>)
{
  chomp($l);
  next if ($l =~ /^#/);
  $j++;
  $l =~ s/^\ //;
  if ($key =~ "GA") {
    @W = split("\t",$l);
    $ID[$j] = $W[0]; # residue id/name
    $RN[$j] = $W[1]; # residue number
    $TT[$j] = $W[2]; # total area of ASA
  }
  if ($key =~ "CNS") {
    @W = split("\ ",$l);
    $ID[$j] = $W[2];
    $RN[$j] = $W[1];
    $TT[$j] = $W[3];  
  }  
}
$nj = $j;

close(I);


# print out
#
my (@ASA) = 0;
for ($j=1;$j<=$nj;$j++) {
  $i = $RN[$j]; $ASA[$i] = $TT[$j];
#  print "$RN[$j]\t$ID[$j]\t$TT[$j]\n";
}

# evaluate ASA for each ES
my ($lam) = 1.0; my($s,$e,$rng,$sq) = ""; my($sum,$sm) = 0;
print "#ES\tRANGE\t RBD\t\tASA\n";
for ($k=1;$k<=23;$k++) {
  ($s,$e) = split("-",$FCSRNGE[$k]);
  $sum = $lam * $ASA[$s];
  for ($i=$s+1;$i<=$e;$i++) {
    $sum = $sum + $lam * $ASA[$i];
  }
  $sm = sprintf("%5d",$sum);
  $rng = $s . "-" . $e;
  $sq = sprintf("%-9s",$FCSSEQ[$k]);
  print "$k\t$rng\t $sq\t$sm\n";  
}
print "\n";


# print the script for pymol visualization
for ($k=1;$k<=23;$k++) {
  ($s,$e) = split("-",$FCSRNGE[$k]);
  my($rng) = $FCSRNGE[$k];
  $rng = $s . ":" . $e;
  $fcs = "ES" . $k;
  $clr = $FCSCLR[$k];
#  print "select $fcs, (chain A and resid $rng)\n";
#  print "color $clr, $fcs\n";
}

# print the script for pymol visualization
for ($k=1;$k<=23;$k++) {
  ($s,$e) = split("-",$FCSRNGE[$k]);
  my($rng) = $FCSRNGE[$k];
  $rng = $s . ":" . $e;
  $fcs = "ES" . $k;
#  print "label A/$s/CA, \"$k\"\n";
}


1; 


