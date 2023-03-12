#!/usr/bin/perl -w
# file: "es-clustering-sel.pl";
# Usage: "es-clustering-sel.pl esf6u simcut selection weight"
# "pdbid.txt" is the original name-pdbid-chains file
# "f6u" is the output file from "epitope-statis.pl" version 7.3 with option 6u
# "simcut" is the similarity requested (>)
# "selection" ES=name|pdbid|classid|4,5,6,..
# 
# Jiansheng Jiang, 12-26-22, v1 
# Jiansheng Jiang, 01-03-23, v2
#
#
use List::MoreUtils qw(uniq);
use String::Similarity;
#
my @EE = (" ", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", 
   "n", "o", "p", "q", "r", "s", "t", "u", "v", "w");
my %EEO = (
	"a" => 1, "b" => 2, "c" => 3, "d" => 4, "e" => 5, "f" => 6, "g" => 7, "h" => 8, 
	"i" => 9, "j" => 10, "k" => 11, "l" => 12, "m" => 13, "n" => 14, "o" => 15, 
	"p" => 16, "q" => 17, "r" => 18, "s" => 19, "t" => 20, "u" => 21, "v" => 22, "w" => 23);
my @CLAS1 = (7,8,9,13,16,18,21,22); my $CLASS1 = "ghimpruv";
my @CLAS2 = (10,11,12,15,17,19,20); my $CLASS2 = "jkloqst";
my @CLAS3 = (1,2,3,14,23); my $CLASS3 = "abcnw";
my @CLAS4 = (4,5,6); my $CLASS4 = "def";
my @RBM = (8,11,12,13,16,18,19,20,21,22); my $RBMS = "hklmprstuv";

#print "dbx: @CLAS1\t@CLAS2\t@CLAS3\t@CLAS4\n";
#print "dbx: $CLASS1\t$CLASS2\t$CLASS3\t$CLASS4\n";
#print "dbx: @RBM\t$RBMS\n";
#
#my $pmode = $ARGV[0];
#my $pmode = $ARGV[0];
my $pdbid = $ARGV[0];
my $f6u = $ARGV[1]; 
my $simcut = $ARGV[2]; if (!$simcut) {$simcut = 0.99;}
#print "dbx: $pdbid\t$f6u\t$simcut\n";

#
my ($l,$name,$pid,$A,$B,$b1,$b2,$bb,$nnam); 
my (@NM,@PD,@CA,@CB);
my (@NAM,@EP,@ES,@EC,@SMI);
my (@W); my(%CTE); 
my ($na,$foo,$nes,$s,$e,$c); 
my (@EST);

my ($i,$nj,$k,$k1,$k2,$k3,$kk,$nk) = 0;
#my (@SUM,@PK,@MARK,@MRK,@CJ);
#
#
# read pdbid list file
my($j) = 0;
open (I,"<$pdbid") || die("ERROR: Could not open the file - $pdbid ?");
while ($l = <I>)
{
  chomp($l);
  next if ($l =~ /^#/);
  ($name,$pid,$A,$B) = split("\t",$l);
  $j = $j + 1;
  $NM[$j] = $name;
  $PD[$j] = $pid; 
  $CA[$j] = $A;
  ($b1,$b2) = split("\,",$B); $bb = substr($b1,0,1);
  $CB[$j] = $bb;
}
close(I);
$nnam = $j;
#print "dbx: nnam=$nnam\n";

#
# read ES f6u file
$NAM[0] = "Start"; $EP[0] = "a";
$j = 0; $c = 0;
open (III,"<$f6u") || die("ERROR: Could not open the file - $f6u ?");
while ($l = <III>)
{
  chomp($l);
  next if ($l =~ /^#/);
  $j = $j + 1;
  @W = split /[\t\s]+/, $l;
  $na = $W[0]; $s = $W[1]; $e = $W[2]; $c = $W[3];
  ($na,$foo) = split("\ ",$W[0]);
  $CTE{$na} = $s;
  $NAM[$j] = $na;
  $EP[$j] = $s; $ES[$j] = $e;  $EC[$j] = $c;
#    print "dbx: na=$na\tEP=$s\tES=$e\tEC=$c\n";
}
$nes = $j;
close (III);
#print "dbx: nes=$nes\n";

#
# parse inputs
# input ES=(),  (name | pdbid | class | ,,,)
my $sele = $ARGV[3];
#print "dbx: $sele\n";
if (!$sele) {
  print "What selections are you looking for?\n";
  print "sel=(name | pdbid | class1 | class2 | class3 | class4 | rbm ) or ES by ,,,\n";
  exit;
}
#
my $weight = $ARGV[4];
if ($weight) {
   $weight =~ s/W=//g;
}

#
my $EPQ = ""; my $find = 0;
if ($sele =~ /ES=/) {
  $sele =~ s/ES\=//g;
}
$sele =~ s/\"+//g; $sele =~ s/\(//g; $sele =~ s/\)//g;
#
if ($sele =~ /,/) {
# ES combination with "," separated
  @EST = split(",",$sele);
  foreach $k (@EST) {
    $s = $EE[$k];
    $EPQ = $EPQ . $s; 
  }
  $find = 1;  
} elsif ($sele =~ "class1") {
  $EPQ = $CLASS1; $find = 1;
} elsif ($sele =~ "class2") {
  $EPQ = $CLASS2; $find = 1;
} elsif ($sele =~ "class3") {
  $EPQ = $CLASS3; $find = 1;
} elsif ($sele =~ "class4") {
  $EPQ = $CLASS4; $find = 1;
} elsif ($sele =~ "rbm") {
  $EPQ = $RBMS; $find = 1;
} else {
# search NAME list and PDBid list
  for ($j=1;$j<=$nnam;$j++) {
    if ($sele eq $NM[$j] or lc($sele) eq $PD[$j]) {
      $EPQ = $EP[$j]; $find = 1;
    }
  }
}
if (!$find) {
  print "ERR: cannot find '$sele' ?\n";
  print "sel=(name | pdbid | class1 | class2 | class3 | class4 | rbm ) or (,,,)\n";
  exit;
} else {
#  my @ESS = convertEPtoES($EPQ);
  my @spl = split("", $EPQ);
  my $leng = scalar(@spl) - 1; my (@epitopes);
  for ($i=0;$i<=$leng;$i++) {
    $epitopes[$i] = $EEO{$spl[$i]};
  }
  $s = join ",", @epitopes;
  if ($sele =~ /,/) {
    print "#ES=($s)\tEP=($EPQ)\n";
  } else {
    print "#$sele\tES=($s)\tEP=($EPQ)\n";
  }
}


# start check if any one is "similar" to $EPQ
my (@NS,@NSL);
my (@NI) = sortnamenleng(@EP);
my $lni = scalar(@NI) - 1;
my $n = 0; my ($s2,$sim);
my $s1 = $EPQ; $k = 1;
for ($i=$lni;$i>1;$i--) {
#  print "i=$i\tk=$k\t$NAM[$k]\t\t$EP[$k]\t$LEN[$k]\n";
  $k1 = $NI[$i];
  $s2 = $EP[$k1];
  $sim = similarity($s1,$s2); $s = sprintf("%4.2f",$sim);
  $SMI[$k1] = $s;
  if ($s >= $simcut) {
    $n = $n + 1;
    $NS[$k][$n] = $k1;
#      print "j=$j\tk1=$k1\t$s2\ts=$s\t";
  }
}
# where $k is the cluster number, alway=1 here for selection
$NSL[$k] = $n; # the size of this cluster
print "#found $n similar ES with similarity=$simcut\n";
#
if (!$weight) {
for ($n=1;$n<=$NSL[$k];$n++) {
  $k1 = $NS[$k][$n];
  $s = sprintf("%15s",$EP[$k1]);
  print "$n\t$NAM[$k1]\t$PD[$k1]\t$SMI[$k1]\t$ES[$k1]\t$EC[$k1]\n";
}
print "\n";
}

# weighting scheme
#
my (@efoo,@cfoo) = ""; my ($n1,$w1) = 0;
my (@KN,@WT) = 0;

if ($weight) {
  for ($n=1;$n<=$NSL[$k];$n++) {
    $k1 = $NS[$k][$n]; $n1 = $n - 1; $KN[$n1] = $k1; $WT[$n1] = 0; 
    $e = $ES[$k1]; $e =~ s/\(//g; $e =~ s/\)//g;
    @efoo = split(",",$e); $li = scalar(@efoo) - 1;
    $c = $EC[$k1]; $c =~ s/\(//g; $c =~ s/\)//g;
    @cfoo = split(",",$c); 
    for ($i=0;$i<=$li;$i++) {
      if ($efoo[$i] eq $weight) {
        $WT[$n1] = $cfoo[$i];
      }
    }
#    print "$n\t$NAM[$k1]\t$PD[$k1]\t$SMI[$k1]\t$ES[$k1]\t$EC[$k1]\t$WT[$n1]\n";
  }
  my (@STI) = sortnumb(@WT);
#  print "WT=@WT\n";
my $len = scalar(@WT) - 1;
for ($i=0;$i<=$len;$i++) {
#  print "dbx: i=$i\tKN=$KN[$i]\tWT=$WT[$i]\t$STI[$i]\n";
  $s1 = $STI[$i]; $w1 = $WT[$s1]; $k1 = $KN[$s1]; $i1 = $i + 1;
#  print "$i1\t$NAM[$k1]\t$PD[$k1]\t$SMI[$k1]\t$ES[$k1]\t$EC[$k1]\t$w1\n";
  print "$i1\t$NAM[$k1]\t$PD[$k1]\t$SMI[$k1]\t$ES[$k1]\t$EC[$k1]\t$w1\n";
}
#  print "STI=@STI\n";
  
}


exit;




# recognize a cluster - auto from top similarity to down
# the size of a cluster should be greater than $sizecl
# pick the first one, then check if the next is inside this cluster
#
# hash array %CLUSTER{'name'}
# index array $CLU[$m][$n],  m is the number of clusters, n is the list member in this cluster
#
my ($ncl) = 0;
my (@FLAG,@NMC,@CLUSTER,@CIDX,@CEP,@CES,@CSZ);
#
# first turn the FLAG on
for ($i=$lni;$i>1;$i--) {
  $k = $NI[$i];
  if ($NSL[$k] > $sizecl) {
    $FLAG[$k] = 1;
  } else {
    $FLAG[$k] = 0;
  }
}
#
#
$c = 0; $CSZ[0] = 0;
for ($i=$lni;$i>1;$i--) {
  $k = $NI[$i];
  if ($FLAG[$k] == 1) {
    $na = $NAM[$k];
    $c = $c + 1; $NMC[$c] = $na;
    $CLUSTER[$c] = $na; $CIDX[$c] = $k; $CSZ[$c] = $NSL[$k];
    $CEP[$c] = "(" . $EP[$k] . ")"; 
    $CES[$c] = $ES[$k]; 
#    print "dbx: c=$c\tk=$k\t$EP[$k]\t$na\n";
    for ($n=1;$n<=$NSL[$k];$n++) {
      $k1 = $NS[$k][$n];
      $CLUSTER[$c] =  $CLUSTER[$c] . "," . $NAM[$k1];
      $CIDX[$c] =  $CIDX[$c] . "," . $k1;
      for ($i1=$i-1;$i1>=1;$i1--) {
        $k2 = $NI[$i1];
        if ($FLAG[$k2]) {
          if ($k2 eq $k1) {
            $FLAG[$k2] = 0;
          }
        }
      }    
    }
  }
}
$ncl = $c;


my (@CNI) = sortnumb(@CSZ);
my ($lcni) = scalar(@CNI) - 1;

print "#Total number of clusters = $ncl (size > $sizecl) with the similarity = $simcut\n";
print "#CLUSTER\tSIZE\tEP\tES\tLIST\n";

my $sum = 0; $k = 0; my ($h) = 0;
for ($i=$lcni;$i>=1;$i--) {
  $h = $CNI[$i];
  $k = $k + 1;
  $sum = $sum + $CSZ[$h];
  print "$k\t$CSZ[$h]\t$CEP[$h]\t$CES[$h]";
  print "\t$CLUSTER[$h]\n";
}
print "#Total\t$sum\n";


exit;


sub convertEPtoES {
  my $string = @_;
  my @spl = split("", $string);
  my $leng = scalar(@spl) - 1; 
  my (@epitopes);
  for ($i=0;$i<=$leng;$i++) {
    $epitopes[$i] = $EEO{$spl[$i]};
  }
  return (@epitopes);
}

sub sortnamenleng {
  my @unsorted = @_;
  my @sorted_indexes = sort{
       length($unsorted[$a]) <=> length($unsorted[$b]) || $unsorted[$a] cmp $unsorted[$b]
    } 0..$#unsorted;
  my @sorted_values = @unsorted[ @sorted_indexes ];
  return (@sorted_indexes);
}


# sorted from $b to $a from large to small
sub sortnumb {
  my @unsorted = @_;
  my @sorted_indexes = sort { $unsorted[$b] <=> $unsorted[$a] } 0..$#unsorted;
  my @sorted_values = @unsorted[ @sorted_indexes ];

  return (@sorted_indexes);
}





exit;