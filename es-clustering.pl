#!/usr/bin/perl -w
# file: "es-clustering-v2.pl";
# Usage: "es-clustering.pl pmode esf6u simcut sizecl"
# "pdbid.txt" is the original name-pdbid-chains file
# "f6u" is the output file from "epitope-statis.pl" version 7.4 with option 6u
# "simcut" is the similarity requested (>)
# "sizecl" is a cluster size cut-off
# 
# Jiansheng Jiang, 12-26-22, v1 
# Jiansheng Jiang, 01-03-23, v2
#
#
use List::MoreUtils qw(uniq);
use String::Similarity;

my $pmode = $ARGV[0];
my $f6u = $ARGV[1];
my $simcut = $ARGV[2]; if (!$simcut) {$simcut = 0.99;}
my $sizecl = $ARGV[3]; if (!$sizecl) {$sizecl = 2;}


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

#
my (@W); my(%PDB,%CTE); 
my ($l,$na,$p,$e,$s,$ies); 
my ($ef,$ee,$foo); 
my ($i,$nj,$k,$k1,$k2,$k3,$kk,$nk) = 0;
my (@NAM,@EP,@LEN,@CL,@NCL,@NM,@ES,@EC,@CLS,@ESS,@CLPDB);
my ($p1,$e1,$ii,$s1,$s2,$si,$str,$nn,$sm,$ssum); 
my (@spl,@unq,@SS); 
my (@SUM,@PK,@MARK,@MRK,@CJ);
#
#
$NAM[0] = "Start"; $EP[0] = "a";
my $j = 0; my $c = 0;
if ($pmode eq "u") {
  open (III,"<$f6u") || die("ERROR: Could not open the file - $f6u ?");
  while ($l = <III>)
  {
    chomp($l);
#  $l =~ s/^\*//;
    next if ($l =~ /^#/);
    $j = $j + 1;
    @W = split /[\t\s]+/, $l;
    $na = $W[0]; $s = $W[1]; $e = $W[2]; $c = $W[3];
    ($na,$foo) = split("\ ",$W[0]);
    $CTE{$na} = $s;
    $NAM[$j] = $na;
    $EP[$j] = $s;  $LEN[$j] = length($s);
    $ES[$j] = $e;  $EC[$j] = $c;
#    print "dbx: na=$na\tEP=$s\tES=$e\tEC=$c\n";
  }
  $nj = $j;
  close (III);
#
}


my (@NI) = sortnamenleng(@EP);
my $lni = scalar(@NI) - 1;
my $n = 0;
for ($i=$lni;$i>1;$i--) {
  $k = $NI[$i];
#  print "i=$i\tk=$k\t$NAM[$k]\t\t$EP[$k]\t$LEN[$k]\n";
  $s1 = $EP[$k]; $n = 0;
  for ($j=$i-1;$j>1;$j--) {
    $k1 = $NI[$j];
    $s2 = $EP[$k1];
    $sim = similarity($s1,$s2); $s = sprintf("%4.2f",$sim);
    if ($s >= $simcut) {
      $n = $n + 1;
      $NS[$k][$n] = $k1;
#      print "j=$j\tk1=$k1\t$s2\ts=$s\t";
    }
  }
  $NSL[$k] = $n; # the size of this cluster
#  print "\n";
}
for ($i=$lni;$i>1;$i--) {
  $k = $NI[$i];
#  print "k=$k\t$EP[$k]\t";
  for ($n=1;$n<=$NSL[$k];$n++) {
    $k1 = $NS[$k][$n];
    $s = sprintf("%13s",$EP[$k1]);
#    print "$s\t";
  }
#  print "\n";
}

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



sub sortnamenleng {
  my @unsorted = @_;
  my @sorted_indexes = sort{
       length($unsorted[$a]) <=> length($unsorted[$b]) || $unsorted[$a] cmp $unsorted[$b]
    } 0..$#unsorted;
  my @sorted_values = @unsorted[ @sorted_indexes ];
  return (@sorted_indexes);
}



sub sortnumb {
  my @unsorted = @_;
  my @sorted_indexes = sort { $unsorted[$a] <=> $unsorted[$b] } 0..$#unsorted;
  my @sorted_values = @unsorted[ @sorted_indexes ];

  return (@sorted_indexes);
}





exit;