#!/usr/bin/perl
# use warnings;
# use strict;
# Usage: "epitope-statis.pl nbs|abs contact-dis.txt pmode cdr selection"
#
# option: nbs or abs, a default namelist file is readin
#
# "contact-dis.txt" is the list of contact data are computed from CNS
# (with added FCS, CDR information)
#
# "pmode" is the printing mode, "a" includes RBD mutations; "b" exclude RBD mutations
# pmode = 1,  STRUCTURAL ALIGNMENT OF PARATOPES AGAINST EPITOPE (RBD)
# pmode = 2,  STRUCTURAL SEQUENCE ALIGNMENT TABLE
# pmode = 3,  DISTRIBUTION and STATISTICS OF CDR LOOPS ON RBD SURFACE
# pmode = 4,  STATISTICS FOR FCS DISTRIBUTION
# pmode = 5,  print out individual structural alignment
#
# CDR-type: [all, cdr3, cdr2, cdr1, cdr0]
#
# "selection": single-name; multiple-names separated by comma; file-name (line-by-line); 
#              widetype * - needs quotations, i.e. "Beta*".
#
# Jiansheng Jiang, 10-06-21, v1
# Jiansheng Jiang, 04-28-22, v5
# Jiansheng Jiang, 05-24-22, v6, improvements
#   (1) change FCS0 to FCS23; 
#   (2) FCS range are defined; 
#   (3) option: "nbs" and "abs" analysis
#   (4) Statistics of CDR loops, percentages.
#   (5) Bits calculation added, not sure correctness or if it is useful. (Crooks,2004; Schneider,1990)
#   (6) Top residues sorting is added.
#   (7) Identify RBD mutations and output to a "mutation-list" file
#   (8) pmode=4, add sub-group analysis: if "4a" includes mutations; if "4b" exclude mutations
#       (i) type of CDR loops [all, cdr3, cdr2, cdr1, cdr0]; 
#       (ii) a "list", is a filename that lists selected antibody-name line by line
#       (iii) widetype - allowed but needs quotations, i.e. "Beta*".
# Jiansheng Jiang, 06-15-22, v6.1, 
#   FCS in relation with Barne's Classification
#   Default: read in the "-namelist.txt" of antibody or nanobody onto @NAM
#   pmode=5a picklist file, reformating print-out
# Jiansheng Jiang, 07-24-22, v6.2, 
#   pmode=3a, CDR statistics re-formating print-out
# 
# Jiansheng Jiang, 10-23-22, v7
# v7.1 Change FCS numbering move FCS23 to FCS1, and push FCS plus 1 
#      change FCS (Frequently Contacted Site) to ES (Epitopic Site) for printing, 
#      although we keep using FCS arraies.
# v7.3 pmode=6
#      Clustering ES, coding ES 1-23 as a-w, pinting out psudo-sequence in fasta form for alignment
#     
#
# presettings
#
# Three letter to one code for aa
my %ONE = (
	"ALA" => "A", "PHE" => "F", "ILE" => "I", "LEU" => "L", 
	"MET" => "M", "PRO" => "P", "VAL" => "V", "TRP" => "W", 
	"ASP" => "D", "GLU" => "E", "LYS" => "K", "ARG" => "R", "HIS" => "H",  
	"CYS" => "C", "GLY" => "G", "ASN" => "N", 
	"GLN" => "Q", "SER" => "S", "THR" => "T", "TYR" => "Y", 
	);
my %AAO = (
	"A" => 1, "R" => 2, "N" => 3, "D" => 4, "C" => 5, "Q" => 6, "E" => 7, "G" => 8, 
	"H" => 9, "I" => 10, "L" => 11, "K" => 12, "M" => 13, "F" => 14, "P" => 15, 
	"S" => 16, "T" => 17, "W" => 18, "Y" => 19, "V" => 20,
	);
my @AA = (" ", "A", "R", "N", "D", "C", "Q", "E", "G", 
	"H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V");
my @EE = (" ", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", 
   "n", "o", "p", "q", "r", "s", "t", "u", "v", "w");
my %EEO = (
	"a" => 1, "b" => 2, "c" => 3, "d" => 4, "e" => 5, "f" => 6, "g" => 7, "h" => 8, 
	"i" => 9, "j" => 10, "k" => 11, "l" => 12, "m" => 13, "n" => 14, "o" => 15, 
	"p" => 16, "q" => 17, "r" => 18, "s" => 19, "t" => 20, "u" => 21, "v" => 22, "w" => 23);
#
# RBD sequence $start-$end
my $rbdseq  = "ITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVS";
   $rbdseq .= "PTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCV";
   $rbdseq .= "IAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGV";
   $rbdseq .= "EGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKK";
my $fcsrbd  = "---------1-----2------3----------------4------5-----6------------------";
   $fcsrbd .= "---7---------8---9-------------------10----11----12---13-----14---15--";
   $fcsrbd .= "--16-----17--18-19--20----21-----22---------23-----------";
my $clsrbd  = "---------3-----3------3----------------4------4-----4------------------";
   $clsrbd .= "---1---------1---1--------------------2-----2-----2----1------3----2--";
   $clsrbd .= "---1------2---1--2---2-----1------1----------3-----------";
# Variants- mutations are defined
my @Alpha = ("E484K","S494P","N501Y");
my @Beta = ("K417N","E484K","N501Y");
my @Delta = ("K417N","L452R","T478K");
my @Omicron = ("G339D","R346K","S371L","S373P","S375F","T376A","D405N","K417N",
   "N440K","G446S","S477N","T478K","E484A","Q493R","G496S","Q498R","N501Y","Y505H");
my @foo = (@Alpha,@Beta,@Omicron);
my @MUT = get_unique(@foo);
my @BetaList =("Beta-6","Beta-22","Beta-24","Beta-26","Beta-27","Beta-29","Beta-32","Beta-38",
   "Beta-40","Beta-44","Beta-47","Beta-49","Beta-50","Beta-53","Beta-54","Beta-55",
   "BD-667","BD-744","BD-771","BD-813","BD-836","P2C-1F11","FI-3A");
my @OmicronList = ("Cong-8D3");
my @VAR = (@BetaList,@OmicroList);

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
#
my($nfcs) = 23;
my ($start) = 332; my ($end) = 529;
# shift 332 i.e. from 332 to 529
for ($i=0;$i<=197;$i++) {
  $i1 = $i + $start;
  $FCSN[$i1] = $FCSN[$i];
}
my (@RBD) = "";
my ($q0) = $start; # shift 332, i.e. from 332 to 529
for ($i=0;$i<=197;$i++) {
  $q = $i + $q0;
  $s = substr($rbdseq,$i,1);
  $sq = $s . sprintf("%3s",$q);
  $RBD[$q] = $sq;
#  print "$q\t$s\t$RBD[$q]\n";
}
#

# read in the key namelist, nbs-namelist.txt, or abs-namelist.txt
$nbsabs = $ARGV[0]; my($namefile);
if ($nbsabs =~ "nbs" || $nbsabs =~ "abs") {
  $namefile = $nbsabs . "-namelist.txt";
} else {
  print "%WRN: what type? opition= nbs or abs?\n";
  print "Usage: epitope-statis.pl nbs|abs contact-dis.txt pmode cdr selection.\n";
  exit;
}
my (@NAM); my($lnam) = 0;
open (I,"<$namefile") || die("ERROR: Could not open the namelist file - $namefile ?");
chomp(@NAM = <I>);
close (I);
$lnam = scalar(@NAM) - 1;
#print "$nbsabs\t$namefile\t$lnam\n";
#print "NAM=@NAM\n";
my (%NAMK); my($ll);
for ($ll=1;$ll<=$lnam;$ll++) {
  $NAMK{$NAM[$ll]} = $ll;
#  print "$ll\t$NAM[$ll]\t$NK{$NAM[$ll]}\n";
}



$file = $ARGV[1]; # input "raw" data file

$pmode = $ARGV[2];
# "4a" with RBD mutation; "4b" without RBD mutations
$pmodeab = substr($pmode,1,1);
$pmodeab = "a" if (!$pmodeab);

$pname = $ARGV[3]; # CDR loop name 
$psele = $ARGV[4]; # selection or a filename of picklist

if ($pmode =~ "5" or $pmode =~ "6") {
  $psele = $pname;
};

my (@PKLS); my($plist) = ""; my($npls) = 0;
if ($psele) {
  if ($psele =~ /\*/) {
    $psele =~ s/\*//;
    (@PKLS) = grep (/$psele/, @NAM);
  } elsif ($psele =~ /\,/) {
    @PKLS = split(/\,/,$psele);
  } elsif (grep { $psele eq $_ } @NAM) {
    $PKLS[0] = $psele;
  } else {
    $plist = $psele;  
    open (IN,"<$plist") || die("%ERR: Could not open or does not exist - $plist ?");
    while ($l = <IN>) {
      chomp($l);
      next if ($l =~ /^#/);
      @W = split("\t",$l);
      next if (! grep { $W[0] eq $_ } @NAM);
      push (@PKLS, $W[0]);
    }
    close (IN);
  }
  $npls = scalar(@PKLS);
  print "Seleted: n=$npls\t@PKLS\n";
  if ($pmode =~ "6p") {
    $pmode = "6"; $pmodeab = "u";
  }
}

# CDR loop type "cdrlp": [-1=all, 0=cdr0, 1=cdr1, 2=cdr2, 3=cdr3]
my ($cdrlp) = -1; 
if ($pmode ne 5) {    
  if ($pname =~ "cdr3") {
    $cdrlp = 3;
  } elsif ($pname =~ "cdr2") {
    $cdrlp = 2;
  } elsif ($pname =~ "cdr1") {
    $cdrlp = 1;
  } elsif ($pname =~ "cdr0") {
    $cdrlp = 0;
  } else {
    $pname =~ "all";
    $cdrlp = -1;
  }
}

my (@PII) = ""; my($npls1) = $npls - 1;
if ($pmode == 5) {
  for ($k=0;$k<=$npls1;$k++) {
    for ($ii=1;$ii<=$lnam;$ii++) {
      if ($PKLS[$k] eq $NAM[$ii]) {
        $PII[$k] = $ii;
      }
    }
  }
}


my(@W); my(@WW);
my($i,$ni,$j,$nj,$mnj) = 0;
my(@NM,@PNM,@PD,@CA,@CB,@EP,@NEP,@FCS,@PAR,@DIS,@PT) = 0;
my(@NBK) = 0;
my(@PAT) = "";
my(@EPI) = "";
my(@RB) = "";

# learn mutation residues
my ($mutfile) = "mutation-list.txt";
open (O,">$mutfile") || die("ERROR: Could not open the file - $mutfile ?");
my(@foo,@MUTLIST) = "";
open (I,"<$file") || die("ERROR: Could not open the file - $file ?");
while ($l = <I>)
{
  chomp($l);
#  $l =~ s/^\*//;
  next if ($l =~ /^#/);
  $j++;
  @W = split("\t",$l);
  $PNM[$j] = $W[0];
  ($PD[$j],$NM[$j]) = split("_",$W[0]);
  $EP[$j] = $W[2];
  $NEP[$j] = $W[3];
  $RB[$j] = $ONE{$EP[$j]} . $NEP[$j]; 
# identify the mutations on RBD
  my $r = $RBD[$NEP[$j]];
  if ($RB[$j] !~ $r) {
    $s = substr($RB[$j],0,1);
    my $t = $r . $s;
    print O "$l\t$t\n";
    push (@foo,$NM[$j]);
  }  
}
$mnj = $j; $j = 0;
close(I);

@MUTLIST = get_unique(@foo);
print O "MUTLIST=@MUTLIST\n";

# check if the selected entries are the RBD mutations
my (@FLG); my ($la); my($npls1) = $npls - 1;
if ($npls) {
  for ($k=0;$k<=$npls1;$k++) {
    $la = $PKLS[$k]; $FLG[$k] = 1;
    if (grep { $la eq $_ } @MUTLIST) {
      print "%WRN: $la is a RBD mutation antibody or nanobody\n";
      if ($pmodeab eq "b") {
        $FLG[$k] = -1;
      }
    }
  }
}


open (II,"<$file") || die("ERROR: Could not open the file - $file ?");
my ($mut,$sel,$k); my($pdb,$name);
while ($l = <II>)
{
  chomp($l);
#  $l =~ s/^\*//;
  next if ($l =~ /^#/);
  @W = split("\t",$l);
  ($pdb,$name) = split("_",$W[0]);

# exclude known mutations
  $mut = -1;
  foreach $k (@MUTLIST) {
    if ($name =~ $k) {
      $mut = 1;
    } 
  }
  next if ($mut == 1 && $pmodeab eq "b");
#
# pickup selected
  if ($pmode == 5 and $npls) {
    $sel = -1;
    for ($k=0;$k<=$npls1;$k++) {
      if ($name =~ $PKLS[$k]) {
        $sel = 1;
      }
#      $sel = -1 if ($FLG[$k] == -1);
    }
    next if ($sel == -1);
  }
#  print "dbx: $name\tmut=$mut\tsel=$sel\n";
  #PDB_NAME	CHAINS	RBD	Rid	FCS	AB	Rid	DIST	CDR  
  $j++;
  ($PD[$j],$NM[$j]) = ($pdb,$name); # pdbid,name
  $PNM[$j] = $W[0]; # pdbid_name
  $CA[$j] = substr($W[1],0,1); # Chain A - normally is RBD chain
  $CB[$j] = substr($W[1],1,1); # Chain B
  $EP[$j] = $W[2]; # RBD residue name
  $NEP[$j] = $W[3]; # RBD residue id
  $FCS[$j] = $W[4]; # ES site id
  $PAR[$j] = $ONE{$W[5]} . $W[6] . "." . $W[8]; # paratope residue name and id with CDR id
  $DIS[$j] = $W[7]; # distance 
  $CDR[$j] = $W[8]; # CDR loop id
#
  $PT[$j] = $NM[$j] . "_" . $PAR[$j];
# packing the contacts to each name based
  for ($ii=0;$ii<=$lnam;$ii++) {
    if ($NM[$j] eq $NAM[$ii]) {
# NBK is a counter that how many contacts are for each name 
      $NBK[$ii] = $NBK[$ii] + 1;
      $k = $NBK[$ii];
      $WW[$ii][$k] = $l;
      $PAT[$ii][$k] = $PAR[$j];
      $EPI[$ii][$NEP[$j]] = $PAR[$j];
      $ni = $ii;
    }
  }
  $RB[$j] = $ONE{$EP[$j]} . $NEP[$j]; 
  
}
$nj = $j;
print "dbx: nj=$nj\n";
close(II);


my ($p,$f) = 0; my($q);
my (@NC,@NFC) = 0;
my (@NA,@FNA) = "";
for ($p=$start;$p<=$end;$p++) {
  $NC[$p] = 0;
}
for ($f=0;$f<=$nfcs;$f++) {
  $NFC[$f] = 0;
}

# count the number of contacts for each epitope
for ($j=1;$j<=$nj;$j++) {
  $p = $NEP[$j];
  $q = $PD[$j] . "_" . $NM[$j] . "_" . $CB[$j] . "_" . $PAR[$j];
  if ($NM[$j] !~ "ACE2") {
    $NC[$p] = $NC[$p] + 1;
    if ($NC[$p] == 1) {
      $NA[$p] = $q;
    } else {
      $NA[$p] = $NA[$p] . "," . $q;
    }
  }
}
my($sum1) = 0;
for ($p=$start;$p<=$end;$p++) {
  $sum1 = $sum1 + $NC[$p];
}


# count the number of contacts for each FCS
if ($cdrlp == -1) {
  for ($j=1;$j<=$nj;$j++) {
    $p = $NEP[$j]; $f = $FCS[$j]; $k = $CDR[$j]; 
    $q = $PD[$j] . "_" . $NM[$j] . "_" . $CB[$j] . "_" . $PAR[$j];
    if ($NM[$j] !~ "ACE2") {
      $NFC[$f] = $NFC[$f] + 1;
      if ($NFC[$f] == 1) {
        $FNA[$f] = $q;
      } else {
        $FNA[$f] = $FNA[$f] . "," . $q;
      }
    }
  }
} else {
  for ($j=1;$j<=$nj;$j++) {
    $p = $NEP[$j]; $f = $FCS[$j]; $k = $CDR[$j]; 
    $q = $PD[$j] . "_" . $NM[$j] . "_" . $CB[$j] . "_" . $PAR[$j];
    if ($NM[$j] !~ "ACE2") {
      if ($k == $cdrlp) {
        $NFC[$f] = $NFC[$f] + 1;
        if ($NFC[$f] == 1) {
          $FNA[$f] = $q;
        } else {
          $FNA[$f] = $FNA[$f] . "," . $q;
        }
      }    
    }
  }
}
my($sum2) = 0;
for ($f=0;$f<=$nfcs;$f++) {
  $sum2 = $sum2 + $NFC[$f];
}


# count the number of contacts for each CDR loop
my(@NCDR,@CCDR);
for ($f=0;$f<=$nfcs;$f++) {
  for ($i=0;$i<=3;$i++) {
    $NCDR[$f][$i] = 0;
    $CCDR[$f][$i] = "";
  }
}
for ($j=1;$j<=$nj;$j++) {
  $f = $FCS[$j];
  $k = $CDR[$j];
  $q = $PD[$j] . "_" . $NM[$j] . "_" . $CB[$j] . "_" . $PAR[$j];
  if ($NM[$j] !~ "ACE2") {
    $NCDR[$f][$k] = $NCDR[$f][$k] + 1;
    if ($NCDR[$f][$k] == 1) {
      $CCDR[$f][$k] = $q;
    } else {
      $CCDR[$f][$k] = $CCDR[$f][$k] . "," . $q;
    }
  }
}

# print out contacts list

if ($pmode == 1) {
  print "\n";
  print "STRUACTURAL ALIGNMENT OF PARATOPES AGAINST EPITOPE (RBD)";
  print "\n\n";

  for ($j=$start;$j<=$end;$j++) {
    print "$RBD[$j]\t$j\t$NC[$j]\t$NA[$j]\n";
  }
  print "\n";

  print "\n";
  print "STATISTICS OF PARATOPE SEQUENCE ON RBD SURFACE BY NUMBERS";
  print "\n";
  print "RBD\tRESID\tES\tHits\t";
  my ($s) = "";
  for ($i=1;$i<=20;$i++) {
    $s = $s . sprintf("%4s",$AA[$i]);
  }
  print "$s\n";
#
  my(@MM,@SM); my($sum1,$sum2) = 0;
  for ($j=$start;$j<=$end;$j++) {
    for ($i=0;$i<=20;$i++) {
      $MM[$j][$i] = 0;
      $SM[$i] = 0;
    }
  }
#
  for ($j=$start;$j<=$end;$j++) {
    (@W) = split(",",$NA[$j]);
    $kn = $NC[$j]; 
    $sum1 = $sum1 + $kn;
    for ($k=0;$k<=$kn-1;$k++) {
      ($p,$n,$c,$q) = split("_",$W[$k]);
      $r = substr($q,0,1);
      $i = $AAO{$r};
      $MM[$j][$i] = $MM[$j][$i] + 1;
      $SM[$i] = $SM[$i] + 1;
    }
  }
  for ($i=1;$i<=20;$i++) {
    $sum2 = $sum2 + $SM[$i];
  }
#
  for ($j=$start;$j<=$end;$j++) {
    print "$RBD[$j]\t$j\t$FCSN[$j]\t$NC[$j]\t";
    for ($i=1;$i<=20;$i++) {
      $s = sprintf("%4s",$MM[$j][$i]);
      print "$s";
    }
    print "\n";
  }
  # print summary
  print "SUM\t---\t---\t$sum1\t";  
  for ($i=1;$i<=20;$i++) {
    $s = sprintf("%4s",$SM[$i]);
    print "$s";
  }
  print "\n";  
#
  if ($sum2 != $sum1) {
    print "%ERR: SOMETHING WRONG! cross-check sum1=$sum1\tsum2=$sum2\n";
  }
#
   
}



if ($pmode == 2) {
  print "\n";
  print "STRUACTURAL SEQUENCE ALIGNMENT TABLE";
  print "\n";
}

my(@SEQQ) = "";
for ($ll=0;$ll<=$lnam;$ll++) {
  next if (!$NBK[$ll]);
  $seq = "";
  for ($q=$start;$q<=520;$q++) {
    $r = substr($EPI[$ll][$q],0,1);
    if (!$r) {
      $r = "-";
    }
    $seq .= $r;
  }
  if ($pmode == 2) {
    my ($nn) = sprintf("%-12s",$NAM[$ll]);
    print "$nn $seq\n";
  }
  $SEQQ[$ll] = $seq;
}
if ($pmode == 2) {
  my ($nn) = sprintf("%-12s","RBD");
  print "$nn $rbdseq\n";
  my ($nn) = sprintf("%-12s","ES");
  print "$nn $fcsrbd\n";
  my ($nn) = sprintf("%-12s","CLASS");
  print "$nn $clsrbd\n";
}

# make a distribution of CDR loops on RBD surface
if ($pmode == 3) {
  print "\n";
  print "DISTRIBUTION OF CDR LOOPS ON RBD SURFACE";
  print "\n";
}

  my(@CDRR) = "";
  for ($ll=0;$ll<=$lnam;$ll++) {
    next if (!$NBK[$ll]);
    $seq = "";
    for ($q=$start;$q<=520;$q++) {
      $r0 = $EPI[$ll][$q]; $len = length($r0);
      $r = substr($r0,$len-1,1);
      if (!$r) {
        $r = "-";
      }
      $seq .= $r;
    }
    $CDRR[$ll] = $seq;
  }

if ($pmode == 3) {
  for ($ll=0;$ll<=$lnam;$ll++) {
    next if (!$NBK[$ll]);
    $seq = $CDRR[$ll];
    my ($nn) = sprintf("%-12s",$NAM[$ll]);
    print "$nn $seq\n";
  }
  my ($nn) = sprintf("%-12s","RBD");
  print "$nn $rbdseq\n";
  my ($nn) = sprintf("%-12s","ES");
  print "$nn $fcsrbd\n";
  my ($nn) = sprintf("%-12s","CLASS");
  print "$nn $clsrbd\n";
  
# Statistics of CDR loops on RBD surface
  print "\n";
  print "STATISTICS OF CDR LOOPS ON RBD SURFACE - Percentage";
  print "\n";
  print "ES\tRANGE\tRBD\tCLS\tHits\t(%)\t   NONE (%)     CDR1 (%)     CDR2 (%)     CDR3 (%)\n";
  my($sum1) = 0;
  for ($f=0;$f<=$nfcs;$f++) {
    $sum = 0;
    for ($i=0;$i<=3;$i++) {
      $sum = $sum + $NCDR[$f][$i];
    }
    $sum1 = $sum1 + $sum;
  }
  for ($f=0;$f<=$nfcs;$f++) {
    $id = sprintf("%7s",$FCSRNGE[$f]) . " " . sprintf("%-8s",$FCSSEQ[$f]);
    $id = $id . sprintf("%3s",$FCSCLS[$f]); 
    $sum = 0;
    for ($i=0;$i<=3;$i++) {
      $sum = $sum + $NCDR[$f][$i];
    }
    if (!$sum) {
      $p = 0;
    } else {
      $p = $sum/$sum1 * 100;
    }
    $p2 = sprintf("%5.2f",$p);
    print "$f\t$id\t$sum\t$p2\t";
    for ($i=0;$i<=3;$i++) {
      $s = sprintf("%6s",$NCDR[$f][$i]);
      if (!$sum) {
        $p = 0;
      } else {
        $p = $NCDR[$f][$i]/$sum1 * 100;
      }
      $p2 = sprintf("%5.2f",$p);
      print "$s  $p2";
    }
    print "\n";
  }
  print "SUM\t------- --------\t$sum1\t100.00\t";  
  for ($i=0;$i<=3;$i++) {
    $ss = 0;
    for ($f=0;$f<=$nfcs;$f++) {
      $ss = $ss + $NCDR[$f][$i];
    }
    $s = sprintf("%6s",$ss);
    $p = $ss/$sum1 * 100;
    $p2 = sprintf("%5.2f",$p);
    print "$s  $p2";
  }
  print "\n";  

}



if ($pmode == 4) {
  print "\n";
  print "STATISTICS OF PARATOPE SEQUENCE ON ES OF RBD SURFACE BY NUMBERS";
  print "\n";
  print "ES\tRANGE\tRBD\tCLS\tHits\t";
  my ($s) = "";
  for ($i=1;$i<=20;$i++) {
    $s = $s . sprintf("%4s",$AA[$i]);
  }
  my ($t) = sprintf("%20s","TOP");
  print "$s\t$t\n";
#
  my(@MM,@SM); my($sum1,$sum2) = 0;
  for ($f=0;$f<=$nfcs;$f++) {
    for ($i=0;$i<=20;$i++) {
      $MM[$f][$i] = 0;
      $SM[$i] = 0;
    }
  }
#
  for ($f=0;$f<=$nfcs;$f++) {
    (@W) = split(",",$FNA[$f]);
    $kn = $NFC[$f]; 
    $sum1 = $sum1 + $kn;
    for ($k=0;$k<=$kn-1;$k++) {
      ($p,$n,$c,$q) = split("_",$W[$k]);
      $r = substr($q,0,1);
      $i = $AAO{$r};
      $MM[$f][$i] = $MM[$f][$i] + 1;
      $SM[$i] = $SM[$i] + 1;
    }
  }
  for ($i=1;$i<=20;$i++) {
    $sum2 = $sum2 + $SM[$i];
  }
# top sorting
  my(@BB,@TOP);
  for ($f=0;$f<=$nfcs;$f++) {
    for ($i=1;$i<=20;$i++) {
       $BB[$i] = $MM[$f][$i];
    }
    my (@DD) = hightolow(@BB);
    my ($tp) = "";
    for ($i=20;$i>=1;$i--) {
      $sq = $AA[$DD[$i]]; $sq = "-" if (!$BB[$DD[$i]]);
      $tp = $tp . $sq;      
    }
    $TOP[$f] = $tp;
  }  
#
  for ($f=0;$f<=$nfcs;$f++) {
    $id = sprintf("%7s",$FCSRNGE[$f]) . " " . sprintf("%-8s",$FCSSEQ[$f]);
    $id = $id . sprintf("%3s",$FCSCLS[$f]); 
    print "$f\t$id\t$NFC[$f]\t";
    for ($i=1;$i<=20;$i++) {
      $s = sprintf("%4s",$MM[$f][$i]);
      print "$s";
    }
    print "\t$TOP[$f]\n";
  }
# print summary
  print "SUM\t------- --------\t$sum1\t";  
  for ($i=1;$i<=20;$i++) {
    $s = sprintf("%4s",$SM[$i]);
    print "$s";
  }
#
  for ($i=1;$i<=20;$i++) {
    $BB[$i] = $SM[$i];
  }
  my (@DD) = hightolow(@BB);
  my ($tp) = "";
  for ($i=20;$i>=1;$i--) {
    $sq = $AA[$DD[$i]]; $sq = "-" if (!$BB[$DD[$i]]);
    $tp = $tp . $sq;      
  }
  print "\t$tp\n";  
# v7.1 add print %
  print "SUM\t------- --------\t (%) \t";  
  for ($i=1;$i<=20;$i++) {
    $s = sprintf("%4s", int(100*$SM[$i]/$sum1+0.5));
    print "$s";
  }
  print "\n";
#
  if ($sum2 != $sum1) {
    print "%ERR: SOMETHING WRONG! cross-check sum1=$sum1\tsum2=$sum2\n";
  }
#
  print "\n";
  print "STATISTICS OF PARATOPE SEQUENCE ON ES of RBD SURFACE BY BITS";
  print "\n";
  print "ES\tRANGE\tRBD\tCLS\tHits\t";
  my ($s) = "";
  for ($i=1;$i<=20;$i++) {
    $s = $s . sprintf("%5s",$AA[$i]);
  }
  print "$s  Bits\tTOP3\n";
#
  my($mbits) = 4.32; my($b,$sb); my(@BB) = 0.0;
  for ($f=0;$f<=$nfcs;$f++) {
    $id = sprintf("%7s",$FCSRNGE[$f]) . " " . sprintf("%-8s",$FCSSEQ[$f]);
    $id = $id . sprintf("%3s",$FCSCLS[$f]); 
    print "$f\t$id\t$NFC[$f]\t";
    $sb = $mbits;
    for ($i=1;$i<=20;$i++) {
      $n = $MM[$f][$i];
      $b = bits($sum2,$n);
      $s = $sb * $b; 
      $BB[$i] = $s;
      $sb = $sb - $b;
    }
    # for low statistical number, apply a correction factor
    if ($NFC[$f] < 5) {
      $sb = 0;
    } elsif ($NFC[$f] < 10) {
      $sb = $sb / 16;
    } elsif ($NFC[$f] < 20) {
      $sb = $sb / 4;
    }
    for ($i=1;$i<=20;$i++) {
      $s = $BB[$i];
      $s = sprintf("%5.2f",$s);
      print "$s";
    }
    $sb = sprintf("%5.2f",$sb);
    print "$sb\t";
    my ($m1,$m2,$m3) = max3(@BB);
    $s3 = $AA[$m3]; $s3 = "-" if (!$m3);
    $s2 = $AA[$m2]; $s2 = "-" if (!$m2);
    $s1 = $AA[$m1]; $s1 = "-" if (!$m1);
    print "$s3$s2$s1\n"; 
  }
  # print summary
  print "SUM\t------- --------\t$sum1\t";  
  $sb = $mbits;
  for ($i=1;$i<=20;$i++) {
    $n = $SM[$i];
    $b = bits($sum2,$n);
    $sb = $sb - $b;
  }
   for ($i=1;$i<=20;$i++) {
    $n = $SM[$i];
    $b = bits($sum2,$n);
    $s = $sb * $b;
    $s = sprintf("%5.2f",$s);
    print "$s";
  }
  $sb = sprintf("%5.2f",$sb);
  print "$sb\t";  
    my ($m1,$m2,$m3) = max3(@SM);
    $s3 = $AA[$m3]; $s3 = "-" if (!$m3);
    $s2 = $AA[$m2]; $s2 = "-" if (!$m2);
    $s1 = $AA[$m1]; $s1 = "-" if (!$m1);
    print "$s3$s2$s1\n"; 

}


if ($pmode == 5 && $npls > 0) {
  print "\n";
  print "STRUCTURAL ALIGNMENT FOR $pname";
  print "\n";

  for ($k=0;$k<=$npls-1;$k++) {

    $ll = $PII[$k];
    my ($nn) = sprintf("%-12s",$NAM[$ll]);
    print "$nn $SEQQ[$ll]\n";
    $seq = "";
    for ($q=$start;$q<=520;$q++) {
      $r = substr($EPI[$ll][$q],0,1);
      if (!$r) {
        $r = "-";
      } elsif ($r) {
        $f = substr($RBD[$q],0,1);
        $r = sprintf("%s",$f);
      }
      $seq .= $r;
    }
    $nn = sprintf("%-12s","RBD");
    print "$nn\t$seq\n";
    $nn = sprintf("%-12s","CDR");
    print "$nn\t$CDRR[$ll]\n";
  }
  my ($nn) = sprintf("%-12s","ES");
  print "$nn\t$fcsrbd\n";
  my ($nn) = sprintf("%-12s","CLASS");
  print "$nn\t$clsrbd\n";
# print the raw data
  print "\n";
  print "#PDB_NAME\tCHAINS\tRBD\tRid\tES\tAbNb\tRid\tDIST\tCDR\n";
  for ($k=0;$k<=$npls-1;$k++) {
    $ll = $PII[$k];
    for ($kk=1;$kk<=$NBK[$ll];$kk++) {
      print "$WW[$ll][$kk]\n";
    }
  }
  print "\n";
}

# 
# 
# clustering ES
# pmode = 6 for printing redundant ES
# pmode = 6u for printing uniq ES
#

my (%EGC,%EGN) = "";
my (@ES,@EP,@EK,@EC); my($e,$p,$ps,$pc); my($i,$j,$k,$f,$f0,$k0) = 0;

if ($pmode == 6) {
# clean
  for ($k=0;$k<=$lnam;$k++) {
    $ES[$k] = ""; $EP[$k] = ""; $EC[$k] = "";
    for ($i=1;$i<=23;$i++){
      $EK[$k][$i] = 0;
    }
  }
# coding EP numbers to alphabeta ES: EP=(1-23); ES=(a-w)
  $k0 = 0;
  for ($j=1;$j<=$nj;$j++) {
    $k = $NAMK{$NM[$j]}; $e = "";
    if ($k != $k0) {
      $f0 = 0; $k0 = $k;
    }
    if ($FCS[$j] >= 1) {
      $f = $FCS[$j];  $EK[$k][$f] = $EK[$k][$f] + 1;
      $e = $EE[$f];
      if ($pmodeab eq "u") {
        if ($f != $f0) {
          $EP[$k] = $EP[$k] . "," . $f;
          $f0 = $f;
          $ES[$k] = $ES[$k] . $e;
        } 
      } else {
          $EP[$k] = $EP[$k] . "," . $f;
          $ES[$k] = $ES[$k] . $e;
      }
    }
  }
#
# frequency EC
  for ($k=1;$k<=$lnam;$k++) {
    for ($i=1;$i<=23;$i++){
      $e1 = $EK[$k][$i];
      if ($e1 >= 1) {
        $EC[$k] = $EC[$k] . "," . $EK[$k][$i];      
      }
    }
  } 
# print out as psudo "sequence" for alignment
#
  if ($pmodeab eq "u") {
    print "#Name\t\t\tES(a-w)\tES(1-23)\tN-contacts\n";
  } else {
    print "#Name\t\tES(a-w)\tES(1-23)\n";
  }
  for ($k=1;$k<=$lnam;$k++) { 
    $name = $NAM[$k];
    $nn = sprintf("%-18s",$NAM[$k]);
    if ($pmodeab eq "u") {
      ($ps = $EP[$k]) =~ s/^,//;
      $ps = "(" . $ps . ")";
      ($pc = $EC[$k]) =~ s/^,//;
      $pc = "(" . $pc . ")";
      print "$nn\t$ES[$k]\t$ps\t$pc\n";
      $EGC{$name} = $name . "\t" . $ES[$k] . "\t" . $pc;
    } else {
      print "$nn\t\t$ES[$k]\t$ps\n";
    }
  }
# Print out (ES) numbers for each 
#  print "#\n";  
#  if ($pmodeab eq "u") {
#    print "#Name\t\t\tES(1-23)\tN-Contacts\n";
#  } else {
#    print "#Name\t\t\tES(1-23)\n";
#  }
  for ($k=1;$k<=$lnam;$k++) {
    $name = $NAM[$k];
    $nn = sprintf("%-18s",$NAM[$k]);
    ($p = $EP[$k]) =~ s/^,//;
    $p = "(" . $p . ")";
    if ($pmodeab eq "u") {
      ($pc = $EC[$k]) =~ s/^,//;
      $pc = "(" . $pc . ")";
#      print "$nn\t$p\t$pc\n";
      $EGN{$name} = $name . "\t" . $p . "\t" . $pc;
    } else {
#      print "$nn\t$p\n";
    }
  }
#
# pickup list
  print "\n";
  print "#picked: n=$npls\t@PKLS\n";
  if ($npls) {
    foreach $name (@PKLS) {
      ($s,$e,$pc) = split("\t",$EGC{$name});
      ($s,$n,$pc) = split("\t",$EGN{$name});
      $nn = sprintf("%-18s",$s);
      print "$nn\t$e\t$n\t$pc\n"; 
    }     
  }

}





exit; 

# below are the subroutine or functions

# bits = - p * log_2(p), p is probability, or frequency
sub bits {
  my ($M,$N) = @_;
  if (!$N) {return $N};
  my $p = $N/$M;
  my $l2p = log($p)/log(2);
  my $bit = -$p * $l2p;
  return $bit;
}


sub maxi {
  my (@BB) = @_;
  my ($max,$mi) = 0;
  for ($i=1;$i<=20;$i++) {
    if ($BB[$i] >= $max) {
      $max = $BB[$i];
      $mi = $i; 
    }
  }
  return ($max,$mi);
}

sub max3 {
  my (@BB) = @_;
  my ($max) = 0.001;
  my ($m1,$m2,$m3) = 0;
  for ($i=1;$i<=20;$i++) {
    if ($BB[$i] >= $max) {
      $max = $BB[$i];
      $m1 = $i;
    }
  }
  $max = 0.001;
  for ($i=1;$i<=20;$i++) {
    next if ($i == $m1);
    if ($BB[$i] >= $max) {
      $max = $BB[$i];
      $m2 = $i;
    }
  }
  $max = 0.001;
  for ($i=1;$i<=20;$i++) {
    next if ($i == $m1 || $i == $m2);
    if ($BB[$i] >= $max) {
      $max = $BB[$i];
      $m3 = $i;
    }
  }
  return ($m1,$m2,$m3);
}

# specific sorting subroutine for ES
# input an array, output an indexed array
sub hightolow {
  my (@BB) = @_; my(@DD); my($nbb) = scalar(@BB) - 1;
  my ($maxi,$mi) = 0;
# index DD ties to the fraction position of BB, assume $nbb < 100;
  for ($i=1;$i<=$nbb;$i++) {
    $DD[$i] = $i; $BB[$i] = sprintf("%5.2f",$BB[$i]+$i/100);
  }  
#  print "bb+dd: @BB\n";
  for ($i=1;$i<=$nbb;$i++) {
    $maxi = $BB[$i];
    for ($j=$i+1;$j<=20;$j++) {
      if ($BB[$j] > $maxi) {
        $maxi = $BB[$j];
        $mi = $j;
      }
    }
    if ($maxi == $BB[$i]) {
    #
    } else {
      $BB[$mi] = $BB[$i];
      $BB[$i] = $maxi;
    }
  }
  for ($i=1;$i<=$nbb;$i++) {
    $b = int($BB[$i]);  $d = ($BB[$i] - $b) * 100 + 0.0001;
    $d = sprintf("%2d",$d);
    $DD[$i] = sprintf("%2s",$d);
    $BB[$i] = sprintf("%2s",$b);
  }
  return (@DD);
}

# usage: my @unique = get_unique(@array);
sub get_unique {
    my %seen;
    grep !$seen{$_}++, @_;
}
