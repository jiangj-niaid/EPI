#!/usr/bin/perl

use Text::Levenshtein qw(distance); 
# if cannot locate "Levenshtein", using the following command to install it
# cpan install Text::Levenshtein


# similarity between different clusters
# for 83nbs 10 clusters
my (@ST) = ("hilmprs",
      "hilmpqrst","himoprst","ghilmprstv","klmqrstu","bjkluv",
      "bcklmnoqst","bckloqstu","hmprstv","bkloqrst","giklmpqrstuv");
my (@SM) = 0.0;
my ($lst) = scalar(@ST) - 1;
my ($i,$j,$s1,$s2,$sim);
my ($simcut) = 0.30;
#my ($limit) = 0.90;
print "@ST\n";
#print "length=$lst\n";

print "   \t\t\t";
for ($j=1;$j<=$lst;$j++) {
  $s2 = sprintf("%-15s",$ST[$j]);
  print "$s2 ";
}
print "\n";

for ($i=0;$i<=$lst;$i++) {
#  $s1 = $ST[$i];
  $s1 = sprintf("%12s",$ST[$i]);
  print "i=$i\t$s1\t";
  for ($j=1;$j<=$lst;$j++) {
    $s2 = sprintf("%12s",$ST[$j]);
    $sim = similarity($s1,$s2); $s = sprintf("%4.2f",$sim);
    if ($s > $simcut) {
      print "$s    \t";
    }
  }
  print "\n";
}


exit;

