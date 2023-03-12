#!/bin/csh
# File:  run-v7-examples.pl $option infile out-name
# run results examples, 06/07/2022
# run results examples, 02/10/2023

# examples
# for abs
# csh ../Scripts/run-v7-examples-abs.csh abs 340absH-contact-dis+fcs+cdr-fix-$date.txt 340abs-$date.
#	(340absH or 340absL or 340absA)
# for nbs
# csh ../Scripts/run-v7-examples-nbs.csh nbs 83nbs-contact-dis+fcs+cdr-0116.txt 83nbs-$date
#

set per="epitope-statis+fcs+cdr-v7d4.pl"
set opi=$1
set inf=$2
set pid=`date "+%m%d"`
set out=$3"-"${pid}

perl ../Scripts/${per} ${opi} ${inf} 1a > results/${out}-v7-1a.txt
perl ../Scripts/${per} ${opi} ${inf} 1b > results/${out}-v7-1b.txt
perl ../Scripts/${per} ${opi} ${inf} 2a > results/${out}-v7-2a.txt
perl ../Scripts/${per} ${opi} ${inf} 2b > results/${out}-v7-2b.txt
perl ../Scripts/${per} ${opi} ${inf} 3a > results/${out}-v7-3a.txt
perl ../Scripts/${per} ${opi} ${inf} 3b > results/${out}-v7-3b.txt

perl ../Scripts/${per} ${opi} ${inf} 4a all > results/${out}-v7-4a-all.txt
perl ../Scripts/${per} ${opi} ${inf} 4a cdr3 > results/${out}-v7-4a-cdr3.txt
perl ../Scripts/${per} ${opi} ${inf} 4a cdr2 > results/${out}-v7-4a-cdr2.txt
perl ../Scripts/${per} ${opi} ${inf} 4a cdr1 > results/${out}-v7-4a-cdr1.txt
perl ../Scripts/${per} ${opi} ${inf} 4a cdr0 > results/${out}-v7-4a-cdr0.txt

perl ../Scripts/${per} ${opi} ${inf} 4b all > results/${out}-v7-4b-all.txt
perl ../Scripts/${per} ${opi} ${inf} 4b cdr3 > results/${out}-v7-4b-cdr3.txt
perl ../Scripts/${per} ${opi} ${inf} 4b cdr2 > results/${out}-v7-4b-cdr2.txt
perl ../Scripts/${per} ${opi} ${inf} 4b cdr1 > results/${out}-v7-4b-cdr1.txt
perl ../Scripts/${per} ${opi} ${inf} 4b cdr0 > results/${out}-v7-4b-cdr0.txt

# for selections from a file that lists antibody names
perl ../Scripts/${per} ${opi} ${inf} 5a  pickup-list.txt > results/${out}-v7-5a-picklist.txt
perl ../Scripts/${per} ${opi} ${inf} 5b  pickup-list.txt > results/${out}-v7-5b-picklist.txt

perl ../Scripts/${per} ${opi} ${inf} 6u > results/${out}-v7-6u-es.txt
perl ../Scripts/${per} ${opi} ${inf} 6 > results/${out}-v7-6-es.txt

# selections 

# by single name
#perl ../Scripts/${per} ${opi} ${inf} 4a cdr3 Beta-50 > results/${out}-v7-4a-cdr3-Beta50.txt
#perl ../Scripts/${per} ${opi} ${inf} 4a cdr3 Sb68 > results/${out}-v7-4a-cdr3-Sb68.txt

# multiple names separate by ","
#perl ../Scripts/${per} ${opi} ${inf} 4a cdr3 Beta-50,Beta-27 > results/${out}-v7-4a-cdr3-Beta50+27.txt
#perl ../Scripts/${per} ${opi} ${inf} 4a cdr3 Sb16,Sb45,Sb14,Sb68 > results/${out}-v7-4a-cdr3-Sbs.txt

# allowed wildtype, "*", needs quotations
perl ../Scripts/${per} ${opi} ${inf} 4a cdr3 "Beta*" > results/${out}-v7-4a-cdr3-Beta-all.txt
perl ../Scripts/${per} ${opi} ${inf} 4b cdr3 "Beta*" > results/${out}-v7-4b-cdr3-Beta-all.txt
#perl ../Scripts/${per} ${opi} ${inf} 4a cdr3 "Sb*" > results/${out}-v7-4a-cdr3-Sb-all.txt
#perl ../Scripts/${per} ${opi} ${inf} 4b cdr3 "Nb*" > results/${out}-v7-4b-cdr3-Nb-all.txt

